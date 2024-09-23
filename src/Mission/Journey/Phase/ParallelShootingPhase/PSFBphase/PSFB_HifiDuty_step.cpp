//parallel shooting finite burn (PSFB) step for EMTGv9
//High-fidelity duty cycle
//Jacob Englander 3-19-2018

#include "doubleType.h"

#include "PSFB_HifiDuty_step.h"

#include "PropagatorFactory.h"
#include "IntegrationSchemeFactory.h"

namespace EMTG
{
    namespace Phases
    {
        PSFB_HifiDuty_step::PSFB_HifiDuty_step() 
        {
            this->StateAfterThrustInertial.resize(10, 1, 0.0);
        };

        PSFB_HifiDuty_step::PSFB_HifiDuty_step(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stepIndex,
            const size_t& stageIndex,
            ParallelShootingStep* previousStep,
            ParallelShootingPhase* myPhase,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) : PSFB_HifiDuty_step()
        {
            this->initialize(name,
                journeyIndex,
                phaseIndex,
                stepIndex,
                stageIndex,
                previousStep,
                myPhase,
                myUniverse,
                mySpacecraft,
                myOptions);
        }//end constructor

        void PSFB_HifiDuty_step::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stepIndex,
            const size_t& stageIndex,
            ParallelShootingStep* previousStep,
            ParallelShootingPhase* myPhase,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            //base class
            this->PSFBstep::initialize(name,
                journeyIndex,
                phaseIndex,
                stepIndex,
                stageIndex,
                previousStep, //typecast to pointer to base, safe because only base methods will be called
                myPhase, //typecast to pointer to base, safe because only base methods will be called
                myUniverse,
                mySpacecraft,
                myOptions);

            //need to shrink the duty cycle slightly otherwise we get no coast at all
            if (this->StepDutyCycle > 1.0 - 1.0e-10)
                this->StepDutyCycle -= 1.0e-10;
            else if (this->StepDutyCycle < 1.0e-10)
                this->StepDutyCycle += 1.0e-10;

            this->dStepTime_dPhaseFlightTime = this->StepDutyCycle / this->myPhase->get_num_steps() / this->num_interior_control_points; //for the thrust arc
            this->dStepTime_dPhaseFlightTime_Coast = (1.0 - this->StepDutyCycle) / this->myPhase->get_num_steps(); //for the coast arc
        }//end initialize
        
        PSFB_HifiDuty_step::~PSFB_HifiDuty_step()
        {
            //acceleration model
            delete this->myCoastSpacecraftAccelerationModel;

            //propagator
            delete this->myCoastPropagator;

            //integration scheme
            delete this->myCoastIntegrationScheme;
        }//end destructor

        //configure propagator
        void PSFB_HifiDuty_step::configure_propagator()
        {
            //base class
            this->PSFBstep::configure_propagator();

            this->ThrustSTM.resize(this->num_interior_control_points, math::Matrix<double>(14, 14, math::identity));
            this->ThrustAugmentedSTM.resize(this->num_interior_control_points, math::Matrix<double>(14, 14, math::identity));
            this->ThrustdPropagatedStatedIndependentVariable.resize(this->num_interior_control_points, math::Matrix<double>(10, 2, 0.0));

            this->CoastSTM.resize(14, 14, math::identity);
            this->CoastAugmentedSTM.resize(14, 14, math::identity);
            this->CoastdPropagatedStatedIndependentVariable.resize(10, 2, 0.0);

            this->CumulativeAugmentedSTM.resize(this->num_interior_control_points + 1, math::Matrix<double>(14, 14, math::identity)); //one extra entry for the coast
            
            //reconfigure the thrust propagator
            this->myPropagators.back().setStateRight(this->StateAfterThrustInertial);
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                this->myPropagators[subStepIndex].setSTM(this->ThrustSTM[subStepIndex]);
                this->myPropagators[subStepIndex].setdStatedIndependentVariable(this->ThrustdPropagatedStatedIndependentVariable[subStepIndex]);
            }
            this->mySpacecraftAccelerationModel->setDutyCycle(1.0);

            //coast acceleration model
            this->myCoastSpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(this->myOptions,
                this->myJourneyOptions,
                this->myUniverse,
                this->Xdescriptions,
                this->mySpacecraft,
                14); //STM size
            this->myCoastSpacecraftAccelerationModel->setDutyCycle(0.0);

            //coast EOM
            this->myCoastEOM.setSpacecraftAccelerationModel(this->myCoastSpacecraftAccelerationModel);

            //coast integration scheme
            this->myCoastIntegrationScheme = CreateIntegrationScheme(&this->myCoastEOM, 10, 14);
            
            //coast propagator
            this->myCoastPropagator = CreatePropagator(this->myOptions,
                this->myUniverse,
                10,
                14,
                this->StateAfterThrustInertial,
                this->StateStepRightInertial,
                this->CoastSTM,
                this->CoastdPropagatedStatedIndependentVariable,
                (Integration::Integrand*) &this->myEOM,
                this->myCoastIntegrationScheme,
                &this->dStepTime_dPhaseFlightTime_Coast,
                    this->myJourneyOptions->override_integration_step_size
                        ? this->myJourneyOptions->integration_step_size
                        : this->myOptions->integration_time_step_size);
        }//end configure propagator

        //master calcbounds
        void PSFB_HifiDuty_step::calcbounds_step()
        {
            this->calcbounds_step_left_state();

            this->calcbounds_step_control();

            this->calcbounds_step_left_match_point_constraints();

            this->calcbounds_step_main();

            this->calcbounds_distance_constraints();

            this->calcbounds_maneuver_constraints();
        }//end calcbounds_step

        //master process
        void PSFB_HifiDuty_step::process_step(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->total_number_of_states_to_integrate = needG
                ? 9 + 12 * 12
                : 9;

            this->process_step_left_state(X, Xindex, F, Findex, G, needG);

            this->process_step_control(X, Xindex, F, Findex, G, needG);

            this->process_step_left_match_point_constraints(X, Xindex, F, Findex, G, needG);

            this->process_step_main(X, Xindex, F, Findex, G, needG);

            this->process_distance_constraints(X, Xindex, F, Findex, G, needG);

            this->process_maneuver_constraints(X, Xindex, F, Findex, G, needG);

            if (needG)
                this->process_derivative_tuples(X, Xindex, F, Findex, G, needG);
        }//end process_step

        void PSFB_HifiDuty_step::process_substep_times(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->StepThrustTime = (this->PhaseFlightTime - this->myPhase->getInitialCoastDuration() - this->myPhase->getTerminalCoastDuration()) * this->dStepTime_dPhaseFlightTime * this->num_interior_control_points;
            this->subStepThrustTime = this->StepThrustTime / this->num_interior_control_points;
            this->StepCoastTime = (this->PhaseFlightTime - this->myPhase->getInitialCoastDuration() - this->myPhase->getTerminalCoastDuration()) * this->dStepTime_dPhaseFlightTime_Coast;
            this->StepFlightTime /= this->StepDutyCycle; //because it's currently scaled by duty cycle
        }//end process_epoch_time()

        //process step main - aka propagate
        void PSFB_HifiDuty_step::process_step_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: get substep times
            this->process_substep_times(X, Xindex, F, Findex, G, needG);

            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                //Step 1: propagate
                if (subStepIndex == 0)
                {
                    this->ThrustdPropagatedStatedIndependentVariable[subStepIndex].assign_zeros();
                    this->myPropagators[subStepIndex].setCurrentEpoch(this->StateStepLeftInertial(7));
                    this->myPropagators[subStepIndex].setCurrentIndependentVariable(this->StateStepLeftInertial(7));
                }
                else
                {
                    this->ThrustdPropagatedStatedIndependentVariable[subStepIndex] = this->ThrustdPropagatedStatedIndependentVariable[subStepIndex - 1];
                    this->myPropagators[subStepIndex].setCurrentEpoch(this->StateAfterSubStepInertial[subStepIndex - 1](7));
                    this->myPropagators[subStepIndex].setCurrentIndependentVariable(this->StateAfterSubStepInertial[subStepIndex - 1](7));
                }
                this->myPropagators[subStepIndex].setIndexOfEpochInStateVec(7);
                this->myPropagators[subStepIndex].propagate(this->subStepThrustTime, this->ControlVector[subStepIndex], needG);

                //Step 2: update the epoch
                if (subStepIndex == 0)
                    this->StateAfterSubStepInertial[subStepIndex](7) = this->StateStepLeftInertial(7) + this->subStepThrustTime;
                else if (subStepIndex == this->num_interior_control_points - 1)
                    this->StateAfterThrustInertial(7) = this->StateAfterSubStepInertial[subStepIndex - 1](7) + this->subStepThrustTime;
                else
                    this->StateAfterSubStepInertial[subStepIndex](7) = this->StateAfterSubStepInertial[subStepIndex - 1](7) + this->subStepThrustTime;

                //Step 3: derivatives
                if (needG)
                {
                    //Step 3.1: form the ThrustAugmented STM
                    //upper right 7x7 is the original STM
                    for (size_t i = 0; i < 7; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                            this->ThrustAugmentedSTM[subStepIndex](i, j) = this->ThrustSTM[subStepIndex](i, j);

                        //Phi_t terms
                        this->ThrustAugmentedSTM[subStepIndex](i, 7) = this->ThrustSTM[subStepIndex](i, 7);
                        this->ThrustAugmentedSTM[subStepIndex](i, 13) = this->ThrustSTM[subStepIndex](i, 13);

                        //control terms
                        this->ThrustAugmentedSTM[subStepIndex](i, 10) = this->ThrustSTM[subStepIndex](i, 10);
                        this->ThrustAugmentedSTM[subStepIndex](i, 11) = this->ThrustSTM[subStepIndex](i, 11);
                        this->ThrustAugmentedSTM[subStepIndex](i, 12) = this->ThrustSTM[subStepIndex](i, 12);
                    }

                    //tanks
                    for (size_t i = 8; i < 10; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                            this->ThrustAugmentedSTM[subStepIndex](i, j) = this->ThrustSTM[subStepIndex](i, j);

                        //Phi_t terms
                        this->ThrustAugmentedSTM[subStepIndex](i, 7) = this->ThrustSTM[subStepIndex](i, 7);
                        this->ThrustAugmentedSTM[subStepIndex](i, 13) = this->ThrustSTM[subStepIndex](i, 13);

                        //control terms
                        this->ThrustAugmentedSTM[subStepIndex](i, 10) = this->ThrustSTM[subStepIndex](i, 10);
                        this->ThrustAugmentedSTM[subStepIndex](i, 11) = this->ThrustSTM[subStepIndex](i, 11);
                        this->ThrustAugmentedSTM[subStepIndex](i, 12) = this->ThrustSTM[subStepIndex](i, 12);
                    }


                    //epoch time with respect to propagation time
                    this->ThrustAugmentedSTM[subStepIndex](7, 13) = this->dStepTime_dPhaseFlightTime;
                }//end STM creation
            }//end loop over substeps

            //Step 4: propagate to the end of the coast arc
            //Step 4.1: propagate
            this->CoastdPropagatedStatedIndependentVariable = this->ThrustdPropagatedStatedIndependentVariable.back();
            this->myCoastPropagator->setCurrentEpoch(this->StateAfterThrustInertial(7));
            this->myCoastPropagator->setIndexOfEpochInStateVec(7);
            this->myCoastPropagator->setCurrentIndependentVariable(this->StateAfterThrustInertial(7));
            this->myCoastPropagator->propagate(this->StepCoastTime, needG);

            //Step 4.2: coast derivatives
            if (needG)
            {
                //upper right 7x7 is the original STM
                for (size_t i = 0; i < 7; ++i)
                {
                    for (size_t j = 0; j < 7; ++j)
                        this->CoastAugmentedSTM(i, j) = this->CoastSTM(i, j);

                    //Phi_t terms
                    this->CoastAugmentedSTM(i, 7) = this->CoastSTM(i, 7);
                    this->CoastAugmentedSTM(i, 13) = this->CoastSTM(i, 13);
                    
                    //control terms
                    this->CoastAugmentedSTM(i, 10) = 0.0;
                    this->CoastAugmentedSTM(i, 11) = 0.0;
                    this->CoastAugmentedSTM(i, 12) = 0.0;
                }

                //tanks
                for (size_t i = 8; i < 10; ++i)
                {
                    for (size_t j = 0; j < 7; ++j)
                        this->CoastAugmentedSTM(i, j) = this->CoastSTM(i, j);

                    //Phi_t terms
                    this->CoastAugmentedSTM(i, 7) = this->CoastSTM(i, 7);
                    this->CoastAugmentedSTM(i, 13) = this->CoastSTM(i, 13);

                    //control terms
                    this->CoastAugmentedSTM(i, 10) = 0.0;
                    this->CoastAugmentedSTM(i, 11) = 0.0;
                    this->CoastAugmentedSTM(i, 12) = 0.0;
                }
                
                //epoch time with respect to propagation time
                this->CoastAugmentedSTM(7, 13) = this->dStepTime_dPhaseFlightTime_Coast;

                //Step 5: construct cumulative STM chains
                //Step 5.1 cumulative thrust STM chains
                this->CumulativeAugmentedSTM.back() = this->CoastAugmentedSTM;
                this->CumulativeAugmentedSTM[this->num_interior_control_points - 1] = this->CoastAugmentedSTM * this->ThrustAugmentedSTM.back();

                for (int subStepIndex = this->num_interior_control_points - 1; subStepIndex >= 0; --subStepIndex)
                {
                    //create stripped version of next step's cumulative STM to remove that substep's control influence
                    math::Matrix<double> StrippedSTM = this->CumulativeAugmentedSTM[subStepIndex + 1];

                    for (size_t i = 0; i < 10; ++i)
                    {
                        for (size_t j = 10; j < 10 + this->num_controls; ++j)
                            StrippedSTM(i, j) = 0.0;
                    }

                    this->CumulativeAugmentedSTM[subStepIndex] = StrippedSTM * this->ThrustAugmentedSTM[subStepIndex];
                }
            }

        }//end process_step_main

        void PSFB_HifiDuty_step::output(std::ofstream& outputfile,
            size_t& eventcount)
        {
            //Step 1: print the thrust arc
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                //Step 1.1: temporarily redirect the propagator to the output state
                this->myPropagators[subStepIndex].setStateRight(this->output_state);

                //Step 1.2: propagate
                if (subStepIndex == 0)
                {
                    this->ThrustdPropagatedStatedIndependentVariable[subStepIndex].assign_zeros();
                    this->myPropagators[subStepIndex].setCurrentEpoch(this->StateStepLeftInertial(7));
                    this->myPropagators[subStepIndex].setCurrentIndependentVariable(this->StateStepLeftInertial(7));
                }
                else
                {
                    this->ThrustdPropagatedStatedIndependentVariable[subStepIndex] = this->dPropagatedStatedIndependentVariable[subStepIndex - 1];
                    this->myPropagators[subStepIndex].setCurrentEpoch(this->StateAfterSubStepInertial[subStepIndex - 1](7));
                    this->myPropagators[subStepIndex].setCurrentIndependentVariable(this->StateAfterSubStepInertial[subStepIndex - 1](7));
                }
                this->myPropagators[subStepIndex].setIndexOfEpochInStateVec(7);
                this->myPropagators[subStepIndex].propagate(this->subStepThrustTime / 2.0, this->ControlVector[subStepIndex], false);

                //Step 1.4: redirect the propagator back to where it is supposed to go
                this->myPropagators[subStepIndex].setStateRight(subStepIndex == this->num_interior_control_points - 1 ? this->StateAfterThrustInertial : this->StateAfterSubStepInertial[subStepIndex]);

                //Step 1.5: figure out spacecrafty things

                //Step 1.5.1: where am I relative to the sun?
                math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                if (this->myUniverse->central_body_SPICE_ID == 10)
                {
                    R_sc_Sun = output_state.getSubMatrix1D(0, 2);
                }
                else
                {
                    //where is the central body relative to the sun?
                    doubleType central_body_state_and_derivatives[12];
                    this->myUniverse->locate_central_body(output_state(7),
                        central_body_state_and_derivatives,
                        *this->myOptions,
                        false);

                    math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                    }

                    R_sc_Sun = output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                }

                //Step 1.5.2: call the power model
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, output_state(7));

                //Step 1.5.3: call the thruster model
                if (this->num_controls == 4)
                    this->mySpacecraft->computeElectricPropulsionPerformance(this->StepDutyCycle, this->ControlVector[subStepIndex](3));
                else
                    this->mySpacecraft->computeElectricPropulsionPerformance(this->StepDutyCycle);

                //Step 1.5.4: store the thruster model outputs
                doubleType max_thrust = this->mySpacecraft->getEPthrust() * 1.0e-3;
                doubleType max_mass_flow_rate = this->mySpacecraft->getEPMassFlowRate();
                doubleType Isp = this->mySpacecraft->getEPIsp();
                doubleType power = this->mySpacecraft->getAvailablePower();
                doubleType active_power = this->mySpacecraft->getEPActivePower();
                size_t number_of_active_engines = this->mySpacecraft->getEPNumberOfActiveThrusters();
                size_t ThrottleLevel = this->mySpacecraft->getEPThrottleLevel();

                std::string event_type;
                if (max_thrust > 5.0e-7 && this->throttle[subStepIndex] > 5.0e-4)
                    event_type = "PSFBthrust";
                else
                    event_type = "coast";

                math::Matrix<doubleType> ThrustVector = this->ControlVector[subStepIndex] * max_thrust * 1000.0 / this->StepDutyCycle;
                math::Matrix<doubleType> deltaV = ThrustVector * this->subStepThrustTime / this->output_state(6);

                this->writey_thing::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    event_type,//event_type
                    "deep-space",//event_location
                    this->subStepThrustTime / 86400.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    output_state,//state
                    deltaV,//dV
                    ThrustVector,//ThrustVector
                    this->throttle[subStepIndex],//throttle
                    max_thrust * 1000.0 / this->StepDutyCycle,//Thrust
                    Isp,//Isp
                    power,//AvailPower
                    max_mass_flow_rate / this->StepDutyCycle,//mdot
                    number_of_active_engines,//number_of_active_engines
                    active_power,
                    this->mySpacecraft->getEPThrottleLevelString());//active_power)
            }//end loop over substeps

            //Step 2: print the coast arc
            {
                //Step 2.1: temporarily redirect the coast propagator to the output state
                this->myCoastPropagator->setStateRight(this->output_state);

                //Step 2.2: propagate
                this->CoastdPropagatedStatedIndependentVariable.assign_zeros();
                this->myCoastPropagator->setCurrentEpoch(this->StateAfterThrustInertial(7));
                this->myCoastPropagator->setIndexOfEpochInStateVec(7);
                this->myCoastPropagator->setCurrentIndependentVariable(this->StateAfterThrustInertial(7));
                this->myCoastPropagator->propagate(this->StepCoastTime / 2.0, false);

                //Step 2.3: redirect the propagator back to where it is supposed to go
                this->myCoastPropagator->setStateRight(this->StateStepRightInertial);

                //Step 2.4: figure out spacecrafty things

                //Step 2.4.1: where am I relative to the sun?
                math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                if (this->myUniverse->central_body_SPICE_ID == 10)
                {
                    R_sc_Sun = output_state.getSubMatrix1D(0, 2);
                }
                else
                {
                    //where is the central body relative to the sun?
                    doubleType central_body_state_and_derivatives[12];
                    this->myUniverse->locate_central_body(output_state(7),
                        central_body_state_and_derivatives,
                        *this->myOptions,
                        false);

                    math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                    }

                    R_sc_Sun = output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                }
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, output_state(7));
                
                doubleType power = this->mySpacecraft->getAvailablePower();
                doubleType active_power = this->mySpacecraft->getEPActivePower();

                //Step 2.5: print
                this->writey_thing::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "nav-coast",//event_type
                    "deep-space",//event_location
                    this->StepCoastTime / 86400.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    output_state,//state
                    math::Matrix<doubleType>(3, 1, 0.0),//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    0.0,//throttle
                    0,//Thrust
                    0,//Isp
                    power,//AvailPower
                    0,//mdot
                    0,//number_of_active_engines
                    active_power,
                    "none");
            }
        }//end output

        void PSFB_HifiDuty_step::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
        {
			size_t ManeuverThrottleLevel = -1;
			double divisor = 1.0;
			
            //print the thrust arc
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                //Step 1: temporarily assign the propagator to the output state
                this->myPropagators[subStepIndex].setStateRight(this->output_state);

                //Step 2: propagate and print the thrust arc, skipping the first entry
                double timeToPropagate = .1;
                double totalPropagationTime = this->subStepThrustTime _GETVALUE;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 2.1: propagate
                    if (subStepIndex == 0)
                    {
                        this->ThrustdPropagatedStatedIndependentVariable[subStepIndex].assign_zeros();
                        this->myPropagators[subStepIndex].setCurrentEpoch(this->StateStepLeftInertial(7));
                        this->myPropagators[subStepIndex].setCurrentIndependentVariable(this->StateStepLeftInertial(7));
                    }
                    else
                    {
                        this->ThrustdPropagatedStatedIndependentVariable[subStepIndex] = this->dPropagatedStatedIndependentVariable[subStepIndex - 1];
                        this->myPropagators[subStepIndex].setCurrentEpoch(this->StateAfterSubStepInertial[subStepIndex - 1](7));
                        this->myPropagators[subStepIndex].setCurrentIndependentVariable(this->StateAfterSubStepInertial[subStepIndex - 1](7));
                    }
                    this->myPropagators[subStepIndex].setIndexOfEpochInStateVec(7);
                    this->myPropagators[subStepIndex].propagate(timeToPropagate, this->ControlVector[subStepIndex], false);

                    this->temp_state = this->output_state;

                    //Step 2.2: convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //Step 2.3: print
                    if ((totalPropagationTime - timeToPropagate) > 0)
                    {
                        //we need an instantaneous power/propulsion state
                        //Step 2.3.1: where am I relative to the sun?
                        math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                        if (this->myUniverse->central_body_SPICE_ID == 10)
                        {
                            R_sc_Sun = this->output_state.getSubMatrix1D(0, 2);
                        }
                        else
                        {
                            //where is the central body relative to the sun?
                            doubleType central_body_state_and_derivatives[12];
                            this->myUniverse->locate_central_body(this->output_state(7),
                                central_body_state_and_derivatives,
                                *this->myOptions,
                                false);

                            math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                            for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                            {
                                R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                            }

                            R_sc_Sun = output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                        }

                        //Step 2.3.2: call the power model
                        doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                        this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                        //Step 2.3.3: call the thruster model
                        if (this->num_controls == 4)
                            this->mySpacecraft->computeElectricPropulsionPerformance(1.0, this->ControlVector[subStepIndex](3));
                        else
                            this->mySpacecraft->computeElectricPropulsionPerformance(1.0);

                        bool ifprint = true;

                        if (timeToPropagate == .1)
                            ManeuverThrottleLevel = this->mySpacecraft->getEPThrottleLevel();

                        if (this->mySpacecraft->getEPThrottleLevel() != ManeuverThrottleLevel)
                        {
                            if (divisor < 86400.0)
                            {
                                timeToPropagate -= this->EphemerisOutputResolution / divisor;
                                divisor *= 2.0;
                                ifprint = false;
                            }
                            else
                            {
                                ManeuverThrottleLevel = this->mySpacecraft->getEPThrottleLevel();
                                divisor = 1.0;
                            }
                        }

                        if (ifprint)
                        {
                            if (this->myOptions->generate_acceleration_model_instrumentation_file)
                            {
                                this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7), this->ControlVector[subStepIndex]);
                            }

                            //Step 2.3.4: print
                            this->write_ephemeris_line(outputfile,
                                output_state,
                                this->ControlVector[subStepIndex],
                                this->mySpacecraft->getEPthrust() * this->throttle[subStepIndex],
                                this->mySpacecraft->getEPMassFlowRate() * this->throttle[subStepIndex],
                                this->mySpacecraft->getEPIsp(),
                                this->mySpacecraft->getEPNumberOfActiveThrusters(),
                                this->mySpacecraft->getEPActivePower(),
                                this->mySpacecraft->getEPThrottleLevelString());
                        }

                    }

                    //Step 2.4: increment propagatedEpoch
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution / divisor
                        ? this->EphemerisOutputResolution / divisor
                        : (totalPropagationTime - timeToPropagate);
                }

                //Step 2.5: reset the propagator to its original state vector
                this->myPropagators[subStepIndex].setStateRight(subStepIndex == this->num_interior_control_points - 1 ? this->StateStepRightInertial : this->StateAfterSubStepInertial[subStepIndex]);
            }//end loop over interior control points

            //Step 3: propagate and print the coast arc
            {
                this->myCoastPropagator->setStateRight(output_state);
                double timeToPropagate = 0.0;
                double totalPropagationTime = this->StepCoastTime _GETVALUE;
				bool nav_coast_end_target = true;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 3.1: propagate
                    this->CoastdPropagatedStatedIndependentVariable.assign_zeros();
                    this->myCoastPropagator->setCurrentEpoch(this->StateAfterThrustInertial(7));
                    this->myCoastPropagator->setIndexOfEpochInStateVec(7);
                    this->myCoastPropagator->setCurrentIndependentVariable(this->StateAfterThrustInertial(7));
                    this->myCoastPropagator->propagate(timeToPropagate, false);
                    output_state(7) = this->StateAfterThrustInertial(7) + timeToPropagate;
                    this->temp_state = this->output_state;

                    //Step 3.2: convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //Step 3.3: print
                    if ((totalPropagationTime - timeToPropagate) > 0)
                    {
                        //we need an instantaneous power/propulsion state
                        //Step 3.2.3.1: where am I relative to the sun?
                        math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                        if (this->myUniverse->central_body_SPICE_ID == 10)
                        {
                            R_sc_Sun = output_state.getSubMatrix1D(0, 2);
                        }
                        else
                        {
                            //where is the central body relative to the sun?
                            doubleType central_body_state_and_derivatives[12];
                            this->myUniverse->locate_central_body(output_state(7),
                                central_body_state_and_derivatives,
                                *this->myOptions,
                                false);

                            math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                            for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                            {
                                R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                            }

                            R_sc_Sun = output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                        }

                        //Step 3.3.2: call the power model
                        doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                        this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                        if (this->myOptions->generate_acceleration_model_instrumentation_file)
                        {
                            this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
                        }

                        math::Matrix<doubleType> zero_vec(3, 1, 0.0);
                        //Step 3.3.3: print
                        this->write_ephemeris_line(outputfile,
                            output_state,
                            zero_vec,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            "none");
                    }

                    //Step 3.4: increment propagatedEpoch                  
					if ((totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution)
					{
						timeToPropagate += this->EphemerisOutputResolution;
					}
					else
					{
						if (nav_coast_end_target)
						{
							timeToPropagate += totalPropagationTime - timeToPropagate - 1.0;
							nav_coast_end_target = false;
						}
						else
						{
							timeToPropagate += totalPropagationTime - timeToPropagate;
						}
					}
                }

                //Step 3.5: reset the propagator to its original state vector
                this->myCoastPropagator->setStateRight(this->StateStepRightInertial);
            }//end output of coast arc
        }//end output_ephemeris()

        void PSFB_HifiDuty_step::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {

            //Step 1: create the target spec
            if (haveManeuverNeedTarget)
            {
                //Step 1.1: initialize target spec object
                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->StateStepLeftInertial(7),
                    this->StateStepLeftInertial);

                //Step 1.2: write target spec object
                myTargetSpecLine.write(target_spec_file);

                //Step 1.3: disable the target flag
                haveManeuverNeedTarget = false;
            }

            //Step 2: create the maneuver spec
            //Step 2.1: initialize maneuver spec object
            maneuver_spec_line myManeuverSpecLine(this->name);

            //Step 2.2: step through the PSFB_HifiDuty step at IntegrationStep intervals, every time throttle level changes save off a thrust step
            size_t ManeuverThrottleLevel;
            doubleType ManeuverStartEpoch;
            doubleType ManeuverStartMass;
            doubleType ManeuverThrustMagnitude;
            doubleType ManeuverMassFlowRate;
            //Step 2.2.1: propulsion characteristics on the left side of the step
            {
                //Step 2.2.1.1: locate the spacecraft relative to the sun
                for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
                    output_state(stateIndex) = this->StateStepLeftInertial(stateIndex);

                math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                if (this->myUniverse->central_body_SPICE_ID == 10)
                {
                    R_sc_Sun = output_state.getSubMatrix1D(0, 2);
                }
                else
                {
                    //where is the central body relative to the sun?
                    doubleType central_body_state_and_derivatives[12];
                    this->myUniverse->locate_central_body(output_state(7),
                        central_body_state_and_derivatives,
                        *this->myOptions,
                        false);

                    math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                    }

                    R_sc_Sun = output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                }

                //Step 2.2.1.2: call the power model
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                //Step 2.2.1.3: call the thruster model
                if (this->num_controls == 4)
                    this->mySpacecraft->computeElectricPropulsionPerformance(1.0, this->ControlVector[0](3));
                else
                    this->mySpacecraft->computeElectricPropulsionPerformance(1.0);

                //Step 2.2.1.3: populate fields
                ManeuverThrottleLevel = this->mySpacecraft->getEPThrottleLevel();
                ManeuverStartEpoch = output_state(7);
                ManeuverStartMass = output_state(6);
                ManeuverThrustMagnitude = this->mySpacecraft->getEPthrust(); //in N
                ManeuverMassFlowRate = this->mySpacecraft->getEPMassFlowRate();
            }

			double divisor = 1.0;
            //Step 2.2.2: integrate through the step and save off new thrust entries as needed
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                //Step 2.2.2.1: hijack the propagator
                this->myPropagators[subStepIndex].setStateRight(this->output_state);

                //Step 2.2.2.2: propagate the thrust arc, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->subStepThrustTime _GETVALUE;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 2.2.2.3: propagate
                    if (subStepIndex == 0)
                    {
                        this->dPropagatedStatedIndependentVariable[subStepIndex].assign_zeros();
                        this->myPropagators[subStepIndex].setCurrentEpoch(this->StateStepLeftInertial(7));
                        this->myPropagators[subStepIndex].setCurrentIndependentVariable(this->StateStepLeftInertial(7));
                    }
                    else
                    {
                        this->dPropagatedStatedIndependentVariable[subStepIndex] = this->dPropagatedStatedIndependentVariable[subStepIndex - 1];
                        this->myPropagators[subStepIndex].setCurrentEpoch(this->StateAfterSubStepInertial[subStepIndex - 1](7));
                        this->myPropagators[subStepIndex].setCurrentIndependentVariable(this->StateAfterSubStepInertial[subStepIndex - 1](7));
                    }
                    this->myPropagators[subStepIndex].setIndexOfEpochInStateVec(7);
                    this->myPropagators[subStepIndex].propagate(timeToPropagate, this->ControlVector[subStepIndex], false);

                    //Step 2.2.2.4: convert to Sun-centered if necessary
                    if (!(boost::to_lower_copy(this->myUniverse->central_body_name) == "sun"))
                    {
                        doubleType body_state[12];

                        this->myUniverse->locate_central_body(output_state(7),
                            body_state,
                            *this->myOptions,
                            false);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            output_state(stateIndex) += body_state[stateIndex];
                    }

                    //we need an instantaneous power/propulsion state
                    //Step 2.2.2.5: where am I relative to the sun?
                    math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                    if (this->myUniverse->central_body_SPICE_ID == 10)
                    {
                        R_sc_Sun = output_state.getSubMatrix1D(0, 2);
                    }
                    else
                    {
                        //where is the central body relative to the sun?
                        doubleType central_body_state_and_derivatives[12];
                        this->myUniverse->locate_central_body(output_state(7),
                            central_body_state_and_derivatives,
                            *this->myOptions,
                            false);

                        math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                        for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                        {
                            R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                        }

                        R_sc_Sun = output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                    }

                    //Step 2.2.2.6: call the power model
                    doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                    this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                    //Step 2.2.2.7: call the thruster model
                    if (this->num_controls == 4)
                        this->mySpacecraft->computeElectricPropulsionPerformance(1.0, this->ControlVector[subStepIndex](3));
                    else
                        this->mySpacecraft->computeElectricPropulsionPerformance(1.0);

                    //Step 2.2.2.8: did the thrust change? if so save off a thrust arc and reset

                    size_t CurrentThrottleLevel = this->mySpacecraft->getEPThrottleLevel();
                    doubleType CurrentThrustMagnitude = this->mySpacecraft->getEPthrust(); //in N
                    doubleType CurrentMassFlowRate = this->mySpacecraft->getEPMassFlowRate();

                    if (CurrentThrottleLevel != ManeuverThrottleLevel
                        || fabs(CurrentThrustMagnitude - ManeuverThrustMagnitude) > 1.0e-3) //if different throttle level or thrust different by more than one mN
                    {
						if (divisor < 86400.0)
						{
							timeToPropagate -= this->EphemerisOutputResolution / divisor;
							divisor *= 2.0;
						}
						else
						{
							divisor = 1.0;
	                        //save off a maneuver spec item
	                        myManeuverSpecLine.append_maneuver_spec_item("EME2000",
	                            ManeuverStartEpoch,
	                            this->ControlVector[subStepIndex],
	                            ManeuverStartMass,
	                            output_state(6),
	                            ManeuverThrustMagnitude,
	                            ManeuverMassFlowRate,
	                            output_state(7) - ManeuverStartEpoch,
	                            1.0);

	                        //reset for next maneuver item
	                        ManeuverThrottleLevel = CurrentThrottleLevel;
	                        ManeuverStartEpoch = output_state(7);
	                        ManeuverStartMass = output_state(6);
	                        ManeuverThrustMagnitude = CurrentThrustMagnitude;
	                        ManeuverMassFlowRate = CurrentMassFlowRate;
						}
                    }
                    	
					//Step 2.2.2.9: increment propagatedEpoch - we deliberately do NOT take the last partial step
                    timeToPropagate += this->EphemerisOutputResolution / divisor; 
                }//end propagation over substeps

                 //Step 2.2.3: reset the propagator
                this->myPropagators[subStepIndex].setStateRight(subStepIndex == this->num_interior_control_points - 1 ? this->StateStepRightInertial : this->StateAfterSubStepInertial[subStepIndex]);
            }

            //Step 2.3: save off the last thrust step
            myManeuverSpecLine.append_maneuver_spec_item("EME2000",
                ManeuverStartEpoch,
                this->ControlVector.back(),
                ManeuverStartMass,
                this->StateAfterThrustInertial(6),
                ManeuverThrustMagnitude,
                ManeuverMassFlowRate,
                StateAfterThrustInertial(7) - ManeuverStartEpoch,
                1.0);

            //Step 2.4: write the maneuver spec line
            myManeuverSpecLine.write(maneuver_spec_file);

            //Step 2.5: set the target flag
            haveManeuverNeedTarget = true;

        }//end output_maneuver_and_target_spec()

        void PSFB_HifiDuty_step::output_STM(size_t& STMindex)
        {
            this->CoastAugmentedSTM.print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_STM_" + std::to_string(STMindex) + "_coast_forward.stm");
            ++STMindex;
        }//end output_STM
    }//close namespace Phases
}//close namespace EMTG