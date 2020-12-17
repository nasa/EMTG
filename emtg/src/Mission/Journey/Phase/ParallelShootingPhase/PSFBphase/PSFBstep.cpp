// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2020 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Other Rights Reserved.

// Licensed under the NASA Open Source License (the "License"); 
// You may not use this file except in compliance with the License. 
// You may obtain a copy of the License at:
// https://opensource.org/licenses/NASA-1.3
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
// express or implied.   See the License for the specific language
// governing permissions and limitations under the License.

//parallel shooting finite burn (PSFB) step for EMTGv9
//Jacob Englander 2-23-2018

#include "doubleType.h"

#include "PSFBstep.h"

#include "PropagatorFactory.h"
#include "IntegrationSchemeFactory.h"

#include "maneuver_spec_line.h"
#include "target_spec_line.h"

namespace EMTG
{
    namespace Phases
    {
        PSFBstep::PSFBstep() {};

        PSFBstep::PSFBstep(const std::string& name,
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

        void PSFBstep::initialize(const std::string& name,
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
            this->ParallelShootingStep::initialize(name,
                journeyIndex,
                phaseIndex,
                stepIndex,
                stageIndex,
                (ParallelShootingStep*) previousStep, //typecast to pointer to base, safe because only base methods will be called
                (ParallelShootingPhase*) myPhase, //typecast to pointer to base, safe because only base methods will be called
                myUniverse,
                mySpacecraft,
                myOptions);

            //configure propagator
            this->configure_propagator();

            //output stuff
            this->output_state.resize(10 + 13 * 13, 1, 0.0);
        }//end initialize

        PSFBstep::~PSFBstep()
        {
            //propagator
            //we no longer need to delete the propagator because boost::ptr_vector has its own garbage cleanup

            //acceleration model
            delete this->mySpacecraftAccelerationModel;

            //integration scheme
            delete this->myIntegrationScheme;
        }//end destructor


        //configure propagator
        void PSFBstep::configure_propagator()
        {
            this->STM.resize(this->num_interior_control_points, math::Matrix<double>(14, 14, math::identity));
            this->AugmentedSTM.resize(this->num_interior_control_points, math::Matrix<double>(14, 14, math::identity));
            this->CumulativeAugmentedSTM.resize(this->num_interior_control_points, math::Matrix<double>(14, 14, math::identity));
            this->dPropagatedStatedIndependentVariable.resize(this->num_interior_control_points, math::Matrix<double>(10, 2, 0.0));

            //acceleration model object
            this->mySpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(this->myOptions,
                this->myJourneyOptions,
                this->myUniverse,
                this->Xdescriptions,
                this->mySpacecraft,
                14); // STM size
            this->mySpacecraftAccelerationModel->setDutyCycle(this->StepDutyCycle);

            //EOM
            this->myEOM.setSpacecraftAccelerationModel(this->mySpacecraftAccelerationModel);

            //integration scheme
            this->myIntegrationScheme = CreateIntegrationScheme(&this->myEOM, 10, 14);

            //propagator
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                this->myPropagators.push_back(CreatePropagator(this->myOptions,
                    this->myUniverse,
                    10,
                    14,
                    subStepIndex == 0 ? this->StateStepLeftInertial : this->StateAfterSubStepInertial[subStepIndex - 1],
                    subStepIndex == this->num_interior_control_points - 1 ? this->StateStepRightInertial : this->StateAfterSubStepInertial[subStepIndex],
                    this->STM[subStepIndex],
                    this->dPropagatedStatedIndependentVariable[subStepIndex],
                    (Integration::Integrand*) &this->myEOM,
                    this->myIntegrationScheme,
                    &this->dStepTime_dPhaseFlightTime,
                    this->myJourneyOptions->override_integration_step_size
                    ? this->myJourneyOptions->integration_step_size
                    : this->myOptions->integration_time_step_size));
            }
        }//end configure propagator

        //master calcbounds
        void PSFBstep::calcbounds_step()
        {
            this->calcbounds_step_left_state();

            this->calcbounds_step_control();

            this->calcbounds_step_left_match_point_constraints();

            this->calcbounds_step_main();

            this->calcbounds_distance_constraints();

            this->calcbounds_maneuver_constraints();
        }//end calcbounds_step

        //master process
        void PSFBstep::process_step(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->total_number_of_states_to_integrate = needG
                ? 10 + 14 * 14
                : 10;

            this->process_step_left_state(X, Xindex, F, Findex, G, needG);

            this->process_step_control(X, Xindex, F, Findex, G, needG);

            this->process_step_left_match_point_constraints(X, Xindex, F, Findex, G, needG);

            this->process_step_main(X, Xindex, F, Findex, G, needG);

            this->process_distance_constraints(X, Xindex, F, Findex, G, needG);

            this->process_maneuver_constraints(X, Xindex, F, Findex, G, needG);

            if (needG)
                this->process_derivative_tuples(X, Xindex, F, Findex, G, needG);
        }//end process_step

        //process step main - aka propagate
        void PSFBstep::process_step_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                //Step 1: propagate
                if (subStepIndex == 0)
                {
                    this->myPropagators[subStepIndex].setCurrentEpoch(this->StateStepLeftInertial(7));
                    this->myPropagators[subStepIndex].setCurrentIndependentVariable(this->StateStepLeftInertial(7));
                }
                else
                {
                    this->myPropagators[subStepIndex].setCurrentEpoch(this->StateAfterSubStepInertial[subStepIndex - 1](7));
                    this->myPropagators[subStepIndex].setCurrentIndependentVariable(this->StateAfterSubStepInertial[subStepIndex - 1](7));
                }
                this->myPropagators[subStepIndex].setIndexOfEpochInStateVec(7);
                this->myPropagators[subStepIndex].propagate(this->StepFlightTime / this->num_interior_control_points, this->ControlVector[subStepIndex], needG);

                //Step 2: update the epoch
                if (subStepIndex == 0)
                    this->StateAfterSubStepInertial[subStepIndex](7) = this->StateStepLeftInertial(7) + this->StepFlightTime / this->num_interior_control_points;
                else if (subStepIndex == this->num_interior_control_points - 1)
                    this->StateStepRightInertial(7) = this->StateAfterSubStepInertial[subStepIndex - 1](7) + this->StepFlightTime / this->num_interior_control_points;
                else
                    this->StateAfterSubStepInertial[subStepIndex](7) = this->StateAfterSubStepInertial[subStepIndex - 1](7) + this->StepFlightTime / this->num_interior_control_points;

                //Step 3: derivatives
                if (needG)
                {
                    //Step 3.1: form the augmented STM
                    //upper right 7x7 is the original STM
                    for (size_t i = 0; i < 7; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                            this->AugmentedSTM[subStepIndex](i, j) = this->STM[subStepIndex](i, j);

                        //Phi_t terms
                        this->AugmentedSTM[subStepIndex](i, 7) = this->STM[subStepIndex](i, 7);
                        this->AugmentedSTM[subStepIndex](i, 13) = this->STM[subStepIndex](i, 13);

                        //control terms
                        this->AugmentedSTM[subStepIndex](i, 10) = this->STM[subStepIndex](i, 10);
                        this->AugmentedSTM[subStepIndex](i, 11) = this->STM[subStepIndex](i, 11);
                        this->AugmentedSTM[subStepIndex](i, 12) = this->STM[subStepIndex](i, 12);
                    }

                    //tanks
                    for (size_t i = 8; i < 10; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                            this->AugmentedSTM[subStepIndex](i, j) = this->STM[subStepIndex](i, j);

                        //Phi_t terms
                        this->AugmentedSTM[subStepIndex](i, 7) = this->STM[subStepIndex](i, 7);
                        this->AugmentedSTM[subStepIndex](i, 13) = this->STM[subStepIndex](i, 13);

                        //control terms
                        this->AugmentedSTM[subStepIndex](i, 10) = this->STM[subStepIndex](i, 10);
                        this->AugmentedSTM[subStepIndex](i, 11) = this->STM[subStepIndex](i, 11);
                        this->AugmentedSTM[subStepIndex](i, 12) = this->STM[subStepIndex](i, 12);
                    }


                    //epoch time with respect to propagation time
                    this->AugmentedSTM[subStepIndex](7, 13) = this->dStepTime_dPhaseFlightTime;
                }//end STM creation
            }//end loop over substeps


            //Step 4: construct cumulative STMs
            this->CumulativeAugmentedSTM.back() = this->AugmentedSTM.back();
            for (int subStepIndex = this->num_interior_control_points - 2; subStepIndex >= 0; --subStepIndex)
            {
                //create stripped version of next step's cumulative STM to remove that substep's control influence
                math::Matrix<double> StrippedSTM = this->CumulativeAugmentedSTM[subStepIndex + 1];

                for (size_t i = 0; i < 10; ++i)
                {
                    for (size_t j = 10; j < 10 + this->num_controls; ++j)
                        StrippedSTM(i, j) = 0.0;
                }

                this->CumulativeAugmentedSTM[subStepIndex] = StrippedSTM * this->AugmentedSTM[subStepIndex];
            }
        }//end process_step_main

        void PSFBstep::output(std::ofstream& outputfile,
            size_t& eventcount)
        {
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                //Step 1: temporarily redirect the propagator to the output state
                this->myPropagators[subStepIndex].setStateRight(this->output_state);

                //Step 2: propagate
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
                this->myPropagators[subStepIndex].propagate(this->StepFlightTime / this->num_interior_control_points / 2.0, this->ControlVector[subStepIndex], false);

                //Step 4: redirect the propagator back to where it is supposed to go
                this->myPropagators[subStepIndex].setStateRight(subStepIndex == this->num_interior_control_points - 1 ? this->StateStepRightInertial : this->StateAfterSubStepInertial[subStepIndex]);

                //Step 5: figure out spacecrafty things

                //Step 5.1: where am I relative to the sun?
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

                //Step 5.2: call the power model
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, output_state(7));

                //Step 5.3: call the thruster model
                if (this->num_controls == 4)
                    this->mySpacecraft->computeElectricPropulsionPerformance(this->StepDutyCycle, this->ControlVector[subStepIndex](3));
                else
                    this->mySpacecraft->computeElectricPropulsionPerformance(this->StepDutyCycle);

                //Step 5.4: store the thruster model outputs
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
                math::Matrix<doubleType> deltaV = ThrustVector * this->StepFlightTime / this->output_state(6);

                this->writey_thing::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    event_type,//event_type
                    "deep-space",//event_location
                    this->StepFlightTime / 86400.0,// timestep_size,
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
        }//end output

        void PSFBstep::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
        {
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                //Step 1: temporarily assign the propagator to the output state
                this->myPropagators[subStepIndex].setStateRight(this->output_state);

                //Step 2: propagate and print the thrust arc, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->StepFlightTime _GETVALUE / this->num_interior_control_points;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 2.1: propagate
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
                    if ((totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution)
                    {
                        //we need an instantaneous power/propulsion state
                        //Step 2.3.1: where am I relative to the sun?
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

                        //Step 2.3.2: call the power model
                        doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                        this->mySpacecraft->computePowerState(r_sc_sun_AU, output_state(7));

                        //Step 2.3.3: call the thruster model
                        if (this->num_controls == 4)
                            this->mySpacecraft->computeElectricPropulsionPerformance(this->StepDutyCycle, this->ControlVector[subStepIndex](3));
                        else
                            this->mySpacecraft->computeElectricPropulsionPerformance(this->StepDutyCycle);

                        //Step 2.3.4: print
                        if (this->myOptions->generate_acceleration_model_instrumentation_file)
                        {
                            this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7), this->ControlVector[subStepIndex]);
                        }

                        this->write_ephemeris_line(outputfile,
                            output_state,
                            this->ControlVector[subStepIndex],
                            this->mySpacecraft->getEPthrust() * 1.0e-3 * this->throttle[subStepIndex],
                            this->mySpacecraft->getEPMassFlowRate() * this->throttle[subStepIndex],
                            this->mySpacecraft->getEPIsp(),
                            this->mySpacecraft->getEPNumberOfActiveThrusters(),
                            this->mySpacecraft->getEPActivePower(),
                            this->mySpacecraft->getEPThrottleLevelString());
                    }

                    //Step 2.4: increment propagatedEpoch
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution
                        ? this->EphemerisOutputResolution
                        : (totalPropagationTime - timeToPropagate);
                }

                //Step 3: reset the propagator to its original state vectors
                this->myPropagators[subStepIndex].setStateRight(subStepIndex == this->num_interior_control_points - 1 ? this->StateStepRightInertial : this->StateAfterSubStepInertial[subStepIndex]);
            }//end loop over substeps
        }//end output_ephemeris()

        void PSFBstep::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {

            //Step 1: create the target spec if a target is needed
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

            //Step 2.2: step through the PSFB step at IntegrationStep intervals, every time throttle level changes save off a thrust step
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
                this->mySpacecraft->computePowerState(r_sc_sun_AU, output_state(7));

                //Step 2.2.1.3: call the thruster model using 100% duty cycle because we want the raw performance information
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

            //Step 2.2.2: integrate through the step and save off new thrust entries as needed
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                //Step 2.2.2.1: hijack the propagator
                this->myPropagators[subStepIndex].setStateRight(this->output_state);

                //Step 2.2.2.2: propagate the thrust arc, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->StepFlightTime _GETVALUE / this->num_interior_control_points;
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
                    this->mySpacecraft->computePowerState(r_sc_sun_AU, output_state(7));

                    //Step 2.2.2.7: call the thruster model using 100% duty cycle because we want the raw performance information
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
                        //save off a maneuver spec item
                        myManeuverSpecLine.append_maneuver_spec_item("EME2000",
                            ManeuverStartEpoch,
                            this->ControlVector[subStepIndex],
                            ManeuverStartMass,
                            output_state(6),
                            ManeuverThrustMagnitude,
                            ManeuverMassFlowRate,
                            output_state(7) - ManeuverStartEpoch,
                            this->StepDutyCycle);

                        //reset for next maneuver item
                        ManeuverThrottleLevel = CurrentThrottleLevel;
                        ManeuverStartEpoch = output_state(7);
                        ManeuverStartMass = output_state(6);
                        ManeuverThrustMagnitude = CurrentThrustMagnitude;
                        ManeuverMassFlowRate = CurrentMassFlowRate;
                    }


                    //Step 2.2.2.9: increment propagatedEpoch - we deliberately do NOT take the last partial step
                    timeToPropagate += this->EphemerisOutputResolution;
                }//end propagation over substeps

                //Step 2.2.3: reset the propagator
                this->myPropagators[subStepIndex].setStateRight(subStepIndex == this->num_interior_control_points - 1 ? this->StateStepRightInertial : this->StateAfterSubStepInertial[subStepIndex]);
            }//end loop over substeps

            //Step 2.3: save off the last thrust step

            myManeuverSpecLine.append_maneuver_spec_item("EME2000",
                ManeuverStartEpoch,
                this->ControlVector.back(),
                ManeuverStartMass,
                this->StateStepRightInertial(6),
                ManeuverThrustMagnitude,
                ManeuverMassFlowRate,
                StateStepRightInertial(7) - ManeuverStartEpoch,
                this->StepDutyCycle);

            //Step 2.4: write the maneuver spec line
            myManeuverSpecLine.write(maneuver_spec_file);

            //Step 2.5: set the target flag
            haveManeuverNeedTarget = true;

        }//end output_maneuver_and_target_spec()

    }//close namespace Phases
}//close namespace EMTG