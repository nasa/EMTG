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

//parallel shooting bounded impulse (PSBI) step for EMTGv9
//Jacob Englander 2/8/2019

#include "doubleType.h"

#include "PSBIstep.h"

#include "PropagatorFactory.h"

#include "maneuver_spec_line.h"
#include "target_spec_line.h"

namespace EMTG
{
    namespace Phases
    {
        PSBIstep::PSBIstep() {};

        PSBIstep::PSBIstep(const std::string& name,
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

        void PSBIstep::initialize(const std::string& name,
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

            //state holders
            this->state_before_maneuver.resize(10, 1, 0.0);
            this->state_after_maneuver.resize(10, 1, 0.0);

            //maneuver holder
            this->DeltaV.resize(3, 1, 0.0);

            //create maneuver object
            this->MTM.resize(11, 11, math::MatrixType::identity);
            this->augmentedMTM.resize(14, 14, math::MatrixType::identity);

            //configure propagator
            this->configure_propagator();


            //acceleration model
                this->mySpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(this->myOptions,
                this->myJourneyOptions,
                this->myUniverse,
                this->Xdescriptions,
                this->mySpacecraft,                                
                14,// STM size
                false);//we do NOT need the central body because Kepler propagation takes care of it
            this->mySpacecraftAccelerationModel->setDutyCycle(this->StepDutyCycle);//we're not actually using the acceleration model for thrust right now

            this->myBoundedImpulseManeuver = ForwardBoundedImpulseManeuver(this->journeyIndex,
                this->phaseIndex,
                this->stageIndex,
                this->myUniverse,
                this->mySpacecraft,
                this->myOptions);

            //epoch derivative
            this->dImpulseEpoch_dPropagationVariable = (this->stepIndex + 1 - 0.5) / this->myPhase->get_num_steps();

            this->myBoundedImpulseManeuver.set_ThrustStepLength(this->StepFlightTime);
            this->myBoundedImpulseManeuver.set_dThrustStepLength_dPropagationVariable(this->dStepTime_dPhaseFlightTime);
            this->myBoundedImpulseManeuver.set_dImpulseEpoch_dPropagationVariable_pointer(this->dImpulseEpoch_dPropagationVariable);
            this->myBoundedImpulseManeuver.set_LaunchDate(this->LaunchDate);
            this->myBoundedImpulseManeuver.set_max_thrust(this->max_thrust);
            this->myBoundedImpulseManeuver.set_max_mass_flow_rate(this->max_mass_flow_rate);
            this->myBoundedImpulseManeuver.set_Isp(this->Isp);
            this->myBoundedImpulseManeuver.set_power(this->power);
            this->myBoundedImpulseManeuver.set_active_power(this->active_power);
            this->myBoundedImpulseManeuver.set_number_of_active_engines(this->number_of_active_engines);
            this->myBoundedImpulseManeuver.set_ThrottleLevel(this->ThrottleLevel);
            this->myBoundedImpulseManeuver.set_ThrottleLevelString(this->ThrottleLevelString);
            this->myBoundedImpulseManeuver.set_DeltaV(this->DeltaV);
            this->myBoundedImpulseManeuver.set_spacecraft_state_minus(this->state_before_maneuver);
            this->myBoundedImpulseManeuver.set_spacecraft_state_plus(this->state_after_maneuver);
            this->myBoundedImpulseManeuver.set_MTM(this->MTM);
            this->myBoundedImpulseManeuver.set_DutyCycle(this->StepDutyCycle);
            this->myBoundedImpulseManeuver.set_Control(this->ControlVector[0]);
            this->myBoundedImpulseManeuver.set_Throttle(this->throttle[0]);
            this->myBoundedImpulseManeuver.set_AccelerationModel(this->mySpacecraftAccelerationModel);

            //output stuff
            this->output_state.resize(10 + 13 * 13, 1, 0.0);
        }//end initialize

        PSBIstep::~PSBIstep()
        {
            delete this->mySpacecraftAccelerationModel;
        }//end destructor


        //configure propagator
        void PSBIstep::configure_propagator()
        {
            //we have only one interior control point by definition 
            this->num_interior_control_points = 1;
            this->STM.resize(this->num_interior_control_points, math::Matrix<double>(14, 14, math::MatrixType::identity));
            this->AugmentedSTM.resize(this->num_interior_control_points, math::Matrix<double>(14, 14, math::MatrixType::identity));
            this->CumulativeAugmentedSTM.resize(this->num_interior_control_points, math::Matrix<double>(14, 14, math::MatrixType::identity));
            this->dPropagatedStatedIndependentVariable.resize(this->num_interior_control_points, math::Matrix<double>(10, 2, 0.0));

            //PSBI-specific stuff
            this->STM1.resize(6, 6, math::MatrixType::identity);
            this->STM2.resize(6, 6, math::MatrixType::identity);
            this->augmentedSTM1.resize(14, 14, math::MatrixType::identity);
            this->augmentedSTM2.resize(14, 14, math::MatrixType::identity);
            this->dStatedIndependentVariable1.resize(6, 1, 0.0);
            this->dStatedIndependentVariable2.resize(6, 1, 0.0);
            this->dHalfStepTime_dPhaseFlightTime = this->dStepTime_dPhaseFlightTime / 2;

            //propagators
            //first half of the phase
            this->myPropagators.push_back(CreatePropagator(this->myOptions,
                this->myUniverse,
                6,
                this->StateStepLeftInertial,
                this->state_before_maneuver,
                this->STM1,
                this->dStatedIndependentVariable1,
                &this->dHalfStepTime_dPhaseFlightTime));

            //second half of the phase
            this->myPropagators.push_back(CreatePropagator(this->myOptions,
                this->myUniverse,
                6,
                this->state_after_maneuver,
                this->StateStepRightInertial,
                this->STM2,
                this->dStatedIndependentVariable2,
                &this->dHalfStepTime_dPhaseFlightTime));
        }//end configure propagator

        //master calcbounds
        void PSBIstep::calcbounds_step()
        {
            this->calcbounds_step_left_state();

            this->calcbounds_step_control();

            this->calcbounds_step_left_match_point_constraints();

            this->calcbounds_step_main();

            this->calcbounds_distance_constraints();

            this->calcbounds_maneuver_constraints();
        }//end calcbounds_step

        //master process
        void PSBIstep::process_step(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
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
        void PSBIstep::process_step_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: propagate to maneuver point

            //Step 2.2: propagate and update epoch and mass
            this->myPropagators[0].propagate(this->StepFlightTime / 2.0, needG);

            this->state_before_maneuver(6) = this->StateStepLeftInertial(6);
            this->state_before_maneuver(7) = this->StateStepLeftInertial(7) + this->StepFlightTime / 2.0;
            this->state_before_maneuver(8) = this->StateStepLeftInertial(8);
            this->state_before_maneuver(9) = this->StateStepLeftInertial(9);

            //Step 2: perform maneuver
            this->myBoundedImpulseManeuver.process_maneuver(needG);

            //Step 3: propagate to end of step
            this->myPropagators[1].propagate(this->StepFlightTime / 2.0, needG);

            this->StateStepRightInertial(6) = this->state_after_maneuver(6);
            this->StateStepRightInertial(7) = this->state_after_maneuver(7) + this->StepFlightTime / 2.0;
            this->StateStepRightInertial(8) = this->state_after_maneuver(8);
            this->StateStepRightInertial(9) = this->state_after_maneuver(9);

            //Step 4: refactor the STMs and MTM into the format that ParallelShootingStep expects
            if (needG)
            {
                //Step 4.1: form the augmented STMs
                //upper right 6x6 is the original STM
                for (size_t i = 0; i < 6; ++i)
                {
                    for (size_t j = 0; j < 6; ++j)
                    {
                        this->augmentedSTM1(i, j) = this->STM1(i, j);
                        this->augmentedSTM2(i, j) = this->STM2(i, j);
                    }

                    //then we need to turn the upper right 6 x 1 into the Phi_t terms developed by Lantoine
                    this->augmentedSTM1(i, 13) = this->dStatedIndependentVariable1(i);
                    this->augmentedSTM2(i, 13) = this->dStatedIndependentVariable2(i);
                }
                this->augmentedSTM1(7, 13) = this->dHalfStepTime_dPhaseFlightTime;
                this->augmentedSTM2(7, 13) = this->dHalfStepTime_dPhaseFlightTime;

                //Step 4.2: form the augmented MTM
                //this is like the original MTM, except that:
                //1. the propagation time entry index is now 13 instead of 10
                //2. The control entries are now derivative of velocity w.r.t. control instead of velocity w.r.t. velocity
                
                //Step 4.2.1: derivatives of velocity
                //leave derivatives of position w.r.t. position and velocity w.r.t. velocity as identity
                //other derivatives of velocity
                double dvmax = this->myBoundedImpulseManeuver.get_dvmax() _GETVALUE;
                for (size_t velocityIndex : {3, 4, 5})
                {
                    //derivative of velocity w.r.t. control is dvmax
                    this->augmentedMTM(velocityIndex, 10 + velocityIndex - 3) = dvmax;

                    //derivative of velocity w.r.t. position is drawn from the original MTM
                    for (size_t positionIndex : {0, 1, 2})
                        this->augmentedMTM(velocityIndex, positionIndex) = this->MTM(velocityIndex, positionIndex);

                    //derivative of velocity w.r.t. mass is drawn from the original MTM
                    this->augmentedMTM(velocityIndex, 6) = this->MTM(velocityIndex, 6);

                    //derivative of velocity w.r.t. previous times is drawn from the original MTM
                    this->augmentedMTM(velocityIndex, 7) = this->MTM(velocityIndex, 7);

                    //derivative of velocity w.r.t. phase time of flight is drawn from the original MTM - index change from 10 to 13
                    this->augmentedMTM(velocityIndex, 13) = this->MTM(velocityIndex, 10);
                }

                //Step 4.2.2: derivatives of mass
                //derivative of mass w.r.t. control
                math::Matrix<double> dMassAfterManeuver_dThrottleComponents = this->myBoundedImpulseManeuver.get_dMassAfterManeuver_dThrottleComponents();
                for (size_t controlIndex : {0, 1, 2})
                    this->augmentedMTM(6, 10 + controlIndex) = dMassAfterManeuver_dThrottleComponents(controlIndex);

                //derivative of mass w.r.t. position is drawn from the original MTM
                for (size_t positionIndex : {0, 1, 2})
                    this->augmentedMTM(6, positionIndex) = this->MTM(6, positionIndex);

                //derivative of mass w.r.t. previous times is drawn from the original MTM
                this->augmentedMTM(6, 7) = this->MTM(6, 7);

                //derivative of mass w.r.t. phase time of flight is drawn from the original MTM - index change from 10 to 13
                this->augmentedMTM(6, 13) = this->MTM(6, 10);

                //Step 4.2.3: derivatives of chemical fuel tank
                //derivative of chemical fuel tank w.r.t. control
                math::Matrix<double> dChemicalFuel_dThrottleComponents = this->myBoundedImpulseManeuver.get_dChemicalFuel_dThrottleComponents();
                for (size_t controlIndex : {0, 1, 2})
                    this->augmentedMTM(8, 10 + controlIndex) = dChemicalFuel_dThrottleComponents(controlIndex);
                
                //derivative of chemical fuel tank w.r.t. phase time of flight is drawn from the original MTM - index change from 10 to 13
                this->augmentedMTM(8, 13) = this->MTM(8, 10);

                //Step 4.2.4: derivatives of electric tank
                //derivative of electric tank w.r.t. control
                math::Matrix<double> dElectricPropellant_dThrottleComponents = this->myBoundedImpulseManeuver.get_dElectricPropellant_dThrottleComponents();
                for (size_t controlIndex : {0, 1, 2})
                    this->augmentedMTM(9, 10 + controlIndex) = dElectricPropellant_dThrottleComponents(controlIndex);

                //derivative of electric tank w.r.t. position is drawn from the original MTM
                for (size_t positionIndex : {0, 1, 2})
                    this->augmentedMTM(9, positionIndex) = this->MTM(9, positionIndex);

                //derivative of electric tank w.r.t. previous times is drawn from the original MTM
                this->augmentedMTM(9, 7) = this->MTM(9, 7);

                //derivative of electric tank w.r.t. phase time of flight is drawn from the original MTM - index change from 10 to 13
                this->augmentedMTM(9, 13) = this->MTM(9, 10);
            }

            //Step 5: construct the cumulative step transition matrix
            this->CumulativeAugmentedSTM[0] = this->augmentedSTM2 * this->augmentedMTM * this->augmentedSTM1;
        }//end process_step_main

        void PSBIstep::output(std::ofstream& outputfile,
            size_t& eventcount)
        {
            //Step 2.1: output the current step
            std::string event_type;
            if (this->max_thrust > 1.0e-7 && this->throttle[0] > 5.0e-4)
                event_type = "PSBIthrust";
            else
                event_type = "coast";

            math::Matrix<doubleType> ThrustVector = this->ControlVector[0] * this->max_thrust * 1000.0 / this->StepDutyCycle;

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
                this->state_before_maneuver,//state
                this->DeltaV,//dV
                ThrustVector,//ThrustVector
                this->DeltaV.norm(),//dVmag
                this->max_thrust * 1000.0 / this->StepDutyCycle,//Thrust
                this->Isp,//Isp
                this->power,//AvailPower
                this->max_mass_flow_rate / this->StepDutyCycle,//mdot
                this->number_of_active_engines,//number_of_active_engines
                this->active_power,//active power
                this->ThrottleLevelString);//throttle level)
        }//end output

        void PSBIstep::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
        {
            //Step 1: propagate before the maneuver
            //Step 1.1: temporarily assign the propagator to the output state
            this->myPropagators[0].setStateRight(this->output_state);

            //Step 1.2: propagate and print the thrust arc, skipping the first entry
            double timeToPropagate = this->EphemerisOutputResolution;
            double totalPropagationTime = this->StepFlightTime _GETVALUE / 2.0;
            while (timeToPropagate < totalPropagationTime)
            {
                //Step 1.2.1: propagate
                this->dStatedIndependentVariable1.assign_zeros();
                this->myPropagators[0].setCurrentEpoch(this->StateStepLeftInertial(7));
                this->myPropagators[0].setCurrentIndependentVariable(this->StateStepLeftInertial(7));
                this->myPropagators[0].setIndexOfEpochInStateVec(7);
                this->myPropagators[0].propagate(timeToPropagate, this->ControlVector[0], false);

                //Step 1.2.2: convert to Sun-centered if necessary
                if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                {
                    double LT_dump;
                    double bodyStateDouble[6];
                    spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                    for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                        output_state(stateIndex) += bodyStateDouble[stateIndex];
                }

                //Step 1.2.3: print
                if ((totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution)
                {
                    //we need an instantaneous power/propulsion state
                    //Step 1.2.3.1: where am I relative to the sun?
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

                    //Step 1.2.3.2: call the power model
                    doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                    this->mySpacecraft->computePowerState(r_sc_sun_AU, output_state(7));

                    //Step 1.2.3.3: call the thruster model
                    if (this->num_controls == 4)
                        this->mySpacecraft->computeElectricPropulsionPerformance(this->StepDutyCycle, this->ControlVector[0](3));
                    else
                        this->mySpacecraft->computeElectricPropulsionPerformance(this->StepDutyCycle);

                    //Step 1.2.3.4: print
                    this->write_ephemeris_line(outputfile,
                        output_state,
                        this->ControlVector[0],
                        this->mySpacecraft->getEPthrust() * 1.0e-3 * this->throttle[0],
                        this->mySpacecraft->getEPMassFlowRate() * this->throttle[0],
                        this->mySpacecraft->getEPIsp(),
                        this->mySpacecraft->getEPNumberOfActiveThrusters(),
                        this->mySpacecraft->getEPActivePower(),
                        this->mySpacecraft->getEPThrottleLevelString());
                }

                //Step 1.2.4: increment propagatedEpoch
                timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution
                    ? this->EphemerisOutputResolution
                    : (totalPropagationTime - timeToPropagate);
            }

            //Step 1.3: reset the propagator to its original state vectors
            this->myPropagators[0].setStateRight(this->state_before_maneuver);

            //Step 2: propagate before the maneuver
            //Step 2.1: temporarily assign the propagator to the output state
            this->myPropagators[1].setStateRight(this->output_state);

            //Step 2.2: propagate and print the thrust arc, skipping the first entry
            timeToPropagate = this->EphemerisOutputResolution;
            totalPropagationTime = this->StepFlightTime _GETVALUE / 2.0;
            while (timeToPropagate < totalPropagationTime)
            {
                //Step 2.2.1: propagate
                this->dStatedIndependentVariable2.assign_zeros();
                this->myPropagators[1].setCurrentEpoch(this->StateStepLeftInertial(7));
                this->myPropagators[1].setCurrentIndependentVariable(this->StateStepLeftInertial(7));
                this->myPropagators[1].setIndexOfEpochInStateVec(7);
                this->myPropagators[1].propagate(timeToPropagate, this->ControlVector[0], false);

                //Step 2.2.2: convert to Sun-centered if necessary
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

                //Step 2.2.3: print
                if ((totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution)
                {
                    //we need an instantaneous power/propulsion state
                    //Step 2.2.3.1: where am I relative to the sun?
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

                    //Step 2.2.3.2: call the power model
                    doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                    this->mySpacecraft->computePowerState(r_sc_sun_AU, output_state(7));

                    //Step 2.2.3.3: call the thruster model
                    if (this->num_controls == 4)
                        this->mySpacecraft->computeElectricPropulsionPerformance(this->StepDutyCycle, this->ControlVector[0](3));
                    else
                        this->mySpacecraft->computeElectricPropulsionPerformance(this->StepDutyCycle);

                    //Step 2.2.3.4: print
                    this->write_ephemeris_line(outputfile,
                        output_state,
                        this->ControlVector[0],
                        this->mySpacecraft->getEPthrust() * 1.0e-3 * this->throttle[0],
                        this->mySpacecraft->getEPMassFlowRate() * this->throttle[0],
                        this->mySpacecraft->getEPIsp(),
                        this->mySpacecraft->getEPNumberOfActiveThrusters(),
                        this->mySpacecraft->getEPActivePower(),
                        this->mySpacecraft->getEPThrottleLevelString());
                }

                //Step 2.2.4: increment propagatedEpoch
                timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution
                    ? this->EphemerisOutputResolution
                    : (totalPropagationTime - timeToPropagate);
            }

            //Step 2.3: reset the propagator to its original state vectors
            this->myPropagators[1].setStateRight(this->StateStepRightInertial);

        }//end output_ephemeris()

        void PSBIstep::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {

            //Step 1: create the target spec if a target is needed
            if (haveManeuverNeedTarget)
            {
                //Step 1.1: initialize target spec object
                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->state_before_maneuver(7),
                    this->state_after_maneuver);

                //Step 1.2: write target spec object
                myTargetSpecLine.write(target_spec_file);

                //Step 1.3: disable the target flag
                haveManeuverNeedTarget = false;
            }

            //Step 2: create the maneuver spec
            //Step 2.1: initialize maneuver spec object
            maneuver_spec_line myManeuverSpecLine(this->name);

            //Step 2.2: only one line, and it uses pre-stored data
            //save off a maneuver spec item
            myManeuverSpecLine.append_maneuver_spec_item("EME2000",
                this->state_before_maneuver(7),
                this->ControlVector[0],
                this->state_before_maneuver(6),
                this->state_after_maneuver(6),
                this->max_thrust,
                this->max_mass_flow_rate,
                this->StepFlightTime * this->throttle[0],
                this->StepDutyCycle);

            //Step 2.3: set the target flag
            haveManeuverNeedTarget = true;

        }//end output_maneuver_and_target_spec()

    }//close namespace Phases
}//close namespace EMTG