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

//EMTGv9 MGALTphase
//Jacob Englander 6-24-2017

#include "MGALTphase.h"

namespace EMTG
{
    namespace Phases
    {
        MGALTphase::MGALTphase(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            phase* previousPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            HardwareModels::LaunchVehicle* myLaunchVehicle,
            missionoptions* myOptions) :
            TwoPointShootingLowThrustPhase::TwoPointShootingLowThrustPhase(name,
                journeyIndex, 
                phaseIndex,
                stageIndex,
                previousPhase, 
                Universe,
                mySpacecraft, 
                myLaunchVehicle,
                myOptions,
                10, //numStatesToPropagate
                9) //numMatchConstraints
        {
            //name our match point constraints
            this->stateVectorNames.push_back("virtual chemical fuel");
            this->stateVectorNames.push_back("virtual electric propellant");
            this->matchPointConstraintNames.push_back("virtual chemical fuel");
            this->matchPointConstraintNames.push_back("virtual electric propellant");
            this->matchPointConstraintStateIndex.push_back(8); //note we skipped the encode of boundary epoch, it does not get a match point constraint
            this->matchPointConstraintStateIndex.push_back(9);
            this->continuity_constraint_scale_factors(7) = 1.0 / this->myJourneyOptions->maximum_mass;
            this->continuity_constraint_scale_factors(8) = 1.0 / this->myJourneyOptions->maximum_mass;
            this->stateIndex_phase_propagation_variable = 10;
            this->stateVectorNames.push_back("phase flight time");

            //turn off derivatives of chemical fuel wrt time unless ACS tracking is enabled
            if (!this->myOptions->trackACS)
            {
                this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[7] = false;
            }
            //chemical fuel has no derivatives at all with respect to anything but itself
            for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[7][stateIndex] = false; //chemical fuel wrt state
            //ACS consumption has no dependency on control since it is just a function of time
            this->TruthTable_MatchConstraints_Derivative_wrt_Control[7] = false;

            //nothing depends on electric propellant except itself
            for (size_t constraintIndex = 0; constraintIndex < 8; ++constraintIndex)
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][9] = false; //chemical fuel wrt electric propellant
            this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates = this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates;


            //thrust step length derivatives
            std::fill(this->dThrustStepLengths_dPropagationVariable.begin(),
                this->dThrustStepLengths_dPropagationVariable.end(),
                1.0 / this->num_timesteps);

            //time vectors
            this->ForwardPropagationStepTimes.resize(this->num_timesteps / 2 + 1);
            this->BackwardPropagationStepTimes.resize(this->num_timesteps / 2 + 1);
            this->dForwardPropagationStepTimes_dPropagationVariable.resize(this->num_timesteps / 2 + 1, 1.0 / this->num_timesteps);
            this->dForwardPropagationStepTimes_dPropagationVariable.front() /= 2.0;
            this->dForwardPropagationStepTimes_dPropagationVariable.back() /= 2.0;
            this->dBackwardPropagationStepTimes_dPropagationVariable.resize(this->num_timesteps / 2 + 1, 1.0 / this->num_timesteps);
            this->dBackwardPropagationStepTimes_dPropagationVariable.front() /= 2.0;
            this->dBackwardPropagationStepTimes_dPropagationVariable.back() /= 2.0;

            this->dImpulseEpoch_dPropagationVariable.resize(this->num_timesteps, 0.0);

            for (size_t stepIndex = 0; stepIndex < this->num_timesteps; ++stepIndex)
            {
                this->dImpulseEpoch_dPropagationVariable[stepIndex] = (stepIndex + 1 - 0.5) / this->num_timesteps;
            }

            //STMs and such
            math::Matrix<double> I6(6, math::identity);

            this->ForwardSTM.resize(this->num_timesteps / 2, I6);
            this->BackwardSTM.resize(this->num_timesteps / 2, I6);
            this->Forward_dStatedIndependentVariable.resize(this->num_timesteps / 2, math::Matrix<double>(6, 1, 0.0));
            this->Backward_dStatedIndependentVariable.resize(this->num_timesteps / 2, math::Matrix<double>(6, 1, 0.0));

            //augmented state order is:
            //x
            //y
            //z
            //xdot (for MGALT this includes u_x)
            //ydot (for MGALT this includes u_y)
            //zdot (for MGALT this includes u_z)
            //mass
            //boundaryTime
            //chemical fuel
            //electric propellant
            //PhaseFlightTime
            //(later) P0
            //(later) u_command
            math::Matrix<double> I(this->numStatesToPropagate + 1, math::identity);

            this->ForwardAugmentedSTM.resize(this->num_timesteps / 2, I);
            this->BackwardAugmentedSTM.resize(this->num_timesteps / 2, I);
            this->MTM_STM_blocks_Forward.resize(this->num_timesteps, I);
            this->MTM_STM_blocks_Backward.resize(this->num_timesteps, I);
            MTM_STM_cumulative_chain_Forward.resize(this->num_timesteps, I);
            MTM_STM_cumulative_chain_Backward.resize(this->num_timesteps, I);

            //initial coast STMs - also used to propage the half-step closest to the boundary both for forward and backward propagation
            this->STM_initial_coast = I6;
            this->STM_terminal_coast = I6;
            this->STM_Augmented_initial_coast = I;
            this->STM_Augmented_terminal_coast = I;
            this->InitialCoast_dStatedIndependentVariable.resize(6, 1, 0.0);
            this->TerminalCoast_dStatedIndependentVariable.resize(6, 1, 0.0);

            //MTMs
            this->ForwardMTM.resize(this->num_timesteps / 2, I);
            this->BackwardMTM.resize(this->num_timesteps / 2, I);
            this->TCMTM = I;

            //acceleration model
            this->mySpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(this->myOptions,
                this->myJourneyOptions,
                this->myUniverse,
                this->Xdescriptions,
                this->mySpacecraft,
                14, // STM size
                false);//we do NOT need the central body because Kepler propagation takes care of it
            this->mySpacecraftAccelerationModel->setDutyCycle(this->PhaseDutyCycle);

            if (this->number_of_distance_constraints > 0)
            {
                this->Forward_cumulative_STM_MTM_triangle.resize(this->num_timesteps);
                this->Backward_cumulative_STM_MTM_triangle.resize(this->num_timesteps);
                this->Forward_boundary_STM_MTM.resize(this->num_timesteps);
                this->Backward_boundary_STM_MTM.resize(this->num_timesteps);
                for (size_t step = 0; step < this->num_timesteps / 2; ++step)
                {
                    for (size_t step2 = 0; step2 < step + 1; ++step2)
                    {
                        this->Forward_cumulative_STM_MTM_triangle[step].push_back(math::Matrix<double>(11, math::identity));
                        this->Backward_cumulative_STM_MTM_triangle[step].push_back(math::Matrix<double>(11, math::identity));
                    }
                }
            }

            //propagator and force model objects
            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                //forward
                
                this->ForwardPropagatorObjects.push_back(Astrodynamics::KeplerPropagatorTimeDomain(6));
                this->ForwardPropagatorObjects.back().setStateLeft(this->spacecraft_state_event_plus[step]);
                this->ForwardPropagatorObjects.back().setSTM(this->ForwardSTM[step]);
                this->ForwardPropagatorObjects.back().setdStatedIndependentVariable(this->Forward_dStatedIndependentVariable[step]);
                this->ForwardPropagatorObjects.back().setCentralBodyGM(this->myUniverse->mu);
                if (step == this->num_timesteps / 2 - 1)
                    this->ForwardPropagatorObjects.back().setStateRight(this->match_point_state_minus);
                else
                    this->ForwardPropagatorObjects.back().setStateRight(this->spacecraft_state_event_minus[step + 1]);
                this->ForwardPropagatorObjects.back().set_dPropagationTime_dIndependentVariable(&this->dForwardPropagationStepTimes_dPropagationVariable[step + 1]);

                this->ForwardManeuverObjects.push_back(ForwardBoundedImpulseManeuver(this->journeyIndex,
                    this->phaseIndex,
                    this->stageIndex,
                    this->myUniverse,
                    this->mySpacecraft,
                    this->myOptions));
                this->ForwardManeuverObjects.back().set_ThrustStepLength(this->ThrustStepLengths[step]);
                this->ForwardManeuverObjects.back().set_dThrustStepLength_dPropagationVariable(this->dThrustStepLengths_dPropagationVariable[step]);
                this->ForwardManeuverObjects.back().set_dImpulseEpoch_dPropagationVariable_pointer(this->dImpulseEpoch_dPropagationVariable[step]);
                this->ForwardManeuverObjects.back().set_LaunchDate(this->LaunchDate);
                this->ForwardManeuverObjects.back().set_max_thrust(this->max_thrust[step]);
                this->ForwardManeuverObjects.back().set_max_mass_flow_rate(this->max_mass_flow_rate[step]);
                this->ForwardManeuverObjects.back().set_Isp(this->Isp[step]);
                this->ForwardManeuverObjects.back().set_power(this->power[step]);
                this->ForwardManeuverObjects.back().set_active_power(this->active_power[step]);
                this->ForwardManeuverObjects.back().set_number_of_active_engines(this->number_of_active_engines[step]);
                this->ForwardManeuverObjects.back().set_ThrottleLevel(this->ThrottleLevel[step]);
                this->ForwardManeuverObjects.back().set_ThrottleLevelString(this->ThrottleLevelString[step]);
                this->ForwardManeuverObjects.back().set_DeltaV(this->stepDeltaV[step]);
                this->ForwardManeuverObjects.back().set_spacecraft_state_minus(this->spacecraft_state_event_minus[step]);
                this->ForwardManeuverObjects.back().set_spacecraft_state_plus(this->spacecraft_state_event_plus[step]);
                this->ForwardManeuverObjects.back().set_MTM(this->ForwardMTM[step]);
                this->ForwardManeuverObjects.back().set_DutyCycle(this->PhaseDutyCycle);
                this->ForwardManeuverObjects.back().set_Control(this->ControlVector[step]);
                this->ForwardManeuverObjects.back().set_Throttle(this->throttle[step]);
                this->ForwardManeuverObjects.back().set_AccelerationModel(this->mySpacecraftAccelerationModel);

                //backward
                size_t backstep = this->num_timesteps - 1 - step;
                this->BackwardPropagatorObjects.push_back(Astrodynamics::KeplerPropagatorTimeDomain(6));
                this->BackwardPropagatorObjects.back().setStateLeft(this->spacecraft_state_event_minus[backstep]);
                this->BackwardPropagatorObjects.back().setSTM(this->BackwardSTM[step]);
                this->BackwardPropagatorObjects.back().setdStatedIndependentVariable(this->Backward_dStatedIndependentVariable[step]);
                this->BackwardPropagatorObjects.back().setCentralBodyGM(this->myUniverse->mu);
                if (step == this->num_timesteps / 2 - 1)
                    this->BackwardPropagatorObjects.back().setStateRight(this->match_point_state_plus);
                else
                    this->BackwardPropagatorObjects.back().setStateRight(this->spacecraft_state_event_plus[backstep - 1]);
                this->BackwardPropagatorObjects.back().set_dPropagationTime_dIndependentVariable(&this->dBackwardPropagationStepTimes_dPropagationVariable[step + 1]);

                this->BackwardManeuverObjects.push_back(BackwardBoundedImpulseManeuver(this->journeyIndex,
                    this->phaseIndex,
                    this->stageIndex,
                    this->myUniverse,
                    this->mySpacecraft,
                    this->myOptions));
                this->BackwardManeuverObjects.back().set_ThrustStepLength(this->ThrustStepLengths[backstep]);
                this->BackwardManeuverObjects.back().set_dThrustStepLength_dPropagationVariable(this->dThrustStepLengths_dPropagationVariable[backstep]);
                this->BackwardManeuverObjects.back().set_dImpulseEpoch_dPropagationVariable_pointer(this->dImpulseEpoch_dPropagationVariable[backstep]);
                this->BackwardManeuverObjects.back().set_LaunchDate(this->LaunchDate);
                this->BackwardManeuverObjects.back().set_max_thrust(this->max_thrust[backstep]);
                this->BackwardManeuverObjects.back().set_max_mass_flow_rate(this->max_mass_flow_rate[backstep]);
                this->BackwardManeuverObjects.back().set_Isp(this->Isp[backstep]);
                this->BackwardManeuverObjects.back().set_power(this->power[backstep]);
                this->BackwardManeuverObjects.back().set_active_power(this->active_power[backstep]);
                this->BackwardManeuverObjects.back().set_number_of_active_engines(this->number_of_active_engines[backstep]);
                this->BackwardManeuverObjects.back().set_ThrottleLevel(this->ThrottleLevel[backstep]);
                this->BackwardManeuverObjects.back().set_ThrottleLevelString(this->ThrottleLevelString[backstep]);
                this->BackwardManeuverObjects.back().set_DeltaV(this->stepDeltaV[backstep]);
                this->BackwardManeuverObjects.back().set_spacecraft_state_minus(this->spacecraft_state_event_minus[backstep]);
                this->BackwardManeuverObjects.back().set_spacecraft_state_plus(this->spacecraft_state_event_plus[backstep]);
                this->BackwardManeuverObjects.back().set_MTM(this->BackwardMTM[step]);
                this->BackwardManeuverObjects.back().set_DutyCycle(this->PhaseDutyCycle);
                this->BackwardManeuverObjects.back().set_Control(this->ControlVector[backstep]);
                this->BackwardManeuverObjects.back().set_Throttle(this->throttle[backstep]);
                this->BackwardManeuverObjects.back().set_AccelerationModel(this->mySpacecraftAccelerationModel);
            }

            //forced coasts
            this->InitialCoastPropagatorObject = Astrodynamics::KeplerPropagatorTimeDomain(6);
            this->InitialCoastPropagatorObject.setStateLeft(this->state_after_initial_TCM);
            this->InitialCoastPropagatorObject.setStateRight(this->spacecraft_state_event_minus[0]);
            this->InitialCoastPropagatorObject.setSTM(this->STM_initial_coast);
            this->InitialCoastPropagatorObject.setdStatedIndependentVariable(this->InitialCoast_dStatedIndependentVariable);
            this->InitialCoastPropagatorObject.set_dPropagationTime_dIndependentVariable(&this->dForwardPropagationStepTimes_dPropagationVariable[0]);
            this->InitialCoastPropagatorObject.setCentralBodyGM(this->myUniverse->mu);

            this->TerminalCoastPropagatorObject = Astrodynamics::KeplerPropagatorTimeDomain(6);
            this->TerminalCoastPropagatorObject.setStateLeft(this->state_at_end_of_phase);
            this->TerminalCoastPropagatorObject.setStateRight(this->spacecraft_state_event_plus.back());
            this->TerminalCoastPropagatorObject.setSTM(this->STM_terminal_coast);
            this->TerminalCoastPropagatorObject.setdStatedIndependentVariable(this->TerminalCoast_dStatedIndependentVariable);
            this->TerminalCoastPropagatorObject.set_dPropagationTime_dIndependentVariable(&this->dBackwardPropagationStepTimes_dPropagationVariable[0]);
            this->TerminalCoastPropagatorObject.setCentralBodyGM(this->myUniverse->mu);

            //output stuff
            this->output_state.resize(this->numStatesToPropagate + 12 * 12, 1, 0.0);
        }//end constructor

        //******************************************calcbounds methods
        MGALTphase::~MGALTphase()
        {
            delete this->mySpacecraftAccelerationModel;
        }
        //end destructor

        //******************************************calcbounds methods
        void MGALTphase::calcbounds()
        {
            this->calcbounds_phase_left_boundary();

            this->calcbounds_phase_flight_time();

            this->calcbounds_phase_right_boundary();

            this->calcbounds_phase_main();
        }//end calcbounds()

        void MGALTphase::calcbounds_phase_main()
        {
            //base class - MGALTphase does not have anything special
            TwoPointShootingLowThrustPhase::calcbounds_phase_main();
        }//end calcbounds_phase_main()

        //******************************************process methods
        void MGALTphase::process_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {

            this->process_phase_left_boundary(X, Xindex, F, Findex, G, needG);

            this->process_phase_flight_time(X, Xindex, F, Findex, G, needG);

            this->process_phase_right_boundary(X, Xindex, F, Findex, G, needG);

            this->process_virtual_propellant_tanks(X, Xindex, F, Findex, G, needG);

            this->process_control(X, Xindex, F, Findex, G, needG);

            this->process_phase_main(X, Xindex, F, Findex, G, needG);

            this->process_match_point_constraints(X, Xindex, F, Findex, G, needG);
            
            if (this->number_of_distance_constraints > 0)
                this->process_distance_constraints(X, Xindex, F, Findex, G, needG);

            this->process_deltav_contribution(X, Xindex, F, Findex, G, needG);
        }

        void MGALTphase::process_phase_flight_time(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //base class
            phase::process_phase_flight_time(X, Xindex, F, Findex, G, needG);

            for (size_t step = 0; step < this->num_timesteps; ++step)
                this->ThrustStepLengths[step] = (this->PhaseFlightTime - this->InitialCoastDuration - this->TerminalCoastDuration) 
                * this->dThrustStepLengths_dPropagationVariable[step];

            for (size_t step = 0; step < this->num_timesteps / 2 + 1; ++step)
            {
                this->ForwardPropagationStepTimes[step] = (this->PhaseFlightTime - this->InitialCoastDuration - this->TerminalCoastDuration)
                    * this->dForwardPropagationStepTimes_dPropagationVariable[step];
                this->BackwardPropagationStepTimes[step] = (this->PhaseFlightTime - this->InitialCoastDuration - this->TerminalCoastDuration)
                    * this->dBackwardPropagationStepTimes_dPropagationVariable[step];
            }
        }//end process_phase_flight_time()

        void MGALTphase::process_forward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 0: initial tank state
            //set the tank states
            this->state_at_beginning_of_phase(8) = this->chemical_fuel_used; //TCM propellant, usually zero
            this->state_at_beginning_of_phase(9) = 0.0; //electric propellant

            //Step 1: initial coast
            this->InitialCoastPropagatorObject.propagate(this->InitialCoastDuration + this->ForwardPropagationStepTimes[0], needG);
            this->spacecraft_state_event_minus.front()(6) = this->state_after_initial_TCM(6);
            this->spacecraft_state_event_minus.front()(7) = this->state_after_initial_TCM(7) + this->InitialCoastDuration + this->ForwardPropagationStepTimes[0];
            this->spacecraft_state_event_minus.front()(8) = this->state_after_initial_TCM(8);
            this->spacecraft_state_event_minus.front()(9) = this->state_after_initial_TCM(9);

            //Step 2: loop over forward steps
            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                //Step 2.1: compute the maneuver
                this->ForwardManeuverObjects[step].process_maneuver(needG);

                //Step 2.2: propagate
                this->ForwardPropagatorObjects[step].propagate(this->ForwardPropagationStepTimes[step + 1], needG);

                //Step 2.3: update the epoch and mass
                if (step == this->num_timesteps / 2 - 1)
                {
                    this->match_point_state_minus(6) = this->spacecraft_state_event_plus[step](6);
                    this->match_point_state_minus(7) = this->spacecraft_state_event_plus[step](7) + this->ForwardPropagationStepTimes[step + 1];
                    this->match_point_state_minus(8) = this->spacecraft_state_event_plus[step](8);
                    this->match_point_state_minus(9) = this->spacecraft_state_event_plus[step](9);
                }
                else
                {

                    this->spacecraft_state_event_minus[step + 1](6) = this->spacecraft_state_event_plus[step](6);
                    this->spacecraft_state_event_minus[step + 1](7) = this->spacecraft_state_event_plus[step](7) + this->ForwardPropagationStepTimes[step + 1];
                    this->spacecraft_state_event_minus[step + 1](8) = this->spacecraft_state_event_plus[step](8);
                    this->spacecraft_state_event_minus[step + 1](9) = this->spacecraft_state_event_plus[step](9);
                }
            }
        }//end process_forward_half_phase

        void MGALTphase::process_backward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 0: tank states
            this->state_at_end_of_phase(8) = this->virtual_chemical_fuel_used;
            this->state_at_end_of_phase(9) = this->virtual_electric_propellant_used;

            //Step 1: terminal coast
            this->TerminalCoastPropagatorObject.propagate(-this->TerminalCoastDuration - this->BackwardPropagationStepTimes[0], needG);
            this->spacecraft_state_event_plus.back()(6) = this->state_at_end_of_phase(6);
            this->spacecraft_state_event_plus.back()(7) = this->state_at_end_of_phase(7) - (this->TerminalCoastDuration + this->BackwardPropagationStepTimes[0]);
            this->spacecraft_state_event_plus.back()(8) = this->state_at_end_of_phase(8);
            this->spacecraft_state_event_plus.back()(9) = this->state_at_end_of_phase(9);

            //Step 2: loop over backward steps
            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                size_t backstep = this->num_timesteps - 1 - step;

                //Step 2.1: compute the maneuver
                this->BackwardManeuverObjects[step].process_maneuver(needG);

                //Step 2.2: propagate
                this->BackwardPropagatorObjects[step].propagate(-this->BackwardPropagationStepTimes[step + 1], needG);

                //Step 2.3: update the epoch, mass, and tanks
                if (step == this->num_timesteps / 2 - 1)
                {
                    this->match_point_state_plus(6) = this->spacecraft_state_event_minus[backstep](6);
                    this->match_point_state_plus(7) = this->spacecraft_state_event_minus[backstep](7) - this->BackwardPropagationStepTimes[step + 1];
                    this->match_point_state_plus(8) = this->spacecraft_state_event_minus[backstep](8);
                    this->match_point_state_plus(9) = this->spacecraft_state_event_minus[backstep](9);
                }
                else
                {

                    this->spacecraft_state_event_plus[backstep - 1](6) = this->spacecraft_state_event_minus[backstep](6);
                    this->spacecraft_state_event_plus[backstep - 1](7) = this->spacecraft_state_event_minus[backstep](7) - this->BackwardPropagationStepTimes[step + 1];
                    this->spacecraft_state_event_plus[backstep - 1](8) = this->spacecraft_state_event_minus[backstep](8);
                    this->spacecraft_state_event_plus[backstep - 1](9) = this->spacecraft_state_event_minus[backstep](9);
                }
            }
        }//end process_backward_half_phase

        void MGALTphase::process_match_point_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: contruct the STM/MTM chains and the forward and backward HPTM
            if (needG)
            {
                //Step 1.1: Forward STM/MTM chains
                {
                    for (int step = this->num_timesteps / 2 - 1; step >= 0; --step)
                    {
                        //Step 1.1.1: Forward augmented STMs
                        //upper right 6x6 is the original STM
                        for (size_t i = 0; i < 6; ++i)
                            for (size_t j = 0; j < 6; ++j)
                                this->ForwardAugmentedSTM[step](i, j) = this->ForwardSTM[step](i, j);
                        //then we need to turn the upper right 6 x 1 into the Phi_t terms developed by Lantoine
                        for (size_t i = 0; i < 6; ++i)
                        {
                            this->ForwardAugmentedSTM[step](i, this->stateIndex_phase_propagation_variable) = this->Forward_dStatedIndependentVariable[step](i);
                        }
                        this->ForwardAugmentedSTM[step](7, this->stateIndex_phase_propagation_variable) = this->dForwardPropagationStepTimes_dPropagationVariable[step + 1];

                        //Step 1.1.2: STM-MTM blocks
                        this->MTM_STM_blocks_Forward[step] = this->ForwardAugmentedSTM[step] * this->ForwardMTM[step];

                        //Step 1.1.3: cumulative STM-MTM chains
                        //this is done by multiplying MTM_STM_blocks_Forward on the right by the NEXT step's MTM_STM_blocks_Forward
                        if (step == this->num_timesteps / 2 - 1)
                            this->MTM_STM_cumulative_chain_Forward[step] = this->MTM_STM_blocks_Forward[step];
                        else
                            this->MTM_STM_cumulative_chain_Forward[step] = this->MTM_STM_cumulative_chain_Forward[step + 1] * this->MTM_STM_blocks_Forward[step];
                    }
                }

                //Step 1.2: Backward STM/MTM chains
                {
                    for (int step = this->num_timesteps / 2 - 1; step >= 0; --step)
                    {
                        //Step 1.1.1: Backward augmented STMs
                        //upper right 6x6 is the original STM
                        for (size_t i = 0; i < 6; ++i)
                            for (size_t j = 0; j < 6; ++j)
                                this->BackwardAugmentedSTM[step](i, j) = this->BackwardSTM[step](i, j);
                        //then we need to turn the upper right 6 x 1 into the Phi_t terms developed by Lantoine
                        for (size_t i = 0; i < 6; ++i)
                        {
                            this->BackwardAugmentedSTM[step](i, this->stateIndex_phase_propagation_variable) = this->Backward_dStatedIndependentVariable[step](i);
                        }
                        this->BackwardAugmentedSTM[step](7, this->stateIndex_phase_propagation_variable) = -this->dBackwardPropagationStepTimes_dPropagationVariable[step + 1];

                        //Step 1.1.2: STM-MTM blocks
                        this->MTM_STM_blocks_Backward[step] = this->BackwardAugmentedSTM[step] * this->BackwardMTM[step];

                        //Step 1.1.3: cumulative STM-MTM chains
                        //this is done by multiplying MTM_STM_blocks_Backward on the right by the NEXT step's MTM_STM_blocks_Backward
                        if (step == this->num_timesteps / 2 - 1)
                            this->MTM_STM_cumulative_chain_Backward[step] = this->MTM_STM_blocks_Backward[step];
                        else
                            this->MTM_STM_cumulative_chain_Backward[step] = this->MTM_STM_cumulative_chain_Backward[step + 1] * this->MTM_STM_blocks_Backward[step];
                    }
                }


                //Step 1.3: do the same for the initial and terminal coasts
                {
                    for (size_t i = 0; i < 6; ++i)
                        for (size_t j = 0; j < 6; ++j)
                        {
                            this->STM_Augmented_initial_coast(i, j) = this->STM_initial_coast(i, j);
                            this->STM_Augmented_terminal_coast(i, j) = this->STM_terminal_coast(i, j);
                        }
                    //then we need to turn the upper right 6 x 1 into the Phi_t terms developed by Lantoine
                    //note since the initial coast does include part of the phase propagation time, we do need to track its time dependence
                    for (size_t i = 0; i < 6; ++i)
                    {
                        this->STM_Augmented_initial_coast(i, this->stateIndex_phase_propagation_variable) = this->InitialCoast_dStatedIndependentVariable(i);
                        this->STM_Augmented_terminal_coast(i, this->stateIndex_phase_propagation_variable) = this->TerminalCoast_dStatedIndependentVariable(i);
                    }
                    this->STM_Augmented_initial_coast(7, this->stateIndex_phase_propagation_variable) = this->dForwardPropagationStepTimes_dPropagationVariable[0];
                    this->STM_Augmented_terminal_coast(7, this->stateIndex_phase_propagation_variable) = -this->dBackwardPropagationStepTimes_dPropagationVariable[0];
                }

                //Step 1.4: construct the HPTMs
                this->ForwardHPTM = this->MTM_STM_cumulative_chain_Forward.front() * this->STM_Augmented_initial_coast * this->TCMTM;
                this->BackwardHPTM = this->MTM_STM_cumulative_chain_Backward.front() * this->STM_Augmented_terminal_coast;
            }

            //Step 2: call the base class
            TwoPointShootingPhase::process_match_point_constraints(X, Xindex, F, Findex, G, needG);

            if (needG)
            {
                //Step 3: derivatives with respect to control
                math::Matrix<double> dStateNow_dDecisionVariable(this->numStatesToPropagate + 1, 1, 0.0);
                math::Matrix<double> dMatchPointState_dDecisionVariable(this->numStatesToPropagate + 1, 1, 0.0);

                for (size_t step = 0; step < this->num_timesteps / 2; ++step)
                {
                    size_t Xindex, Gindex;
                    //Step 3.1: forward control
                    double dvmax = this->ForwardManeuverObjects[step].get_dvmax()_GETVALUE;
                    math::Matrix<double> dMassAfterManeuver_dThrottleComponents = this->ForwardManeuverObjects[step].get_dMassAfterManeuver_dThrottleComponents();
                    math::Matrix<double> dChemicalFuel_dThrottleComponents = this->ForwardManeuverObjects[step].get_dChemicalFuel_dThrottleComponents();
                    math::Matrix<double> dElectricPropellant_dThrottleComponents = this->ForwardManeuverObjects[step].get_dElectricPropellant_dThrottleComponents();

                    for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
                    {
                        dStateNow_dDecisionVariable.assign_zeros();
                        dStateNow_dDecisionVariable(3 + controlIndex) = dvmax;
                        dStateNow_dDecisionVariable(6) = dMassAfterManeuver_dThrottleComponents(controlIndex);
                        dStateNow_dDecisionVariable(8) = dChemicalFuel_dThrottleComponents(controlIndex);
                        dStateNow_dDecisionVariable(9) = dElectricPropellant_dThrottleComponents(controlIndex);

                        if (step == this->num_timesteps / 2 - 1)
                            dMatchPointState_dDecisionVariable = this->ForwardAugmentedSTM[step] * dStateNow_dDecisionVariable;
                        else
                            dMatchPointState_dDecisionVariable = this->MTM_STM_cumulative_chain_Forward[step + 1]
                            * (this->ForwardAugmentedSTM[step] * dStateNow_dDecisionVariable);


                        //match-point constraint
                        for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                        {
                            if (this->TruthTable_MatchConstraints_Derivative_wrt_Control[constraintIndex])
                            {
                                size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                                Gindex = this->G_indices_match_point_constraints_wrt_Control[constraintIndex][step][controlIndex];
                                Xindex = this->jGvar->operator[](Gindex);

                                G[Gindex] = -this->X_scale_factors->operator[](Xindex)
                                    * dMatchPointState_dDecisionVariable(stateIndex)
                                    * this->continuity_constraint_scale_factors(constraintIndex);
                            }
                        }
                    }

                    //Step 3.2: backward control
                    size_t backstep = this->num_timesteps - step - 1;
                    dvmax = this->BackwardManeuverObjects[step].get_dvmax()_GETVALUE;
                    math::Matrix<double> dMassBeforeManeuver_dThrottleComponents = this->BackwardManeuverObjects[step].get_dMassBeforeManeuver_dThrottleComponents();
                    dChemicalFuel_dThrottleComponents = this->BackwardManeuverObjects[step].get_dChemicalFuel_dThrottleComponents();
                    dElectricPropellant_dThrottleComponents = this->BackwardManeuverObjects[step].get_dElectricPropellant_dThrottleComponents();

                    for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
                    {
                        dStateNow_dDecisionVariable.assign_zeros();
                        dStateNow_dDecisionVariable(3 + controlIndex) = dvmax;
                        dStateNow_dDecisionVariable(6) = dMassBeforeManeuver_dThrottleComponents(controlIndex);
                        dStateNow_dDecisionVariable(8) = dChemicalFuel_dThrottleComponents(controlIndex);
                        dStateNow_dDecisionVariable(9) = dElectricPropellant_dThrottleComponents(controlIndex);

                        if (step == this->num_timesteps / 2 - 1)
                            dMatchPointState_dDecisionVariable = this->MTM_STM_cumulative_chain_Backward[step] * dStateNow_dDecisionVariable;
                        else
                            dMatchPointState_dDecisionVariable = this->MTM_STM_cumulative_chain_Backward[step]
                            * dStateNow_dDecisionVariable;


                        //match-point constraint
                        for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                        {
                            if (this->TruthTable_MatchConstraints_Derivative_wrt_Control[constraintIndex])
                            {
                                size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                                Gindex = this->G_indices_match_point_constraints_wrt_Control[constraintIndex][backstep][controlIndex];
                                Xindex = this->jGvar->operator[](Gindex);

                                G[Gindex] = -this->X_scale_factors->operator[](Xindex)
                                    * dMatchPointState_dDecisionVariable(stateIndex)
                                    * this->continuity_constraint_scale_factors(constraintIndex);
                            }
                        }
                    }
                }


                //Step 4: derivatives with respect to propagation time
                {
                    dStateNow_dDecisionVariable.assign_zeros();
                    dStateNow_dDecisionVariable(this->stateIndex_phase_propagation_variable) = 1.0;


                    //Forward
                    dMatchPointState_dDecisionVariable = this->ForwardHPTM * dStateNow_dDecisionVariable;

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                    {
                        size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                        size_t Gindex = this->G_indices_match_point_constraints_wrt_PhaseFlightTime[constraintIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += -this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }

                    //Backward
                    dMatchPointState_dDecisionVariable = this->BackwardHPTM * dStateNow_dDecisionVariable;

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                    {
                        size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                        size_t Gindex = this->G_indices_match_point_constraints_wrt_PhaseFlightTime[constraintIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }
                }

                //Step 5: virtual chemical fuel
                {
                    size_t Gindex = this->Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel;
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->continuity_constraint_scale_factors(7);
                }

                //Step 6: virtual electric propellant
                {
                    size_t Gindex = this->Gindex_dVirtualElectricPropellantConstraint_dVirtualElectricPropellant;
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->continuity_constraint_scale_factors(8);
                }
            }//end derivatives
        }//end process_match_point_constraints



        void MGALTphase::process_distance_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: apply the distance constraints
            for (size_t step = 0; step < this->num_timesteps; ++step)
            {
                for (size_t body = 0; body < this->number_of_distance_constraints; ++body)
                {
                    //Step 1.1: get the position vector of the spacecraft relative to the body
                    if (std::get<0>(this->distance_constraint_definitions[body]) == -2)
                    {
                        for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                        {
                            this->distance_constraint_relative_position[step][body](stateIndex) = this->spacecraft_state_event_minus[step](stateIndex);
                            this->distance_constraint_body_position_time_derivatives[step][body](stateIndex) = 0.0;
                        }
                    }
                    else
                    {
                        doubleType body_state[12];
                        this->myUniverse->bodies[std::get<0>(this->distance_constraint_definitions[body]) - 1].locate_body(
                            this->spacecraft_state_event_minus[step](7),
                            body_state,
                            (needG),
                            *this->myOptions);
                        for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                        {
                            this->distance_constraint_relative_position[step][body](stateIndex) = this->spacecraft_state_event_minus[step](stateIndex) - body_state[stateIndex];
                            this->distance_constraint_body_position_time_derivatives[step][body](stateIndex) = body_state[6 + stateIndex] _GETVALUE;
                        }
                    }

                    //Step 1.2: apply the constraint
                    this->distance_from_body[step][body] = this->distance_constraint_relative_position[step][body].norm();
                    F[Findex] = (this->distance_from_body[step][body] - std::get<2>(this->distance_constraint_definitions[body])) / this->myUniverse->LU;
                    ++Findex;
                }
            }

            //Step 2: derivatives, if appropriate
            if (needG)
            {
                //Step 2.0: create temporary matrices
                math::Matrix<double> dStateNow_dDecisionVariable(this->numStatesToPropagate + 1, 1, 0.0);
                math::Matrix<double> dConstraintPointState_dDecisionVariable(this->numStatesToPropagate + 1, 1, 0.0);

                //Step 2.1: build up the forward STM/MTM cumulative triangle
                for (size_t distancestep = 1; distancestep < this->num_timesteps / 2; ++distancestep)
                {
                    this->Forward_cumulative_STM_MTM_triangle[distancestep][distancestep - 1] = this->MTM_STM_blocks_Forward[distancestep - 1];
                    this->Backward_cumulative_STM_MTM_triangle[distancestep][distancestep - 1] = this->MTM_STM_blocks_Backward[distancestep - 1];
                    for (int cstep = distancestep - 1; cstep > 0; --cstep)
                    {
                        this->Forward_cumulative_STM_MTM_triangle[distancestep][cstep - 1] = this->Forward_cumulative_STM_MTM_triangle[distancestep][cstep] * this->MTM_STM_blocks_Forward[cstep - 1];
                        this->Backward_cumulative_STM_MTM_triangle[distancestep][cstep - 1] = this->Backward_cumulative_STM_MTM_triangle[distancestep][cstep] * this->MTM_STM_blocks_Backward[cstep - 1];
                    }
                    this->Forward_boundary_STM_MTM[distancestep] = this->Forward_cumulative_STM_MTM_triangle[distancestep][0] * this->STM_Augmented_initial_coast;
                    this->Backward_boundary_STM_MTM[distancestep] = this->Backward_cumulative_STM_MTM_triangle[distancestep][0] * this->STM_Augmented_terminal_coast;
                }
                this->Forward_boundary_STM_MTM[0] = this->STM_Augmented_initial_coast;
                this->Backward_boundary_STM_MTM[0] = this->STM_Augmented_terminal_coast;

                //Step 2.2: with respect to left boundary
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase_wrt_Time = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();//Xindex, stateIndex, derivative value
                for (size_t distancestep = 0; distancestep < this->num_timesteps / 2; ++distancestep)
                {                    
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfLeftBoundaryByVariable.size(); ++varIndex)
                    {
                        dStateNow_dDecisionVariable.assign_zeros();

                        for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByVariable[varIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->DerivativesOfLeftBoundaryByVariable[varIndex][entryIndex];
                            size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase[dIndex]);

                            dStateNow_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateBeforePhase[dIndex]);
                        }

                        dConstraintPointState_dDecisionVariable = this->Forward_boundary_STM_MTM[distancestep] * dStateNow_dDecisionVariable;

                        for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                        {
                            size_t Gindex = this->G_indices_distance_constraints_wrt_LeftBoundaryState[distancestep][bodyIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                / this->distance_from_body[distancestep][bodyIndex] _GETVALUE
                                * (   this->distance_constraint_relative_position[distancestep][bodyIndex](0) * dConstraintPointState_dDecisionVariable(0)
                                    + this->distance_constraint_relative_position[distancestep][bodyIndex](1) * dConstraintPointState_dDecisionVariable(1)
                                    + this->distance_constraint_relative_position[distancestep][bodyIndex](2) * dConstraintPointState_dDecisionVariable(2)) _GETVALUE
                                / this->myUniverse->LU;
                        }
                    }

                    //with respect to left boundary time variables
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfLeftBoundaryByTimeVariable.size(); ++varIndex)
                    {
                        dStateNow_dDecisionVariable.assign_zeros();

                        for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByTimeVariable[varIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->DerivativesOfLeftBoundaryByTimeVariable[varIndex][entryIndex];
                            size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);

                            dStateNow_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);
                        }

                        dConstraintPointState_dDecisionVariable = this->Forward_boundary_STM_MTM[distancestep] * dStateNow_dDecisionVariable;

                        for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                        {
                            size_t Gindex = this->G_indices_distance_constraints_wrt_LeftBoundaryTime[distancestep][bodyIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                / this->distance_from_body[distancestep][bodyIndex] _GETVALUE
                                * (this->distance_constraint_relative_position[distancestep][bodyIndex](0) * dConstraintPointState_dDecisionVariable(0)
                                    + this->distance_constraint_relative_position[distancestep][bodyIndex](1) * dConstraintPointState_dDecisionVariable(1)
                                    + this->distance_constraint_relative_position[distancestep][bodyIndex](2) * dConstraintPointState_dDecisionVariable(2)) _GETVALUE
                                / this->myUniverse->LU;
                        }
                    }
                }

                //Step 2.3: with respect to right boundary
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterPhase = this->myArrivalEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterPhase_wrt_Time = this->myArrivalEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value
                for (size_t distancestep = 0; distancestep < this->num_timesteps / 2; ++distancestep)
                {
                    size_t distancebackstep = this->num_timesteps - distancestep - 1;
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfRightBoundaryByVariable.size(); ++varIndex)
                    {
                        dStateNow_dDecisionVariable.assign_zeros();

                        for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByVariable[varIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->DerivativesOfRightBoundaryByVariable[varIndex][entryIndex];
                            size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase[dIndex]);

                            dStateNow_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateAfterPhase[dIndex]);
                        }

                        dConstraintPointState_dDecisionVariable = this->Backward_boundary_STM_MTM[distancestep] * dStateNow_dDecisionVariable;

                        for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                        {
                            size_t Gindex = this->G_indices_distance_constraints_wrt_RightBoundaryState[distancebackstep][bodyIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                / this->distance_from_body[distancebackstep][bodyIndex] _GETVALUE
                                * (this->distance_constraint_relative_position[distancebackstep][bodyIndex](0) * dConstraintPointState_dDecisionVariable(0)
                                    + this->distance_constraint_relative_position[distancebackstep][bodyIndex](1) * dConstraintPointState_dDecisionVariable(1)
                                    + this->distance_constraint_relative_position[distancebackstep][bodyIndex](2) * dConstraintPointState_dDecisionVariable(2)) _GETVALUE
                                / this->myUniverse->LU;
                        }
                    }

                    //with respect to right boundary time variables
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfRightBoundaryByTimeVariable.size(); ++varIndex)
                    {
                        dStateNow_dDecisionVariable.assign_zeros();

                        for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByTimeVariable[varIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->DerivativesOfRightBoundaryByTimeVariable[varIndex][entryIndex];
                            size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);

                            dStateNow_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);
                        }

                        dConstraintPointState_dDecisionVariable = this->Backward_boundary_STM_MTM[distancestep] * dStateNow_dDecisionVariable;

                        for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                        {
                            size_t Gindex = this->G_indices_distance_constraints_wrt_RightBoundaryTime[distancebackstep][bodyIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                / this->distance_from_body[distancebackstep][bodyIndex] _GETVALUE
                                * (this->distance_constraint_relative_position[distancebackstep][bodyIndex](0) * dConstraintPointState_dDecisionVariable(0)
                                    + this->distance_constraint_relative_position[distancebackstep][bodyIndex](1) * dConstraintPointState_dDecisionVariable(1)
                                    + this->distance_constraint_relative_position[distancebackstep][bodyIndex](2) * dConstraintPointState_dDecisionVariable(2)) _GETVALUE
                                / this->myUniverse->LU;
                        }
                    }
                }

                for (size_t distancestep = 1; distancestep < this->num_timesteps / 2; ++distancestep)
                {
                    size_t distancebackstep = this->num_timesteps - distancestep - 1;

                    for (size_t controlstep = 0; controlstep < distancestep; ++controlstep)
                    {
                        //Step 2.4: with respect to forward control
                        double dvmax = this->ForwardManeuverObjects[controlstep].get_dvmax()_GETVALUE;
                        math::Matrix<double> dMassAfterManeuver_dThrottleComponents = this->ForwardManeuverObjects[controlstep].get_dMassAfterManeuver_dThrottleComponents();
                        math::Matrix<double> dChemicalFuel_dThrottleComponents = this->ForwardManeuverObjects[controlstep].get_dChemicalFuel_dThrottleComponents();
                        math::Matrix<double> dElectricPropellant_dThrottleComponents = this->ForwardManeuverObjects[controlstep].get_dElectricPropellant_dThrottleComponents();

                        for (size_t controlIndex = 0; controlIndex < this->num_controls; ++controlIndex)
                        {
                            dStateNow_dDecisionVariable.assign_zeros();
                            dStateNow_dDecisionVariable(3 + controlIndex) = dvmax;
                            dStateNow_dDecisionVariable(6) = dMassAfterManeuver_dThrottleComponents(controlIndex);
                            dStateNow_dDecisionVariable(8) = dChemicalFuel_dThrottleComponents(controlIndex);
                            dStateNow_dDecisionVariable(9) = dElectricPropellant_dThrottleComponents(controlIndex);

                            if (distancestep - controlstep == 1)
                            {
                                dConstraintPointState_dDecisionVariable = this->ForwardAugmentedSTM[controlstep]
                                    * dStateNow_dDecisionVariable;
                            }
                            else
                            {
                                dConstraintPointState_dDecisionVariable = this->Forward_cumulative_STM_MTM_triangle[distancestep][controlstep + 1] * this->ForwardAugmentedSTM[controlstep]
                                    * dStateNow_dDecisionVariable;
                            }

                            for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                            {
                                size_t Gindex = this->G_indices_distance_constraints_wrt_ForwardControl[distancestep][bodyIndex][controlstep][controlIndex];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                    / this->distance_from_body[distancestep][bodyIndex] _GETVALUE
                                    * (this->distance_constraint_relative_position[distancestep][bodyIndex](0) * dConstraintPointState_dDecisionVariable(0)
                                        + this->distance_constraint_relative_position[distancestep][bodyIndex](1) * dConstraintPointState_dDecisionVariable(1)
                                        + this->distance_constraint_relative_position[distancestep][bodyIndex](2) * dConstraintPointState_dDecisionVariable(2)) _GETVALUE
                                    / this->myUniverse->LU;
                            }
                        }//end loop over control index

                        //Step 2.5: with respect to backward control
                        size_t controlbackstep = this->num_timesteps - controlstep - 1;
                        dvmax = this->BackwardManeuverObjects[controlstep].get_dvmax()_GETVALUE;
                        math::Matrix<double> dMassBeforeManeuver_dThrottleComponents = this->BackwardManeuverObjects[controlstep].get_dMassBeforeManeuver_dThrottleComponents();
                        dChemicalFuel_dThrottleComponents = this->BackwardManeuverObjects[controlstep].get_dChemicalFuel_dThrottleComponents();
                        dElectricPropellant_dThrottleComponents = this->BackwardManeuverObjects[controlstep].get_dElectricPropellant_dThrottleComponents();

                        for (size_t controlIndex = 0; controlIndex < this->num_controls; ++controlIndex)
                        {
                            dStateNow_dDecisionVariable.assign_zeros();

                            dStateNow_dDecisionVariable(3 + controlIndex) = dvmax;

                            for (size_t ccIndex = 0; ccIndex < 3; ++ccIndex)
                            {
                                dStateNow_dDecisionVariable(3 + ccIndex) += (this->stepDeltaV[controlbackstep](ccIndex) / this->spacecraft_state_event_minus[controlbackstep](6)
                                    * dMassBeforeManeuver_dThrottleComponents(controlIndex))_GETVALUE;
                            }

                            dStateNow_dDecisionVariable(6) = dMassBeforeManeuver_dThrottleComponents(controlIndex);
                            dStateNow_dDecisionVariable(8) = dChemicalFuel_dThrottleComponents(controlIndex);
                            dStateNow_dDecisionVariable(9) = dElectricPropellant_dThrottleComponents(controlIndex);

                            if (distancestep - controlstep == 1)
                            {
                                dConstraintPointState_dDecisionVariable = this->BackwardAugmentedSTM[controlstep]
                                    * dStateNow_dDecisionVariable;
                            }
                            else
                            {
                                dConstraintPointState_dDecisionVariable = this->Backward_cumulative_STM_MTM_triangle[distancestep][controlstep + 1] * this->BackwardAugmentedSTM[controlstep]
                                    * dStateNow_dDecisionVariable;
                            }

                            for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                            {
                                size_t Gindex = this->G_indices_distance_constraints_wrt_BackwardControl[distancestep][bodyIndex][controlstep][controlIndex];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                G[Gindex] = -this->X_scale_factors->operator[](Xindex)
                                    / this->distance_from_body[distancebackstep][bodyIndex] _GETVALUE
                                    * (this->distance_constraint_relative_position[distancebackstep][bodyIndex](0) * dConstraintPointState_dDecisionVariable(0)
                                        + this->distance_constraint_relative_position[distancebackstep][bodyIndex](1) * dConstraintPointState_dDecisionVariable(1)
                                        + this->distance_constraint_relative_position[distancebackstep][bodyIndex](2) * dConstraintPointState_dDecisionVariable(2)) _GETVALUE
                                    / this->myUniverse->LU;
                            }
                        }//end loop over control index
                    }//end loop over control steps
                }//end loop over distance steps

                //Step 2.6: with respect to phase flight time
                for (size_t distancestep = 0; distancestep < this->num_timesteps / 2; ++distancestep)
                {
                    //forward
                    dStateNow_dDecisionVariable.assign_zeros();
                    dStateNow_dDecisionVariable(this->stateIndex_phase_propagation_variable) = 1.0;

                    dConstraintPointState_dDecisionVariable = this->Forward_boundary_STM_MTM[distancestep] * dStateNow_dDecisionVariable;

                    for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                    {
                        size_t Gindex = this->G_indices_distance_constraints_wrt_PhaseFlightTime[distancestep][bodyIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            / this->distance_from_body[distancestep][bodyIndex] _GETVALUE
                            * (this->distance_constraint_relative_position[distancestep][bodyIndex](0) * (dConstraintPointState_dDecisionVariable(0) - this->distance_constraint_body_position_time_derivatives[distancestep][bodyIndex](0))
                                + this->distance_constraint_relative_position[distancestep][bodyIndex](1) * (dConstraintPointState_dDecisionVariable(1) - this->distance_constraint_body_position_time_derivatives[distancestep][bodyIndex](1))
                                + this->distance_constraint_relative_position[distancestep][bodyIndex](2) * (dConstraintPointState_dDecisionVariable(2) - this->distance_constraint_body_position_time_derivatives[distancestep][bodyIndex](2))) _GETVALUE
                            / this->myUniverse->LU;
                    }

                    //backward
                    size_t distancebackstep = this->num_timesteps - distancestep - 1;

                    dConstraintPointState_dDecisionVariable = this->Backward_boundary_STM_MTM[distancestep] * dStateNow_dDecisionVariable;
                    for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                    {
                        size_t Gindex = this->G_indices_distance_constraints_wrt_PhaseFlightTime[distancebackstep][bodyIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            / this->distance_from_body[distancebackstep][bodyIndex] _GETVALUE
                            * (this->distance_constraint_relative_position[distancebackstep][bodyIndex](0) * (dConstraintPointState_dDecisionVariable(0) - this->distance_constraint_body_position_time_derivatives[distancebackstep][bodyIndex](0))
                                + this->distance_constraint_relative_position[distancebackstep][bodyIndex](1) * (dConstraintPointState_dDecisionVariable(1) - this->distance_constraint_body_position_time_derivatives[distancebackstep][bodyIndex](1))
                                + this->distance_constraint_relative_position[distancebackstep][bodyIndex](2) * (dConstraintPointState_dDecisionVariable(2) - this->distance_constraint_body_position_time_derivatives[distancebackstep][bodyIndex](2))) _GETVALUE
                            / this->myUniverse->LU;
                    }
                }
            }//end derivatives
        }//end process_distance_constraints

        void MGALTphase::process_deltav_contribution(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //compute the total delta-v
            this->PhaseTotalDeterministicDeltav = this->myDepartureEvent->get_DeterministicDeltav()
                + this->myArrivalEvent->get_DeterministicDeltav();

            for (size_t step = 0; step < this->num_timesteps; ++step)
                this->PhaseTotalDeterministicDeltav += this->stepDeltaV[step].norm();

            this->PhaseTotalStatisticalDeltav = this->myDepartureEvent->get_StatisticalDeltav()
                + this->myArrivalEvent->get_StatisticalDeltav()
                + this->initial_TCM_magnitude;

            //TODO: derivatives
            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV
                && needG)
            {
                //do stuff
            }
        }

        void MGALTphase::output(std::ofstream& outputfile,
            size_t& eventcount)
        {
            //Step 1: output the departure event
            this->myDepartureEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //Step 2: output the initial TCM if applicable
            this->output_initial_TCM(outputfile, eventcount);

            //Step 3: output the state at the middle of the initial coast, if applicable
            if (this->hasInitialCoast)
            {
                //Step 3.1: create a temporary output state and assign the initial coast propagator to it
                this->InitialCoastPropagatorObject.setStateRight(output_state);

                //Step 3.2: propagate
                this->InitialCoastPropagatorObject.propagate(this->InitialCoastDuration / 2.0, false);
                output_state(6) = this->state_after_initial_TCM(6);
                output_state(7) = this->state_after_initial_TCM(7) + this->InitialCoastDuration / 2.0;

                //Step 3.3: assign the initial coast propagator back to where it is supposed to go
                this->InitialCoastPropagatorObject.setStateRight(this->spacecraft_state_event_minus[0]);

                //Step 3.4: print
                phase::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "force-coast",//event_type
                    "deep-space",//event_location
                    this->InitialCoastDuration / 86400.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    output_state,//state
                    math::Matrix<doubleType>(3, 1, 0.0),//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    0.0,//dVmag
                    0.0,//Thrust
                    0.0,//Isp
                    0.0,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    0.0,
                    "none");//active_power)
            }

            //Step 4: output the maneuver events
            for (size_t step = 0; step < this->num_timesteps; ++step)
            {
                //Step 2.1: output the current step
                std::string event_type;
                if (this->max_thrust[step] > 1.0e-7 && this->throttle[step] > 5.0e-4)
                    event_type = "SFthrust";
                else
                    event_type = "coast";

                math::Matrix<doubleType> ThrustVector = this->ControlVector[step] * this->max_thrust[step] * 1000.0 / this->PhaseDutyCycle;

                phase::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    event_type,//event_type
                    "deep-space",//event_location
                    this->ThrustStepLengths[step] / 86400.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    this->spacecraft_state_event_minus[step],//state
                    this->stepDeltaV[step],//dV
                    ThrustVector,//ThrustVector
                    this->stepDeltaV[step].norm(),//dVmag
                    this->max_thrust[step] * 1000.0 / this->PhaseDutyCycle,//Thrust
                    this->Isp[step],//Isp
                    this->power[step],//AvailPower
                    this->max_mass_flow_rate[step] / this->PhaseDutyCycle,//mdot
                    this->number_of_active_engines[step],//number_of_active_engines
                    this->active_power[step],
                    this->ThrottleLevelString[step]);//active_power)

                //Step 2.2: if this is the middle step, output the match point
                if (step == this->num_timesteps / 2 - 1)
                {

                    phase::write_output_line(outputfile,//outputfile
                        eventcount,//eventcount
                        "match_point",//event_type
                        "deep-space",//event_location
                        0.0,// timestep_size,
                        -1,//flyby_altitude,
                        0,//BdotR
                        0,//BdotT
                        0,//angle1
                        0,//angle2
                        0,//C3
                        this->match_point_state_minus,//state
                        math::Matrix<doubleType>(3, 1, 0.0),//dV
                        math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                        0.0,//dVmag
                        0.0,//Thrust
                        0.0,//Isp
                        0.0,//AvailPower
                        0.0,//mdot
                        0,//number_of_active_engines
                        0.0,
                        this->ThrottleLevelString[step]);//active_power)
                }
            }

            //Step 5: output the state at the middle of the terminal coast, if applicable
            if (this->hasTerminalCoast)
            {
                //Step 5.1: create a temporary output state and assign the Terminal coast propagator to it
                this->TerminalCoastPropagatorObject.setStateRight(output_state);

                //Step 5.2: propagate
                this->TerminalCoastPropagatorObject.propagate(-this->TerminalCoastDuration / 2.0, false);
                output_state(6) = this->state_at_end_of_phase(6);
                output_state(7) = this->state_at_end_of_phase(7) - this->TerminalCoastDuration / 2.0;

                //Step 5.3: assign the Terminal coast propagator back to where it is supposed to go
                this->TerminalCoastPropagatorObject.setStateRight(this->spacecraft_state_event_plus.back());

                //Step 5.4: print
                phase::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "force-coast",//event_type
                    "deep-space",//event_location
                    this->TerminalCoastDuration / 86400.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    output_state,//state
                    math::Matrix<doubleType>(3, 1, 0.0),//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    0.0,//dVmag
                    0.0,//Thrust
                    0.0,//Isp
                    0.0,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    0.0,
                    "none");//active_power)
            }

            //Step 6: output the arrival event
            this->myArrivalEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
        }

        void MGALTphase::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
        {
            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //Step 1: output the departure event
            this->myDepartureEvent->output_ephemeris(outputfile);

            //Step 2: output the initial coast
            {
                //Step 2.1: temporarily assign the initial coast propagator to the output state
                this->InitialCoastPropagatorObject.setStateRight(output_state);

                //Step 2.2: propagate and print, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->InitialCoastDuration + this->ForwardPropagationStepTimes[0]_GETVALUE;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 2.2.1: propagate
                    this->InitialCoastPropagatorObject.propagate(timeToPropagate, false);
                    output_state(6) = this->state_after_initial_TCM(6);
                    output_state(7) = this->state_after_initial_TCM(7) + timeToPropagate;

                    //Step 2.2.2: convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //Step 2.2.3: print
                    this->write_ephemeris_line(outputfile,
                        output_state,
                        math::Matrix<doubleType>(3, 1, 0.0),//control vector
                        0.0,
                        0.0,
                        0.0,
                        0,
                        0.0,
                        "none");

                    //Step 2.2.4: increment propagatedEpoch
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution
                        ? this->EphemerisOutputResolution
                        : (totalPropagationTime - timeToPropagate);
                }

                //Step 2.3: assign the initial coast propagator back to where it is supposed to go
                this->InitialCoastPropagatorObject.setStateRight(this->spacecraft_state_event_minus[0]);
            }

            //Step 3: output the thrust arcs
            //forward step
            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                //Step 3.1: temporarily assign the propagator to the output state
                this->ForwardPropagatorObjects[step].setStateRight(output_state);

                //Step 3.2: propagate and print, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->ForwardPropagationStepTimes[step + 1]_GETVALUE;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 3.2.1: propagate
                    this->ForwardPropagatorObjects[step].propagate(timeToPropagate, false);
                    output_state(6) = this->spacecraft_state_event_plus[step](6);
                    output_state(7) = this->spacecraft_state_event_plus[step](7) + timeToPropagate;

                    //Step 3.2.2: convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //Step 3.2.3: print
                    this->write_ephemeris_line(outputfile,
                        output_state,
                        this->ControlVector[step],
                        this->max_thrust[step] * this->throttle[step] * 1000.0 / this->PhaseDutyCycle,
                        this->max_mass_flow_rate[step] * this->throttle[step] / this->PhaseDutyCycle,
                        this->Isp[step],
                        this->number_of_active_engines[step],
                        this->active_power[step],
                        this->ThrottleLevelString[step]);

                    //Step 3.2.4: increment propagatedEpoch
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution
                        ? this->EphemerisOutputResolution
                        : (totalPropagationTime - timeToPropagate);
                }

                //Step 3.3: reset the propagator to its original state vectors
                if (step == this->num_timesteps / 2 - 1)
                    this->ForwardPropagatorObjects[step].setStateRight(this->match_point_state_minus);
                else
                    this->ForwardPropagatorObjects[step].setStateRight(this->spacecraft_state_event_minus[step + 1]);
            }//end forward step SPK output

            //backward step
            for (int step = this->num_timesteps / 2 - 1; step >= 0; --step)
            {
                size_t backstep = this->num_timesteps - 1 - step;

                //Step 3.1: temporarily assign the propagator to the output state
                if (step == this->num_timesteps / 2 - 1)
                    this->BackwardPropagatorObjects[step].setStateLeft(this->match_point_state_plus);
                else
                    this->BackwardPropagatorObjects[step].setStateLeft(this->spacecraft_state_event_plus[backstep - 1]);
                this->BackwardPropagatorObjects[step].setStateRight(output_state);

                //Step 3.2: propagate and print, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->BackwardPropagationStepTimes[step + 1]_GETVALUE;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 3.2.1: propagate
                    this->BackwardPropagatorObjects[step].propagate(timeToPropagate, false);
                    if (step == this->num_timesteps / 2 - 1)
                    {
                        output_state(6) = this->match_point_state_plus(6);
                        output_state(7) = this->match_point_state_plus(7) + timeToPropagate;
                    }
                    else
                    {
                        output_state(6) = this->spacecraft_state_event_plus[backstep - 1](6);
                        output_state(7) = this->spacecraft_state_event_plus[backstep - 1](7) + timeToPropagate;
                    }

                    //Step 3.2.2: convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //Step 3.2.3: print
                    this->write_ephemeris_line(outputfile,
                        output_state,
                        this->ControlVector[backstep],
                        this->max_thrust[backstep] * this->throttle[backstep] * 1000.0 / this->PhaseDutyCycle,
                        this->max_mass_flow_rate[backstep] * this->throttle[backstep] / this->PhaseDutyCycle,
                        this->Isp[backstep],
                        this->number_of_active_engines[backstep],
                        this->active_power[backstep],
                        this->ThrottleLevelString[backstep]);

                    //Step 3.2.4: increment propagatedEpoch
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution
                        ? this->EphemerisOutputResolution
                        : (totalPropagationTime - timeToPropagate);
                }

                //Step 3.3: reset the propagator to its original state vectors
                this->BackwardPropagatorObjects[step].setStateLeft(this->spacecraft_state_event_minus[backstep]);
                if (step == this->num_timesteps / 2 - 1)
                    this->BackwardPropagatorObjects[step].setStateRight(this->match_point_state_plus);
                else
                    this->BackwardPropagatorObjects[step].setStateRight(this->spacecraft_state_event_plus[backstep - 1]);
            }//end backward step SPK output

            //Step 4: output the terminal coast
            {
                //Step 4.1: temporarily assign the initial coast propagator to the output state
                math::Matrix<doubleType> output_state(8, 1, 0.0);
                this->TerminalCoastPropagatorObject.setStateLeft(this->spacecraft_state_event_plus.back());
                this->TerminalCoastPropagatorObject.setStateRight(output_state);

                //Step 4.2: propagate and print, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->TerminalCoastDuration + this->BackwardPropagationStepTimes[0]_GETVALUE;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 4.2.1: propagate
                    this->TerminalCoastPropagatorObject.propagate(timeToPropagate, false);
                    output_state(6) = this->spacecraft_state_event_plus.back()(6);
                    output_state(7) = this->spacecraft_state_event_plus.back()(7) + timeToPropagate;

                    //Step 4.2.2: convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //Step 4.2.3: print
                    this->write_ephemeris_line(outputfile,
                        output_state,
                        math::Matrix<doubleType>(3, 1, 0.0),//control vector
                        0.0,
                        0.0,
                        0.0,
                        0,
                        0.0,
                        "none");

                    //Step 4.2.4: increment propagatedEpoch
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution
                        ? this->EphemerisOutputResolution
                        : (totalPropagationTime - timeToPropagate);
                }

                //Step 4.3: assign the terminal coast propagator back to where it is supposed to go
                this->TerminalCoastPropagatorObject.setStateLeft(this->state_at_end_of_phase);
                this->TerminalCoastPropagatorObject.setStateRight(this->spacecraft_state_event_plus.back());
            }

            //Step 5: output the arrival event
            this->myArrivalEvent->output_ephemeris(outputfile);
        }//end output_ephemeris()

        void MGALTphase::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //should write out a maneuver/target spec for the departure maneuver if there is one

            maneuver_spec_file << "MGALTphase does not yet write maneuver spec files" << std::endl;

            target_spec_file << "MGALTphase does not yet write target spec files" << std::endl;

            //should write out a maneuver/target spec for the arrival maneuver if there is one
        }//end output_maneuver_and_target_spec()
    }//end namespace Phases
}//end namespace EMTG