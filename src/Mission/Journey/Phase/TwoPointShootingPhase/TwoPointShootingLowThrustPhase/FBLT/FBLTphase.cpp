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

//EMTGv9 FBLTphase
//Jacob Englander 1-2-2018

#include "FBLTphase.h"

#include "PropagatorFactory.h"
#include "IntegrationSchemeFactory.h"

namespace EMTG
{
    namespace Phases
    {
        FBLTphase::FBLTphase(const std::string& name,
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
            this->stateIndex_phase_propagation_variable = 13;
            this->stateVectorNames.push_back("u_x");
            this->stateVectorNames.push_back("u_y");
            this->stateVectorNames.push_back("u_z");
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
            this->ForcedCoast_dStepSize_dPropagationVariable = 0.0;

            //time vectors
            this->ForwardPropagationStepTimes.resize(this->num_timesteps / 2);
            this->BackwardPropagationStepTimes.resize(this->num_timesteps / 2);
            this->dForwardPropagationStepTimes_dPropagationVariable.resize(this->num_timesteps / 2, 1.0 / this->num_timesteps);
            this->dBackwardPropagationStepTimes_dPropagationVariable.resize(this->num_timesteps / 2, 1.0 / this->num_timesteps);

            //STMs and such
            math::Matrix<double> I14(14, math::identity);

            this->ForwardSTM.resize(this->num_timesteps / 2, I14);
            this->BackwardSTM.resize(this->num_timesteps / 2, I14);
            this->Forward_dStatedIndependentVariable.resize(this->num_timesteps / 2, math::Matrix<double>(10, 2, 0.0));
            this->Backward_dStatedIndependentVariable.resize(this->num_timesteps / 2, math::Matrix<double>(10, 2, 0.0));

            //augmented state order is:
            //x
            //y
            //z
            //xdot
            //ydot
            //zdot
            //mass
            //boundaryTime
            //chemical fuel
            //electric propellant
            //PhaseFlightTime
            //u_x
            //u_y
            //u_z
            //(later) u_command
            //(later) P0

            this->ForwardAugmentedSTM.resize(this->num_timesteps / 2, I14);
            this->BackwardAugmentedSTM.resize(this->num_timesteps / 2, I14);
            this->ForwardStrippedAugmentedSTM.resize(this->num_timesteps / 2, I14);
            this->BackwardStrippedAugmentedSTM.resize(this->num_timesteps / 2, I14);
            this->STM_cumulative_chain_Forward.resize(this->num_timesteps, I14);
            this->STM_cumulative_chain_Backward.resize(this->num_timesteps, I14);
            this->STM_stripped_cumulative_chain_Forward.resize(this->num_timesteps, I14);
            this->STM_stripped_cumulative_chain_Backward.resize(this->num_timesteps, I14);

            //acceleration model object
            this->mySpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(this->myOptions,
                                                                                                 this->myJourneyOptions,
                                                                                                 this->myUniverse,
                                                                                                 this->Xdescriptions,
                                                                                                 this->mySpacecraft,
                                                                                                 14); // STM size
            this->mySpacecraftAccelerationModel->setDutyCycle(this->PhaseDutyCycle);

            //EOM
            this->myEOM.setSpacecraftAccelerationModel(this->mySpacecraftAccelerationModel);

            //integration scheme
            this->myIntegrationScheme = CreateIntegrationScheme(&this->myEOM, this->numStatesToPropagate, 14);

            //propagator objects
            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                //forward
                this->ForwardPropagatorObjects.push_back(CreatePropagator(this->myOptions,
                    this->myUniverse,                     
                    this->numStatesToPropagate,
                    14,
                    this->spacecraft_state_event_minus[step],
                    this->spacecraft_state_event_plus[step],
                    this->ForwardSTM[step],
                    this->Forward_dStatedIndependentVariable[step],
                    (Integration::Integrand*) &this->myEOM,
                    this->myIntegrationScheme,
                    &dThrustStepLengths_dPropagationVariable[step],
                    this->myJourneyOptions->override_integration_step_size
                        ? this->myJourneyOptions->integration_step_size
                        : this->myOptions->integration_time_step_size));

                //backward
                size_t backstep = this->num_timesteps - 1 - step;
                this->BackwardPropagatorObjects.push_back(CreatePropagator(this->myOptions, 
                    this->myUniverse,                      
                    this->numStatesToPropagate,
                    14,
                    this->spacecraft_state_event_plus[backstep],
                    this->spacecraft_state_event_minus[backstep],
                    this->BackwardSTM[step],
                    this->Backward_dStatedIndependentVariable[step],
                    (Integration::Integrand*) &this->myEOM,
                    this->myIntegrationScheme,
                    &dThrustStepLengths_dPropagationVariable[backstep],
                    this->myJourneyOptions->override_integration_step_size
                        ? this->myJourneyOptions->integration_step_size
                        : this->myOptions->integration_time_step_size));
            }

            //forced coast STMs and such
            //we need dStatedIndependentVariable *even if we don't have forced coasts*
            this->InitialCoast_dStatedIndependentVariable.resize(10, 2, 0.0);
            this->TerminalCoast_dStatedIndependentVariable.resize(10, 2, 0.0);

            if (this->hasInitialCoast)
            {
                this->STM_initial_coast = I14;
                this->STM_Augmented_initial_coast = I14;

                this->InitialCoastPropagatorObject = CreatePropagator(this->myOptions,
                    this->myUniverse,                     
                    this->numStatesToPropagate,
                    14,
                    this->state_after_initial_TCM,
                    this->spacecraft_state_event_minus[0],
                    this->STM_initial_coast,
                    this->InitialCoast_dStatedIndependentVariable,
                    (Integration::Integrand*) &this->myEOM,
                    this->myIntegrationScheme,
                    &this->ForcedCoast_dStepSize_dPropagationVariable,
                    this->myJourneyOptions->override_integration_step_size
                        ? this->myJourneyOptions->integration_step_size
                        : this->myOptions->integration_time_step_size);
            }

            if (this->hasTerminalCoast)
            {
                this->STM_terminal_coast = I14;
                this->STM_Augmented_terminal_coast = I14;

                this->TerminalCoastPropagatorObject = CreatePropagator(this->myOptions,
                    this->myUniverse,                    
                    this->numStatesToPropagate,
                    14,
                    this->state_at_end_of_phase,
                    this->spacecraft_state_event_plus.back(),
                    this->STM_terminal_coast,
                    this->TerminalCoast_dStatedIndependentVariable,
                    (Integration::Integrand*) &this->myEOM,
                    this->myIntegrationScheme,
                    &this->ForcedCoast_dStepSize_dPropagationVariable,
                    this->myJourneyOptions->override_integration_step_size
                        ? this->myJourneyOptions->integration_step_size
                        : this->myOptions->integration_time_step_size);
            }

            //MTMs
            this->TCMTM = I14;

            //distance constraint stuff
            if (this->number_of_distance_constraints > 0)
            {
                this->Forward_cumulative_STM_triangle.resize(this->num_timesteps);
                this->Backward_cumulative_STM_triangle.resize(this->num_timesteps);
                this->Forward_boundary_STM.resize(this->num_timesteps);
                this->Backward_boundary_STM.resize(this->num_timesteps);
                for (size_t step = 0; step < this->num_timesteps / 2; ++step)
                {
                    for (size_t step2 = 0; step2 < step + 1; ++step2)
                    {
                        this->Forward_cumulative_STM_triangle[step].push_back(math::Matrix<double>(14, math::identity));
                        this->Backward_cumulative_STM_triangle[step].push_back(math::Matrix<double>(14, math::identity));
                    }
                }
                this->Forward_cumulative_stripped_STM_triangle = this->Forward_cumulative_STM_triangle;
                this->Backward_cumulative_stripped_STM_triangle = this->Backward_cumulative_STM_triangle;
            }

            //output stuff
            this->output_state.resize(this->numStatesToPropagate + 14 * 14, 1, 0.0);
        }//end constructor

        FBLTphase::~FBLTphase()
        {
            //acceleration model
            delete this->mySpacecraftAccelerationModel;

            //integration scheme
            delete this->myIntegrationScheme;

            //step propagators delete themselves thanks to boost::ptr_vector

            //forced coast STMs and such
            if (this->hasInitialCoast)
            {
                delete this->InitialCoastPropagatorObject;
            }

            if (this->hasTerminalCoast)
            {
                delete this->TerminalCoastPropagatorObject;
            }
        }//end destructor

        //******************************************calcbounds methods
        void FBLTphase::calcbounds()
        {
            this->calcbounds_phase_left_boundary();

            this->calcbounds_phase_flight_time();

            this->calcbounds_phase_right_boundary();

            this->calcbounds_phase_main();
        }//end calcbounds()

        void FBLTphase::calcbounds_phase_main()
        {
            //base class - MGALTphase does not have anything special
            TwoPointShootingLowThrustPhase::calcbounds_phase_main();
        }//end calcbounds_phase_main()

         //******************************************process methods
        void FBLTphase::process_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //do we want to propagate STMs or not?
            this->total_number_of_states_to_integrate = needG
                ? this->numStatesToPropagate - 1 + 12 * 12
                : this->numStatesToPropagate - 1;

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

        void FBLTphase::process_phase_flight_time(const std::vector<doubleType>& X,
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

            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                this->ForwardPropagationStepTimes[step] = (this->PhaseFlightTime - this->InitialCoastDuration - this->TerminalCoastDuration)
                    * this->dForwardPropagationStepTimes_dPropagationVariable[step];
                this->BackwardPropagationStepTimes[step] = (this->PhaseFlightTime - this->InitialCoastDuration - this->TerminalCoastDuration)
                    * this->dBackwardPropagationStepTimes_dPropagationVariable[step];
            }
            
            //acceleration model needs to know the launch epoch
            this->mySpacecraftAccelerationModel->setLaunchEpoch(this->LaunchDate);

        }//end process_phase_flight_time()

        void FBLTphase::process_forward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 0: initial tank state
            //set the tank states
            this->state_after_initial_TCM(8) = this->chemical_fuel_used; //TCM propellant, usually zero
            this->state_after_initial_TCM(9) = 0.0; //electric propellant

            //Step 1: initial coast
            if (this->hasInitialCoast)
            {
                this->InitialCoastPropagatorObject->setCurrentEpoch(this->state_after_initial_TCM(7));
                this->InitialCoastPropagatorObject->setIndexOfEpochInStateVec(7);
                this->InitialCoastPropagatorObject->setCurrentIndependentVariable(this->state_after_initial_TCM(7));
                this->InitialCoastPropagatorObject->propagate(this->InitialCoastDuration, needG);
                this->spacecraft_state_event_minus.front()(7) = this->state_after_initial_TCM(7) + this->InitialCoastDuration;
            }
            else
                this->spacecraft_state_event_minus.front() = this->state_after_initial_TCM;

            //Step 2: loop over forward steps
            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                //Step 2.1: state copy, if applicable
                if (step > 0)
                    this->spacecraft_state_event_minus[step] = this->spacecraft_state_event_plus[step - 1];

                //Step 2.2: propagate
                this->ForwardPropagatorObjects[step].setCurrentEpoch(this->spacecraft_state_event_minus[step](7));
                this->ForwardPropagatorObjects[step].setIndexOfEpochInStateVec(7);
                this->ForwardPropagatorObjects[step].setCurrentIndependentVariable(this->spacecraft_state_event_minus[step](7));
                this->ForwardPropagatorObjects[step].propagate(this->ForwardPropagationStepTimes[step], this->ControlVector[step], needG);

                //Step 2.3: update the epoch
                this->spacecraft_state_event_plus[step](7) = this->spacecraft_state_event_minus[step](7) + this->ForwardPropagationStepTimes[step];
            }

            //Step 3: match point state
            this->match_point_state_minus = this->spacecraft_state_event_plus[this->num_timesteps / 2 - 1];
        }//end process_forward_half_phase

        void FBLTphase::process_backward_half_phase(const std::vector<doubleType>& X,
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
            if (this->hasTerminalCoast)
            {
                this->TerminalCoastPropagatorObject->setCurrentEpoch(this->state_at_end_of_phase(7));
                this->TerminalCoastPropagatorObject->setIndexOfEpochInStateVec(7);
                this->TerminalCoastPropagatorObject->setCurrentIndependentVariable(this->state_at_end_of_phase(7));
                this->TerminalCoastPropagatorObject->propagate(-this->TerminalCoastDuration, needG);
                this->spacecraft_state_event_plus.back()(7) = this->state_at_end_of_phase(7) - (this->TerminalCoastDuration);
            }
            else
                this->spacecraft_state_event_plus.back() = this->state_at_end_of_phase;

            //Step 2: loop over backward steps
            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                size_t backstep = this->num_timesteps - 1 - step;

                //Step 2.1: state copy, if applicable
                if (step > 0)
                    this->spacecraft_state_event_plus[backstep] = this->spacecraft_state_event_minus[backstep + 1];

                //Step 2.2: propagate
                this->BackwardPropagatorObjects[step].setCurrentEpoch(this->spacecraft_state_event_plus[backstep](7));
                this->BackwardPropagatorObjects[step].setIndexOfEpochInStateVec(7);
                this->BackwardPropagatorObjects[step].setCurrentIndependentVariable(this->spacecraft_state_event_plus[backstep](7));
                this->BackwardPropagatorObjects[step].propagate(-this->BackwardPropagationStepTimes[step], this->ControlVector[backstep], needG);

                //Step 2.3: update the epoch
                this->spacecraft_state_event_minus[backstep](7) = this->spacecraft_state_event_plus[backstep](7) - this->BackwardPropagationStepTimes[step];
            }

            //Step 3: match point state
            this->match_point_state_plus = this->spacecraft_state_event_minus[this->num_timesteps / 2];
        }//end process_backward_half_phase



        void FBLTphase::process_match_point_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: contruct the STM chains and the forward and backward HPTM
            if (needG)
            {
                //Step 1.1: Forward STM chains
                {
                    for (int step = this->num_timesteps / 2 - 1; step >= 0; --step)
                    {
                        //Step 1.1.1: Forward augmented STMs
                        //upper right 7x7 is the original STM
                        for (size_t i = 0; i < 7; ++i)
                        {
                            for (size_t j = 0; j < 7; ++j)
                                this->ForwardAugmentedSTM[step](i, j) = this->ForwardSTM[step](i, j);

                            //Phi_t terms
                            this->ForwardAugmentedSTM[step](i, 7) = this->ForwardSTM[step](i, 7);
                            this->ForwardAugmentedSTM[step](i, this->stateIndex_phase_propagation_variable) = this->ForwardSTM[step](i, 13);

                            //control terms
                            this->ForwardAugmentedSTM[step](i, 10) = this->ForwardSTM[step](i, 10);
                            this->ForwardAugmentedSTM[step](i, 11) = this->ForwardSTM[step](i, 11);
                            this->ForwardAugmentedSTM[step](i, 12) = this->ForwardSTM[step](i, 12);
                        }

                        //tanks
                        for (size_t i = 8; i < 10; ++i)
                        {
                            for (size_t j = 0; j < 7; ++j)
                                this->ForwardAugmentedSTM[step](i, j) = this->ForwardSTM[step](i, j);

                            //Phi_t terms
                            this->ForwardAugmentedSTM[step](i, 7) = this->ForwardSTM[step](i, 7);
                            this->ForwardAugmentedSTM[step](i, this->stateIndex_phase_propagation_variable) = this->ForwardSTM[step](i, 13);

                            //control terms
                            this->ForwardAugmentedSTM[step](i, 10) = this->ForwardSTM[step](i, 10);
                            this->ForwardAugmentedSTM[step](i, 11) = this->ForwardSTM[step](i, 11);
                            this->ForwardAugmentedSTM[step](i, 12) = this->ForwardSTM[step](i, 12);                            
                        }
                        

                        //epoch time with respect to propagation time
                        this->ForwardAugmentedSTM[step](7, this->stateIndex_phase_propagation_variable) = this->dForwardPropagationStepTimes_dPropagationVariable[step];

                        //Step 1.1.2: cumulative STM chains
                        //this is done by multiplying ForwardAugmentedSTM on the right by the NEXT step's ForwardAugmentedSTM
                        if (step == this->num_timesteps / 2 - 1)
                            this->STM_cumulative_chain_Forward[step] = this->ForwardAugmentedSTM[step];
                        else
                            this->STM_cumulative_chain_Forward[step] = this->STM_stripped_cumulative_chain_Forward[step + 1] * this->ForwardAugmentedSTM[step];

                        //create the stripped version of the cumulative chain for this step
                        this->STM_stripped_cumulative_chain_Forward[step].shallow_copy(this->STM_cumulative_chain_Forward[step]);
                        this->ForwardStrippedAugmentedSTM[step].shallow_copy(this->ForwardAugmentedSTM[step]);
                        for (size_t i = 0; i < 14; ++i)
                        {
                            this->STM_stripped_cumulative_chain_Forward[step](i, 10) = 0.0;
                            this->STM_stripped_cumulative_chain_Forward[step](i, 11) = 0.0;
                            this->STM_stripped_cumulative_chain_Forward[step](i, 12) = 0.0;
                            this->ForwardStrippedAugmentedSTM[step](i, 10) = 0.0;
                            this->ForwardStrippedAugmentedSTM[step](i, 11) = 0.0;
                            this->ForwardStrippedAugmentedSTM[step](i, 12) = 0.0;
                        }
                    }
                }

                //Step 2.2: Backward STM chains
                {
                    for (int step = this->num_timesteps / 2 - 1; step >= 0; --step)
                    {
                        //Step 2.1.1: Backward augmented STMs
                        //upper right 7x7 is the original STM
                        for (size_t i = 0; i < 7; ++i)
                        {
                            for (size_t j = 0; j < 7; ++j)
                                this->BackwardAugmentedSTM[step](i, j) = this->BackwardSTM[step](i, j);

                            //Phi_t terms
                            this->BackwardAugmentedSTM[step](i, 7) = this->BackwardSTM[step](i, 7);
                            this->BackwardAugmentedSTM[step](i, this->stateIndex_phase_propagation_variable) = this->BackwardSTM[step](i, 13);

                            //control terms
                            this->BackwardAugmentedSTM[step](i, 10) = this->BackwardSTM[step](i, 10);
                            this->BackwardAugmentedSTM[step](i, 11) = this->BackwardSTM[step](i, 11);
                            this->BackwardAugmentedSTM[step](i, 12) = this->BackwardSTM[step](i, 12);
                        }

                        //tanks
                        for (size_t i = 8; i < 10; ++i)
                        {
                            for (size_t j = 0; j < 7; ++j)
                                this->BackwardAugmentedSTM[step](i, j) = this->BackwardSTM[step](i, j);

                            //Phi_t terms
                            this->BackwardAugmentedSTM[step](i, 7) = this->BackwardSTM[step](i, 7);
                            this->BackwardAugmentedSTM[step](i, this->stateIndex_phase_propagation_variable) = this->BackwardSTM[step](i, 13);

                            //control terms
                            this->BackwardAugmentedSTM[step](i, 10) = this->BackwardSTM[step](i, 10);
                            this->BackwardAugmentedSTM[step](i, 11) = this->BackwardSTM[step](i, 11);
                            this->BackwardAugmentedSTM[step](i, 12) = this->BackwardSTM[step](i, 12);
                        }


                        //epoch time with respect to propagation time
                        this->BackwardAugmentedSTM[step](7, this->stateIndex_phase_propagation_variable) = -this->dBackwardPropagationStepTimes_dPropagationVariable[step];

                        //Step 2.1.2: cumulative STM chains
                        //this is done by multiplying BackwardAugmentedSTM on the right by the NEXT step's BackwardAugmentedSTM
                        if (step == this->num_timesteps / 2 - 1)
                            this->STM_cumulative_chain_Backward[step] = this->BackwardAugmentedSTM[step];
                        else
                            this->STM_cumulative_chain_Backward[step] = this->STM_stripped_cumulative_chain_Backward[step + 1] * this->BackwardAugmentedSTM[step];

                        //create the stripped version of the cumulative chain for this step
                        this->STM_stripped_cumulative_chain_Backward[step].shallow_copy(this->STM_cumulative_chain_Backward[step]);
                        this->BackwardStrippedAugmentedSTM[step].shallow_copy(this->BackwardAugmentedSTM[step]);

                        for (size_t i = 0; i < 14; ++i)
                        {
                            this->STM_stripped_cumulative_chain_Backward[step](i, 10) = 0.0;
                            this->STM_stripped_cumulative_chain_Backward[step](i, 11) = 0.0;
                            this->STM_stripped_cumulative_chain_Backward[step](i, 12) = 0.0;
                            this->BackwardStrippedAugmentedSTM[step](i, 10) = 0.0;
                            this->BackwardStrippedAugmentedSTM[step](i, 11) = 0.0;
                            this->BackwardStrippedAugmentedSTM[step](i, 12) = 0.0;
                        }
                    }
                }


                //Step 1.3: do the same for the initial coast
                if (this->hasInitialCoast)
                {
                    for (size_t i = 0; i < 7; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                            this->STM_Augmented_initial_coast(i, j) = this->STM_initial_coast(i, j);

                        //Phi_t terms
                        this->STM_Augmented_initial_coast(i, 7) = this->STM_initial_coast(i, 7);
                        this->STM_Augmented_initial_coast(i, this->stateIndex_phase_propagation_variable) = this->STM_initial_coast(i, 13);
                    }

                    //tanks
                    for (size_t i = 8; i < 10; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                            this->STM_Augmented_initial_coast(i, j) = this->STM_initial_coast(i, j);

                        //Phi_t terms
                        this->STM_Augmented_initial_coast(i, 7) = this->STM_initial_coast(i, 7);
                        this->STM_Augmented_initial_coast(i, this->stateIndex_phase_propagation_variable) = this->STM_initial_coast(i, 13);
                    }

                    //forward HPTM
                    this->ForwardHPTM = this->STM_stripped_cumulative_chain_Forward.front() * this->STM_Augmented_initial_coast * this->TCMTM;
                }
                else
                    this->ForwardHPTM = this->STM_stripped_cumulative_chain_Forward.front() * this->TCMTM;

                //Step 1.4: terminal coast
                if (this->hasTerminalCoast)
                {
                    for (size_t i = 0; i < 7; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                            this->STM_Augmented_terminal_coast(i, j) = this->STM_terminal_coast(i, j);

                        //Phi_t terms
                        this->STM_Augmented_terminal_coast(i, 7) = this->STM_terminal_coast(i, 7);
                        this->STM_Augmented_terminal_coast(i, this->stateIndex_phase_propagation_variable) = this->STM_terminal_coast(i, 13);
                    }

                    //tanks
                    for (size_t i = 8; i < 10; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                            this->STM_Augmented_terminal_coast(i, j) = this->STM_terminal_coast(i, j);

                        //Phi_t terms
                        this->STM_Augmented_terminal_coast(i, 7) = this->STM_terminal_coast(i, 7);
                        this->STM_Augmented_terminal_coast(i, this->stateIndex_phase_propagation_variable) = this->STM_terminal_coast(i, 13);
                    }

                    //backward HPTM
                    this->BackwardHPTM = this->STM_stripped_cumulative_chain_Backward.front() * this->STM_Augmented_terminal_coast;
                }
                else
                    this->BackwardHPTM = this->STM_stripped_cumulative_chain_Backward.front();
                 
            }

            //Step 2: call the base class
            TwoPointShootingPhase::process_match_point_constraints(X, Xindex, F, Findex, G, needG);

            if (needG)
            {
                //Step 3: derivatives with respect to control
                math::Matrix<double> dStateNow_dDecisionVariable(this->ForwardHPTM.get_n(), 1, 0.0);
                math::Matrix<double> dMatchPointState_dDecisionVariable(this->ForwardHPTM.get_n(), 1, 0.0);

                for (size_t step = 0; step < this->num_timesteps / 2; ++step)
                {
                    size_t Xindex, Gindex;
                    //Step 3.1: forward control

                    for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
                    {
                        dStateNow_dDecisionVariable.assign_zeros();
                        dStateNow_dDecisionVariable(10 + controlIndex) = 1.0;

                        dMatchPointState_dDecisionVariable = this->STM_cumulative_chain_Forward[step] * dStateNow_dDecisionVariable;


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

                    for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
                    {
                        dStateNow_dDecisionVariable.assign_zeros();
                        dStateNow_dDecisionVariable(10 + controlIndex) = -1.0;

                        dMatchPointState_dDecisionVariable = this->STM_cumulative_chain_Backward[step] * dStateNow_dDecisionVariable;


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



        void FBLTphase::process_distance_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: apply the distance constraints
            for (size_t step = 0; step < this->num_timesteps; ++step)
            {
                math::Matrix<doubleType>& SpacecraftState = step < this->num_timesteps / 2
                    ? this->spacecraft_state_event_plus[step]
                    : this->spacecraft_state_event_minus[step];

                for (size_t body = 0; body < this->number_of_distance_constraints; ++body)
                {
                    //Step 1.1: get the position vector of the spacecraft relative to the body
                    if (std::get<0>(this->distance_constraint_definitions[body]) == -2)
                    {
                        for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                        {
                            this->distance_constraint_relative_position[step][body](stateIndex) = SpacecraftState(stateIndex);
                            this->distance_constraint_body_position_time_derivatives[step][body](stateIndex) = 0.0;
                        }
                    }
                    else
                    {
                        doubleType body_state[12];
                        this->myUniverse->bodies[std::get<0>(this->distance_constraint_definitions[body]) - 1].locate_body(
                            SpacecraftState(7),
                            body_state,
                            (needG),
                            *this->myOptions);
                        for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                        {
                            this->distance_constraint_relative_position[step][body](stateIndex) = SpacecraftState(stateIndex) - body_state[stateIndex];
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
                math::Matrix<double> dStateNow_dDecisionVariable(this->ForwardHPTM.get_n(), 1, 0.0);
                math::Matrix<double> dConstraintPointState_dDecisionVariable(this->ForwardHPTM.get_n(), 1, 0.0);

                //Step 2.1: build up the STM cumulative triangle
                for (size_t distancestep = 0; distancestep < this->num_timesteps / 2; ++distancestep)
                {
                    for (int controlstep = distancestep; controlstep >= 0; --controlstep)
                    {
                        if (controlstep == distancestep)
                        {
                            this->Forward_cumulative_STM_triangle[distancestep][controlstep] = this->ForwardAugmentedSTM[controlstep];
                            this->Backward_cumulative_STM_triangle[distancestep][controlstep] = this->BackwardAugmentedSTM[controlstep];
                        }
                        else if (controlstep == distancestep - 1)
                        {
                            this->Forward_cumulative_stripped_STM_triangle[distancestep][controlstep + 1] = this->ForwardStrippedAugmentedSTM[controlstep + 1];
                            this->Forward_cumulative_STM_triangle[distancestep][controlstep] = this->Forward_cumulative_stripped_STM_triangle[distancestep][controlstep + 1] * this->ForwardAugmentedSTM[controlstep];
                            this->Backward_cumulative_stripped_STM_triangle[distancestep][controlstep + 1] = this->BackwardStrippedAugmentedSTM[controlstep + 1];
                            this->Backward_cumulative_STM_triangle[distancestep][controlstep] = this->Backward_cumulative_stripped_STM_triangle[distancestep][controlstep + 1] * this->BackwardAugmentedSTM[controlstep];
                        }
                        else
                        {
                            this->Forward_cumulative_stripped_STM_triangle[distancestep][controlstep + 1] = this->Forward_cumulative_stripped_STM_triangle[distancestep][controlstep + 2] * this->ForwardStrippedAugmentedSTM[controlstep + 1];
                            this->Forward_cumulative_STM_triangle[distancestep][controlstep] = this->Forward_cumulative_stripped_STM_triangle[distancestep][controlstep + 1] * this->ForwardAugmentedSTM[controlstep];
                            this->Backward_cumulative_stripped_STM_triangle[distancestep][controlstep + 1] = this->Backward_cumulative_stripped_STM_triangle[distancestep][controlstep + 2] * this->BackwardStrippedAugmentedSTM[controlstep + 1];
                            this->Backward_cumulative_STM_triangle[distancestep][controlstep] = this->Backward_cumulative_stripped_STM_triangle[distancestep][controlstep + 1] * this->BackwardAugmentedSTM[controlstep];
                        }
                    }

                    if (distancestep > 0)
                    {
                        this->Forward_cumulative_stripped_STM_triangle[distancestep][0] = this->Forward_cumulative_stripped_STM_triangle[distancestep][1] * this->ForwardStrippedAugmentedSTM[0];
                        this->Backward_cumulative_stripped_STM_triangle[distancestep][0] = this->Backward_cumulative_stripped_STM_triangle[distancestep][1] * this->BackwardStrippedAugmentedSTM[0];
                    }
                    else
                    {
                        this->Forward_cumulative_stripped_STM_triangle[distancestep][0] = this->ForwardStrippedAugmentedSTM[0];
                        this->Backward_cumulative_stripped_STM_triangle[distancestep][0] = this->BackwardStrippedAugmentedSTM[0];
                    }

                    this->Forward_boundary_STM[distancestep] = this->Forward_cumulative_stripped_STM_triangle[distancestep][0];
                    this->Backward_boundary_STM[distancestep] = this->Backward_cumulative_stripped_STM_triangle[distancestep][0];

                    if (this->hasInitialCoast)
                    {
                        this->Forward_boundary_STM[distancestep] *= this->STM_Augmented_initial_coast;
                    }
                    if (this->hasTerminalCoast)
                    {
                        this->Backward_boundary_STM[distancestep] *= this->STM_Augmented_terminal_coast;
                    }
                }

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

                        dConstraintPointState_dDecisionVariable = this->Forward_boundary_STM[distancestep] * dStateNow_dDecisionVariable;

                        for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                        {
                            size_t Gindex = this->G_indices_distance_constraints_wrt_LeftBoundaryState[distancestep][bodyIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                / this->distance_from_body[distancestep][bodyIndex] _GETVALUE
                                * (this->distance_constraint_relative_position[distancestep][bodyIndex](0) * dConstraintPointState_dDecisionVariable(0)
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

                        dConstraintPointState_dDecisionVariable = this->Forward_boundary_STM[distancestep] * dStateNow_dDecisionVariable;

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

                        dConstraintPointState_dDecisionVariable = this->Backward_boundary_STM[distancestep] * dStateNow_dDecisionVariable;

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

                        dConstraintPointState_dDecisionVariable = this->Backward_boundary_STM[distancestep] * dStateNow_dDecisionVariable;

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
                        for (size_t controlIndex = 0; controlIndex < this->num_controls; ++controlIndex)
                        {
                            dStateNow_dDecisionVariable.assign_zeros();
                            dStateNow_dDecisionVariable(10 + controlIndex) = 1.0;

                            dConstraintPointState_dDecisionVariable = this->Forward_cumulative_STM_triangle[distancestep][controlstep] * dStateNow_dDecisionVariable;

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
                        
                        for (size_t controlIndex = 0; controlIndex < this->num_controls; ++controlIndex)
                        {
                            dStateNow_dDecisionVariable.assign_zeros();
                            dStateNow_dDecisionVariable(10 + controlIndex) = 1.0;

                            dConstraintPointState_dDecisionVariable = this->Backward_cumulative_STM_triangle[distancestep][controlstep] * dStateNow_dDecisionVariable;

                            for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                            {
                                size_t Gindex = this->G_indices_distance_constraints_wrt_BackwardControl[distancestep][bodyIndex][controlstep][controlIndex];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                G[Gindex] = this->X_scale_factors->operator[](Xindex)
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

                    dConstraintPointState_dDecisionVariable = this->Forward_boundary_STM[distancestep] * dStateNow_dDecisionVariable;

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

                    dConstraintPointState_dDecisionVariable = this->Backward_boundary_STM[distancestep] * dStateNow_dDecisionVariable;
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

        void FBLTphase::process_deltav_contribution(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //compute the total delta-v
            this->PhaseTotalDeterministicDeltav = this->myDepartureEvent->get_DeterministicDeltav()
                + this->myArrivalEvent->get_DeterministicDeltav();

            //***we do not compute the delta-v from the finite burn. in fact we *can't* compute it because we don't know what it is without integrating the entire thrust history
            //max_thrust is not populated until the call to output(). And even then it's not correct because it only tells you the thrust at the center of each control step
            //if you want to do this right, write out an ephemeris file at a tight (~1 day) resolution and apply the rocket equation between each line.
            //Since EMTG does not optimize on delta-v for the FBLT transcription and the computational cost of doing this "correctly" is high because
            //we would have to augment the state vector, I opted not to put this into the code for now.
            //for (size_t step = 0; step < this->num_timesteps; ++step)
            //{
            //    this->PhaseTotalDeterministicDeltav += this->throttle[step] * this->max_thrust[step] * this->PhaseDutyCycle / this->spacecraft_state_event_minus[step](6) * this->ThrustStepLengths[step];
            //}

            this->PhaseTotalStatisticalDeltav = this->myDepartureEvent->get_StatisticalDeltav()
                + this->myArrivalEvent->get_StatisticalDeltav()
                + this->initial_TCM_magnitude;

            //TODO: derivatives
            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV
                && needG)
            {
                //do stuff
            }
        }//end process_deltav_contribution

        void FBLTphase::output(std::ofstream& outputfile,
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
                this->InitialCoastPropagatorObject->setStateRight(this->output_state);

                //Step 3.2: propagate
                this->InitialCoastPropagatorObject->setCurrentEpoch(this->state_after_initial_TCM(7));
                this->InitialCoastPropagatorObject->setIndexOfEpochInStateVec(7);
                this->InitialCoastPropagatorObject->setCurrentIndependentVariable(this->state_after_initial_TCM(7));
                this->InitialCoastPropagatorObject->propagate(this->InitialCoastDuration / 2.0, false);
                this->output_state(6) = this->state_after_initial_TCM(6);
                this->output_state(7) = this->state_after_initial_TCM(7) + this->InitialCoastDuration / 2.0;

                //Step 3.3: assign the initial coast propagator back to where it is supposed to go
                this->InitialCoastPropagatorObject->setStateRight(this->spacecraft_state_event_minus[0]);

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

                    R_sc_Sun = this->output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                }
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                doubleType power = this->mySpacecraft->getAvailablePower();
                doubleType active_power = this->mySpacecraft->getEPActivePower();

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
                    this->output_state,//state
                    math::Matrix<doubleType>(3, 1, 0.0),//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    0.0,//dVmag
                    0.0,//Thrust
                    0.0,//Isp
                    power,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    active_power,
                    "none");//active_power)
            }

            //Step 4: output the maneuver events
            for (size_t step = 0; step < this->num_timesteps; ++step)
            {
                //Step 4.1: output the current step
                std::string event_type;

                //first we need to propagate to the step mid-point

                if (step < this->num_timesteps / 2) //forward timestep
                {
                    //Step 4.1.1: temporarily redirect the propagator to the output state
                    this->ForwardPropagatorObjects[step].setStateRight(this->output_state);

                    //Step 4.1.2: propagate
                    this->ForwardPropagatorObjects[step].setCurrentEpoch(this->spacecraft_state_event_minus[step](7));
                    this->ForwardPropagatorObjects[step].setIndexOfEpochInStateVec(7);
                    this->ForwardPropagatorObjects[step].setCurrentIndependentVariable(this->spacecraft_state_event_minus[step](7));
                    this->ForwardPropagatorObjects[step].propagate(this->ForwardPropagationStepTimes[step] / 2.0, this->ControlVector[step], false);

                    //Step 4.1.3: update the epoch
                    this->output_state(7) = this->spacecraft_state_event_minus[step](7) + this->ForwardPropagationStepTimes[step] / 2.0;
                    
                    //Step 4.1.4: redirect the propagator back to where it is supposed to go
                    this->ForwardPropagatorObjects[step].setStateRight(this->spacecraft_state_event_plus[step]);

                    //Step 4.1.5: figure out spacecrafty things

                    //Step 4.1.5.1: where am I relative to the sun?
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

                        R_sc_Sun = this->output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                    }

                    //Step 4.1.5.2: call the power model
                    doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                    this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                    //Step 4.1.5.3: call the thruster model
                    if (this->ControlVector[step].get_n() == 4)
                        this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle, this->ControlVector[step](3));
                    else
                        this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle);

                    //Step 4.1.5.4: store the thruster model outputs
                    this->max_thrust[step] = this->mySpacecraft->getEPthrust() * 1.0e-3;
                    this->max_mass_flow_rate[step] = this->mySpacecraft->getEPMassFlowRate();
                    this->Isp[step] = this->mySpacecraft->getEPIsp();
                    this->power[step] = this->mySpacecraft->getAvailablePower();
                    this->active_power[step] = this->mySpacecraft->getEPActivePower();
                    this->number_of_active_engines[step] = this->mySpacecraft->getEPNumberOfActiveThrusters();
                    this->ThrottleLevel[step] = this->mySpacecraft->getEPThrottleLevel();
                    this->ThrottleLevelString[step] = this->mySpacecraft->getEPThrottleLevelString();
                }
                else //backward timestep
                {
                    size_t backstep = step;
                    size_t sstep = this->num_timesteps - 1 - step;
                    //Step 4.1.1: temporarily redirect the propagator to the output state
                    this->BackwardPropagatorObjects[sstep].setStateRight(this->output_state);

                    //Step 4.1.2: propagate
                    this->BackwardPropagatorObjects[sstep].setCurrentEpoch(this->spacecraft_state_event_plus[backstep](7));
                    this->BackwardPropagatorObjects[sstep].setIndexOfEpochInStateVec(7);
                    this->BackwardPropagatorObjects[sstep].setCurrentIndependentVariable(this->spacecraft_state_event_plus[backstep](7));
                    this->BackwardPropagatorObjects[sstep].propagate(-this->BackwardPropagationStepTimes[sstep] / 2.0, this->ControlVector[backstep], false);

                    //Step 4.1.3: update the epoch
                    this->output_state(7) = this->spacecraft_state_event_plus[backstep](7) - this->BackwardPropagationStepTimes[sstep] / 2.0;

                    //Step 4.1.4: redirect the propagator back to where it is supposed to go
                    this->BackwardPropagatorObjects[sstep].setStateRight(this->spacecraft_state_event_minus[backstep]);

                    //Step 4.1.5: figure out spacecrafty things

                    //Step 4.1.5.1: where am I relative to the sun?
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

                        R_sc_Sun = this->output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                    }

                    //Step 4.1.5.2: call the power model
                    doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                    this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                    //Step 4.1.5.3: call the thruster model
                    if (this->ControlVector[backstep].get_n() == 4)
                        this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle, this->ControlVector[backstep](3));
                    else
                        this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle);

                    //Step 4.1.5.4: store the thruster model outputs
                    this->max_thrust[backstep] = this->mySpacecraft->getEPthrust() * 1.0e-3;
                    this->max_mass_flow_rate[backstep] = this->mySpacecraft->getEPMassFlowRate();
                    this->Isp[backstep] = this->mySpacecraft->getEPIsp();
                    this->power[backstep] = this->mySpacecraft->getAvailablePower();
                    this->active_power[backstep] = this->mySpacecraft->getEPActivePower();
                    this->number_of_active_engines[backstep] = this->mySpacecraft->getEPNumberOfActiveThrusters();
                    this->ThrottleLevel[backstep] = this->mySpacecraft->getEPThrottleLevel();
                    this->ThrottleLevelString[backstep] = this->mySpacecraft->getEPThrottleLevelString();
                }

                if (this->max_thrust[step] > 5.0e-7 && this->throttle[step] > 5.0e-4)
                    event_type = "FBLTthrust";
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
                    this->output_state,//state
                    this->stepDeltaV[step],//dV
                    ThrustVector,//ThrustVector
                    this->throttle[step],//throttle
                    this->max_thrust[step] * 1000.0 / this->PhaseDutyCycle,//Thrust
                    this->Isp[step],//Isp
                    this->power[step],//AvailPower
                    this->max_mass_flow_rate[step] / this->PhaseDutyCycle,//mdot
                    this->number_of_active_engines[step],//number_of_active_engines
                    this->active_power[step],
                    this->ThrottleLevelString[step]);//active_power)

                //Step 4.2: if this is the middle step, output the match point
                if (step == this->num_timesteps / 2 - 1)
                {

                    math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                    if (this->myUniverse->central_body_SPICE_ID == 10)
                    {
                        R_sc_Sun = this->match_point_state_minus.getSubMatrix1D(0, 2);
                    }
                    else
                    {
                        //where is the central body relative to the sun?
                        doubleType central_body_state_and_derivatives[12];
                        this->myUniverse->locate_central_body(this->match_point_state_minus(7),
                            central_body_state_and_derivatives,
                            *this->myOptions,
                            false);

                        math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                        for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                        {
                            R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                        }

                        R_sc_Sun = this->match_point_state_minus.getSubMatrix1D(0, 2) + R_CB_Sun;
                    }
                    doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                    this->mySpacecraft->computePowerState(r_sc_sun_AU, this->match_point_state_minus(7));

                    doubleType power = this->mySpacecraft->getAvailablePower();
                    doubleType active_power = this->mySpacecraft->getEPActivePower();

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
                        "none");//active_power)
                }
            }

            //Step 5: output the state at the middle of the terminal coast, if applicable
            if (this->hasTerminalCoast)
            {
                //Step 5.1: create a temporary output state and assign the Terminal coast propagator to it
                this->TerminalCoastPropagatorObject->setStateRight(this->output_state);

                //Step 5.2: propagate
                this->TerminalCoastPropagatorObject->setCurrentEpoch(this->state_at_end_of_phase(7));
                this->TerminalCoastPropagatorObject->setIndexOfEpochInStateVec(7);
                this->TerminalCoastPropagatorObject->setCurrentIndependentVariable(this->state_at_end_of_phase(7));
                this->TerminalCoastPropagatorObject->propagate(-this->TerminalCoastDuration / 2.0, false);
                this->output_state(6) = this->state_at_end_of_phase(6);
                this->output_state(7) = this->state_at_end_of_phase(7) - this->TerminalCoastDuration / 2.0;

                //Step 5.3: assign the Terminal coast propagator back to where it is supposed to go
                this->TerminalCoastPropagatorObject->setStateRight(this->spacecraft_state_event_plus.back());

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

                    R_sc_Sun = this->output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                }
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                doubleType power = this->mySpacecraft->getAvailablePower();
                doubleType active_power = this->mySpacecraft->getEPActivePower();

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
                    this->output_state,//state
                    math::Matrix<doubleType>(3, 1, 0.0),//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    0.0,//dVmag
                    0.0,//Thrust
                    0.0,//Isp
                    power,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    active_power,
                    "none");//active_power)
            }

            //Step 6: output the arrival event
            this->myArrivalEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
        }

        void FBLTphase::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
        {
            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            if (this->myJourneyOptions->override_integration_step_size)
                this->EphemerisOutputResolution = this->myJourneyOptions->integration_step_size;
            else
                this->EphemerisOutputResolution = this->myOptions->integration_time_step_size;

            //Step 1: output the departure event
            this->temp_state = this->myDepartureEvent->get_state_before_event();
            this->myDepartureEvent->output_ephemeris(outputfile);

            if (this->myOptions->generate_acceleration_model_instrumentation_file)
            {
                this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
            }

            //Step 2: output the initial coast
            if (this->hasInitialCoast)
            {
                //Step 2.1: temporarily assign the initial coast propagator to the output state
                this->InitialCoastPropagatorObject->setStateRight(this->output_state);
                this->InitialCoastPropagatorObject->setCurrentEpoch(this->state_after_initial_TCM(7));
                this->InitialCoastPropagatorObject->setIndexOfEpochInStateVec(7);
                this->InitialCoastPropagatorObject->setCurrentIndependentVariable(this->state_after_initial_TCM(7));

                //Step 2.2: propagate and print, skipping the first and last entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->InitialCoastDuration;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 2.2.1: propagate
					this->InitialCoastPropagatorObject->setCurrentEpoch(this->state_after_initial_TCM(7));
					this->InitialCoastPropagatorObject->setIndexOfEpochInStateVec(7);
					this->InitialCoastPropagatorObject->setCurrentIndependentVariable(this->state_after_initial_TCM(7));
                    this->InitialCoastPropagatorObject->propagate(timeToPropagate, false);
                    this->output_state(7) = this->state_after_initial_TCM(7) + timeToPropagate;
                    this->temp_state = this->output_state;
                    //Step 2.2.2: convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, this->output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            this->output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //Step 2.2.3: print
                    if (this->myOptions->generate_acceleration_model_instrumentation_file)
                    {
                        this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
                    }

                    if ((totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution)
                        this->write_ephemeris_line(outputfile,
                            this->output_state,
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
                this->InitialCoastPropagatorObject->setStateRight(this->spacecraft_state_event_minus[0]);
            }

            //Step 3: output the thrust arcs
            //forward step
            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                //Step 3.1: temporarily assign the propagator to the output state
                this->ForwardPropagatorObjects[step].setStateRight(this->output_state);
                this->ForwardPropagatorObjects[step].setCurrentEpoch(this->spacecraft_state_event_minus[step](7));
                this->ForwardPropagatorObjects[step].setIndexOfEpochInStateVec(7);
                this->ForwardPropagatorObjects[step].setCurrentIndependentVariable(this->spacecraft_state_event_minus[step](7));

                //Step 3.2: propagate and print, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->ForwardPropagationStepTimes[step]_GETVALUE;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 3.2.1: propagate
                    this->ForwardPropagatorObjects[step].setIndexOfEpochInStateVec(7);
                    this->ForwardPropagatorObjects[step].setCurrentIndependentVariable(this->spacecraft_state_event_minus[step](7));
                    this->ForwardPropagatorObjects[step].setCurrentEpoch(this->spacecraft_state_event_minus[step](7));
                    this->ForwardPropagatorObjects[step].propagate(timeToPropagate, this->ControlVector[step], false);
                    this->output_state(7) = this->spacecraft_state_event_minus[step](7) + timeToPropagate;
                    this->temp_state = this->output_state;
                    //Step 3.2.2: convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, this->output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            this->output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //Step 3.2.3: print
                    if ((totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution)
                    {
                        //we need an instantaneous power/propulsion state
                        //Step 3.2.3.1: where am I relative to the sun?
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

                            R_sc_Sun = this->output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                        }

                        //Step 3.2.3.2: call the power model
                        doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                        this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                        //Step 3.2.3.3: call the thruster model
                        if (this->ControlVector[step].get_n() == 4)
                            this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle, this->ControlVector[step](3));
                        else
                            this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle);

                        //Step 3.2.3.4: print
                        if (this->myOptions->generate_acceleration_model_instrumentation_file)
                        {
                            this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7), this->ControlVector[step]);
                        }

                        this->write_ephemeris_line(outputfile,
                            this->output_state,
                            this->ControlVector[step],
                            this->mySpacecraft->getEPthrust() * 1.0e-3 * this->throttle[step],
                            this->mySpacecraft->getEPMassFlowRate() * this->throttle[step],
                            this->mySpacecraft->getEPIsp(),
                            this->mySpacecraft->getEPNumberOfActiveThrusters(),
                            this->mySpacecraft->getEPActivePower(),
                            this->mySpacecraft->getEPThrottleLevelString());
                    }

                    //Step 3.2.4: increment propagatedEpoch
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution
                        ? this->EphemerisOutputResolution
                        : (totalPropagationTime - timeToPropagate);
                }

                //Step 3.3: reset the propagator to its original state vectors
                if (step == this->num_timesteps / 2 - 1)
                    this->ForwardPropagatorObjects[step].setStateRight(this->match_point_state_minus);
                else
                    this->ForwardPropagatorObjects[step].setStateRight(this->spacecraft_state_event_plus[step]);
            }//end forward step SPK output

            //backward step
            for (int step = this->num_timesteps / 2 - 1; step >= 0; --step)
            {
                size_t backstep = this->num_timesteps - 1 - step;

                //Step 3.1: temporarily assign the propagator to the output state
                if (backstep == this->num_timesteps / 2)
                    this->BackwardPropagatorObjects[step].setStateLeft(this->match_point_state_plus);
                else
                    this->BackwardPropagatorObjects[step].setStateLeft(this->spacecraft_state_event_minus[backstep]);

                this->BackwardPropagatorObjects[step].setStateRight(this->output_state);

                if (step == this->num_timesteps / 2 - 1)
                {
                    this->BackwardPropagatorObjects[step].setCurrentEpoch(this->match_point_state_plus(7));
                    this->BackwardPropagatorObjects[step].setIndexOfEpochInStateVec(7);
                    this->BackwardPropagatorObjects[step].setCurrentIndependentVariable(this->match_point_state_plus(7));
                }
                else
                {
                    this->BackwardPropagatorObjects[step].setCurrentEpoch(this->spacecraft_state_event_minus[backstep](7));
                    this->BackwardPropagatorObjects[step].setIndexOfEpochInStateVec(7);
                    this->BackwardPropagatorObjects[step].setCurrentIndependentVariable(this->spacecraft_state_event_minus[backstep](7));
                }

                //Step 3.2: propagate and print, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->BackwardPropagationStepTimes[step]_GETVALUE;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 3.2.1: propagate
                    if (step == this->num_timesteps / 2 - 1)
                    {
                        this->BackwardPropagatorObjects[step].setCurrentEpoch(this->match_point_state_plus(7));
                        this->BackwardPropagatorObjects[step].setIndexOfEpochInStateVec(7);
                        this->BackwardPropagatorObjects[step].setCurrentIndependentVariable(this->match_point_state_plus(7));
                    }
                    else
                    {
                        this->BackwardPropagatorObjects[step].setCurrentEpoch(this->spacecraft_state_event_minus[backstep](7));
                        this->BackwardPropagatorObjects[step].setIndexOfEpochInStateVec(7);
                        this->BackwardPropagatorObjects[step].setCurrentIndependentVariable(this->spacecraft_state_event_minus[backstep](7));
                    }
                    this->BackwardPropagatorObjects[step].propagate(timeToPropagate, this->ControlVector[backstep], false);
                    if (step == this->num_timesteps / 2 - 1)
                    {
                        this->output_state(7) = this->match_point_state_plus(7) + timeToPropagate;
                    }
                    else
                    {
                        this->output_state(7) = this->spacecraft_state_event_minus[backstep](7) + timeToPropagate;
                    }
                    this->temp_state = this->output_state;
                    //Step 3.2.2: convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, this->output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            this->output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //Step 3.2.3: print
                    if ((totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution)
                    {
                        //we need an instantaneous power/propulsion state
                        //Step 3.2.3.1: where am I relative to the sun?
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

                            R_sc_Sun = this->output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                        }

                        //Step 3.2.3.2: call the power model
                        doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                        this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                        //Step 3.2.3.3: call the thruster model
                        if (this->ControlVector[step].get_n() == 4)
                            this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle, this->ControlVector[backstep](3));
                        else
                            this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle);

                        //Step 3.2.3.4: print
                        if (this->myOptions->generate_acceleration_model_instrumentation_file)
                        {
                            this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7), this->ControlVector[backstep]);
                        }

                        this->write_ephemeris_line(outputfile,
                            this->output_state,
                            this->ControlVector[backstep],
                            this->mySpacecraft->getEPthrust() * 1.0e-3 * this->throttle[backstep],
                            this->mySpacecraft->getEPMassFlowRate() * this->throttle[backstep],
                            this->mySpacecraft->getEPIsp(),
                            this->mySpacecraft->getEPNumberOfActiveThrusters(),
                            this->mySpacecraft->getEPActivePower(),
                            this->mySpacecraft->getEPThrottleLevelString());
                    }

                    //Step 3.2.4: increment propagatedEpoch
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution
                        ? this->EphemerisOutputResolution
                        : (totalPropagationTime - timeToPropagate);
                }

                //Step 3.3: reset the propagator to its original state vectors
                this->BackwardPropagatorObjects[step].setStateLeft(this->spacecraft_state_event_plus[backstep]);
                if (step == this->num_timesteps / 2 - 1)
                    this->BackwardPropagatorObjects[step].setStateRight(this->match_point_state_plus);
                else
                    this->BackwardPropagatorObjects[step].setStateRight(this->spacecraft_state_event_minus[backstep]);
            }//end backward step SPK output

             //Step 4: output the terminal coast
            if (this->hasTerminalCoast)
            {
                //Step 4.1: temporarily assign the terminal coast propagator to the output state
                this->TerminalCoastPropagatorObject->setStateLeft(this->spacecraft_state_event_plus.back());
                this->TerminalCoastPropagatorObject->setStateRight(this->output_state);
                this->TerminalCoastPropagatorObject->setCurrentEpoch(this->state_at_end_of_phase(7));
                this->TerminalCoastPropagatorObject->setIndexOfEpochInStateVec(7);
                this->TerminalCoastPropagatorObject->setCurrentIndependentVariable(this->state_at_end_of_phase(7));

                //Step 4.2: propagate and print, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->TerminalCoastDuration;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 4.2.1: propagate
					this->TerminalCoastPropagatorObject->setCurrentEpoch(this->state_at_end_of_phase(7));
					this->TerminalCoastPropagatorObject->setIndexOfEpochInStateVec(7);
					this->TerminalCoastPropagatorObject->setCurrentIndependentVariable(this->state_at_end_of_phase(7));
                    this->TerminalCoastPropagatorObject->propagate(timeToPropagate, false);
                    this->output_state(7) = this->spacecraft_state_event_plus.back()(7) + timeToPropagate;
                    this->temp_state = this->output_state;
                    //Step 4.2.2: convert to Sun-centered if necessary
                    if (!(boost::to_lower_copy(this->myUniverse->central_body_name) == "sun"))
                    {
                        doubleType body_state[12];

                        this->myUniverse->locate_central_body(this->output_state(7),
                            body_state,
                            *this->myOptions,
                            false);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            this->output_state(stateIndex) += body_state[stateIndex];
                    }

                    //Step 4.2.3: print
                    if (this->myOptions->generate_acceleration_model_instrumentation_file)
                    {
                        this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
                    }

                    if ((totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution)
                        this->write_ephemeris_line(outputfile,
                            this->output_state,
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
                this->TerminalCoastPropagatorObject->setStateLeft(this->state_at_end_of_phase);
                this->TerminalCoastPropagatorObject->setStateRight(this->spacecraft_state_event_plus.back());
            }

            //Step 5: output the arrival event
            this->myArrivalEvent->output_ephemeris(outputfile);
            this->temp_state = this->myArrivalEvent->get_state_before_event();
            if (this->myOptions->generate_acceleration_model_instrumentation_file)
            {
                this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
            }

        }//end output_ephemeris()

        void FBLTphase::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //should write out a maneuver/target spec for the departure maneuver if there is one
            this->myDepartureEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);

            // The DepartureEvent may have a maneuver, in which case, we may need a target
            if (haveManeuverNeedTarget)
            {
                // initialize target spec object
                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->spacecraft_state_event_minus[0](7),
                    this->spacecraft_state_event_minus[0]); // this is the state immediately before the first thrust segment
            
                // write target spec object
                myTargetSpecLine.write(target_spec_file);
            
                // disable the target flag
                haveManeuverNeedTarget = false;
            }

            // Loop through the forward control segments
            for (size_t segment = 0; segment < this->num_timesteps / 2; ++segment)
            {
                std::string segment_name = this->name + "_Segment" + std::to_string(segment);

                // Initialize a maneuver spec line object
                maneuver_spec_line myManeuverSpecLine(segment_name);

                // Step through the FBLT step at IntegrationStep intervals
                // Every time the throttle level changes save off a thrust step
                size_t ManeuverThrottleLevel;
                doubleType ManeuverStartEpoch;
                doubleType ManeuverStartMass;
                doubleType ManeuverThrustMagnitude;
                doubleType ManeuverMassFlowRate;
                // Compute the propulsion characteristics on the left side of the segment
                {
                    // Locate the spacecraft relative to the sun
                    for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
                    {
                        this->output_state(stateIndex) = this->spacecraft_state_event_minus[segment](stateIndex);
                    }

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

                        R_sc_Sun = this->output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                    }

                    //Step 2.2.1.2: call the power model
                    doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                    this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                    //Step 2.2.1.3: call the thruster model using 100% duty cycle because we want the raw performance information
                    if (this->num_controls == 4)
                    {
                        this->mySpacecraft->computeElectricPropulsionPerformance(1.0, this->ControlVector[segment](3));
                    }
                    else
                    {
                        this->mySpacecraft->computeElectricPropulsionPerformance(1.0);
                    }

                    //Step 2.2.1.3: populate fields
                    ManeuverThrottleLevel = this->mySpacecraft->getEPThrottleLevel();
                    ManeuverStartEpoch = this->output_state(7);
                    ManeuverStartMass = this->output_state(6);
                    ManeuverThrustMagnitude = this->mySpacecraft->getEPthrust(); //in N
                    ManeuverMassFlowRate = this->mySpacecraft->getEPMassFlowRate();

                } // end segment left side propulsion characteristics 

                //Step 2.2.2.1: hijack the propagator
                this->ForwardPropagatorObjects[segment].setStateRight(this->output_state);

                //Step 2.2.2.2: propagate the control segment, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->ForwardPropagationStepTimes[segment] _GETVALUE;
                double divisor = 1.0;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 2.2.2.3: propagate
                    if (segment == 0)
                    {
                        this->Forward_dStatedIndependentVariable[segment].assign_zeros();
                        this->ForwardPropagatorObjects[segment].setCurrentEpoch(this->spacecraft_state_event_minus[0](7));
                        this->ForwardPropagatorObjects[segment].setCurrentIndependentVariable(this->spacecraft_state_event_minus[0](7));
                    }
                    else
                    {
                        this->Forward_dStatedIndependentVariable[segment] = this->Forward_dStatedIndependentVariable[segment - 1];
                        this->ForwardPropagatorObjects[segment].setCurrentEpoch(this->spacecraft_state_event_plus[segment - 1](7));
                        this->ForwardPropagatorObjects[segment].setCurrentIndependentVariable(this->spacecraft_state_event_plus[segment - 1](7));
                    }
                    this->ForwardPropagatorObjects[segment].setIndexOfEpochInStateVec(7);
                    this->ForwardPropagatorObjects[segment].propagate(timeToPropagate, this->ControlVector[segment], false);

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

                    //Step 2.2.2.6: call the power model
                    doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                    this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                    //Step 2.2.2.7: call the thruster model
                    if (this->num_controls == 4)
                        this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle, this->ControlVector[segment](3));
                    else
                        this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle);

                    //Step 2.2.2.8: did the thrust change? if so save off a thrust arc and reset

                    size_t CurrentThrottleLevel = this->mySpacecraft->getEPThrottleLevel();
                    doubleType CurrentThrustMagnitude = this->mySpacecraft->getEPthrust(); //in N
                    doubleType CurrentMassFlowRate = this->mySpacecraft->getEPMassFlowRate();

                    if (CurrentThrottleLevel != ManeuverThrottleLevel
                        || fabs(CurrentThrustMagnitude - ManeuverThrustMagnitude) > 1.0e-3) //if different throttle level or thrust different by more than one mN
                    {
                        //save off a maneuver spec item
                        if (divisor < this->EphemerisOutputResolution)
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
                                this->ControlVector[segment],
                                ManeuverStartMass,
                                output_state(6),
                                ManeuverThrustMagnitude,
                                ManeuverMassFlowRate,
                                output_state(7) - ManeuverStartEpoch,
                                this->PhaseDutyCycle);

                            //reset for next maneuver item
                            ManeuverThrottleLevel = CurrentThrottleLevel;
                            ManeuverStartEpoch = output_state(7);
                            ManeuverStartMass = output_state(6);
                            ManeuverThrustMagnitude = CurrentThrustMagnitude;
                            ManeuverMassFlowRate = CurrentMassFlowRate;
                        }
                    }

                    //Step 2.2.2.9: increment propagatedEpoch - we deliberately do NOT take the last partial step
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution / divisor
                        ? this->EphemerisOutputResolution / divisor
                        : (totalPropagationTime - timeToPropagate);

                } // end loop through current forward control segment

                // reset the propagator right-hand state container
                if (segment == this->num_timesteps / 2 - 1)
                {
                    this->ForwardPropagatorObjects[segment].setStateRight(this->match_point_state_minus);
                }
                else
                {
                    this->ForwardPropagatorObjects[segment].setStateRight(this->spacecraft_state_event_plus[segment]);
                }

                //Step 2.3: save off the last thrust step
                myManeuverSpecLine.append_maneuver_spec_item("EME2000",
                    ManeuverStartEpoch,
                    this->ControlVector[segment],
                    ManeuverStartMass,
                    this->spacecraft_state_event_plus[segment](6),
                    ManeuverThrustMagnitude,
                    ManeuverMassFlowRate,
                    this->spacecraft_state_event_plus[segment](7) - ManeuverStartEpoch,
                    this->PhaseDutyCycle);

                // write the maneuver spec line
                myManeuverSpecLine.write(maneuver_spec_file);

                target_spec_line segmentEndTargetLine(segment_name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->spacecraft_state_event_plus[segment](7),
                    this->spacecraft_state_event_plus[segment]); // this is the state at the end of the current thrust segment

                // write target spec object
                segmentEndTargetLine.write(target_spec_file);

            } // end loop through forward control segments

            // Loop over backward control segments
            // Start at the control segment nearest the match point, but integrate forward in time
            // Integrate from spacecraft_state_event_minus[n/2 - 1] to spacecraft_state_event_plus[n/2 - 1]
            for (int segment = this->num_timesteps / 2 - 1; segment >= 0; --segment)
            {

                size_t backstep = this->num_timesteps - 1 - segment;

                // Initialize a maneuver spec line object
                std::string segment_name = this->name + "_Segment" + std::to_string(backstep);
                maneuver_spec_line myManeuverSpecLine(segment_name);

                // Step through the FBLT segment at IntegrationStep intervals
                // Every time the throttle level changes save off a new maneuver
                size_t ManeuverThrottleLevel;
                doubleType ManeuverStartEpoch;
                doubleType ManeuverStartMass;
                doubleType ManeuverThrustMagnitude;
                doubleType ManeuverMassFlowRate;
                // Compute the propulsion characteristics on the left side of the segment
                {
                    // Locate the spacecraft relative to the sun
                    for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
                    {
                        if (backstep == this->num_timesteps / 2)
                        {
                            this->output_state(stateIndex) = this->match_point_state_plus(stateIndex);
                        }
                        else
                        {
                            this->output_state(stateIndex) = this->spacecraft_state_event_minus[backstep](stateIndex);
                        }
                    }

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

                        R_sc_Sun = this->output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                    }

                    //Step 2.2.1.2: call the power model
                    doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                    this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                    //Step 2.2.1.3: call the thruster model using 100% duty cycle because we want the raw performance information
                    if (this->num_controls == 4)
                    {
                        this->mySpacecraft->computeElectricPropulsionPerformance(1.0, this->ControlVector[backstep](3));
                    }
                    else
                    {
                        this->mySpacecraft->computeElectricPropulsionPerformance(1.0);
                    }

                    //Step 2.2.1.3: populate fields
                    ManeuverThrottleLevel = this->mySpacecraft->getEPThrottleLevel();
                    ManeuverStartEpoch = this->output_state(7);
                    ManeuverStartMass = this->output_state(6);
                    ManeuverThrustMagnitude = this->mySpacecraft->getEPthrust(); //in N
                    ManeuverMassFlowRate = this->mySpacecraft->getEPMassFlowRate();

                } // end segment left side propulsion characteristics 

                //Step 2.2.2.1: hijack the propagator
                //Step 3.1: temporarily assign the propagator to the output state
                if (backstep == this->num_timesteps / 2)
                {
                    this->BackwardPropagatorObjects[segment].setStateLeft(this->match_point_state_plus);
                }
                else
                {
                    this->BackwardPropagatorObjects[segment].setStateLeft(this->spacecraft_state_event_minus[backstep]);
                }
                this->BackwardPropagatorObjects[segment].setStateRight(this->output_state);

                //Step 2.2.2.2: propagate the control segment, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->BackwardPropagationStepTimes[segment] _GETVALUE;

                double divisor = 1.0;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 2.2.2.3: propagate
                    if (segment == this->num_timesteps / 2 - 1)
                    {
                        this->BackwardPropagatorObjects[segment].setCurrentEpoch(this->match_point_state_plus(7));
                        this->BackwardPropagatorObjects[segment].setCurrentIndependentVariable(this->match_point_state_plus(7));
                    }
                    else
                    {
                        this->BackwardPropagatorObjects[segment].setCurrentEpoch(this->spacecraft_state_event_minus[backstep](7));
                        this->BackwardPropagatorObjects[segment].setCurrentIndependentVariable(this->spacecraft_state_event_minus[backstep](7));
                    }
                    this->BackwardPropagatorObjects[segment].setIndexOfEpochInStateVec(7);

                    // propagate up to the current location in the thrust segment
                    this->BackwardPropagatorObjects[segment].propagate(timeToPropagate, this->ControlVector[backstep], false);

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

                    //Step 2.2.2.6: call the power model
                    doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                    this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                    //Step 2.2.2.7: call the thruster model using 100% duty cycle because we want the raw performance information
                    if (this->num_controls == 4)
                        this->mySpacecraft->computeElectricPropulsionPerformance(1.0, this->ControlVector[backstep](3));
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
                        if (divisor < this->EphemerisOutputResolution)
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
                                this->ControlVector[backstep],
                                ManeuverStartMass,
                                output_state(6),
                                ManeuverThrustMagnitude,
                                ManeuverMassFlowRate,
                                output_state(7) - ManeuverStartEpoch,
                                this->PhaseDutyCycle);

                            //reset for next maneuver item
                            ManeuverThrottleLevel = CurrentThrottleLevel;
                            ManeuverStartEpoch = output_state(7);
                            ManeuverStartMass = output_state(6);
                            ManeuverThrustMagnitude = CurrentThrustMagnitude;
                            ManeuverMassFlowRate = CurrentMassFlowRate;
                        }
                                                                                        
                    }

                    //Step 2.2.2.9: increment propagatedEpoch - we deliberately do NOT take the last partial step
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution / divisor
                        ? this->EphemerisOutputResolution / divisor
                        : (totalPropagationTime - timeToPropagate);
                    

                } // end loop through current backward control segment

                // reset the propagator left and right-hand state containers       
                this->BackwardPropagatorObjects[segment].setStateLeft(this->spacecraft_state_event_plus[backstep]);
                if (segment == this->num_timesteps / 2 - 1)
                {
                    this->BackwardPropagatorObjects[segment].setStateRight(this->match_point_state_plus);
                }
                else
                {
                    this->BackwardPropagatorObjects[segment].setStateRight(this->spacecraft_state_event_minus[backstep]);
                }

                //Step 2.3: save off the last thrust step
                myManeuverSpecLine.append_maneuver_spec_item("EME2000",
                    ManeuverStartEpoch,
                    this->ControlVector[backstep],
                    ManeuverStartMass,
                    this->spacecraft_state_event_plus[backstep](6),
                    ManeuverThrustMagnitude,
                    ManeuverMassFlowRate,
                    this->spacecraft_state_event_plus[backstep](7) - ManeuverStartEpoch,
                    this->PhaseDutyCycle);

                // write the maneuver spec line
                myManeuverSpecLine.write(maneuver_spec_file);

                target_spec_line segmentEndTargetLine(segment_name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->spacecraft_state_event_plus[backstep](7),
                    this->spacecraft_state_event_plus[backstep]); // this is the state at the end of the current thrust segment

                // write target spec object
                // if we are the last segment before the right hand phase boundary, then we will skip writing a target
                // and just use the right boundary event target to avoid possible back-to-back targets
                if (segment != 0)
                {
                    segmentEndTargetLine.write(target_spec_file);
                }

            } // end loop through backward control segments


            //set the target flag
            haveManeuverNeedTarget = true;

            //should write out a maneuver/target spec for the arrival maneuver if there is one
            this->myArrivalEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);



        }//end output_maneuver_and_target_spec()

    }//end namespace phases
}//end namespace EMTG