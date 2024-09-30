// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2024 United States Government as represented by the
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

//EMTGv9 ControlLawThrustPhase
//Jacob Englander 9-5-2020

#include "ControlLawThrustPhase.h"

#include "PropagatorFactory.h"
#include "IntegrationSchemeFactory.h"


namespace EMTG
{
    namespace Phases
    {
        ControlLawThrustPhase::ControlLawThrustPhase(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            phase* previousPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            HardwareModels::LaunchVehicle* myLaunchVehicle,
            missionoptions* myOptions) :
            ControlLawThrustPhase::CoastPhase(name,
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
            this->initialize();
        }//end constructor

        void ControlLawThrustPhase::initialize()
        {
            this->num_timesteps = this->myJourneyOptions->override_num_steps
                ? this->myJourneyOptions->number_of_steps
                : this->myOptions->num_timesteps;

            this->stateVectorNames.push_back("virtual chemical fuel");
            this->stateVectorNames.push_back("virtual electric propellant");

            //name our match point constraints
            this->matchPointConstraintNames.push_back("virtual chemical fuel");
            this->matchPointConstraintStateIndex.push_back(8); //note we skipped the encode of boundary epoch, it does not get a match point constraint
            this->continuity_constraint_scale_factors(7) = 1.0 / (1.0 * this->myJourneyOptions->maximum_mass);
            this->matchPointConstraintNames.push_back("virtual electric propellant");
            this->matchPointConstraintStateIndex.push_back(9); //note we skipped the encode of boundary epoch, it does not get a match point constraint
            this->continuity_constraint_scale_factors(8) = 1.0 / (1.0 * this->myJourneyOptions->maximum_mass);
            this->stateIndex_phase_propagation_variable = this->numStatesToPropagate;
            this->stateVectorNames.push_back("phase flight time");

            //electric maneuver - eventually we need to allow other types of maneuver
            this->hasBipropManeuver = false;
            this->hasMonopropManeuver = false;
            this->hasElectricManeuver = true;

            //size the TCM transition matrix
            this->TCMTM = math::Matrix<double>(this->numStatesToPropagate + 1, math::identity);

            //configure propagators
            this->MatchPointFraction = this->myJourneyOptions->CoastPhaseMatchPointFraction;
            this->dForwardStepSize_dPropagationVariable = this->MatchPointFraction;
            this->dBackwardStepSize_dPropagationVariable = 1.0 - this->MatchPointFraction;
            this->ForwardIntegrationStepLength = this->myJourneyOptions->CoastPhaseForwardIntegrationStepLength;
            this->BackwardIntegrationStepLength = this->myJourneyOptions->CoastPhaseBackwardIntegrationStepLength;

            //assume integrated propagator
            {
                this->isKeplerian = false;
                this->ForwardSTM = math::Matrix<double>(11, math::identity);
                this->BackwardSTM = math::Matrix<double>(11, math::identity);
                this->Forward_dPropagatedStatedIndependentVariable.resize(this->numStatesToPropagate + 1, 2, 0.0);
                this->Backward_dPropagatedStatedIndependentVariable.resize(this->numStatesToPropagate + 1, 2, 0.0);

                //acceleration model object
                this->mySpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(this->myOptions,
                    this->myJourneyOptions,
                    this->myUniverse,
                    this->Xdescriptions,
                    this->mySpacecraft,
                    11); // STM size
                this->mySpacecraftAccelerationModel->setDutyCycle(this->PhaseDutyCycle);
                this->mySpacecraftAccelerationModel->setThrustControlLaw(this->myJourneyOptions->thrust_control_law);

                //EOM
                this->myEOM.setSpacecraftAccelerationModel(this->mySpacecraftAccelerationModel);

                //integration scheme
                this->myIntegrationScheme = CreateIntegrationScheme(&this->myEOM, this->numStatesToPropagate, 11);

                this->ForwardHalfPhasePropagator = CreatePropagator(this->myOptions,
                    this->myUniverse,
                    this->numStatesToPropagate,
                    11,
                    this->state_after_initial_TCM,
                    this->match_point_state_minus,
                    this->ForwardSTM,
                    this->Forward_dPropagatedStatedIndependentVariable,
                    (Integration::Integrand*) &this->myEOM,
                    this->myIntegrationScheme,
                    &this->dForwardStepSize_dPropagationVariable,
                    this->ForwardIntegrationStepLength);

                this->BackwardHalfPhasePropagator = CreatePropagator(this->myOptions,
                    this->myUniverse,                    
                    this->numStatesToPropagate,
                    11,
                    this->state_at_end_of_phase,
                    this->match_point_state_plus,
                    this->BackwardSTM,
                    this->Backward_dPropagatedStatedIndependentVariable,
                    (Integration::Integrand*) &this->myEOM,
                    this->myIntegrationScheme,
                    &this->dBackwardStepSize_dPropagationVariable,
                    this->BackwardIntegrationStepLength);
            }

            //SPTMs
            this->ForwardSPTM = math::Matrix<double>(this->numStatesToPropagate + 1, math::identity);
            this->BackwardSPTM = math::Matrix<double>(this->numStatesToPropagate + 1, math::identity);

            //derivative truth table

            //unless SRP is enabled, only mass affects mass
            if (!this->myOptions->perturb_SRP)
            {
                //only mass variables affect mass constraints
                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                {
                    this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[6][stateIndex] = false;
                    this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[6][stateIndex] = false;
                }
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[6][7] = false; //epoch
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[6][7] = false; //epoch

                //mass variables do not affect other constraints if the phase is Keplerian (which it can't be)
                if (this->isKeplerian)
                {
                    for (size_t constraintIndex = 0; constraintIndex < 6; ++constraintIndex)
                    {
                        this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][6] = false; //mass variables
                        this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][6] = false; //mass variables
                    }
                }
            }

            //always true - propellant tanks are only affected by themselves and mass
            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
            {
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[7][stateIndex] = false; //fuel mass wrt state
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[7][stateIndex] = false; //fuel mass wrt state
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[8][stateIndex] = false; //electric propellant/oxidizer mass wrt state
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[8][stateIndex] = false; //electric propellant/oxidizer mass wrt state
            }
            this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[7][7] = false; //fuel mass wrt epoch
            this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[7][9] = false; //fuel mass wrt electric propellant/oxidizer mass
            this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[7][7] = false; //fuel mass wrt epoch
            this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[7][9] = false; //fuel mass wrt electric propellant/oxidizer mass
            //always true - propellant tank states do not affect anything but themselves
            for (size_t constraintIndex = 0; constraintIndex < 7; ++constraintIndex)
            {
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][8] = false; //other states wrt fuel mass
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][8] = false; //other states wrt fuel mass
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][9] = false; //other states wrt electric propellant/oxidizer mass
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][9] = false; //other states wrt electric propellant/oxidizer mass
            }

            //output stuff
            this->output_state.resize(this->numStatesToPropagate + 11 * 11, 1, 0.0);
        }//end initialize

        void ControlLawThrustPhase::calcbounds_match_point_constraints()
        {
            //base class
            this->TwoPointShootingPhase::calcbounds_match_point_constraints();

            //virtual chemical fuel
            this->create_sparsity_entry(this->Findices_match_point_constraints[7],
                this->X_index_of_first_decision_variable_in_this_phase,
                true,
                this->prefix + "virtual chemical fuel",
                this->Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel);

            //virtual electric propellant
            this->create_sparsity_entry(this->Findices_match_point_constraints[8],
                this->X_index_of_first_decision_variable_in_this_phase,
                true,
                this->prefix + "virtual electric propellant",
                this->Gindex_dVirtualElectricPropellantConstraint_dVirtualElectricPropellant);
        }//end calcbounds_match_point_constraints()
        
        void ControlLawThrustPhase::calcbounds_virtual_propellant_tanks()
        {
            //encode virtual propellant tank variables

            //Step 1: virtual chemical fuel
            {
                this->Xlowerbounds->push_back(0.0);
                this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
                this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
                this->Xdescriptions->push_back(prefix + "virtual chemical fuel");
                size_t Xindex = this->Xdescriptions->size() - 1;

                //tell the spacecraft about this variable
                //global constraint
                this->mySpacecraft->appendGlobalChemicalFuelTank_Xindices(Xindex);
                this->mySpacecraft->appendGlobalChemicalFuelTank_Xscaleranges(this->Xupperbounds->back(), this->Xlowerbounds->back());
                //stage constraint
                this->mySpacecraft->appendChemicalFuelTank_Xindices(this->stageIndex, Xindex);
                this->mySpacecraft->appendChemicalFuelTank_Xscaleranges(this->stageIndex, this->Xupperbounds->back(), this->Xlowerbounds->back());
            }

            //Step 2: virtual electric propellant
            {
                this->Xlowerbounds->push_back(0.0);
                this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
                this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
                this->Xdescriptions->push_back(prefix + "virtual electric propellant");
                size_t Xindex = this->Xdescriptions->size() - 1;

                //tell the spacecraft about this variable
                //global constraint
                this->mySpacecraft->appendGlobalElectricPropellantTank_Xindices(Xindex);
                this->mySpacecraft->appendGlobalElectricPropellantTank_Xscaleranges(this->Xupperbounds->back(), this->Xlowerbounds->back());
                //stage constraint
                this->mySpacecraft->appendElectricPropellantTank_Xindices(this->stageIndex, Xindex);
                this->mySpacecraft->appendElectricPropellantTank_Xscaleranges(this->stageIndex, this->Xupperbounds->back(), this->Xlowerbounds->back());
            }//end virtual electric propellant tank
        }//end calcbounds_virtual_propellant_tanks()

        //************************************process methods
        void ControlLawThrustPhase::process_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //do we want to propagate STMs or not?
            this->total_number_of_states_to_integrate = needG
                ? this->numStatesToPropagate + 11 * 11
                : this->numStatesToPropagate;

            this->process_phase_left_boundary(X, Xindex, F, Findex, G, needG);

            this->process_phase_flight_time(X, Xindex, F, Findex, G, needG);

            this->process_phase_right_boundary(X, Xindex, F, Findex, G, needG);

            this->process_phase_main(X, Xindex, F, Findex, G, needG);
        }//end process_phase()

        void ControlLawThrustPhase::process_forward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: set the tank states
            this->state_after_initial_TCM(8) = this->chemical_fuel_used; //TCM propellant, usually zero
            this->state_after_initial_TCM(9) = 0.0; //electric propellant

            //Step 2: call the forward propagator
            this->Forward_dPropagatedStatedIndependentVariable.assign_zeros();
            this->ForwardHalfPhasePropagator->setCurrentEpoch(this->state_after_initial_TCM(7));
            this->ForwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
            this->ForwardHalfPhasePropagator->setCurrentIndependentVariable(this->state_after_initial_TCM(7));

            //propagate with an empty control vector, just to get the acceleration model to activate the thrust term. It will be VNB
            this->ForwardHalfPhasePropagator->propagate(this->PhaseFlightTime * this->MatchPointFraction, needG);

            //Step 3: form the SPTM
            if (needG)
            {
                //Step 3.1: upper left 6x6 is the regular STM
                size_t stateMax = 11;
                for (size_t i = 0; i < stateMax; ++i)
                    for (size_t j = 0; j < stateMax; ++j)
                        this->ForwardSPTM(i, j) = this->ForwardSTM(i, j);

                //Step 3.2: turn the upper right 10 x 2 into the Phi_t terms developed by Lantoine
                for (size_t i = 0; i < stateMax; ++i)
                {
                    this->ForwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->ForwardSTM(i, 10);
                }
            }
        }//end process_forward_half_phase()

        void ControlLawThrustPhase::process_backward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: set the tank states
            this->state_at_end_of_phase(8) = this->virtual_chemical_fuel_used;
            this->state_at_end_of_phase(9) = this->virtual_electric_propellant_used;

            //Step 2: call the backward propagator
            this->Backward_dPropagatedStatedIndependentVariable.assign_zeros();
            this->BackwardHalfPhasePropagator->setCurrentEpoch(this->state_at_end_of_phase(7));
            this->BackwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
            this->BackwardHalfPhasePropagator->setCurrentIndependentVariable(this->state_at_end_of_phase(7));

            //propagate with an empty control vector, just to get the acceleration model to activate the thrust term. It will be VNB
            this->BackwardHalfPhasePropagator->propagate(-this->PhaseFlightTime * (1.0 - this->MatchPointFraction), needG);

            //Step 3: form the SPTM                                                                                     
            if (needG)
            {
                //Step 3.1: upper left 6x6 is the regular STM
                size_t stateMax = 11;
                for (size_t i = 0; i < stateMax; ++i)
                    for (size_t j = 0; j < stateMax; ++j)
                        this->BackwardSPTM(i, j) = this->BackwardSTM(i, j);
                //Step 3.2: turn the upper right 6 x 2 into the Phi_t terms developed by Lantoine
                for (size_t i = 0; i < 6; ++i)
                {
                    this->BackwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->BackwardSTM(i, 10);
                }
            }
        }//end process_backward_half_phase()

        void ControlLawThrustPhase::process_match_point_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: contruct the SPTM chains and the forward and backward HPTM
            if (needG)
            {
                this->ForwardHPTM = this->ForwardSPTM * this->TCMTM;
                this->BackwardHPTM = this->BackwardSPTM;
            }

            //Step 2: call the base class
            TwoPointShootingPhase::process_match_point_constraints(X, Xindex, F, Findex, G, needG);


            if (needG)
            {
                math::Matrix<double> dStateNow_dDecisionVariable(this->numStatesToPropagate + 1, 1, 0.0);
                math::Matrix<double> dMatchPointState_dDecisionVariable(this->numStatesToPropagate + 1, 1, 0.0);

                //Step 3: Derivatives with respect to propagation time
                dStateNow_dDecisionVariable.assign_zeros();
                dStateNow_dDecisionVariable(this->stateIndex_phase_propagation_variable) = 1.0;


                //Forward
                dMatchPointState_dDecisionVariable = this->ForwardHPTM * dStateNow_dDecisionVariable;

                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {
                    if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                    {
                        size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                        size_t Gindex = this->G_indices_match_point_constraints_wrt_PhaseFlightTime[constraintIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += -this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }
                }

                //Backward
                dMatchPointState_dDecisionVariable = this->BackwardHPTM * dStateNow_dDecisionVariable;

                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {
                    if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                    {
                        size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                        size_t Gindex = this->G_indices_match_point_constraints_wrt_PhaseFlightTime[constraintIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }
                }

                //virtual chemical fuel
                {
                    size_t Gindex = this->Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel;
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->continuity_constraint_scale_factors(7);
                }

                //virtual electric propellant
                {
                    size_t Gindex = this->Gindex_dVirtualElectricPropellantConstraint_dVirtualElectricPropellant;
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->continuity_constraint_scale_factors(7);
                }
            }//end derivatives
        }//end process_match_point_constraints()

        void ControlLawThrustPhase::process_virtual_propellant_tanks(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->virtual_chemical_fuel_used = X[Xindex++];

            this->virtual_electric_propellant_used = X[Xindex++];
        }//end process_virtual_propellant_tanks()

        void ControlLawThrustPhase::output(std::ofstream& outputfile,
            size_t& eventcount)
        {
            doubleType ForwardOutputTimestep = this->PhaseFlightTime * this->MatchPointFraction / (this->num_timesteps / 2);
            doubleType BackwardOutputTimestep = this->PhaseFlightTime * (1.0 - this->MatchPointFraction) / (this->num_timesteps / 2);
            math::Matrix<doubleType> empty3(3, 1, 0.0);
            math::Matrix<doubleType> ThrustVector;

            //Step 1: output the departure event
            this->myDepartureEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
            this->mySpacecraft->setActiveStage(this->stageIndex);


            //Step 2: output the initial TCM if applicable
            this->output_initial_TCM(outputfile, eventcount);

            //Step 3: print the forward half phase
            //set the propagator temporarily to dump into the output state
            this->ForwardHalfPhasePropagator->setStateRight(output_state);

            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                //propagate to the halfway point of this step
                this->Forward_dPropagatedStatedIndependentVariable.assign_zeros();
                this->ForwardHalfPhasePropagator->setCurrentEpoch(this->state_after_initial_TCM(7));
                this->ForwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
                this->ForwardHalfPhasePropagator->setCurrentIndependentVariable(this->state_after_initial_TCM(7));
                this->ForwardHalfPhasePropagator->propagate(ForwardOutputTimestep * (step + 0.5), false);

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

                //call the power model
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                //call the propulsion model
                this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle);

                //store the thruster model outputs
                doubleType max_thrust = this->mySpacecraft->getEPthrust() * 1.0e-3;
                doubleType max_mass_flow_rate = this->mySpacecraft->getEPMassFlowRate();
                doubleType Isp = this->mySpacecraft->getEPIsp();
                doubleType power = this->mySpacecraft->getAvailablePower();
                doubleType active_power = this->mySpacecraft->getEPActivePower();
                size_t number_of_active_engines = this->mySpacecraft->getEPNumberOfActiveThrusters();
                size_t ThrottleLevel = this->mySpacecraft->getEPThrottleLevel();
                std::string ThrottleLevelString = this->mySpacecraft->getEPThrottleLevelString();

                math::Matrix<doubleType> velocityUnitVector = output_state.getSubMatrix1D(3, 5).unitize();
                
                if (this->myJourneyOptions->thrust_control_law == ThrustControlLaw::Velocity)
                {
                    ThrustVector = velocityUnitVector * max_thrust * 1000.0 / this->PhaseDutyCycle;
                }
                else if (this->myJourneyOptions->thrust_control_law == ThrustControlLaw::AntiVelocity)
                {
                    ThrustVector = -velocityUnitVector * max_thrust * 1000.0 / this->PhaseDutyCycle;
                }

                //print
                this->write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "FBLTthrust",//event_type
                    "deep-space",//event_location
                    ForwardOutputTimestep / 86400.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    this->output_state,//state
                    empty3,//dV
                    ThrustVector,//ThrustVector
                    1.0,//throttle
                    max_thrust * 1000.0 / this->PhaseDutyCycle,//Thrust
                    Isp,//Isp
                    power,//AvailPower
                    max_mass_flow_rate / this->PhaseDutyCycle,//mdot
                    number_of_active_engines,//number_of_active_engines
                    active_power,
                    ThrottleLevelString);//active_power)
            }

            //reset the propagator to go where it is supposed to
            this->ForwardHalfPhasePropagator->setStateRight(this->match_point_state_minus);

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

            //Step 4: print the match point
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
                power,//AvailPower
                0.0,//mdot
                0,//number_of_active_engines
                active_power,
                "none");//active_power

            //Step 5: print the backward half phase
            //set the propagator temporarily to dump into the output state
            this->BackwardHalfPhasePropagator->setStateRight(output_state);

            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                size_t backstep = this->num_timesteps / 2 - step - 1;

                //propagate to the halfway point of this step
                this->Backward_dPropagatedStatedIndependentVariable.assign_zeros();
                this->BackwardHalfPhasePropagator->setCurrentEpoch(this->state_at_end_of_phase(7));
                this->BackwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
                this->BackwardHalfPhasePropagator->setCurrentIndependentVariable(this->state_at_end_of_phase(7));
                this->BackwardHalfPhasePropagator->propagate(-BackwardOutputTimestep * (backstep + 0.5), false);

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

                //call the power model
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                //call the propulsion model
                this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle);

                //store the thruster model outputs
                doubleType max_thrust = this->mySpacecraft->getEPthrust() * 1.0e-3;
                doubleType max_mass_flow_rate = this->mySpacecraft->getEPMassFlowRate();
                doubleType Isp = this->mySpacecraft->getEPIsp();
                doubleType power = this->mySpacecraft->getAvailablePower();
                doubleType active_power = this->mySpacecraft->getEPActivePower();
                size_t number_of_active_engines = this->mySpacecraft->getEPNumberOfActiveThrusters();
                size_t ThrottleLevel = this->mySpacecraft->getEPThrottleLevel();
                std::string ThrottleLevelString = this->mySpacecraft->getEPThrottleLevelString();

                math::Matrix<doubleType> velocityUnitVector = output_state.getSubMatrix1D(3, 5).unitize();

                if (this->myJourneyOptions->thrust_control_law == ThrustControlLaw::Velocity)
                {
                    ThrustVector = velocityUnitVector * max_thrust * 1000.0 / this->PhaseDutyCycle;
                }
                else if (this->myJourneyOptions->thrust_control_law == ThrustControlLaw::AntiVelocity)
                {
                    ThrustVector = -velocityUnitVector * max_thrust * 1000.0 / this->PhaseDutyCycle;
                }

                //print
                this->write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "FBLTthrust",//event_type
                    "deep-space",//event_location
                    BackwardOutputTimestep / 86400.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    this->output_state,//state
                    empty3,//dV
                    ThrustVector,//ThrustVector
                    1.0,//throttle
                    max_thrust * 1000.0 / this->PhaseDutyCycle,//Thrust
                    Isp,//Isp
                    power,//AvailPower
                    max_mass_flow_rate / this->PhaseDutyCycle,//mdot
                    number_of_active_engines,//number_of_active_engines
                    active_power,
                    ThrottleLevelString);//active_power)
            }

            //reset the propagator to go where it is supposed to
            this->BackwardHalfPhasePropagator->setStateRight(this->match_point_state_plus);

            //Step 6: print the arrival event
            this->myArrivalEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
        }


        void ControlLawThrustPhase::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
        {
            math::Matrix<doubleType> ThrustVector;

            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //Step 1: output the departure event
            this->myDepartureEvent->output_ephemeris(outputfile);

            if (this->myOptions->generate_acceleration_model_instrumentation_file && !this->isKeplerian)
            {
                this->temp_state = this->myDepartureEvent->get_state_before_event();
                this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
            }

            //Step 2: print the phase
            {
                //Step 2.1: we'll need an output vector
                this->temp_state = this->output_state;
                //Step 1: temporarily assign the propagator to the output state
                this->ForwardHalfPhasePropagator->setStateLeft(this->temp_state);
                this->ForwardHalfPhasePropagator->setStateRight(output_state);
                this->ForwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
                this->temp_state = this->state_after_initial_TCM;

                //Step 3.2: propagate and print, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
				double propagated_time = 0;
                while (propagated_time < this->PhaseFlightTime)
                {
                    //Step 3.2.1: propagate
                    this->Forward_dPropagatedStatedIndependentVariable.assign_zeros();
                    this->ForwardHalfPhasePropagator->setCurrentEpoch(this->temp_state(7));
                    this->ForwardHalfPhasePropagator->setCurrentIndependentVariable(this->temp_state(7));
                    this->ForwardHalfPhasePropagator->propagate(timeToPropagate, false);
                    this->temp_state = output_state;
					propagated_time += timeToPropagate;

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
                    if (this->myOptions->generate_acceleration_model_instrumentation_file && !this->isKeplerian)
                    {
                        this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
                    }

                    //where am I relative to the Sun?
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

                    //call the power model
                    doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                    this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                    //call the propulsion model
                    this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle);

                    //store the thruster model outputs
                    doubleType max_thrust = this->mySpacecraft->getEPthrust() * 1.0e-3;
                    doubleType max_mass_flow_rate = this->mySpacecraft->getEPMassFlowRate();
                    doubleType Isp = this->mySpacecraft->getEPIsp();
                    doubleType power = this->mySpacecraft->getAvailablePower();
                    doubleType active_power = this->mySpacecraft->getEPActivePower();
                    doubleType number_of_active_engines = this->mySpacecraft->getEPNumberOfActiveThrusters();
                    size_t ThrottleLevel = this->mySpacecraft->getEPThrottleLevel();
                    std::string ThrottleLevelString = this->mySpacecraft->getEPThrottleLevelString();

                    math::Matrix<doubleType> velocityUnitVector = output_state.getSubMatrix1D(3, 5).unitize();

                    if (this->myJourneyOptions->thrust_control_law == ThrustControlLaw::Velocity)
                    {
                        ThrustVector = velocityUnitVector * max_thrust * 1000.0 / this->PhaseDutyCycle;
                    }
                    else if (this->myJourneyOptions->thrust_control_law == ThrustControlLaw::AntiVelocity)
                    {
                        ThrustVector = -velocityUnitVector * max_thrust * 1000.0 / this->PhaseDutyCycle;
                    }

                    this->write_ephemeris_line(outputfile,
                        output_state, 
                        ThrustVector,
                        this->mySpacecraft->getEPthrust() * 1.0e-3,
                        this->mySpacecraft->getEPMassFlowRate(),
                        this->mySpacecraft->getEPIsp(),
                        this->mySpacecraft->getEPNumberOfActiveThrusters(),
                        this->mySpacecraft->getEPActivePower(),
                        this->mySpacecraft->getEPThrottleLevelString());

                    //Step 3.2.4: increment propagatedEpoch
                    timeToPropagate = (this->PhaseFlightTime _GETVALUE - propagated_time) > this->EphemerisOutputResolution
                        ? this->EphemerisOutputResolution
                        : (this->PhaseFlightTime _GETVALUE - propagated_time);
                }

                //Step 3.3: reset the propagator to its original state vectors
                this->ForwardHalfPhasePropagator->setStateLeft(this->state_after_initial_TCM);
                this->ForwardHalfPhasePropagator->setStateRight(this->match_point_state_minus);
            }

            //Step 3: print the arrival event
            this->myArrivalEvent->output_ephemeris(outputfile);
            this->temp_state = this->myArrivalEvent->get_state_before_event();
            if (this->myOptions->generate_acceleration_model_instrumentation_file && !this->isKeplerian)
            {
                this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
            }
        }//end output_ephemeris()
        
        void ControlLawThrustPhase::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
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
                    this->state_after_initial_TCM(7),
                    this->state_after_initial_TCM); // this is the state immediately before we start thrusting

                // write target spec object
                myTargetSpecLine.write(target_spec_file);

                // disable the target flag
                haveManeuverNeedTarget = false;
            }

            //write a VNB maneuver every time the throttle level changes
            maneuver_spec_line myManeuverSpecLine(this->name);
            math::Matrix<doubleType> ControlVector;

            if (this->myJourneyOptions->thrust_control_law == ThrustControlLaw::Velocity)
            {
                ControlVector = math::Matrix<doubleType>(3, 1, { 1.0, 0.0, 0.0 });
            }
            else if (this->myJourneyOptions->thrust_control_law == ThrustControlLaw::AntiVelocity)
            {
                ControlVector = math::Matrix<doubleType>(3, 1, { -1.0, 0.0, 0.0 });
            }

            size_t ManeuverThrottleLevel;
            doubleType ManeuverStartEpoch;
            doubleType ManeuverStartMass;
            doubleType ManeuverThrustMagnitude;
            doubleType ManeuverMassFlowRate;

            //forward half-phase
            {
                // Step through the forward half-phase at IntegrationStep intervals
                // Every time the throttle level changes save off a thrust step
                // Compute the propulsion characteristics on the left side of the phase
                {
                    // Locate the spacecraft relative to the sun
                    for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
                    {
                        this->output_state(stateIndex) = this->state_after_initial_TCM(stateIndex);
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

                    //Step 2.2.1.3: call the thruster model
                    this->mySpacecraft->computeElectricPropulsionPerformance(this->PhaseDutyCycle);

                    //Step 2.2.1.3: populate fields
                    ManeuverThrottleLevel = this->mySpacecraft->getEPThrottleLevel();
                    ManeuverStartEpoch = this->output_state(7);
                    ManeuverStartMass = this->output_state(6);
                    ManeuverThrustMagnitude = this->mySpacecraft->getEPthrust(); //in N
                    ManeuverMassFlowRate = this->mySpacecraft->getEPMassFlowRate();

                } // end segment left side propulsion characteristics 

                //Step 2.2.2.1: hijack the propagator
                this->ForwardHalfPhasePropagator->setStateRight(this->output_state);

                //Step 2.2.2.2: propagate the control segment, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = (this->PhaseFlightTime * this->MatchPointFraction) _GETVALUE;
                double divisor = 1.0;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 2.2.2.3: propagate
                    this->ForwardHalfPhasePropagator->setCurrentEpoch(this->state_after_initial_TCM(7));
                    this->ForwardHalfPhasePropagator->propagate(timeToPropagate, false);

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
                    this->mySpacecraft->computeElectricPropulsionPerformance(1.0); //compute with 100% duty cycle because we want to write out actual thrust/mdot

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
                            myManeuverSpecLine.append_maneuver_spec_item("VNB",
                                ManeuverStartEpoch,
                                ControlVector,
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

                }

                // reset the propagator right-hand state container
                this->ForwardHalfPhasePropagator->setStateRight(this->match_point_state_minus);
            } // end forward propagation

            //backward propagation
            {
                //Step 2.2.2.1: hijack the propagator
                //Step 3.1: temporarily assign the propagator to the output state
                this->BackwardHalfPhasePropagator->setStateRight(this->output_state);

                //Step 2.2.2.2: propagate the control segment, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = (this->PhaseFlightTime * (1.0 - this->MatchPointFraction)) _GETVALUE;

                double divisor = 1.0;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 2.2.2.3: propagate
                    // propagate up to the current location in the thrust segment
                    this->BackwardHalfPhasePropagator->setCurrentEpoch(this->match_point_state_plus(7));
                    this->BackwardHalfPhasePropagator->propagate(timeToPropagate, false);

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
                    this->mySpacecraft->computeElectricPropulsionPerformance(1.0); //compute with 100% duty cycle because we want to write out actual thrust/mdot

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
                                ControlVector,
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
                this->BackwardHalfPhasePropagator->setStateRight(this->match_point_state_plus);
            } //end backward propagation

            //save off the last thrust step
            myManeuverSpecLine.append_maneuver_spec_item("VNB",
                ManeuverStartEpoch,
                ControlVector,
                ManeuverStartMass,
                this->state_at_end_of_phase(6),
                ManeuverThrustMagnitude,
                ManeuverMassFlowRate,
                this->state_at_end_of_phase(7) - ManeuverStartEpoch,
                this->PhaseDutyCycle);

            // write the maneuver spec line
            myManeuverSpecLine.write(maneuver_spec_file);

            //set the have maneuver, need target flag
            haveManeuverNeedTarget = true;

            //should write out a maneuver/target spec for the arrival maneuver if there is one
            this->myArrivalEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
        }//end output_maneuver_and_target_spec()

    }//end namespace Phases
}//end namespace EMTG