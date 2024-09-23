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

//EMTGv9 CoastPhase
//Jacob Englander 11-20-2017

#include "CoastPhase.h"

#include "PropagatorFactory.h"
#include "IntegrationSchemeFactory.h"


namespace EMTG
{
    namespace Phases
    {
        CoastPhase::CoastPhase(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            phase* previousPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            HardwareModels::LaunchVehicle* myLaunchVehicle,
            missionoptions* myOptions) :
            TwoPointShootingPhase::TwoPointShootingPhase(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                previousPhase,
                Universe,
                mySpacecraft,
                myLaunchVehicle,
                myOptions,
                9, //numStatesToPropagate
                8) //numMatchConstraints
        {
            this->initialize();
        }//end constructor

        CoastPhase::CoastPhase(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            phase* previousPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            HardwareModels::LaunchVehicle* myLaunchVehicle,
            missionoptions* myOptions,
            const size_t& numStatesToPropagate,
            const size_t& numMatchConstraints) :
            TwoPointShootingPhase::TwoPointShootingPhase(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                previousPhase,
                Universe,
                mySpacecraft,
                myLaunchVehicle,
                myOptions,
                numStatesToPropagate,
                numMatchConstraints)
        {
        }

        void CoastPhase::initialize()
        {
            this->num_timesteps = this->myJourneyOptions->override_num_steps
                ? this->myJourneyOptions->number_of_steps
                : this->myOptions->num_timesteps;

            this->stateVectorNames.push_back("virtual chemical fuel");

            //name our match point constraints
            this->matchPointConstraintNames.push_back("virtual chemical fuel");
            this->matchPointConstraintStateIndex.push_back(8); //note we skipped the encode of boundary epoch, it does not get a match point constraint
            this->continuity_constraint_scale_factors(7) = 1.0 / (1.0 * this->myJourneyOptions->maximum_mass);
            this->stateIndex_phase_propagation_variable = this->numStatesToPropagate;
            this->stateVectorNames.push_back("phase flight time");

            //no maneuvers
            this->hasBipropManeuver = false;
            this->hasMonopropManeuver = false;
            this->hasElectricManeuver = false;

            //size the TCM transition matrix
            this->TCMTM = math::Matrix<double>(this->numStatesToPropagate + 1, math::identity);

            //configure propagators
            this->MatchPointFraction = this->myJourneyOptions->CoastPhaseMatchPointFraction;
            this->dForwardStepSize_dPropagationVariable = this->MatchPointFraction;
            this->dBackwardStepSize_dPropagationVariable = 1.0 - this->MatchPointFraction;
            this->ForwardIntegrationStepLength = this->myJourneyOptions->CoastPhaseForwardIntegrationStepLength;
            this->BackwardIntegrationStepLength = this->myJourneyOptions->CoastPhaseBackwardIntegrationStepLength;

            //which propagator do we want?
            if (this->myPropagatorType == PropagatorType::KeplerPropagator)
            {
                this->isKeplerian = true;
                this->ForwardSTM = math::Matrix<double>(6, math::identity);
                this->BackwardSTM = math::Matrix<double>(6, math::identity);
                this->Forward_dPropagatedStatedIndependentVariable = math::Matrix<double>(6, 1, 0.0);
                this->Backward_dPropagatedStatedIndependentVariable = math::Matrix<double>(6, 1, 0.0);

                this->ForwardHalfPhasePropagator = Astrodynamics::CreatePropagator(this->myOptions,
                    this->myUniverse,
                    6,
                    this->state_after_initial_TCM,
                    this->match_point_state_minus,
                    this->ForwardSTM,
                    this->Forward_dPropagatedStatedIndependentVariable,
                    &this->dForwardStepSize_dPropagationVariable);

                this->BackwardHalfPhasePropagator = Astrodynamics::CreatePropagator(this->myOptions,
                    this->myUniverse,
                    6,
                    this->state_at_end_of_phase,
                    this->match_point_state_plus,
                    this->BackwardSTM,
                    this->Backward_dPropagatedStatedIndependentVariable,
                    &this->dBackwardStepSize_dPropagationVariable);
            }
            else //integrated propagator
            {
                this->isKeplerian = false;
                this->ForwardSTM = math::Matrix<double>(10, math::identity);
                this->BackwardSTM = math::Matrix<double>(10, math::identity);
                this->Forward_dPropagatedStatedIndependentVariable.resize(this->numStatesToPropagate + 1, 2, 0.0);
                this->Backward_dPropagatedStatedIndependentVariable.resize(this->numStatesToPropagate + 1, 2, 0.0);

                //acceleration model object
                this->mySpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(this->myOptions,
                    this->myJourneyOptions,
                    this->myUniverse,
                    this->Xdescriptions,
                    this->mySpacecraft,
                    10); // STM size
                this->mySpacecraftAccelerationModel->setDutyCycle(this->PhaseDutyCycle);

                //EOM
                this->myEOM.setSpacecraftAccelerationModel(this->mySpacecraftAccelerationModel);

                //integration scheme
                this->myIntegrationScheme = CreateIntegrationScheme(&this->myEOM, this->numStatesToPropagate, 10);

                this->ForwardHalfPhasePropagator = CreatePropagator(this->myOptions,
                    this->myUniverse,
                    this->numStatesToPropagate,
                    10,
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
                    10,
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

                //mass variables do not affect other constraints
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
            }
            this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[7][7] = false; //fuel mass wrt epoch
            //this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[7][9] = false; //fuel mass wrt oxidizer mass
            this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[7][7] = false; //fuel mass wrt epoch
            //this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[7][9] = false; //fuel mass wrt oxidizer mass
            //always true - propellant tank states do not affect anything but themselves
            for (size_t constraintIndex = 0; constraintIndex < 7; ++constraintIndex)
            {
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][8] = false; //other states wrt fuel mass
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][8] = false; //other states wrt fuel mass
            }

            //output stuff
            this->output_state.resize(this->numStatesToPropagate + 10 * 10, 1, 0.0);
        }//end initialize

        CoastPhase::~CoastPhase()
        {
            delete this->ForwardHalfPhasePropagator;
            delete this->BackwardHalfPhasePropagator;

            if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
            {
                delete this->mySpacecraftAccelerationModel;
                delete this->myIntegrationScheme;
            }
        }//end destructor

        //************************************calcbounds methods
        void CoastPhase::setup_calcbounds(std::vector<double>* Xupperbounds,
            std::vector<double>* Xlowerbounds,
            std::vector<double>* X_scale_factors,
            std::vector<double>* Fupperbounds,
            std::vector<double>* Flowerbounds,
			std::vector<double>* F_scale_factors,
            std::vector<std::string>* Xdescriptions,
            std::vector<std::string>* Fdescriptions,
            std::vector<size_t>* iGfun,
            std::vector<size_t>* jGvar,
            std::vector<std::string>* Gdescriptions,
            std::vector<size_t>* iAfun,
            std::vector<size_t>* jAvar,
            std::vector<std::string>* Adescriptions,
            std::vector<double>* A)
        {
            //base class
            this->phase::setup_calcbounds(Xupperbounds,
                Xlowerbounds,
                X_scale_factors,
                Fupperbounds,
                Flowerbounds,
				F_scale_factors,
                Xdescriptions,
                Fdescriptions,
                iGfun,
                jGvar,
                Gdescriptions,
                iAfun,
                jAvar,
                Adescriptions,
                A);
        }//end setup_calcbounds()

        void CoastPhase::calcbounds()
        {
            this->calcbounds_phase_left_boundary();

            this->calcbounds_phase_flight_time();

            this->calcbounds_phase_right_boundary();

            this->calcbounds_phase_main();
        }//end calcbounds()

        void CoastPhase::calcbounds_phase_main()
        {
            this->calcbounds_virtual_propellant_tanks();

            this->calcbounds_match_point_constraints();

            this->calcbounds_deltav_contribution();
        }//end calcbounds_phase_main()

        void CoastPhase::calcbounds_match_point_constraints()
        {
            //base class
            this->TwoPointShootingPhase::calcbounds_match_point_constraints();

            //virtual chemical fuel
            this->create_sparsity_entry(this->Findices_match_point_constraints[7],
                this->X_index_of_first_decision_variable_in_this_phase,
                true,
                this->prefix + "virtual chemical fuel",
                this->Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel);
        }//end calcbounds_match_point_constraints()
        
        void CoastPhase::calcbounds_virtual_propellant_tanks()
        {
            //encode virtual propellant tank variables

            //virtual chemical fuel
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
        }//end calcbounds_virtual_propellant_tanks()



        void CoastPhase::calcbounds_deltav_contribution()
        {
            //nothing right now
        }//end calcbounds_deltav_contribution()

        //************************************process methods
        void CoastPhase::process_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //do we want to propagate STMs or not?
            this->total_number_of_states_to_integrate = needG
                ? this->numStatesToPropagate + 10 * 10
                : this->numStatesToPropagate;

            this->process_phase_left_boundary(X, Xindex, F, Findex, G, needG);

            this->process_phase_flight_time(X, Xindex, F, Findex, G, needG);

            this->process_phase_right_boundary(X, Xindex, F, Findex, G, needG);

            this->process_phase_main(X, Xindex, F, Findex, G, needG);
        }//end process_phase()

        void CoastPhase::process_phase_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_virtual_propellant_tanks(X, Xindex, F, Findex, G, needG);

            this->process_forward_half_phase(X, Xindex, F, Findex, G, needG);

            this->process_backward_half_phase(X, Xindex, F, Findex, G, needG);

            this->process_match_point_constraints(X, Xindex, F, Findex, G, needG);

            this->process_deltav_contribution(X, Xindex, F, Findex, G, needG);
        }//end process_phase_main()

        void CoastPhase::process_forward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: set the tank states
            this->state_after_initial_TCM(8) = this->chemical_fuel_used; //TCM propellant, usually zero

            //Step 2: call the forward propagator
            this->Forward_dPropagatedStatedIndependentVariable.assign_zeros();
            this->ForwardHalfPhasePropagator->setCurrentEpoch(this->state_after_initial_TCM(7));
            this->ForwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
            this->ForwardHalfPhasePropagator->setCurrentIndependentVariable(this->state_after_initial_TCM(7));
            this->ForwardHalfPhasePropagator->propagate(this->PhaseFlightTime * this->MatchPointFraction, needG);

            if (this->isKeplerian) //otherwise this happens in the propagator
            {
                doubleType ACS_fuel_used = (this->myOptions->trackACS ? this->myOptions->ACS_kg_per_day * this->PhaseFlightTime  * this->MatchPointFraction / 86400.0 : 0.0);

                this->match_point_state_minus(6) = this->state_after_initial_TCM(6) - ACS_fuel_used;
                this->match_point_state_minus(7) = this->state_after_initial_TCM(7) + this->PhaseFlightTime * this->MatchPointFraction;
                this->match_point_state_minus(8) = this->state_after_initial_TCM(8) + ACS_fuel_used; //virtual fuel

                if (this->myOptions->trackACS)
                {
                    //mass
                    this->ForwardSPTM(6, 6) = (this->match_point_state_minus(6) / this->state_after_initial_TCM(6)) _GETVALUE;
                    this->ForwardSPTM(6, this->stateIndex_phase_propagation_variable) = -this->MatchPointFraction * this->myOptions->ACS_kg_per_day / 86400.0 * this->ForwardSPTM(6, 6);

                    //tanks
                    this->ForwardSPTM(8, 6) = ((this->chemical_fuel_used - ACS_fuel_used) / this->match_point_state_minus(6)) _GETVALUE;
                    this->ForwardSPTM(8, this->stateIndex_phase_propagation_variable) = this->MatchPointFraction * this->myOptions->ACS_kg_per_day / 86400.0;
                }
            }

            //MTM entries for ACS
            if (needG)
            {
                //Step 3: form the SPTM
                //Step 3.1: upper left 6x6 is the regular STM
                size_t stateMax = this->isKeplerian ? 6 : 10;
                for (size_t i = 0; i < stateMax; ++i)
                    for (size_t j = 0; j < stateMax; ++j)
                        this->ForwardSPTM(i, j) = this->ForwardSTM(i, j);

                //Step 3.2: turn the upper right 10 x 2 into the Phi_t terms developed by Lantoine
                for (size_t i = 0; i < stateMax; ++i)
                {
                    //this really ought to be cleaned up...
                    if (this->isKeplerian)
                    {
                        this->ForwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->Forward_dPropagatedStatedIndependentVariable(i);
                    }
                    else //integrated propagator
                    {
                        this->ForwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->ForwardSTM(i, 9);
                    }
                }
            }
        }//end process_forward_half_phase()

        void CoastPhase::process_backward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: set the tank states
            this->state_at_end_of_phase(8) = this->virtual_chemical_fuel_used;

            //Step 2: call the backward propagator
            this->Backward_dPropagatedStatedIndependentVariable.assign_zeros();
            this->BackwardHalfPhasePropagator->setCurrentEpoch(this->state_at_end_of_phase(7));
            this->BackwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
            this->BackwardHalfPhasePropagator->setCurrentIndependentVariable(this->state_at_end_of_phase(7));
            this->BackwardHalfPhasePropagator->propagate(-this->PhaseFlightTime * (1.0 - this->MatchPointFraction), needG);


            if (this->isKeplerian) //otherwise this happens in the propagator
            {
                doubleType ACS_fuel_used = (this->myOptions->trackACS ? this->myOptions->ACS_kg_per_day * this->PhaseFlightTime * (1.0 - this->MatchPointFraction) / 86400.0 : 0.0);

                this->match_point_state_plus(6) = this->state_at_end_of_phase(6) + ACS_fuel_used;
                this->match_point_state_plus(7) = this->state_at_end_of_phase(7) - this->PhaseFlightTime * (1.0 - this->MatchPointFraction);
                this->match_point_state_plus(8) = this->state_at_end_of_phase(8) - ACS_fuel_used; //virtual fuel

                if (this->myOptions->trackACS)
                {
                    //mass
                    this->BackwardSPTM(6, 6) = (this->match_point_state_plus(6) / this->state_at_end_of_phase(6)) _GETVALUE;
                    this->BackwardSPTM(6, this->stateIndex_phase_propagation_variable) = (1.0 - this->MatchPointFraction) * this->myOptions->ACS_kg_per_day / 86400.0 * this->BackwardSPTM(6, 6);

                    //tanks
                    this->BackwardSPTM(8, 6) = -((this->chemical_fuel_used - ACS_fuel_used) / this->state_at_end_of_phase(6)) _GETVALUE;
                    this->BackwardSPTM(8, this->stateIndex_phase_propagation_variable) = -(1.0 - this->MatchPointFraction) * this->myOptions->ACS_kg_per_day / 86400.0;
                }
            }

                                                                                                 
            if (needG)
            {
                //Step 3: form the SPTM
                //Step 3.1: upper left 6x6 is the regular STM
                size_t stateMax = this->isKeplerian ? 6 : 10;
                for (size_t i = 0; i < stateMax; ++i)
                    for (size_t j = 0; j < stateMax; ++j)
                        this->BackwardSPTM(i, j) = this->BackwardSTM(i, j);
                //Step 3.2: turn the upper right 6 x 2 into the Phi_t terms developed by Lantoine
                for (size_t i = 0; i < 6; ++i)
                {
                    //this really ought to be cleaned up...
                    if (this->isKeplerian)
                    {
                        this->BackwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->Backward_dPropagatedStatedIndependentVariable(i);
                    }
                    else //integrated propagator
                    {
                        this->BackwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->BackwardSTM(i, 9);
                    }
                    
                }
            }
        }//end process_backward_half_phase()

        void CoastPhase::process_match_point_constraints(const std::vector<doubleType>& X,
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
            }//end derivatives
        }//end process_match_point_constraints()

        void CoastPhase::process_virtual_propellant_tanks(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->virtual_chemical_fuel_used = X[Xindex++];
        }//end process_virtual_propellant_tanks()

        void CoastPhase::process_deltav_contribution(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->PhaseTotalDeterministicDeltav = this->myDepartureEvent->get_DeterministicDeltav() + this->myArrivalEvent->get_DeterministicDeltav();
            this->PhaseTotalStatisticalDeltav = this->initial_TCM_magnitude;
        }


        void CoastPhase::output(std::ofstream& outputfile,
            size_t& eventcount)
        {
            doubleType ForwardOutputTimestep = this->PhaseFlightTime * this->MatchPointFraction / (this->num_timesteps / 2);
            doubleType BackwardOutputTimestep = this->PhaseFlightTime * (1.0 - this->MatchPointFraction) / (this->num_timesteps / 2);
            math::Matrix<doubleType> empty3(3, 1, 0.0);

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

                if (this->isKeplerian)
                {
                    doubleType ACS_fuel_used = (this->myOptions->trackACS ? this->myOptions->ACS_kg_per_day * ForwardOutputTimestep * (step + 0.5) / 86400.0 : 0.0);
                    output_state(6) = this->state_after_initial_TCM(6) - ACS_fuel_used;
                    output_state(7) = this->state_after_initial_TCM(7) + ForwardOutputTimestep * (step + 0.5);
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
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                doubleType power = this->mySpacecraft->getAvailablePower();
                doubleType active_power = this->mySpacecraft->getEPActivePower();

                //print
                this->write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "coast",//event_type
                    "deep-space",//event_location
                    ForwardOutputTimestep / 86400.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    output_state,//state
                    empty3,//delta-v
                    empty3,//ThrustVector
                    0.0,//dVmag
                    0.0,//Thrust
                    0.0,//Isp
                    power,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    active_power, //active_power
                    "none");
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

                if (this->isKeplerian)
                {
                    doubleType ACS_fuel_used = (this->myOptions->trackACS ? this->myOptions->ACS_kg_per_day * BackwardOutputTimestep * (step + 0.5) / 86400.0 : 0.0);
                    output_state(6) = this->state_at_end_of_phase(6) + ACS_fuel_used;
                    output_state(7) = this->state_at_end_of_phase(7) - BackwardOutputTimestep * (backstep + 0.5);
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
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                doubleType power = this->mySpacecraft->getAvailablePower();
                doubleType active_power = this->mySpacecraft->getEPActivePower();

                //print
                this->write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "coast",//event_type
                    "deep-space",//event_location
                    BackwardOutputTimestep / 86400.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    output_state,//state
                    empty3,//delta-v
                    empty3,//ThrustVector
                    0.0,//dVmag
                    0.0,//Thrust
                    0.0,//Isp
                    power,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    active_power,
                    "none");//active_power
            }

            //reset the propagator to go where it is supposed to
            this->BackwardHalfPhasePropagator->setStateRight(this->match_point_state_plus);

            //Step 6: print the arrival event
            this->myArrivalEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
        }


        void CoastPhase::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
        {
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

                    if (this->isKeplerian)
                    {
                        doubleType ACS_fuel_used = (this->myOptions->trackACS ? this->myOptions->ACS_kg_per_day * propagated_time / 86400.0 : 0.0);
                        output_state(6) = this->state_after_initial_TCM(6) - ACS_fuel_used;
                        output_state(7) = this->state_after_initial_TCM(7) + propagated_time;
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
                    if (this->myOptions->generate_acceleration_model_instrumentation_file && !this->isKeplerian)
                    {
                        this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
                    }

                    this->write_ephemeris_line(outputfile,
                        output_state,
                        math::Matrix<doubleType>(3, 1, 0.0),//control vector
                        0.0,
                        0.0,
                        0.0,
                        0,
                        0.0,
                        "none");

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

        void CoastPhase::output_STMs()
        {
            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::ofstream summaryfile(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + ".stm_description");
            for (size_t stateIndex = 0; stateIndex < this->stateVectorNames.size(); ++stateIndex)
                summaryfile << this->stateVectorNames[stateIndex] << std::endl;
            summaryfile.close();

            this->ForwardSPTM.print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_STM_0_forward.stm");
            this->BackwardSPTM.print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_STM_1_backward.stm");
        }//end output_STMs()
        
        void CoastPhase::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //should write out a maneuver/target spec for the departure maneuver if there is one
            this->myDepartureEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);

            //coast phases don't write maneuver and target spec files anyway, so nothing for the phase body

            //should write out a maneuver/target spec for the arrival maneuver if there is one
            this->myArrivalEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
        }//end output_maneuver_and_target_spec()

    }//end namespace Phases
}//end namespace EMTG