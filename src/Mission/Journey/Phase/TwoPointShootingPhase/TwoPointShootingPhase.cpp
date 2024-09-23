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

//EMTGv9 TwoPointShootingPhase
//Jacob Englander 6-22-2017

#include "TwoPointShootingPhase.h"

namespace EMTG
{
    namespace Phases
    {
        TwoPointShootingPhase::TwoPointShootingPhase(const std::string& name,
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
            phase::phase(name,
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
            this->match_point_state_minus.resize(this->numStatesToPropagate, 1, 0.0);
            this->match_point_state_plus.resize(this->numStatesToPropagate, 1, 0.0);

            this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables.resize(this->numMatchConstraints);
            this->G_indices_match_point_constraints_wrt_RightBoundaryVariables.resize(this->numMatchConstraints);
            this->G_indices_match_point_constraints_wrt_LeftBoundaryTimeVariables.resize(this->numMatchConstraints);
            this->G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables.resize(this->numMatchConstraints);
            this->G_indices_match_point_constraints_wrt_PhaseFlightTime.resize(this->numMatchConstraints);

            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
            {
                this->matchPointConstraintNames.push_back(this->stateVectorNames[stateIndex]);
                this->matchPointConstraintStateIndex.push_back(stateIndex);
            }
        }//end constructor

        //******************************************calcbounds methods
        void TwoPointShootingPhase::calcbounds_match_point_constraints()
        {

            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterPhase = this->myArrivalEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase_wrt_Time = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterPhase_wrt_Time = this->myArrivalEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value

            //Step 1: make a list of variables affecting both boundaries
            //Step 1.1: Left Boundary
            for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforePhase.size(); ++dIndex)
            {
                bool alreadyListed = false;
                size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase[dIndex]);
                for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingLeftBoundary.size(); ++listIndex)
                {
                    if (this->ListOfVariablesAffectingLeftBoundary[listIndex] == Xindex)
                    {
                        alreadyListed = true;
                        break;
                    }
                }
                if (!alreadyListed)
                {
                    this->ListOfVariablesAffectingLeftBoundary.push_back(Xindex);
                }
            }

            this->DerivativesOfLeftBoundaryByVariable.resize(ListOfVariablesAffectingLeftBoundary.size());

            for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingLeftBoundary.size(); ++listIndex)
            {
                for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforePhase.size(); ++dIndex)
                {
                    size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase[dIndex]);

                    if (Xindex == this->ListOfVariablesAffectingLeftBoundary[listIndex])
                    {
                        this->DerivativesOfLeftBoundaryByVariable[listIndex].push_back(dIndex);
                    }
                }
            }

            //Step 1.2: Right Boundary
            for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterPhase.size(); ++dIndex)
            {
                bool alreadyListed = false;
                size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase[dIndex]);
                for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingRightBoundary.size(); ++listIndex)
                {
                    if (this->ListOfVariablesAffectingRightBoundary[listIndex] == Xindex)
                    {
                        alreadyListed = true;
                        break;
                    }
                }
                if (!alreadyListed)
                {
                    this->ListOfVariablesAffectingRightBoundary.push_back(Xindex);
                }
            }

            this->DerivativesOfRightBoundaryByVariable.resize(ListOfVariablesAffectingRightBoundary.size());

            for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingRightBoundary.size(); ++listIndex)
            {
                for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterPhase.size(); ++dIndex)
                {
                    size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase[dIndex]);

                    if (Xindex == this->ListOfVariablesAffectingRightBoundary[listIndex])
                    {
                        this->DerivativesOfRightBoundaryByVariable[listIndex].push_back(dIndex);
                    }
                }
            }

            //Step 2: make a list of time variables affecting both boundaries
            {
                //Step 2.1: Left Boundary
                for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforePhase_wrt_Time.size(); ++dIndex)
                {
                    bool alreadyListed = false;
                    size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);
                    for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingLeftBoundary.size(); ++listIndex)
                    {
                        if (this->ListOfTimeVariablesAffectingLeftBoundary[listIndex] == Xindex)
                        {
                            alreadyListed = true;
                            break;
                        }
                    }
                    if (!alreadyListed)
                    {
                        this->ListOfTimeVariablesAffectingLeftBoundary.push_back(Xindex);
                    }
                }

                this->DerivativesOfLeftBoundaryByTimeVariable.resize(ListOfTimeVariablesAffectingLeftBoundary.size());

                for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingLeftBoundary.size(); ++listIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforePhase_wrt_Time.size(); ++dIndex)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);

                        if (Xindex == this->ListOfTimeVariablesAffectingLeftBoundary[listIndex])
                        {
                            this->DerivativesOfLeftBoundaryByTimeVariable[listIndex].push_back(dIndex);
                        }
                    }
                }

                //Step 2.2: Right Boundary
                for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterPhase_wrt_Time.size(); ++dIndex)
                {
                    bool alreadyListed = false;
                    size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);
                    for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingRightBoundary.size(); ++listIndex)
                    {
                        if (this->ListOfTimeVariablesAffectingRightBoundary[listIndex] == Xindex)
                        {
                            alreadyListed = true;
                            break;
                        }
                    }
                    if (!alreadyListed)
                    {
                        this->ListOfTimeVariablesAffectingRightBoundary.push_back(Xindex);
                    }
                }

                this->DerivativesOfRightBoundaryByTimeVariable.resize(ListOfTimeVariablesAffectingRightBoundary.size());

                for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingRightBoundary.size(); ++listIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterPhase_wrt_Time.size(); ++dIndex)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);

                        if (Xindex == this->ListOfTimeVariablesAffectingRightBoundary[listIndex])
                        {
                            this->DerivativesOfRightBoundaryByTimeVariable[listIndex].push_back(dIndex);
                        }
                    }
                }
            }

            //Step 3: construct the match point constraints
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
            {
                Flowerbounds->push_back(-math::SMALL);
                Fupperbounds->push_back(math::SMALL);
                Fdescriptions->push_back(prefix + "match point " + this->matchPointConstraintNames[constraintIndex]);
                this->Findices_match_point_constraints.push_back(Fdescriptions->size() - 1);

                //derivatives with respect to phase flight time
                if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                {
                    this->create_sparsity_entry(this->Findices_match_point_constraints[constraintIndex],
                        this->Xindex_PhaseFlightTime,
                        this->G_indices_match_point_constraints_wrt_PhaseFlightTime[constraintIndex]);
                }


                //derivatives with respect to decision variables that affect boundary states
                //with respect to left boundary
                std::vector<bool> constraint_LeftBoundaryVariableAffectsMatchConstraint(this->ListOfVariablesAffectingLeftBoundary.size(), false);
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingLeftBoundary.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingLeftBoundary[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->DerivativesOfLeftBoundaryByVariable[varIndex][entryIndex];
                        size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase[dIndex]);

                        //hooray, this variable affects the constraint. Make the Jacobian entry!
                        if (this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][stateIndex])
                        {
                            constraint_LeftBoundaryVariableAffectsMatchConstraint[varIndex] = true;
                            break;
                        }
                    }

                    if (constraint_LeftBoundaryVariableAffectsMatchConstraint[varIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints[constraintIndex],
                            Xindex,
                            this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables[constraintIndex]);
                    }
                    else
                        this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables[constraintIndex].push_back(0);
                }
                this->LeftBoundaryVariableAffectsMatchConstraint.push_back(constraint_LeftBoundaryVariableAffectsMatchConstraint);
                
                //with respect to right boundary
                std::vector<bool> constraint_RightBoundaryVariableAffectsMatchConstraint(this->ListOfVariablesAffectingRightBoundary.size(), false);
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingRightBoundary.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingRightBoundary[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->DerivativesOfRightBoundaryByVariable[varIndex][entryIndex];
                        size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase[dIndex]);

                        //hooray, this variable affects the constraint. Make the Jacobian entry!
                        if (this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][stateIndex])
                        {
                            constraint_RightBoundaryVariableAffectsMatchConstraint[varIndex] = true;
                            break;
                        }
                    }

                    if (constraint_RightBoundaryVariableAffectsMatchConstraint[varIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints[constraintIndex],
                            Xindex,
                            this->G_indices_match_point_constraints_wrt_RightBoundaryVariables[constraintIndex]);
                    }
                    else
                        this->G_indices_match_point_constraints_wrt_RightBoundaryVariables[constraintIndex].push_back(0);
                }
                this->RightBoundaryVariableAffectsMatchConstraint.push_back(constraint_RightBoundaryVariableAffectsMatchConstraint);

                //derivatives with respect to time variables that affect boundary states
                //with respect to left boundary
                std::vector<bool> constraint_LeftBoundaryTimeVariableAffectsMatchConstraint(this->ListOfTimeVariablesAffectingLeftBoundary.size(), false);
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingLeftBoundary.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingLeftBoundary[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByTimeVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->DerivativesOfLeftBoundaryByTimeVariable[varIndex][entryIndex];
                        size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);

                        //hooray, this variable affects the constraint. Make the Jacobian entry!
                        if (this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][stateIndex])
                        {
                            constraint_LeftBoundaryTimeVariableAffectsMatchConstraint[varIndex] = true;
                            break;
                        }
                    }

                    if (constraint_LeftBoundaryTimeVariableAffectsMatchConstraint[varIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints[constraintIndex],
                            Xindex,
                            this->G_indices_match_point_constraints_wrt_LeftBoundaryTimeVariables[constraintIndex]);
                    }
                    else
                        this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables[constraintIndex].push_back(0);
                }
                this->LeftBoundaryTimeVariableAffectsMatchConstraint.push_back(constraint_LeftBoundaryTimeVariableAffectsMatchConstraint);

                //with respect to right boundary
                std::vector<bool> constraint_RightBoundaryTimeVariableAffectsMatchConstraint(this->ListOfTimeVariablesAffectingRightBoundary.size(), false);
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingRightBoundary.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingRightBoundary[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByTimeVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->DerivativesOfRightBoundaryByTimeVariable[varIndex][entryIndex];
                        size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);

                        //hooray, this variable affects the constraint. Make the Jacobian entry!
                        if (this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][stateIndex])
                        {
                            constraint_RightBoundaryTimeVariableAffectsMatchConstraint[varIndex] = true;
                            break;
                        }
                    }

                    if (constraint_RightBoundaryTimeVariableAffectsMatchConstraint[varIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints[constraintIndex],
                            Xindex,
                            this->G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables[constraintIndex]);
                    }
                    else
                        this->G_indices_match_point_constraints_wrt_RightBoundaryVariables[constraintIndex].push_back(0);
                }
                this->RightBoundaryTimeVariableAffectsMatchConstraint.push_back(constraint_RightBoundaryTimeVariableAffectsMatchConstraint);
            }
        }//end calcbounds_match_point_constraints()

        //******************************************process methods
        void TwoPointShootingPhase::process_phase_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_forward_half_phase(X, Xindex, F, Findex, G, needG);
            
            this->process_backward_half_phase(X, Xindex, F, Findex, G, needG);
        }

        void TwoPointShootingPhase::process_match_point_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: set the match point constraints

            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
            {
                size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                F[Findex++] = (this->match_point_state_plus(stateIndex) - this->match_point_state_minus(stateIndex))
                    * this->continuity_constraint_scale_factors(constraintIndex);
            }

            //Step 2: derivatives of the match point constraints with respect to the boundary variables
            if (needG)
            {
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterPhase = this->myArrivalEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase_wrt_Time = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterPhase_wrt_Time = this->myArrivalEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value
                math::Matrix<double> dBoundaryState_dDecisionVariable(this->ForwardHPTM.get_n(), 1, 0.0);
                math::Matrix<double> dMatchPointState_dDecisionVariable(this->ForwardHPTM.get_n(), 1, 0.0);

                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {
                    size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                    //Step 2.1: left boundary
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfLeftBoundaryByVariable.size(); ++varIndex)
                    {
                        if (this->LeftBoundaryVariableAffectsMatchConstraint[constraintIndex][varIndex])
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByVariable[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfLeftBoundaryByVariable[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateBeforePhase[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->ForwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables[constraintIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Step 2.2: right boundary
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfRightBoundaryByVariable.size(); ++varIndex)
                    {
                        if (this->RightBoundaryVariableAffectsMatchConstraint[constraintIndex][varIndex])
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByVariable[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfRightBoundaryByVariable[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateAfterPhase[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->BackwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_RightBoundaryVariables[constraintIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Step 2.3: left boundary epoch
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfLeftBoundaryByTimeVariable.size(); ++varIndex)
                    {
                        if (this->LeftBoundaryTimeVariableAffectsMatchConstraint[constraintIndex][varIndex])
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByTimeVariable[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfLeftBoundaryByTimeVariable[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->ForwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_LeftBoundaryTimeVariables[constraintIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Step 2.4: right boundary epoch
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfRightBoundaryByTimeVariable.size() - 1; ++varIndex)
                    {
                        if (this->RightBoundaryTimeVariableAffectsMatchConstraint[constraintIndex][varIndex])
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByTimeVariable[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfRightBoundaryByTimeVariable[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->BackwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables[constraintIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Step 2.5: right boundary epoch due to current phase flight time
                    {
                        if (this->RightBoundaryTimeVariableAffectsMatchConstraint[constraintIndex].back())
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByTimeVariable.back().size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfRightBoundaryByTimeVariable.back()[entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->BackwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables[constraintIndex].back();
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }//end loop over match point constraints
            }//end match point constraint derivatives
        }//end process_match_point_constraints
    }//close namespace Phases
}//close namespace EMTG