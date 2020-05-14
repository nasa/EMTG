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

//base parallel shooting "first step" for EMTGv9
//Jacob Englander 2-21-2018

#include "ParallelShootingLastStep.h"

namespace EMTG
{
    namespace Phases
    {
        ParallelShootingLastStep::ParallelShootingLastStep() 
        {
            this->Gindices_StepRightMatchPoint_wrt_CurrentStepRightStateVariables.resize(9);//one for each constraint
            this->Gindices_StepRightMatchPoint_wrt_CurrentStepRightStateTimeVariables.resize(9);//one for each constraint
        }

        ParallelShootingLastStep::ParallelShootingLastStep(const std::string& name,
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
            this->ParallelShootingStep::initialize(name,
                journeyIndex,
                phaseIndex,
                stepIndex,
                stageIndex,
                previousStep,
                myPhase,
                myUniverse,
                mySpacecraft,
                myOptions);
        }

        void ParallelShootingLastStep::calcbounds_step_right_match_point_constraints()
        {
            //the left match point constraint is defined as state_before_terminal_coast - current_left_state

            //Step 1: derivatives with respect to state after terminal coast
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeTerminalCoast = this->myPhase->get_Derivatives_of_state_before_terminal_coast();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeTerminalCoast_wrt_Time = this->myPhase->get_Derivatives_of_state_before_terminal_coast_wrt_Time();//Xindex, stateIndex, derivative value

                                                                                                                                                                                           //Step 1.1: non-time                                                                                                                                                                        ////Step 1.1: state variables
            {
                for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforeTerminalCoast.size(); ++dIndex) //loop over derivative entries for this state
                    {
                        size_t previousStateIndex = std::get<1>(Derivatives_of_StateBeforeTerminalCoast[dIndex]);

                        if (stateIndex == previousStateIndex)
                        {
                            size_t Xindex = std::get<0>(Derivatives_of_StateBeforeTerminalCoast[dIndex]);

                            //if this variable is not yet in the list of variables that affect the previous step's right state, add it to the list
                            if (std::find(this->ListOfVariablesAffectingPhaseTerminalCoastLeftState.begin(), this->ListOfVariablesAffectingPhaseTerminalCoastLeftState.end(), Xindex)
                                == this->ListOfVariablesAffectingPhaseTerminalCoastLeftState.end())
                                this->ListOfVariablesAffectingPhaseTerminalCoastLeftState.push_back(Xindex);
                        }
                    }
                }

                this->DerivativesOfPhaseTerminalCoastLeftStateByVariable.resize(ListOfVariablesAffectingPhaseTerminalCoastLeftState.size());

                for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingPhaseTerminalCoastLeftState.size(); ++listIndex)
                {
                    for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                    {
                        for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforeTerminalCoast.size(); ++dIndex) //loop over derivative entries for this state
                        {
                            size_t previousStateIndex = std::get<1>(Derivatives_of_StateBeforeTerminalCoast[dIndex]);

                            if (stateIndex == previousStateIndex)
                            {
                                size_t Xindex = std::get<0>(Derivatives_of_StateBeforeTerminalCoast[dIndex]);

                                if (Xindex == this->ListOfVariablesAffectingPhaseTerminalCoastLeftState[listIndex])
                                {
                                    this->DerivativesOfPhaseTerminalCoastLeftStateByVariable[listIndex].push_back({ stateIndex, dIndex });
                                }
                            }
                        }
                    }
                }
            }//end non-time

             //Step 1.2: time variables
            {
                for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforeTerminalCoast_wrt_Time.size(); ++dIndex) //loop over derivative entries for this state
                    {
                        size_t previousStateIndex = std::get<1>(Derivatives_of_StateBeforeTerminalCoast_wrt_Time[dIndex]);

                        if (stateIndex == previousStateIndex)
                        {
                            size_t Xindex = std::get<0>(Derivatives_of_StateBeforeTerminalCoast_wrt_Time[dIndex]);

                            //if this variable is not yet in the list of variables that affect the previous step's right state, add it to the list
                            if (std::find(this->ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState.begin(), this->ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState.end(), Xindex)
                                == this->ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState.end())
                                this->ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState.push_back(Xindex);
                        }
                    }
                }

                this->DerivativesOfPhaseTerminalCoastLeftStateByTimeVariable.resize(ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState.size());

                for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState.size(); ++listIndex)
                {
                    for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                    {
                        for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforeTerminalCoast_wrt_Time.size(); ++dIndex) //loop over derivative entries for this state
                        {
                            size_t previousStateIndex = std::get<1>(Derivatives_of_StateBeforeTerminalCoast_wrt_Time[dIndex]);

                            if (stateIndex == previousStateIndex)
                            {
                                size_t Xindex = std::get<0>(Derivatives_of_StateBeforeTerminalCoast_wrt_Time[dIndex]);

                                if (Xindex == this->ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState[listIndex])
                                {
                                    this->DerivativesOfPhaseTerminalCoastLeftStateByTimeVariable[listIndex].push_back({ stateIndex, dIndex });
                                }
                            }
                        }
                    }
                }
            }//end time variables
            
            //Step 2: (current step variables already done in calcbounds_step_main())

            //Step 3: construct the match point constraints
            std::vector<std::string>& matchPointConstraintNames = this->myPhase->get_matchPointConstraintNames();
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
            {
                this->Flowerbounds->push_back(-math::SMALL);
                this->Fupperbounds->push_back(math::SMALL);
                this->Fdescriptions->push_back(prefix + "right match point " + matchPointConstraintNames[constraintIndex]);
                this->Findices_right_match_point_constraints.push_back(Fdescriptions->size() - 1);
            }

            //Step 4: derivatives with respect to the terminal coast state
            {
                //Step 4.1: non-time
                this->Gindices_RightMatchPoint_wrt_PhaseTerminalCoastLeftState.resize(this->ListOfVariablesAffectingPhaseTerminalCoastLeftState.size());
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingPhaseTerminalCoastLeftState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingPhaseTerminalCoastLeftState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPhaseTerminalCoastLeftStateByVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPhaseTerminalCoastLeftStateByVariable[varIndex][entryIndex]);

                        size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                        this->create_sparsity_entry(this->Findices_right_match_point_constraints[constraintIndex],
                            Xindex,
                            this->Gindices_RightMatchPoint_wrt_PhaseTerminalCoastLeftState[varIndex]);
                    }
                }

                //Step 4.2: time
                this->Gindices_RightMatchPoint_wrt_PhaseTerminalCoastLeftStateTime.resize(this->ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState.size());
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPhaseTerminalCoastLeftStateByTimeVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPhaseTerminalCoastLeftStateByTimeVariable[varIndex][entryIndex]);

                        size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                        this->create_sparsity_entry(this->Findices_right_match_point_constraints[constraintIndex],
                            Xindex,
                            this->Gindices_RightMatchPoint_wrt_PhaseTerminalCoastLeftStateTime[varIndex]);
                    }
                }
            }//end derivatives with respect to phase terminal coast

            //Step 5: derivatives with respect to current step right boundary
            {
                //Step 5.1: non-time
                this->Gindices_StepRightMatchPoint_wrt_CurrentStepRightStateVariables.resize(this->ListOfVariablesAffectingCurrentStepLeftState.size());
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingCurrentStepLeftState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingCurrentStepLeftState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfCurrentStepRightStateByVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfCurrentStepRightStateByVariable[varIndex][entryIndex]);

                        size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                        this->create_sparsity_entry(this->Findices_right_match_point_constraints[constraintIndex],
                            Xindex,
                            this->Gindices_StepRightMatchPoint_wrt_CurrentStepRightStateVariables[varIndex]);
                    }
                }

                //Step 5.2: time
                this->Gindices_StepRightMatchPoint_wrt_CurrentStepRightStateTimeVariables.resize(this->ListOfTimeVariablesAffectingCurrentStepLeftState.size());
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingCurrentStepLeftState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingCurrentStepLeftState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfCurrentStepRightStateByTimeVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfCurrentStepRightStateByTimeVariable[varIndex][entryIndex]);

                        size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                        this->create_sparsity_entry(this->Findices_right_match_point_constraints[constraintIndex],
                            Xindex,
                            this->Gindices_StepRightMatchPoint_wrt_CurrentStepRightStateTimeVariables[varIndex]);
                    }
                }
            }//end derivatives with respect to previous step right boundary

            //Step 6: derivatives with respect to control
            this->Gindices_StepRightMatchPoint_wrt_CurrentStepControlVariables.resize(this->num_interior_control_points, std::vector< std::vector<size_t> >(3));
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                for (size_t stateIndex : {0, 1, 2, 3, 4, 5, 6, 9})//no affect on epoch or chemical fuel
                {
                    for (size_t controlIndex = 0; controlIndex < this->dIndex_right_state_wrt_control[subStepIndex][stateIndex].size(); ++controlIndex)
                    {
                        size_t dIndex = this->dIndex_right_state_wrt_control[subStepIndex][stateIndex][controlIndex];
                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepRightInertial[dIndex]);

                        size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                        this->create_sparsity_entry(this->Findices_right_match_point_constraints[constraintIndex],
                            Xindex,
                            this->Gindices_StepRightMatchPoint_wrt_CurrentStepControlVariables[subStepIndex][controlIndex]);
                    }
                }
            }
        }//end calcbounds_step_right_match_point_constraints

        void ParallelShootingLastStep::process_step_right_match_point_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: compute the match point constraints
            std::vector<size_t>& matchPointConstraintStateIndex = this->myPhase->get_matchPointConstraintStateIndex();
            math::Matrix<double>& continuity_constraint_scale_factors = this->myPhase->get_continuity_constraint_scale_factors();
            math::Matrix<doubleType>& StateBeforeTerminalCoast = this->myPhase->get_StateBeforeTerminalCoast();
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
            {
                size_t stateIndex = matchPointConstraintStateIndex[constraintIndex];
                F[Findex++] = (StateBeforeTerminalCoast(stateIndex) - this->StateStepRightInertial(stateIndex))
                    * continuity_constraint_scale_factors(constraintIndex);
            }

            //Step 2: derivatives of the match point constraints
            if (needG)
            {
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeTerminalCoast = this->myPhase->get_Derivatives_of_state_before_terminal_coast();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeTerminalCoast_wrt_Time = this->myPhase->get_Derivatives_of_state_before_terminal_coast_wrt_Time();//Xindex, stateIndex, derivative value

                math::Matrix<double> dMatchState_dDecisionVariable(this->StateStepRightInertial.get_n(), 1, 0.0);

                //Step 2.1: non-time terminal coast
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingPhaseTerminalCoastLeftState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingPhaseTerminalCoastLeftState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPhaseTerminalCoastLeftStateByVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPhaseTerminalCoastLeftStateByVariable[varIndex][entryIndex]);

                        if (stateIndex != 7) //there is no constraint on epoch
                        {
                            size_t dIndex = std::get<1>(this->DerivativesOfPhaseTerminalCoastLeftStateByVariable[varIndex][entryIndex]);

                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            double TheDerivative = std::get<2>(Derivatives_of_StateBeforeTerminalCoast[dIndex]);

                            size_t Gindex = this->Gindices_RightMatchPoint_wrt_PhaseTerminalCoastLeftState[varIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * TheDerivative
                                * continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }

                //Step 2.2: time terminal coast
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPhaseTerminalCoastLeftStateByTimeVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPhaseTerminalCoastLeftStateByTimeVariable[varIndex][entryIndex]);

                        if (stateIndex != 7) //there is no constraint on epoch
                        {
                            size_t dIndex = std::get<1>(this->DerivativesOfPhaseTerminalCoastLeftStateByTimeVariable[varIndex][entryIndex]);

                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            double TheDerivative = std::get<2>(Derivatives_of_StateBeforeTerminalCoast_wrt_Time[dIndex]);

                            size_t Gindex = this->Gindices_RightMatchPoint_wrt_PhaseTerminalCoastLeftStateTime[varIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * TheDerivative
                                * continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }
                
                //Step 2.3: non-time current step right side
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingCurrentStepLeftState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingCurrentStepLeftState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfCurrentStepRightStateByVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfCurrentStepRightStateByVariable[varIndex][entryIndex]);

                        if (stateIndex != 7) //there is no constraint on epoch
                        {
                            size_t dIndex = std::get<1>(this->DerivativesOfCurrentStepRightStateByVariable[varIndex][entryIndex]);

                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            double TheDerivative = std::get<2>(this->Derivatives_of_StateStepRightInertial[dIndex]);

                            size_t Gindex = this->Gindices_StepRightMatchPoint_wrt_CurrentStepRightStateVariables[varIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -TheDerivative
                                * continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }

                //Step 2.4: time current step right side
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingCurrentStepLeftState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingCurrentStepLeftState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfCurrentStepRightStateByTimeVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfCurrentStepRightStateByTimeVariable[varIndex][entryIndex]);

                        if (stateIndex != 7) //there is no constraint on epoch
                        {
                            size_t dIndex = std::get<1>(this->DerivativesOfCurrentStepRightStateByTimeVariable[varIndex][entryIndex]);

                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            double TheDerivative = std::get<2>(Derivatives_of_StateStepRightInertial_wrt_Time[dIndex]);

                            size_t Gindex = this->Gindices_StepRightMatchPoint_wrt_CurrentStepRightStateTimeVariables[varIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            
                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * -TheDerivative
                                * continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }


                //Step 2.5: derivatives with respect to control
                const size_t stateIndices[] = { 0, 1, 2, 3, 4, 5, 6, 9 };//no affect on epoch or chemical fuel
                for (size_t entryIndex = 0; entryIndex < 8; ++entryIndex)
                {
                    size_t stateIndex = stateIndices[entryIndex];
                    for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
                    {
                        for (size_t controlIndex = 0; controlIndex < this->dIndex_right_state_wrt_control[subStepIndex][stateIndex].size(); ++controlIndex)
                        {
                            size_t dIndex = this->dIndex_right_state_wrt_control[subStepIndex][stateIndex][controlIndex];
                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepRightInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepRightMatchPoint_wrt_CurrentStepControlVariables[subStepIndex][controlIndex][entryIndex];

                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            double TheDerivative = std::get<2>(this->Derivatives_of_StateStepRightInertial[dIndex]);


                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -TheDerivative
                                * continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }
                
            }//end match point constraint derivatives
        }//end process_step_right_match_point_constraints
    }//end namespace Phases
}//end namespace EMTG