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

#include "ParallelShootingFirstStep.h"

namespace EMTG
{
    namespace Phases
    {
        ParallelShootingFirstStep::ParallelShootingFirstStep() {};

        ParallelShootingFirstStep::ParallelShootingFirstStep(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stepIndex,
            const size_t& stageIndex,
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
                NULL, //first step doesn't have a previous step!
                myPhase,
                myUniverse,
                mySpacecraft,
                myOptions);
        }

        void ParallelShootingFirstStep::calcbounds_step_left_match_point_constraints()
        {
            //the left match point constraint is defined as current_left_state - previous_right_state
            //this is the specialized version for the first step, in which the "previous step" is the state after the phase's initial coast

            //Step 1: derivatives with respect to previous step right state
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_PreviousStepRightState = this->myPhase->get_Derivatives_of_state_after_initial_coast();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_PreviousStepRightState_wrt_Time = this->myPhase->get_Derivatives_of_state_after_initial_coast_wrt_Time();//Xindex, stateIndex, derivative value

            
            //Step 1.1: non-time                                                                                                                                                                        ////Step 1.1: state variables
            {
                for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_PreviousStepRightState.size(); ++dIndex) //loop over derivative entries for this state
                    {
                        size_t previousStateIndex = std::get<1>(Derivatives_of_PreviousStepRightState[dIndex]);

                        if (stateIndex == previousStateIndex)
                        {
                            size_t Xindex = std::get<0>(Derivatives_of_PreviousStepRightState[dIndex]);

                            //if this variable is not yet in the list of variables that affect the previous step's right state, add it to the list
                            if (std::find(this->ListOfVariablesAffectingPreviousStepRightState.begin(), this->ListOfVariablesAffectingPreviousStepRightState.end(), Xindex)
                                == this->ListOfVariablesAffectingPreviousStepRightState.end())
                                this->ListOfVariablesAffectingPreviousStepRightState.push_back(Xindex);
                        }
                    }
                }

                this->DerivativesOfPreviousStepRightStateByVariable.resize(ListOfVariablesAffectingPreviousStepRightState.size());

                for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingPreviousStepRightState.size(); ++listIndex)
                {
                    for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                    {
                        for (size_t dIndex = 0; dIndex < Derivatives_of_PreviousStepRightState.size(); ++dIndex) //loop over derivative entries for this state
                        {
                            size_t previousStateIndex = std::get<1>(Derivatives_of_PreviousStepRightState[dIndex]);

                            if (stateIndex == previousStateIndex)
                            {
                                size_t Xindex = std::get<0>(Derivatives_of_PreviousStepRightState[dIndex]);

                                if (Xindex == this->ListOfVariablesAffectingPreviousStepRightState[listIndex])
                                {
                                    this->DerivativesOfPreviousStepRightStateByVariable[listIndex].push_back({ stateIndex, dIndex });
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
                    for (size_t dIndex = 0; dIndex < Derivatives_of_PreviousStepRightState_wrt_Time.size(); ++dIndex) //loop over derivative entries for this state
                    {
                        size_t previousStateIndex = std::get<1>(Derivatives_of_PreviousStepRightState_wrt_Time[dIndex]);

                        if (stateIndex == previousStateIndex)
                        {
                            size_t Xindex = std::get<0>(Derivatives_of_PreviousStepRightState_wrt_Time[dIndex]);

                            //if this variable is not yet in the list of variables that affect the previous step's right state, add it to the list
                            if (std::find(this->ListOfTimeVariablesAffectingPreviousStepRightState.begin(), this->ListOfTimeVariablesAffectingPreviousStepRightState.end(), Xindex)
                                == this->ListOfTimeVariablesAffectingPreviousStepRightState.end())
                                this->ListOfTimeVariablesAffectingPreviousStepRightState.push_back(Xindex);
                        }
                    }
                }

                this->DerivativesOfPreviousStepRightStateByTimeVariable.resize(ListOfTimeVariablesAffectingPreviousStepRightState.size());

                for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingPreviousStepRightState.size(); ++listIndex)
                {
                    for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                    {
                        for (size_t dIndex = 0; dIndex < Derivatives_of_PreviousStepRightState_wrt_Time.size(); ++dIndex) //loop over derivative entries for this state
                        {
                            size_t previousStateIndex = std::get<1>(Derivatives_of_PreviousStepRightState_wrt_Time[dIndex]);

                            if (stateIndex == previousStateIndex)
                            {
                                size_t Xindex = std::get<0>(Derivatives_of_PreviousStepRightState_wrt_Time[dIndex]);

                                if (Xindex == this->ListOfTimeVariablesAffectingPreviousStepRightState[listIndex])
                                {
                                    this->DerivativesOfPreviousStepRightStateByTimeVariable[listIndex].push_back({ stateIndex, dIndex });
                                }
                            }
                        }
                    }
                }
            }//end time variables

             //Step 2: derivatives with respect to current step left state
             //Step 2.1: non-time variables
            {
                for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                {
                    for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateStepLeftInertial.size(); ++dIndex) //loop over derivative entries for this state
                    {
                        size_t previousStateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        if (stateIndex == previousStateIndex)
                        {
                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                            //if this variable is not yet in the list of variables that affect the current step's left state, add it to the list
                            if (std::find(this->ListOfVariablesAffectingCurrentStepLeftState.begin(), this->ListOfVariablesAffectingCurrentStepLeftState.end(), Xindex)
                                == this->ListOfVariablesAffectingCurrentStepLeftState.end())
                            {
                                this->ListOfVariablesAffectingCurrentStepLeftState.push_back(Xindex);
                                if (stateIndex < 3)
                                    this->ListOfVariablesAffectingCurrentStepLeftPosition.push_back(Xindex);
                                else if (stateIndex < 6)
                                    this->ListOfVariablesAffectingCurrentStepLeftVelocity.push_back(Xindex);
                            }
                        }
                    }
                }

                this->DerivativesOfCurrentStepLeftStateByVariable.resize(this->ListOfVariablesAffectingCurrentStepLeftState.size());

                for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingCurrentStepLeftState.size(); ++listIndex)
                {
                    for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                    {
                        for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateStepLeftInertial.size(); ++dIndex) //loop over derivative entries for this state
                        {
                            size_t previousStateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                            if (stateIndex == previousStateIndex)
                            {
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                if (Xindex == this->ListOfVariablesAffectingCurrentStepLeftState[listIndex])
                                {
                                    this->DerivativesOfCurrentStepLeftStateByVariable[listIndex].push_back({ stateIndex, dIndex });
                                }
                            }
                        }
                    }
                }
            }//end non-time

             //Step 2.2: time variables
            {
                for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                {
                    for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateStepLeftInertial_wrt_Time.size(); ++dIndex) //loop over derivative entries for this state
                    {
                        size_t previousStateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial_wrt_Time[dIndex]);

                        if (stateIndex == previousStateIndex)
                        {
                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial_wrt_Time[dIndex]);

                            //if this variable is not yet in the list of variables that affect the previous step's right state, add it to the list
                            if (std::find(this->ListOfTimeVariablesAffectingCurrentStepLeftState.begin(), this->ListOfTimeVariablesAffectingCurrentStepLeftState.end(), Xindex)
                                == this->ListOfTimeVariablesAffectingCurrentStepLeftState.end())
                                this->ListOfTimeVariablesAffectingCurrentStepLeftState.push_back(Xindex);
                        }
                    }
                }

                this->DerivativesOfCurrentStepLeftStateByTimeVariable.resize(this->ListOfTimeVariablesAffectingCurrentStepLeftState.size());

                for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingCurrentStepLeftState.size(); ++listIndex)
                {
                    for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                    {
                        for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateStepLeftInertial_wrt_Time.size(); ++dIndex) //loop over derivative entries for this state
                        {
                            size_t previousStateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial_wrt_Time[dIndex]);

                            if (stateIndex == previousStateIndex)
                            {
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial_wrt_Time[dIndex]);

                                if (Xindex == this->ListOfTimeVariablesAffectingCurrentStepLeftState[listIndex])
                                {
                                    this->DerivativesOfCurrentStepLeftStateByTimeVariable[listIndex].push_back({ stateIndex, dIndex });
                                }
                            }
                        }
                    }
                }
            }//end time

            //Step 3: construct the match point constraints
            std::vector<std::string>& matchPointConstraintNames = this->myPhase->get_matchPointConstraintNames();
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
            {
                Flowerbounds->push_back(-math::SMALL);
                Fupperbounds->push_back(math::SMALL);
                Fdescriptions->push_back(prefix + "left match point " + this->myPhase->get_matchPointConstraintNames()[constraintIndex]);
                this->Findices_left_match_point_constraints.push_back(Fdescriptions->size() - 1);
            }//end construction of match point constraints

            //Step 4: derivatives with respect to previous step right boundary
            {
                //Step 4.1: non-time
                this->Gindices_LeftMatchPoint_wrt_PreviousStepRightState.resize(this->ListOfVariablesAffectingPreviousStepRightState.size());
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingPreviousStepRightState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingPreviousStepRightState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPreviousStepRightStateByVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPreviousStepRightStateByVariable[varIndex][entryIndex]);

                        size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[constraintIndex],
                            Xindex,
                            this->Gindices_LeftMatchPoint_wrt_PreviousStepRightState[varIndex]);
                    }
                }

                //Step 4.2: time
                this->Gindices_LeftMatchPoint_wrt_PreviousStepRightStateTime.resize(this->ListOfTimeVariablesAffectingPreviousStepRightState.size());
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingPreviousStepRightState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingPreviousStepRightState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex][entryIndex]);

                        size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[constraintIndex],
                            Xindex,
                            this->Gindices_LeftMatchPoint_wrt_PreviousStepRightStateTime[varIndex]);
                    }
                }
            }//end derivatives with respect to previous step right boundary

            //Step 5: derivatives with respect to current step left state
            for (size_t encodedStateIndex : {0, 1, 2, 3, 4, 5})
            {
                for (size_t cartesianStateIndex : {0, 1, 2, 3, 4, 5})
                {
                    if (this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[cartesianStateIndex][encodedStateIndex])
                    {
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_StateElements[encodedStateIndex][cartesianStateIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[cartesianStateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_StepLeftEncodedState);
                    }
                }
            }

            //mass
            this->create_sparsity_entry(this->Findices_left_match_point_constraints[6],
                this->Xindex_mass,
                this->Gindex_StepLeftMatchPoint_wrt_Mass);
            //virtual chemical fuel
            this->create_sparsity_entry(this->Findices_left_match_point_constraints[7],
                this->Xindex_chemical_fuel,
                this->Gindex_StepLeftMatchPoint_wrt_ChemicalFuel);
            //virtual electric propellant
            this->create_sparsity_entry(this->Findices_left_match_point_constraints[8],
                this->Xindex_electric_propellant,
                this->Gindex_StepLeftMatchPoint_wrt_ElectricPropellant);
        }//end calcbounds_step_left_match_point_constraints

        void ParallelShootingFirstStep::process_step_left_match_point_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //specialized version in which "previous step" is actually the phase's state after initial coast
            //Step 1: compute the match point constraints
            std::vector<size_t>& matchPointConstraintStateIndex = this->myPhase->get_matchPointConstraintStateIndex();
            math::Matrix<double>& continuity_constraint_scale_factors = this->myPhase->get_continuity_constraint_scale_factors();
            math::Matrix<doubleType>& PreviousStepRightInertial = this->myPhase->get_StateAfterInitialCoast();
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
            {
                size_t stateIndex = matchPointConstraintStateIndex[constraintIndex];
                F[Findex++] = (this->StateStepLeftInertial(stateIndex) - PreviousStepRightInertial(stateIndex))
                    * continuity_constraint_scale_factors(constraintIndex);
            }

            //Step 2: derivatives of the match point constraints
            if (needG)
            {
                math::Matrix<double> dMatchState_dDecisionVariable(this->StateStepLeftInertial.get_n(), 1, 0.0);

                //Step 2.1: with respect to the current step//first clear the entries
                for (size_t& Gindex : this->Gindices_StepLeftMatchPoint_wrt_StepLeftEncodedState)
                {
                    G[Gindex] = 0.0;
                }

                //now populate them
                for (size_t encodedStateIndex : {0, 1, 2, 3, 4, 5})
                {
                    for (size_t cartesianStateIndex : {0, 1, 2, 3, 4, 5})
                    {

                        if (this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[cartesianStateIndex][encodedStateIndex])
                        {
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_StateElements[encodedStateIndex][cartesianStateIndex];

                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_StepLeftEncodedState[dIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(cartesianStateIndex);
                        }
                    }
                }

                //mass
                {
                    size_t Gindex = this->Gindex_StepLeftMatchPoint_wrt_Mass;
                    size_t Xindex = this->Xindex_mass;

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * continuity_constraint_scale_factors(6);
                }
                //virtual chemical fuel
                {
                    size_t Gindex = this->Gindex_StepLeftMatchPoint_wrt_ChemicalFuel;
                    size_t Xindex = this->Xindex_chemical_fuel;

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * continuity_constraint_scale_factors(7);
                }
                //virtual electric propellant
                {
                    size_t Gindex = this->Gindex_StepLeftMatchPoint_wrt_ElectricPropellant;
                    size_t Xindex = this->Xindex_electric_propellant;

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * continuity_constraint_scale_factors(8);
                }


                //Step 2.2: with respect to the previous step
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_PreviousStepRightState = this->myPhase->get_Derivatives_of_state_after_initial_coast();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_PreviousStepRightState_wrt_Time = this->myPhase->get_Derivatives_of_state_after_initial_coast_wrt_Time();//Xindex, stateIndex, derivative value

                                                                                                                                                                                           //Step 4.1: non-time
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingPreviousStepRightState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingPreviousStepRightState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPreviousStepRightStateByVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPreviousStepRightStateByVariable[varIndex][entryIndex]);

                        if (stateIndex != 7) //there is no constraint on epoch
                        {
                            size_t dIndex = std::get<1>(this->DerivativesOfPreviousStepRightStateByVariable[varIndex][entryIndex]);

                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            double TheDerivative = std::get<2>(Derivatives_of_PreviousStepRightState[dIndex]);

                            size_t Gindex = this->Gindices_LeftMatchPoint_wrt_PreviousStepRightState[varIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -TheDerivative
                                * continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }

                //Step 4.2: time
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingPreviousStepRightState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingPreviousStepRightState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex][entryIndex]);

                        if (stateIndex != 7) //there is no constraint on epoch
                        {
                            size_t dIndex = std::get<1>(this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex][entryIndex]);

                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            double TheDerivative = std::get<2>(Derivatives_of_PreviousStepRightState_wrt_Time[dIndex]);

                            size_t Gindex = this->Gindices_LeftMatchPoint_wrt_PreviousStepRightStateTime[varIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -TheDerivative
                                * continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }
            }//end match point constraint derivatives
        }//end process_step_left_match_point_constraints
    }//close namespace Phases
}//close namespace EMTG