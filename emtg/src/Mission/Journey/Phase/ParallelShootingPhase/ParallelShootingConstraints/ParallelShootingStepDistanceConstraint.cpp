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

//phase distance constraint factory
//for TwoPointShootingLowThrustPhase
//Jacob Englander 1/18/2018

#include "ParallelShootingStepDistanceConstraint.h"

#include <vector>

#include "boost/algorithm/string.hpp"

namespace EMTG
{
    namespace Phases
    {
        ParallelShootingStepDistanceConstraint::ParallelShootingStepDistanceConstraint() :
            distance_constraint_relative_position(3, 1, 0.0),
            distance_constraint_body_position_time_derivatives(3, 1, 0.0),
            distance_from_body(0.0)
        {}

        ParallelShootingStepDistanceConstraint::ParallelShootingStepDistanceConstraint(const std::string& ConstraintDefinition,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stepIndex,
            ParallelShootingStep* myStep,
            missionoptions* myOptions,
            Astrodynamics::universe* myUniverse) :
            ParallelShootingStepDistanceConstraint()
        {
            this->myStep = myStep;
            this->journeyIndex = journeyIndex;
            this->phaseIndex = phaseIndex;
            this->stepIndex = stepIndex;
            this->myOptions = myOptions;
            this->myUniverse = myUniverse;
            this->myJourneyOptions = &this->myOptions->Journeys[this->journeyIndex];

            std::vector<std::string> ConstraintDefinitionCell;
            boost::split(ConstraintDefinitionCell, ConstraintDefinition, boost::is_any_of("_"), boost::token_compress_on);

            if (boost::to_lower_copy(ConstraintDefinitionCell[1]) == "cb")
            {
                this->isDefinedRelativeToCentralBody = true;
                this->myBody = &this->myUniverse->central_body;
            }
            else
            {
                this->isDefinedRelativeToCentralBody = false;
                int bodyIndex = std::stoi(ConstraintDefinitionCell[1]);
                this->myBody = &this->myUniverse->bodies[bodyIndex - 1];
            }

            this->name = this->myStep->getName() + " distance from " + this->myBody->name;

            this->parse_constraint_definition(ConstraintDefinition);
        }//end constructor

        void ParallelShootingStepDistanceConstraint::parse_constraint_definition(const std::string& ConstraintDefinition)
        {
            std::vector<std::string> ConstraintDefinitionCell;
            boost::split(ConstraintDefinitionCell, ConstraintDefinition, boost::is_any_of("_"), boost::token_compress_on);

            int bodyIndex;
            if (boost::to_lower_copy(ConstraintDefinitionCell[1]) == "cb")
                bodyIndex = -2;
            else
                bodyIndex = std::stoi(ConstraintDefinitionCell[1]);
            
            if (ConstraintDefinitionCell[2].find("km") < 1024)
            {
                this->lowerBound = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[2], "km"));
            }
            else if (ConstraintDefinitionCell[2].find("au") < 1024)
            {
                this->lowerBound = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[2], "au")) * this->myOptions->AU;
            }
            else if (ConstraintDefinitionCell[2].find("lu") < 1024)
            {
                this->lowerBound = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[2], "lu")) * this->myUniverse->LU;
            }
            else
            {
                this->lowerBound = std::stod(ConstraintDefinitionCell[2]);
            }

            if (ConstraintDefinitionCell[3].find("km") < 1024)
            {
                this->upperBound = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[3], "km"));
            }
            else if (ConstraintDefinitionCell[3].find("au") < 1024)
            {
                this->upperBound = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[3], "au")) * this->myOptions->AU;
            }
            else if (ConstraintDefinitionCell[3].find("lu") < 1024)
            {
                this->upperBound = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[3], "lu")) * this->myUniverse->LU;
            }
            else
            {
                this->upperBound = std::stod(ConstraintDefinitionCell[3]);
            }

            this->distance_constraint_definition = std::tuple<int, double, double>({ bodyIndex, lowerBound, upperBound });
        }//end parse_constraint_definition()

        void ParallelShootingStepDistanceConstraint::calcbounds()
        {
            //Step 1: make the constraint
            Flowerbounds->push_back((this->lowerBound - this->upperBound) / this->myUniverse->LU);
            Fupperbounds->push_back(0.0);
            Fdescriptions->push_back(this->name);

            //Step 2: sparsity pattern
            if (this->isDefinedRelativeToCentralBody
                && (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalRADEC
                    || this->myOptions->ParallelShootingConstraintStateRepresentation == StateRepresentation::SphericalAZFPA))
            {
                //lucky us, there's only one entry! and no time dependence!
                size_t Xindex = this->myStep->getXindex_state_elements()[0];

                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    Xindex,
                    this->G_indices_distance_constraints_wrt_StepLeftPosition);
            }
            else
            {
                //Step 2.1: non-time variables
                std::vector<size_t>& ListOfVariablesAffectingCurrentStepLeftPosition = this->myStep->getListOfVariablesAffectingCurrentStepLeftPosition();

                std::vector< std::vector< std::tuple<size_t, size_t> > >& DerivativesOfCurrentStepLeftStateByVariable = this->myStep->getDerivativesOfCurrentStepLeftStateByVariable();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StepLeftState = this->myStep->get_Derivatives_of_StateStepLeftInertial();//Xindex, stateIndex, derivative value

                for (size_t varIndex = 0; varIndex < ListOfVariablesAffectingCurrentStepLeftPosition.size(); ++varIndex)
                {
                    size_t Xindex = ListOfVariablesAffectingCurrentStepLeftPosition[varIndex];
                    std::vector<size_t> dIndex_distance_constraints_wrt_StepLeftPositionVariable;

                    std::vector< std::tuple<size_t, size_t> > DerivativesOfCurrentStepLeftStateThisVariable = DerivativesOfCurrentStepLeftStateByVariable[varIndex];

                    bool madeEntryFlag = false;

                    for (std::tuple<size_t, size_t> thisDerivativeEntry : DerivativesOfCurrentStepLeftStateThisVariable)
                    {
                        size_t stateIndex = std::get<0>(thisDerivativeEntry);


                        if (stateIndex < 3)
                        {
                            size_t dIndex = std::get<1>(thisDerivativeEntry);
                            dIndex_distance_constraints_wrt_StepLeftPositionVariable.push_back(dIndex);


                            if (!madeEntryFlag)
                            {
                                this->create_sparsity_entry(Fdescriptions->size() - 1,
                                    Xindex,
                                    this->G_indices_distance_constraints_wrt_StepLeftPosition);

                                madeEntryFlag = true;
                            }
                        }
                    }

                    this->dIndex_distance_constraints_wrt_StepLeftPosition.push_back(dIndex_distance_constraints_wrt_StepLeftPositionVariable);
                }

                if (!this->isDefinedRelativeToCentralBody)
                {
                    //Step 2.2: time variables - only relevant if reference body is NOT the central body
                    std::vector<size_t>& ListOfTimeVariablesAffectingCurrentStepLeftState = this->myStep->getListOfTimeVariablesAffectingCurrentStepLeftState();
                    for (size_t varIndex = 0; varIndex < ListOfTimeVariablesAffectingCurrentStepLeftState.size(); ++varIndex)
                    {
                        size_t Xindex = ListOfTimeVariablesAffectingCurrentStepLeftState[varIndex];

                        this->create_sparsity_entry(Fdescriptions->size() - 1,
                            Xindex,
                            this->G_indices_distance_constraints_wrt_StepLeftTime);
                    }
                }
            }//end derivatives for the non-trivial case
        }//end calcbounds

        void ParallelShootingStepDistanceConstraint::process(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: evaluate the constraint
            if (this->isDefinedRelativeToCentralBody
                && (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalRADEC
                    || this->myOptions->ParallelShootingConstraintStateRepresentation == StateRepresentation::SphericalAZFPA))
            {
                this->distance_from_body = X[this->myStep->getXindex_state_elements()[0]];
            }
            else
            {
                math::Matrix<doubleType>& SpacecraftState = this->myStep->get_StateStepLeftInertial();

                //Step 1.1: get the position vector of the spacecraft relative to the body
                doubleType body_state[12];
                this->myBody->locate_body(
                    SpacecraftState(7),
                    body_state,
                    (needG),
                    *this->myOptions);
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->distance_constraint_relative_position(stateIndex) = SpacecraftState(stateIndex) - body_state[stateIndex];
                    this->distance_constraint_body_position_time_derivatives(stateIndex) = body_state[6 + stateIndex] _GETVALUE;
                }

                //Step 1.2: apply the constraint
                this->distance_from_body = this->distance_constraint_relative_position.norm();
            }

            //Step 1.3: put the constraint into the constraint vector
            F[Findex++] = (this->distance_from_body - this->upperBound) / this->myUniverse->LU;

            //Step 2: derivatives
            if (needG)
            {
                if (this->isDefinedRelativeToCentralBody
                    && (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalRADEC
                        || this->myOptions->ParallelShootingConstraintStateRepresentation == StateRepresentation::SphericalAZFPA))
                {
                    //lucky us, there's only one entry! and no time dependence!

                    size_t Gindex = this->G_indices_distance_constraints_wrt_StepLeftPosition[0];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        / this->myUniverse->LU;
                }
                else
                {
                    math::Matrix<double> dConstraintPointState_dDecisionVariable(10, 1, 0.0);

                    //Step 2.1 non-time variables
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StepLeftState = this->myStep->get_Derivatives_of_StateStepLeftInertial();//Xindex, stateIndex, derivative value

                    for (size_t varIndex = 0; varIndex < this->dIndex_distance_constraints_wrt_StepLeftPosition.size(); ++varIndex)
                    {
                        dConstraintPointState_dDecisionVariable.assign_zeros();

                        for (size_t dIndex : this->dIndex_distance_constraints_wrt_StepLeftPosition[varIndex])
                        {
                            size_t stateIndex = std::get<1>(Derivatives_of_StepLeftState[dIndex]);

                            dConstraintPointState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StepLeftState[dIndex]);
                        }

                        size_t Gindex = this->G_indices_distance_constraints_wrt_StepLeftPosition[varIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            / this->distance_from_body _GETVALUE
                            * (this->distance_constraint_relative_position(0) * dConstraintPointState_dDecisionVariable(0)
                                + this->distance_constraint_relative_position(1) * dConstraintPointState_dDecisionVariable(1)
                                + this->distance_constraint_relative_position(2) * dConstraintPointState_dDecisionVariable(2)) _GETVALUE
                            / this->myUniverse->LU;

                        if (!this->isDefinedRelativeToCentralBody)
                        {
                            //Step 2.2: time variables
                            std::vector< std::vector< std::tuple<size_t, size_t> > >& DerivativesOfCurrentStepLeftStateByTimeVariable = this->myStep->getDerivativesOfCurrentStepLeftStateByTimeVariable();
                            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StepLeftState_wrt_Time = this->myStep->get_Derivatives_of_StateStepLeftInertial_wrt_Time();//Xindex, stateIndex, derivative value


                            for (size_t varIndex = 0; varIndex < DerivativesOfCurrentStepLeftStateByTimeVariable.size(); ++varIndex)
                            {
                                dConstraintPointState_dDecisionVariable.assign_zeros();

                                for (size_t entryIndex = 0; entryIndex < DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex].size(); ++entryIndex)
                                {
                                    size_t stateIndex = std::get<0>(DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex][entryIndex]);
                                    size_t dIndex = std::get<1>(DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex][entryIndex]);

                                    dConstraintPointState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StepLeftState_wrt_Time[dIndex]);
                                }

                                size_t Gindex = this->G_indices_distance_constraints_wrt_StepLeftTime[varIndex];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                size_t num_steps = this->myJourneyOptions->override_num_steps ? this->myJourneyOptions->number_of_steps : this->myOptions->num_timesteps;
                                double timeScale = varIndex < DerivativesOfCurrentStepLeftStateByTimeVariable.size() - 1 ?
                                    1.0 : (double)(this->stepIndex) / num_steps;

                                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                    / this->distance_from_body _GETVALUE
                                    * (this->distance_constraint_relative_position(0) * dConstraintPointState_dDecisionVariable(0)
                                        + this->distance_constraint_relative_position(1) * dConstraintPointState_dDecisionVariable(1)
                                        + this->distance_constraint_relative_position(2) * dConstraintPointState_dDecisionVariable(2)
                                        - this->distance_constraint_relative_position(0) * distance_constraint_body_position_time_derivatives(0) * timeScale
                                        - this->distance_constraint_relative_position(1) * distance_constraint_body_position_time_derivatives(1) * timeScale
                                        - this->distance_constraint_relative_position(2) * distance_constraint_body_position_time_derivatives(2) * timeScale
                                        ) _GETVALUE
                                    / this->myUniverse->LU;
                            }
                        }//end time derivatives, only relevant when the constraint is defined to a body other than the central body
                    }//end non-central body constraint derivatives
                }//end derivatives for the non-trivial case
            }//end derivatives
        }//end process
    }//end namespace Phases
}//end namespace EMTG