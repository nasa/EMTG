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

#include "BoundaryDistanceConstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            BoundaryDistanceConstraint::BoundaryDistanceConstraint(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                BoundaryEventBase* myBoundaryEvent,
                const std::string& constraintDefinition) :
                SpecializedBoundaryConstraintBase::SpecializedBoundaryConstraintBase(name,
                    journeyIndex,
                    phaseIndex,
                    stageIndex,
                    Universe,
                    mySpacecraft,
                    myOptions,
                    myBoundaryEvent,
                    constraintDefinition)
            {
                this->r_sc_CB = math::Matrix<doubleType>(3, 1, 0.0);
                this->r_body_CB = math::Matrix<doubleType>(3, 1, 0.0);
                this->dr_sc_CB_dt = math::Matrix<double>(3, 1, 0.0);
                this->dr_body_CB_dt = math::Matrix<double>(3, 1, 0.0);

                this->myReferenceFrame = ReferenceFrame::ICRF;
            }

            void BoundaryDistanceConstraint::calcbounds()
            {
                // Step 1: parse the constraint definition
				// constraintDefinition is like j1p0_departure_distanceconstraint_cb_0.9au_10.0au
                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                //Step 2: which body are we doing this with respect to?
                if (ConstraintDefinitionCell[3] == "cb")
                {
                    this->isCentralBodyConstraint = true;
                    this->bodyName = this->myUniverse->central_body_name;
                }
                else
                {
                    this->isCentralBodyConstraint = false;

                    int bodyIndex = std::stoi(ConstraintDefinitionCell[3]) - 1;
                    this->BoundaryDistanceConstraint::myBody = &this->myUniverse->bodies[bodyIndex];
                    this->bodyName = this->BoundaryDistanceConstraint::myBody->name;
                }

                //Step 3: create the constraint
                //figure out the lower and upper bounds - depends on units
                if (ConstraintDefinitionCell[4].find("km") < 1024)
                {
                    this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[4], "km")) / this->myUniverse->LU);
                }
                else if (ConstraintDefinitionCell[4].find("au") < 1024)
                {
                    this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[4], "au")) * this->myOptions->AU / this->myUniverse->LU);
                }
                else if (ConstraintDefinitionCell[4].find("lu") < 1024)
                {
                    this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[4], "lu")) * this->myUniverse->LU / this->myUniverse->LU);
                }
                else
                {
                    this->Flowerbounds->push_back(std::stod(ConstraintDefinitionCell[4]) / this->myUniverse->LU);
                }

                if (ConstraintDefinitionCell[5].find("km") < 1024)
                {
                    this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[5], "km")) / this->myUniverse->LU);
                }
                else if (ConstraintDefinitionCell[5].find("au") < 1024)
                {
                    this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[5], "au")) * this->myOptions->AU / this->myUniverse->LU);
                }
                else if (ConstraintDefinitionCell[5].find("lu") < 1024)
                {
                    this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[5], "lu")) * this->myUniverse->LU / this->myUniverse->LU);
                }
                else
                {
                    this->Fupperbounds->push_back(std::stod(ConstraintDefinitionCell[5]) / this->myUniverse->LU);
                }
                this->Fdescriptions->push_back(prefix + "distance constraint from " + bodyName);

                //Step 5: sparsity pattern
                //derivatives with respect to anything influencing the boundary event's right-hand position vector
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    //non-time variables
                    std::vector<size_t> state_dIndex_distance_with_respect_to_StateAfterEvent;
                    std::vector<size_t> state_Gindex_distance_constraint_wrt_StateAfterEvent_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent[dIndex]) == stateIndex)
                        {
                            state_dIndex_distance_with_respect_to_StateAfterEvent.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent[dIndex]),
                                state_Gindex_distance_constraint_wrt_StateAfterEvent_variables);
                        }
                    }
                    this->dIndex_distance_with_respect_to_StateAfterEvent.push_back(state_dIndex_distance_with_respect_to_StateAfterEvent);
                    this->Gindex_distance_constraint_wrt_StateAfterEvent_variables.push_back(state_Gindex_distance_constraint_wrt_StateAfterEvent_variables);

                    //time variables
                    std::vector<size_t> state_dIndex_distance_with_respect_to_StateAfterEvent_wrt_Time;
                    std::vector<size_t> state_Gindex_distance_constraint_wrt_StateAfterEvent_time_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]) == stateIndex)
                        {
                            state_dIndex_distance_with_respect_to_StateAfterEvent_wrt_Time.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]),
                                state_Gindex_distance_constraint_wrt_StateAfterEvent_time_variables);
                        }
                    }
                    this->dIndex_distance_with_respect_to_StateAfterEvent_wrt_Time.push_back(state_dIndex_distance_with_respect_to_StateAfterEvent_wrt_Time);
                    this->Gindex_distance_constraint_wrt_StateAfterEvent_time_variables.push_back(state_Gindex_distance_constraint_wrt_StateAfterEvent_time_variables);
                }

                //derivatives with respect to time variables that affect body position
                if (!this->isCentralBodyConstraint)
                {
                    std::vector<size_t> timeVariables = this->myBoundaryEvent->get_Xindices_EventRightEpoch();
                    for (size_t Xindex : timeVariables)
                    {
                        this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                            Xindex,
                            this->Gindex_distance_constraint_wrt_time_variables);
                    }
                }
            }//end calcbounds()

            void BoundaryDistanceConstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {
                //Step 1: get the state relative to the central body
                math::Matrix<doubleType>& StateRelativeToCentralBody = this->myBoundaryEvent->get_state_after_event();
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    this->r_sc_CB(stateIndex) = StateRelativeToCentralBody(stateIndex);

                //Step 2: where am I relative to the body of interest?
                if (this->isCentralBodyConstraint)
                {
                    this->r_sc_body = this->r_sc_CB;
                }//end distance from central body constraint
                else //distance from any body other than the central body
                {
                    //Step 2.1: where is the body of interest relative to the central body?
                    doubleType temp_body_state[12];
                    this->BoundaryDistanceConstraint::myBody->locate_body(StateRelativeToCentralBody(7),
                        temp_body_state,
                        needG,
                        *this->myOptions);

                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        this->r_body_CB(stateIndex) = temp_body_state[stateIndex];
                        this->dr_body_CB_dt(stateIndex) = temp_body_state[stateIndex + 6] _GETVALUE;
                    }

                    this->r_sc_body = this->r_sc_CB - this->r_body_CB;
                }//end distance from some other body constraint

                    //Step 3: compute the distance
                this->Distance = this->r_sc_body.norm();
                F[Findex++] = this->Distance / this->myUniverse->LU;

                //Step 4: derivatives
                if (needG)
                {
                    //Step 4.1: we are going to be adding components on top of components, so let's start by zeroing out all of our G entries
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        for (size_t entryIndex = 0; entryIndex < this->Gindex_distance_constraint_wrt_StateAfterEvent_variables[stateIndex].size(); ++entryIndex)
                        {
                            G[this->Gindex_distance_constraint_wrt_StateAfterEvent_variables[stateIndex][entryIndex]] = 0.0;
                        }
                        for (size_t entryIndex = 0; entryIndex < this->Gindex_distance_constraint_wrt_StateAfterEvent_time_variables[stateIndex].size(); ++entryIndex)
                        {
                            G[this->Gindex_distance_constraint_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex]] = 0.0;
                        }
                    }
                    for (size_t entryIndex = 0; entryIndex < this->Gindex_distance_constraint_wrt_time_variables.size(); ++entryIndex)
                    {
                        G[this->Gindex_distance_constraint_wrt_time_variables[entryIndex]] = 0.0;
                    }

                    //Step 4.2: derivatives with respect to non-time variables affecting state after boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();

                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        for (size_t entryIndex = 0; entryIndex < this->dIndex_distance_with_respect_to_StateAfterEvent[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_distance_with_respect_to_StateAfterEvent[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_distance_constraint_wrt_StateAfterEvent_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (this->r_sc_body(stateIndex) / this->Distance * std::get<2>(Derivatives_of_StateAfterEvent[dIndex])) _GETVALUE;

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry
                                / this->myUniverse->LU;
                        }
                    }

                    //Step 4.3: derivatives with respect to time variables affecting state after boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        for (size_t entryIndex = 0; entryIndex < this->dIndex_distance_with_respect_to_StateAfterEvent_wrt_Time[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_distance_with_respect_to_StateAfterEvent_wrt_Time[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_distance_constraint_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (this->r_sc_body(stateIndex) / this->Distance * std::get<2>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex])) _GETVALUE;

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry
                                / this->myUniverse->LU;
                        }
                    }


                    //Step 4.4: derivatives of the body-of-interest's position with respect to time variables
                    if (!this->isCentralBodyConstraint)
                    {
                        double DistanceTimeDerivative = ((this->r_sc_body(0) * (this->dr_sc_CB_dt(0) - this->dr_body_CB_dt(0))
                            + this->r_sc_body(1) * (this->dr_sc_CB_dt(1) - this->dr_body_CB_dt(1))
                            + this->r_sc_body(2) * (this->dr_sc_CB_dt(2) - this->dr_body_CB_dt(2))) / this->Distance)_GETVALUE;

                        for (size_t entryIndex = 0; entryIndex < this->Gindex_distance_constraint_wrt_time_variables.size(); ++entryIndex)
                        {
                            size_t Gindex = this->Gindex_distance_constraint_wrt_time_variables[entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * DistanceTimeDerivative
                                / this->myUniverse->LU;
                        }
                    }
                }//end derivatives
            }//end process_constraint()

            void BoundaryDistanceConstraint::output(std::ofstream& outputfile)
            {
                std::string framestring = this->myUniverse->central_body_name + " " + this->ReferenceFrameStrings[this->myReferenceFrame];

                outputfile << this->myBoundaryEvent->getName() << " distance from " << this->bodyName << " (km, " << framestring << "): " << this->Distance _GETVALUE << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG