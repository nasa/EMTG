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

#include "BoundaryMassConstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            BoundaryMassConstraint::BoundaryMassConstraint(const std::string& name,
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
            }

            void BoundaryMassConstraint::calcbounds()
            {
                //Step 1: parse the constraint definition
                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                //Step 2: what is the constraint value?
                this->BoundaryMassLowerBound = std::stod(ConstraintDefinitionCell[3]);
                this->BoundaryMassUpperBound = std::stod(ConstraintDefinitionCell[4]);

                //Step 3: create the constraint
                //figure out the lower and upper bounds
                this->Flowerbounds->push_back(this->BoundaryMassLowerBound / this->BoundaryMassUpperBound);
                this->Fupperbounds->push_back(1.0);

                this->Fdescriptions->push_back(prefix + "boundary mass constraint");

                //Step 5: sparsity pattern
                //derivatives with respect to anything influencing the boundary event's right-hand position vector
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                {
                    size_t stateIndex = 6;
                    //non-time variables
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent[dIndex]) == stateIndex)
                        {
                            this->dIndex_BoundaryMassConstraint_wrt_StateAfterEvent_variables.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent[dIndex]),
                                this->Gindex_BoundaryMassConstraint_wrt_StateAfterEvent_variables);
                        }
                    }

                    //time variables
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]) == stateIndex)
                        {
                            this->dIndex_BoundaryMassConstraint_wrt_StateAfterEvent_time_variables.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]),
                                this->Gindex_BoundaryMassConstraint_wrt_StateAfterEvent_time_variables);
                        }
                    }
                }
            }//end calcbounds()

            void BoundaryMassConstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {
                //Step 1: get the spacecraft mass
                math::Matrix<doubleType>& SpacecraftState = this->myBoundaryEvent->get_state_after_event();
                this->ActualMass = SpacecraftState(6);

                //Step 2: apply the constraint
                F[Findex++] = this->ActualMass / this->BoundaryMassUpperBound;

                //Step 4: derivatives
                if (needG)
                {
                    //Step 4.1: we are going to be adding components on top of components, so let's start by zeroing out all of our G entries
                    for (size_t Gindex : this->Gindex_BoundaryMassConstraint_wrt_StateAfterEvent_variables)
                        G[Gindex] = 0.0;
                    for (size_t Gindex : this->Gindex_BoundaryMassConstraint_wrt_StateAfterEvent_time_variables)
                        G[Gindex] = 0.0;
                    
                    //Step 4.2: derivatives with respect to non-time variables affecting state after boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();
                    size_t stateIndex = 6;

                    for (size_t entryIndex = 0; entryIndex < this->dIndex_BoundaryMassConstraint_wrt_StateAfterEvent_variables.size(); ++entryIndex)
                    {
                        size_t dIndex = this->dIndex_BoundaryMassConstraint_wrt_StateAfterEvent_variables[entryIndex];
                        size_t Gindex = this->Gindex_BoundaryMassConstraint_wrt_StateAfterEvent_variables[entryIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        double Gentry = std::get<2>(Derivatives_of_StateAfterEvent[dIndex]);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * Gentry
                            / this->BoundaryMassUpperBound;
                    }

                    //Step 4.3: derivatives with respect to time variables affecting state after boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();
                    
                    for (size_t entryIndex = 0; entryIndex < this->dIndex_BoundaryMassConstraint_wrt_StateAfterEvent_time_variables.size(); ++entryIndex)
                    {
                        size_t dIndex = this->dIndex_BoundaryMassConstraint_wrt_StateAfterEvent_time_variables[entryIndex];
                        size_t Gindex = this->Gindex_BoundaryMassConstraint_wrt_StateAfterEvent_time_variables[entryIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        double Gentry = std::get<2>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * Gentry
                            / this->BoundaryMassUpperBound;
                    }
                }//end derivatives
            }//end process_constraint()

            void BoundaryMassConstraint::output(std::ofstream& outputfile)
            {
                outputfile << this->myBoundaryEvent->getName() << " mass " << this->ActualMass _GETVALUE << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG