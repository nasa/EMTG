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

#include "BoundaryBCFlatitudeConstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            BoundaryBCFlatitudeConstraint::BoundaryBCFlatitudeConstraint(const std::string& name,
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
                this->myReferenceFrame = ReferenceFrame::J2000_BCF;
            }//end constructor

            void BoundaryBCFlatitudeConstraint::calcbounds()
            {
                // Step 1: parse the constraint definition
				// constraintDefinition is like p0_departure_longitude_84.0_85.0
                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                //Step 2: lower bound
                this->Flowerbounds->push_back(std::stod(ConstraintDefinitionCell[3]) * math::deg2rad);

                //Step 3: upper bound
                this->Fupperbounds->push_back(std::stod(ConstraintDefinitionCell[4]) * math::deg2rad);
                this->Fdescriptions->push_back(prefix + "BCF latitude");

                //Step 5: sparsity pattern
                //derivatives with respect to anything influencing the boundary event's right-hand position vector
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();

                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    //non-time variables
                    std::vector<size_t> state_dIndex_with_respect_to_StateBeforeEvent;
                    std::vector<size_t> state_Gindex_constraint_wrt_StateBeforeEvent_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforeEvent.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateBeforeEvent[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_StateBeforeEvent.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateBeforeEvent[dIndex]),
                                state_Gindex_constraint_wrt_StateBeforeEvent_variables);
                        }
                    }
                    this->dIndex_with_respect_to_StateBeforeEvent.push_back(state_dIndex_with_respect_to_StateBeforeEvent);
                    this->Gindex_constraint_wrt_StateBeforeEvent_variables.push_back(state_Gindex_constraint_wrt_StateBeforeEvent_variables);

                    //time variables
                    std::vector<size_t> state_dIndex_with_respect_to_StateBeforeEvent_wrt_Time;
                    std::vector<size_t> state_Gindex_constraint_wrt_StateBeforeEvent_time_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_StateBeforeEvent_wrt_Time.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]),
                                state_Gindex_constraint_wrt_StateBeforeEvent_time_variables);
                        }
                    }
                    this->dIndex_with_respect_to_StateBeforeEvent_wrt_Time.push_back(state_dIndex_with_respect_to_StateBeforeEvent_wrt_Time);
                    this->Gindex_constraint_wrt_StateBeforeEvent_time_variables.push_back(state_Gindex_constraint_wrt_StateBeforeEvent_time_variables);
                }

                //derivatives with respect to time variables that affect current epoch, because this is a BCF constraint
                std::vector<size_t> timeVariables = this->myBoundaryEvent->get_Xindices_EventRightEpoch();
                for (size_t Xindex : timeVariables)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        Xindex,
                        this->Gindex_constraint_wrt_time_variables);
                }
            }//end calcbounds()

            void BoundaryBCFlatitudeConstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {
                //Step 1: get the state in ICRF
                math::Matrix<doubleType>& StateBeforeEventICRF = this->myBoundaryEvent->get_state_before_event();

                //Step 2: convert to BCF
                this->myUniverse->LocalFrame.construct_rotation_matrices(StateBeforeEventICRF(7), needG);
                math::Matrix<doubleType> R_from_ICRF_to_J2000BCF = this->myUniverse->LocalFrame.get_R_from_ICRF_to_J2000BCF();
                this->PositionBeforeEventBCF = R_from_ICRF_to_J2000BCF * StateBeforeEventICRF.getSubMatrix1D(0, 2);


                //Step 3: compute longitude
                doubleType rxy_J2000BCF = sqrt(this->PositionBeforeEventBCF(0) * this->PositionBeforeEventBCF(0) + this->PositionBeforeEventBCF(1) * this->PositionBeforeEventBCF(1));
                this->Latitude = atan2(this->PositionBeforeEventBCF(2), rxy_J2000BCF);

                F[Findex++] = this->Latitude;

                //Step 4: derivatives
                if (needG)
                {
                    //Step 4.1: we are going to be adding components on top of components, so let's start by zeroing out all of our G entries
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_StateBeforeEvent_variables[stateIndex].size(); ++entryIndex)
                        {
                            G[this->Gindex_constraint_wrt_StateBeforeEvent_variables[stateIndex][entryIndex]] = 0.0;
                        }
                        for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_StateBeforeEvent_time_variables[stateIndex].size(); ++entryIndex)
                        {
                            G[this->Gindex_constraint_wrt_StateBeforeEvent_time_variables[stateIndex][entryIndex]] = 0.0;
                        }
                    }
                    for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_time_variables.size(); ++entryIndex)
                    {
                        G[this->Gindex_constraint_wrt_time_variables[entryIndex]] = 0.0;
                    }

                    //Step 4.2: derivatives of longitude w.r.t. BCF states
                    double dLatitude_dzJ2000BCF = (rxy_J2000BCF / (this->PositionBeforeEventBCF(2) * this->PositionBeforeEventBCF(2) + rxy_J2000BCF * rxy_J2000BCF))_GETVALUE;
                    double dLatitude_drxyJ2000BCF = (-this->PositionBeforeEventBCF(2) / (this->PositionBeforeEventBCF(2) * this->PositionBeforeEventBCF(2) + rxy_J2000BCF * rxy_J2000BCF))_GETVALUE;
                    double dLatitude_dxJ2000BCF = dLatitude_drxyJ2000BCF * (this->PositionBeforeEventBCF(0) / rxy_J2000BCF) _GETVALUE;
                    double dLatitude_dyJ2000BCF = dLatitude_drxyJ2000BCF * (this->PositionBeforeEventBCF(1) / rxy_J2000BCF) _GETVALUE;

                    //Step 4.3: derivatives with respect to non-time variables affecting state Before boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent();

                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        double dxJ2000BCF_dstateICRF = R_from_ICRF_to_J2000BCF(0, stateIndex) _GETVALUE;
                        double dyJ2000BCF_dstateICRF = R_from_ICRF_to_J2000BCF(1, stateIndex) _GETVALUE;
                        double dzJ2000BCF_dstateICRF = R_from_ICRF_to_J2000BCF(2, stateIndex) _GETVALUE;

                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_StateBeforeEvent[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_StateBeforeEvent[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_StateBeforeEvent_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (dLatitude_dxJ2000BCF * dxJ2000BCF_dstateICRF + dLatitude_dyJ2000BCF * dyJ2000BCF_dstateICRF + dLatitude_dzJ2000BCF * dzJ2000BCF_dstateICRF)
                                * std::get<2>(Derivatives_of_StateBeforeEvent[dIndex]);
                            
                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }

                    //Step 4.3: derivatives with respect to time variables affecting state Before boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();

                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        double dxJ2000BCF_dstateICRF = R_from_ICRF_to_J2000BCF(0, stateIndex) _GETVALUE;
                        double dyJ2000BCF_dstateICRF = R_from_ICRF_to_J2000BCF(1, stateIndex) _GETVALUE;
                        double dzJ2000BCF_dstateICRF = R_from_ICRF_to_J2000BCF(2, stateIndex) _GETVALUE;

                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_StateBeforeEvent_wrt_Time[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_StateBeforeEvent_wrt_Time[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_StateBeforeEvent_time_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (dLatitude_dxJ2000BCF * dxJ2000BCF_dstateICRF + dLatitude_dyJ2000BCF * dyJ2000BCF_dstateICRF + dLatitude_dzJ2000BCF * dzJ2000BCF_dstateICRF)
                                * std::get<2>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }


                    //Step 4.4: derivatives of longitude with respect to frame rotation due to change in epoch
                    math::Matrix<doubleType> dR_from_ICRF_to_J2000BCF_dt = this->myUniverse->LocalFrame.get_dR_from_ICRF_to_J2000BCF_dt();
                    this->dPositionBeforeEventBCF_dt = dR_from_ICRF_to_J2000BCF_dt * StateBeforeEventICRF.getSubMatrix1D(0, 2);

                    double timeDerivative = (dLatitude_dxJ2000BCF * this->dPositionBeforeEventBCF_dt(0)
                        + dLatitude_dyJ2000BCF * this->dPositionBeforeEventBCF_dt(1)
                        + dLatitude_dyJ2000BCF * this->dPositionBeforeEventBCF_dt(2)) _GETVALUE;

                    for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_time_variables.size(); ++entryIndex)
                    {
                        size_t Gindex = this->Gindex_constraint_wrt_time_variables[entryIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * timeDerivative;
                    }
                }//end derivatives
            }//end process_constraint()

            void BoundaryBCFlatitudeConstraint::output(std::ofstream& outputfile)
            {
                std::string framestring = this->myUniverse->central_body_name + " " + this->ReferenceFrameStrings[this->myReferenceFrame];

                outputfile << this->myBoundaryEvent->getName() << " BCF latitude (degrees, " << framestring << "): " << this->Latitude _GETVALUE / math::deg2rad  << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG