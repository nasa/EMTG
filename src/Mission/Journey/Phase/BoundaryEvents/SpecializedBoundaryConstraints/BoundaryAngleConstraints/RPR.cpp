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

#include "RPR.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {
            RPRconstraint::RPRconstraint(const std::string& name,
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
                this->R_cb_ref1 = math::Matrix<doubleType>(3, 1, 0.0);
                this->R_cb_ref2 = math::Matrix<doubleType>(3, 1, 0.0);
                this->R_probe_ref1 = math::Matrix<doubleType>(3, 1, 0.0);
                this->R_probe_ref2 = math::Matrix<doubleType>(3, 1, 0.0);
                this->dcosRPR_dR_probe_ref1 = math::Matrix<double>(3, 1, 0.0);
                this->dcosRPR_dR_probe_ref2 = math::Matrix<double>(3, 1, 0.0);
                this->dR_cb_ref1_dt = math::Matrix<double>(3, 1, 0.0);
                this->dR_cb_ref2_dt = math::Matrix<double>(3, 1, 0.0);
                this->dR_probe_ref1_dt = math::Matrix<double>(3, 1, 0.0);
                this->dR_probe_ref2_dt = math::Matrix<double>(3, 1, 0.0);

                //derived class version of myBoundaryEvent has a different type than SpecializedBoundaryConstraintBase uses
                this->myBoundaryEvent = myBoundaryEvent;
            }

            void RPRconstraint::calcbounds()
            {
                // Step 1: parse the constraint definition
                // constraintDefinition is like p0_arrival_RRP_Ref1ID_Ref2ID_0.0_10.0 (numbers are in degrees)
                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                //Step 2: which bodies are we doing this with respect to?
                //Step 2.1: reference 1
                if (ConstraintDefinitionCell[3] == "cb")
                {
                    this->ref1IsCentralBody = true;
                    this->ref1Name = this->myUniverse->central_body_name;
                }
                else
                {
                    this->ref1IsCentralBody = false;

                    int bodyIndex = std::stoi(ConstraintDefinitionCell[3]) - 1;
                    this->myReference1 = &this->myUniverse->bodies[bodyIndex];
                    this->ref1Name = this->myReference1->name;
                }

                //Step 2.2: reference 2
                if (ConstraintDefinitionCell[4] == "cb")
                {
                    this->ref2IsCentralBody = true;
                    this->ref2Name = this->myUniverse->central_body_name;
                }
                else
                {
                    this->ref2IsCentralBody = false;

                    int bodyIndex = std::stoi(ConstraintDefinitionCell[4]) - 1;
                    this->myReference2 = &this->myUniverse->bodies[bodyIndex];
                    this->ref2Name = this->myReference2->name;
                }

                //Step 3: if the two references are the same, throw an error
                if (this->ref1Name == this->ref2Name)
                {
                    throw std::invalid_argument("RPR constraint must have two different reference body. Constraint \"" + constraintDefinition + "\" in Journey " + std::to_string(journeyIndex) + " uses the same body twice.");
                }

                //Step 4: create the constraint
                //figure out the lower and upper bounds - convert from degrees to radians
                //the bigger number has the smaller cosine
                this->Flowerbounds->push_back(cos(std::stod(ConstraintDefinitionCell[6]) * math::PI / 180.0));

                this->Fupperbounds->push_back(cos(std::stod(ConstraintDefinitionCell[5]) * math::PI / 180.0));

                this->Fdescriptions->push_back(prefix + this->ref1Name + "-Probe-" + this->ref2Name + " constraint");

                //Step 5: sparsity pattern
                //derivatives with respect to anything influencing the boundary event's right-hand position and velocity vectors
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();
                std::vector<size_t> dIndex_VbeforeEvent_dVinfinity_in = this->myBoundaryEvent->getdIndex_VbeforeEvent_dVinfinity_in();

                for (size_t stateIndex : {0, 1, 2})
                {
                    //non-time variables
                    std::vector<size_t> state_dIndex_with_respect_to_PositionAfterEvent;
                    std::vector<size_t> state_Gindex_constraint_wrt_PositionAfterEvent_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_PositionAfterEvent.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent[dIndex]),
                                state_Gindex_constraint_wrt_PositionAfterEvent_variables);
                        }
                    }
                    this->dIndex_with_respect_to_PositionAfterEvent.push_back(state_dIndex_with_respect_to_PositionAfterEvent);
                    this->Gindex_constraint_wrt_PositionAfterEvent_variables.push_back(state_Gindex_constraint_wrt_PositionAfterEvent_variables);

                    //time variables
                    std::vector<size_t> state_dIndex_with_respect_to_PositionAfterEvent_wrt_Time;
                    std::vector<size_t> state_Gindex_constraint_wrt_PositionAfterEvent_time_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_PositionAfterEvent_wrt_Time.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]),
                                state_Gindex_constraint_wrt_PositionAfterEvent_time_variables);
                        }
                    }
                    this->dIndex_with_respect_to_PositionAfterEvent_wrt_Time.push_back(state_dIndex_with_respect_to_PositionAfterEvent_wrt_Time);
                    this->Gindex_constraint_wrt_PositionAfterEvent_time_variables.push_back(state_Gindex_constraint_wrt_PositionAfterEvent_time_variables);

                }

                //derivatives with respect to time variables that affect reference body position
                if (!(this->ref1IsCentralBody || this->ref2IsCentralBody))
                {
                    std::vector<size_t> timeVariables = this->myBoundaryEvent->get_Xindices_EventRightEpoch();
                    for (size_t Xindex : timeVariables)
                    {
                        this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                            Xindex,
                            this->Gindex_constraint_wrt_time_variables);
                    }
                }
            }//end calcbounds()

            void RPRconstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {

                //Step 1: get the state relative to the central body
                math::Matrix<doubleType>& SpacecraftStateRelativeToCentralBody = this->myBoundaryEvent->get_state_after_event();

                //Step 2: get the vector from the central body to the reference bodies
                //Step 2.1: reference 1
                if (this->ref1IsCentralBody) //is the reference body the central body?
                {
                    for (size_t stateIndex : {0, 1, 2})
                        this->R_cb_ref1(stateIndex) = 0.0;
                }
                else //reference body is some other body in the Universe
                {
                    //Step 2.1: find the reference body relative to the central body
                    doubleType temp_body_state[12];
                    this->myReference1->locate_body(SpacecraftStateRelativeToCentralBody(7),
                        temp_body_state,
                        needG,
                        *this->myOptions);

                    //Step 2.2: compute the vector from the body to the reference body
                    for (size_t stateIndex : {0, 1, 2})
                    {
                        this->R_cb_ref1(stateIndex) = temp_body_state[stateIndex];
                        this->dR_cb_ref1_dt(stateIndex) = temp_body_state[stateIndex + 6] _GETVALUE;
                    }
                }
                //Step 2.2: reference 2
                if (this->ref2IsCentralBody) //is the reference body the central body?
                {
                    for (size_t stateIndex : {0, 1, 2})
                        this->R_cb_ref2(stateIndex) = 0.0;
                }
                else //reference body is some other body in the Universe
                {
                    //Step 2.1: find the reference body relative to the central body
                    doubleType temp_body_state[12];
                    this->myReference2->locate_body(SpacecraftStateRelativeToCentralBody(7),
                        temp_body_state,
                        needG,
                        *this->myOptions);

                    //Step 2.2: compute the vector from the body to the reference body
                    for (size_t stateIndex : {0, 1, 2})
                    {
                        this->R_cb_ref2(stateIndex) = temp_body_state[stateIndex];
                        this->dR_cb_ref2_dt(stateIndex) = temp_body_state[stateIndex + 6] _GETVALUE;
                    }
                }

                //Step 3: get the vector from reference 1 to reference 1
                for (size_t stateIndex : {0, 1, 2})
                {
                    this->R_probe_ref1(stateIndex) = SpacecraftStateRelativeToCentralBody(stateIndex) - this->R_cb_ref1(stateIndex);
                }

                //Step 4: get the vector from reference 2 to the probe
                for (size_t stateIndex : {0, 1, 2})
                {
                    this->R_probe_ref2(stateIndex) = SpacecraftStateRelativeToCentralBody(stateIndex) - this->R_cb_ref2(stateIndex);
                }

                //Step 5: compute the angle between the vectors
                this->cosRPR = this->R_probe_ref1.dot(this->R_probe_ref2) / (this->R_probe_ref1.norm() * this->R_probe_ref2.norm());
                this->Angle = math::safe_acos(this->cosRPR);

                //Step 6: apply the constraint - we will constrain the cos(Angle) instead of the Angle itself
                F[Findex++] = this->cosRPR;

                //Step 7: derivatives
                if (needG)
                {
                    //Step 6.1: we are going to be adding components on top of components, so let's start by zeroing out all of our G entries
                    for (size_t stateIndex : {0, 1, 2})
                    {
                        for (size_t Gindex : this->Gindex_constraint_wrt_PositionAfterEvent_variables[stateIndex])
                        {
                            G[Gindex] = 0.0;
                        }
                        for (size_t Gindex : this->Gindex_constraint_wrt_PositionAfterEvent_time_variables[stateIndex])
                        {
                            G[Gindex] = 0.0;
                        }
                    }
                    for (size_t Gindex : this->Gindex_constraint_wrt_time_variables)
                    {
                        G[Gindex] = 0.0;
                    }

                    //Step 6.2: derivatives with respect to non-time variables affecting state after boundary event
                    doubleType r_probe_ref1 = this->R_probe_ref1.norm();
                    double dr_probe_ref1_drx = (this->R_probe_ref1(0) / r_probe_ref1)_GETVALUE;
                    double dr_probe_ref1_dry = (this->R_probe_ref1(1) / r_probe_ref1)_GETVALUE;
                    double dr_probe_ref1_drz = (this->R_probe_ref1(2) / r_probe_ref1)_GETVALUE;
                    doubleType r_probe_ref2 = this->R_probe_ref2.norm();
                    double dr_probe_ref2_drx = (this->R_probe_ref2(0) / r_probe_ref2)_GETVALUE;
                    double dr_probe_ref2_dry = (this->R_probe_ref2(1) / r_probe_ref2)_GETVALUE;
                    double dr_probe_ref2_drz = (this->R_probe_ref2(2) / r_probe_ref2)_GETVALUE;

                    double r1x = this->R_probe_ref1(0)_GETVALUE;
                    double r1y = this->R_probe_ref1(1)_GETVALUE;
                    double r1z = this->R_probe_ref1(2)_GETVALUE;
                    double r2x = this->R_probe_ref2(0)_GETVALUE;
                    double r2y = this->R_probe_ref2(1)_GETVALUE;
                    double r2z = this->R_probe_ref2(2)_GETVALUE;
                    double r1 = r_probe_ref1 _GETVALUE;
                    double r2 = r_probe_ref2 _GETVALUE;

                    this->dcosRPR_dR_probe_ref1(0) = r2x / (r1 * r2) - (r1x * r2x + r1y * r2y + r1z * r2z) * dr_probe_ref1_drx / (r1 * r1 * r2);
                    this->dcosRPR_dR_probe_ref1(1) = r2y / (r1 * r2) - (r1x * r2x + r1y * r2y + r1z * r2z) * dr_probe_ref1_dry / (r1 * r1 * r2);
                    this->dcosRPR_dR_probe_ref1(2) = r2z / (r1 * r2) - (r1x * r2x + r1y * r2y + r1z * r2z) * dr_probe_ref1_drz / (r1 * r1 * r2);

                    this->dcosRPR_dR_probe_ref2(0) = r1x / (r1 * r2) - (r1x * r2x + r1y * r2y + r1z * r2z) * dr_probe_ref2_drx / (r1 * r2 * r2);
                    this->dcosRPR_dR_probe_ref2(1) = r1y / (r1 * r2) - (r1x * r2x + r1y * r2y + r1z * r2z) * dr_probe_ref2_dry / (r1 * r2 * r2);
                    this->dcosRPR_dR_probe_ref2(2) = r1z / (r1 * r2) - (r1x * r2x + r1y * r2y + r1z * r2z) * dr_probe_ref2_drz / (r1 * r2 * r2);

                    //Step 6.3: derivatives with respect to non-time variables affecting ONLY the probe's position
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();

                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_PositionAfterEvent[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_PositionAfterEvent[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_PositionAfterEvent_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (this->dcosRPR_dR_probe_ref1(stateIndex) + this->dcosRPR_dR_probe_ref2(stateIndex))
                                * std::get<2>(Derivatives_of_StateAfterEvent[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;

                        }
                    }//end loop over states

                    //Step 6.4: derivative with respect to time
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_PositionAfterEvent_wrt_Time[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_PositionAfterEvent_wrt_Time[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_PositionAfterEvent_time_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = this->dcosRPR_dR_probe_ref1(stateIndex) * (this->dR_cb_ref1_dt(stateIndex) - std::get<2>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]))
                                + this->dcosRPR_dR_probe_ref2(stateIndex) * (this->dR_cb_ref2_dt(stateIndex) - std::get<2>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]));
                                
                            G[Gindex] -= this->X_scale_factors->operator[](Xindex)
                                * Gentry;

                        }
                    }//end loop over states
                }//end derivatives

            }//end process_constraint()

            void RPRconstraint::output(std::ofstream& outputfile)
            {
                outputfile << this->myBoundaryEvent->getName() << " " << this->ref1Name + "-Probe-" + this->ref2Name + " angle (deg): " << this->Angle _GETVALUE * 180.0 / math::PI << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG