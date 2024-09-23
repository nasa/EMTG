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

#include "BoundaryIsApseConstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            BoundaryIsApseConstraint::BoundaryIsApseConstraint(const std::string& name,
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
                this->R = math::Matrix<doubleType>(3, 1, 0.0);
                this->V = math::Matrix<doubleType>(3, 1, 0.0);

                this->myReferenceFrame = ReferenceFrame::ICRF;
            }

            void BoundaryIsApseConstraint::calcbounds()
            {
                //Step 1: we don't have to parse the constraint because if we got here, it's all the same. So go straight to creating the constraint
                this->Flowerbounds->push_back(-math::SMALL);
                this->Fupperbounds->push_back(math::SMALL);
                this->Fdescriptions->push_back(prefix + "RdotV = 0");

                //Step 2: sparsity pattern
                //derivatives with respect to anything influencing the boundary event's right-hand position vector
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    //POSITION
                    //non-time variables
                    std::vector<size_t> state_dIndex_periapse_position_with_respect_to_StateAfterEvent;
                    std::vector<size_t> state_Gindex_periapse_position_wrt_StateAfterEvent_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent[dIndex]) == stateIndex)
                        {
                            state_dIndex_periapse_position_with_respect_to_StateAfterEvent.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent[dIndex]),
                                state_Gindex_periapse_position_wrt_StateAfterEvent_variables);
                        }
                    }
                    this->dIndex_periapse_position_with_respect_to_StateAfterEvent.push_back(state_dIndex_periapse_position_with_respect_to_StateAfterEvent);
                    this->Gindex_periapse_position_wrt_StateAfterEvent_variables.push_back(state_Gindex_periapse_position_wrt_StateAfterEvent_variables);

                    //time variables
                    std::vector<size_t> state_dIndex_periapse_position_with_respect_to_StateAfterEvent_wrt_Time;
                    std::vector<size_t> state_Gindex_periapse_position_wrt_StateAfterEvent_time_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]) == stateIndex)
                        {
                            state_dIndex_periapse_position_with_respect_to_StateAfterEvent_wrt_Time.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]),
                                state_Gindex_periapse_position_wrt_StateAfterEvent_time_variables);
                        }
                    }
                    this->dIndex_periapse_position_with_respect_to_StateAfterEvent_wrt_Time.push_back(state_dIndex_periapse_position_with_respect_to_StateAfterEvent_wrt_Time);
                    this->Gindex_periapse_position_wrt_StateAfterEvent_time_variables.push_back(state_Gindex_periapse_position_wrt_StateAfterEvent_time_variables);

                    //VELOCITY
                    //non-time variables
                    std::vector<size_t> state_dIndex_periapse_velocity_with_respect_to_StateAfterEvent;
                    std::vector<size_t> state_Gindex_periapse_velocity_wrt_StateAfterEvent_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent[dIndex]) == stateIndex + 3)
                        {
                            state_dIndex_periapse_velocity_with_respect_to_StateAfterEvent.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent[dIndex]),
                                state_Gindex_periapse_velocity_wrt_StateAfterEvent_variables);
                        }
                    }
                    this->dIndex_periapse_velocity_with_respect_to_StateAfterEvent.push_back(state_dIndex_periapse_velocity_with_respect_to_StateAfterEvent);
                    this->Gindex_periapse_velocity_wrt_StateAfterEvent_variables.push_back(state_Gindex_periapse_velocity_wrt_StateAfterEvent_variables);

                    //time variables
                    std::vector<size_t> state_dIndex_periapse_velocity_with_respect_to_StateAfterEvent_wrt_Time;
                    std::vector<size_t> state_Gindex_periapse_velocity_wrt_StateAfterEvent_time_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]) == stateIndex + 3)
                        {
                            state_dIndex_periapse_velocity_with_respect_to_StateAfterEvent_wrt_Time.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]),
                                state_Gindex_periapse_velocity_wrt_StateAfterEvent_time_variables);
                        }
                    }
                    this->dIndex_periapse_velocity_with_respect_to_StateAfterEvent_wrt_Time.push_back(state_dIndex_periapse_velocity_with_respect_to_StateAfterEvent_wrt_Time);
                    this->Gindex_periapse_velocity_wrt_StateAfterEvent_time_variables.push_back(state_Gindex_periapse_velocity_wrt_StateAfterEvent_time_variables);
                }
            }//end calcbounds()

            void BoundaryIsApseConstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {
                //Step 1: get the state
                math::Matrix<doubleType>& StateAfterEvent = this->myBoundaryEvent->get_state_after_event();
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->R(stateIndex) = StateAfterEvent(stateIndex);
                    this->V(stateIndex) = StateAfterEvent(stateIndex + 3);
                }

                //Step 2: compute R dot V
                doubleType r = this->R.norm();
                doubleType v = this->V.norm();
                this->RdotV = this->R.dot(this->V) / (r * v);

                //Step 3: apply the constraint
                F[Findex++] = this->RdotV;

                //Step 4: derivatives
                if (needG)
                {
                    //Step 4.1: components of the derivative
                    double x = StateAfterEvent(0) _GETVALUE;
                    double y = StateAfterEvent(1) _GETVALUE;
                    double z = StateAfterEvent(2) _GETVALUE;
                    double xdot = StateAfterEvent(3) _GETVALUE;
                    double ydot = StateAfterEvent(4) _GETVALUE;
                    double zdot = StateAfterEvent(5) _GETVALUE;
                    double dRdotV_dR[3], dRdotV_dV[3];
                    dRdotV_dR[0] = (xdot*y*y - x*ydot*y + xdot*z*z - x*zdot*z) / (r*r*r*v) _GETVALUE;
                    dRdotV_dR[1] = (ydot*x*x - xdot*y*x + ydot*z*z - y*zdot*z) / (r*r*r*v) _GETVALUE;
                    dRdotV_dR[2] = (zdot*x*x - xdot*z*x + zdot*y*y - ydot*z*y) / (r*r*r*v) _GETVALUE;
                    dRdotV_dV[0] = (x*ydot*ydot - xdot*y*ydot + x*zdot*zdot - xdot*z*zdot) / (r*v*v*v) _GETVALUE;
                    dRdotV_dV[1] = (y*xdot*xdot - x*ydot*xdot + y*zdot*zdot - ydot*z*zdot) / (r*v*v*v) _GETVALUE;
                    dRdotV_dV[2] = (z*xdot*xdot - x*zdot*xdot + z*ydot*ydot - y*zdot*ydot) / (r*v*v*v) _GETVALUE;

                    //Step 4.2: chain the components with the derivatives of state variables

                    //Step 4.2.1: non-time variables first
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();

                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        for (size_t entryIndex = 0; entryIndex < this->dIndex_periapse_position_with_respect_to_StateAfterEvent[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_periapse_position_with_respect_to_StateAfterEvent[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_periapse_position_wrt_StateAfterEvent_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = dRdotV_dR[stateIndex] * std::get<2>(Derivatives_of_StateAfterEvent[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }

                        for (size_t entryIndex = 0; entryIndex < this->dIndex_periapse_velocity_with_respect_to_StateAfterEvent[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_periapse_velocity_with_respect_to_StateAfterEvent[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_periapse_velocity_wrt_StateAfterEvent_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = dRdotV_dV[stateIndex] * std::get<2>(Derivatives_of_StateAfterEvent[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }

                    //Step 4.2.1: time variables
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        for (size_t entryIndex = 0; entryIndex < this->dIndex_periapse_position_with_respect_to_StateAfterEvent_wrt_Time[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_periapse_position_with_respect_to_StateAfterEvent_wrt_Time[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_periapse_position_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = dRdotV_dR[stateIndex] * std::get<2>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }

                        for (size_t entryIndex = 0; entryIndex < this->dIndex_periapse_velocity_with_respect_to_StateAfterEvent_wrt_Time[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_periapse_velocity_with_respect_to_StateAfterEvent_wrt_Time[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_periapse_velocity_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = dRdotV_dV[stateIndex] * std::get<2>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }
                }
            }//end process_constraint()

            void BoundaryIsApseConstraint::output(std::ofstream& outputfile)
            {
                std::string framestring = this->myUniverse->central_body_name + " " + this->ReferenceFrameStrings[this->myReferenceFrame];

                outputfile << this->myBoundaryEvent->getName() << " R dot V (" << framestring << "): " << this->RdotV _GETVALUE << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG