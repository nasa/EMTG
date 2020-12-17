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

#include "BoundaryVelocityVFPAconstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            BoundaryVelocityVFPAconstraint::BoundaryVelocityVFPAconstraint(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                ReferenceFrame constraint_reference_frame,
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
                this->myReferenceFrame = constraint_reference_frame;
            }//end constructor

            void BoundaryVelocityVFPAconstraint::calcbounds()
            {
                // Step 1: parse the constraint definition
				// constraintDefinition is like p0_departure_VFPA_84.0_85.0_J2000BCF
                //remember that we are actually constraining cos(VFPA)
                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                //Step 2: lower bound
				//Note that the lower bound of cosine(VFPA) is the upper bound of VFPA, hence the switch of the indices
                this->Flowerbounds->push_back(cos(std::stod(ConstraintDefinitionCell[4]) * math::deg2rad)); // upper bound of angle, lower bound of cos(angle)

                //Step 3: upper bound
                this->Fupperbounds->push_back(cos(std::stod(ConstraintDefinitionCell[3]) * math::deg2rad)); // lower bound of angle, upper bound of cos(angle)
                this->Fdescriptions->push_back(prefix + "cosine of vertical flight path angle");

                //Step 5: sparsity pattern
                //derivatives with respect to anything influencing the boundary event's right-hand velocity vector
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();

                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
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

                //derivatives with respect to time variables that affect current epoch
                std::vector<size_t> timeVariables = this->myBoundaryEvent->get_Xindices_EventRightEpoch();
                for (size_t Xindex : timeVariables)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        Xindex,
                        this->Gindex_constraint_wrt_time_variables);
                }
            }//end calcbounds()

            void BoundaryVelocityVFPAconstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {
                //Step 1: get the state in ICRF
                math::Matrix<doubleType>& StateBeforeEventICRF = this->myBoundaryEvent->get_state_before_event();

                //Step 2: convert to ConstraintFrame
                this->myUniverse->LocalFrame.construct_rotation_matrices(StateBeforeEventICRF(7), needG);
                math::Matrix<doubleType> R_from_ICRF_to_ConstraintFrame = this->myUniverse->LocalFrame.get_R(ReferenceFrame::ICRF, this->myReferenceFrame);
                this->PositionBeforeEventConstraintFrame = R_from_ICRF_to_ConstraintFrame * StateBeforeEventICRF.getSubMatrix1D(0, 2);
                this->VelocityBeforeEventConstraintFrame = R_from_ICRF_to_ConstraintFrame * StateBeforeEventICRF.getSubMatrix1D(3, 5);


                //Step 3: compute VFPA
                doubleType rdotv_ConstraintFrame = this->PositionBeforeEventConstraintFrame.dot(this->VelocityBeforeEventConstraintFrame);
                doubleType cosVFPA = rdotv_ConstraintFrame / (this->PositionBeforeEventConstraintFrame.norm() * this->VelocityBeforeEventConstraintFrame.norm() + math::SMALL);
                this->VFPA = math::safe_acos(cosVFPA);

                F[Findex++] = cosVFPA;

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

                    //Step 4.2: derivatives of longitude w.r.t. CF states
                    double r = this->PositionBeforeEventConstraintFrame.norm() _GETVALUE;
                    double v = this->VelocityBeforeEventConstraintFrame.norm() _GETVALUE;
                    double rdotv = rdotv_ConstraintFrame _GETVALUE;
                    double rx = this->PositionBeforeEventConstraintFrame(0) _GETVALUE;
                    double ry = this->PositionBeforeEventConstraintFrame(1) _GETVALUE;
                    double rz = this->PositionBeforeEventConstraintFrame(2) _GETVALUE;
                    double vx = this->VelocityBeforeEventConstraintFrame(0) _GETVALUE;
                    double vy = this->VelocityBeforeEventConstraintFrame(1) _GETVALUE;
                    double vz = this->VelocityBeforeEventConstraintFrame(2) _GETVALUE;

                    //component derivatives in the constraint frame
                    double dr_drx = rx / r;
                    double dr_dry = ry / r;
                    double dr_drz = rz / r;
                    double dv_dvx = vx / v;
                    double dv_dvy = vy / v;
                    double dv_dvz = vz / v;
                    double dcosVFPA_drx = vx / (r*v) - rdotv * dr_drx / (r*r * v);
                    double dcosVFPA_dry = vy / (r*v) - rdotv * dr_dry / (r*r * v);
                    double dcosVFPA_drz = vz / (r*v) - rdotv * dr_drz / (r*r * v);
                    double dcosVFPA_dvx = rx / (r*v) - rdotv * dv_dvx / (r*v*v);
                    double dcosVFPA_dvy = ry / (r*v) - rdotv * dv_dvy / (r*v*v);
                    double dcosVFPA_dvz = rz / (r*v) - rdotv * dv_dvz / (r*v*v);

                    //Step 4.3: derivatives with respect to non-time variables affecting state Before boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent();

                    //position
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        double drxConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(0, stateIndex) _GETVALUE;
                        double dryConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(1, stateIndex) _GETVALUE;
                        double drzConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(2, stateIndex) _GETVALUE;

                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_StateBeforeEvent[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_StateBeforeEvent[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_StateBeforeEvent_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (  dcosVFPA_drx * drxConstraintFrame_dstateICRF
                                             + dcosVFPA_dry * dryConstraintFrame_dstateICRF
                                             + dcosVFPA_drz * drzConstraintFrame_dstateICRF)
                                * std::get<2>(Derivatives_of_StateBeforeEvent[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }
                    //velocity
                    for (size_t stateIndex = 3; stateIndex < 6; ++stateIndex)
                    {
                        double dvxConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(0, stateIndex - 3) _GETVALUE;
                        double dvyConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(1, stateIndex - 3) _GETVALUE;
                        double dvzConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(2, stateIndex - 3) _GETVALUE;

                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_StateBeforeEvent[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_StateBeforeEvent[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_StateBeforeEvent_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (  dcosVFPA_dvx * dvxConstraintFrame_dstateICRF
                                             + dcosVFPA_dvy * dvyConstraintFrame_dstateICRF
                                             + dcosVFPA_dvz * dvzConstraintFrame_dstateICRF)
                                * std::get<2>(Derivatives_of_StateBeforeEvent[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }

                    //Step 4.3: derivatives with respect to time variables affecting state Before boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();

                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        double drxConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(0, stateIndex) _GETVALUE;
                        double dryConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(1, stateIndex) _GETVALUE;
                        double drzConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(2, stateIndex) _GETVALUE;
                        double dvxConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(0, stateIndex) _GETVALUE;
                        double dvyConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(1, stateIndex) _GETVALUE;
                        double dvzConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(2, stateIndex) _GETVALUE;

                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_StateBeforeEvent_wrt_Time[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_StateBeforeEvent_wrt_Time[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_StateBeforeEvent_time_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (  dcosVFPA_drx * drxConstraintFrame_dstateICRF
                                             + dcosVFPA_dry * dryConstraintFrame_dstateICRF
                                             + dcosVFPA_drz * drzConstraintFrame_dstateICRF
                                             + dcosVFPA_dvx * dvxConstraintFrame_dstateICRF
                                             + dcosVFPA_dvy * dvyConstraintFrame_dstateICRF
                                             + dcosVFPA_dvz * dvzConstraintFrame_dstateICRF)
                                * std::get<2>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }


                    //Step 4.4: derivatives of VFPA with respect to frame rotation due to change in epoch
                    math::Matrix<doubleType> dR_from_ICRF_to_ConstraintFrame_dt = this->myUniverse->LocalFrame.get_dRdt(ReferenceFrame::ICRF, this->myReferenceFrame);
                    this->dPositionBeforeEventConstraintFrame_dt = dR_from_ICRF_to_ConstraintFrame_dt * StateBeforeEventICRF.getSubMatrix1D(0, 2);
                    this->dVelocityBeforeEventConstraintFrame_dt = dR_from_ICRF_to_ConstraintFrame_dt * StateBeforeEventICRF.getSubMatrix1D(3, 5);

                    double timeDerivative = (  dcosVFPA_drx * this->dPositionBeforeEventConstraintFrame_dt(0)
                                             + dcosVFPA_dry * this->dPositionBeforeEventConstraintFrame_dt(1)
                                             + dcosVFPA_drz * this->dPositionBeforeEventConstraintFrame_dt(2) 
                                             + dcosVFPA_dvx * this->dVelocityBeforeEventConstraintFrame_dt(0)
                                             + dcosVFPA_dvy * this->dVelocityBeforeEventConstraintFrame_dt(1)
                                             + dcosVFPA_dvz * this->dVelocityBeforeEventConstraintFrame_dt(2)) _GETVALUE;

                    for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_time_variables.size(); ++entryIndex)
                    {
                        size_t Gindex = this->Gindex_constraint_wrt_time_variables[entryIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * timeDerivative;
                    }
                }//end derivatives
            }//end process_constraint()

            void BoundaryVelocityVFPAconstraint::output(std::ofstream& outputfile)
            {
                std::string framestring = this->myUniverse->central_body_name + " " + this->ReferenceFrameStrings[this->myReferenceFrame];

                outputfile << this->myBoundaryEvent->getName() << " vertical flight path angle (degrees, " << framestring << "): " << this->VFPA _GETVALUE / math::deg2rad  << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG