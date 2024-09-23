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

#include <algorithm>
#include "BoundaryTargetBodyDeticElevationConstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            BoundaryTargetBodyDeticElevationConstraint::BoundaryTargetBodyDeticElevationConstraint(const std::string& name,
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
                this->r_cb2target_ICRF = math::Matrix<doubleType>(3, 1, 0.0);
                this->dr_cb2target_dt_ICRF = math::Matrix<double>(3, 1, 0.0);
                this->r_cb2target_constraint_frame = math::Matrix<doubleType>(3, 1, 0.0);
                this->r_sc2target_ICRF = math::Matrix<doubleType>(3, 1, 0.0);
                this->r_sc2target_constraint_frame = math::Matrix<doubleType>(3, 1, 0.0);
                this->surface_normal_vec = math::Matrix<doubleType> (3, 1, 0.0);

                this->myReferenceFrame = constraint_reference_frame;
            }//end constructor

            void BoundaryTargetBodyDeticElevationConstraint::calcbounds()
            {
                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                if (std::stoi(ConstraintDefinitionCell[3]) > this->myUniverse->bodies.size())
                {
                    throw std::invalid_argument("The target body universe file index (" + ConstraintDefinitionCell[3]  + ") for this detic elevation constraint is invalid. Constraint \"" + constraintDefinition + "\" in Journey " + std::to_string(journeyIndex));
                }

                std::string out1 = boost::algorithm::to_lower_copy(this->target_body_name);
                std::string out2 = boost::algorithm::to_lower_copy(this->myUniverse->central_body_name);

                if (this->myUniverse->central_body_name == this->target_body_name)
                {
                    throw std::invalid_argument("The target body name: " + this->target_body_name + " is the same as this Journey's central body: " + this->myUniverse->central_body_name + ". Constraint \"" + constraintDefinition + "\" in Journey " + std::to_string(journeyIndex));
                }

                int bodyIndex = std::stoi(ConstraintDefinitionCell[3]) - 1;
                this->target_body = &this->myUniverse->bodies[bodyIndex];
                this->target_body_name = this->target_body->name;

                this->Flowerbounds->push_back(sin(std::stod(ConstraintDefinitionCell[4]) * math::PI / 180.0));
                this->Fupperbounds->push_back(sin(std::stod(ConstraintDefinitionCell[5]) * math::PI / 180.0));
                this->Fdescriptions->push_back(prefix + this->target_body_name + " detic elevation constraint");

                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();

                // create sparsity pattern
                for (size_t stateIndex : {0, 1, 2})
                {
                    //non-time variables
                    std::vector<size_t> state_dIndex_with_respect_to_PositionBeforeEvent;
                    std::vector<size_t> state_Gindex_constraint_wrt_PositionBeforeEvent_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforeEvent.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateBeforeEvent[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_PositionBeforeEvent.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateBeforeEvent[dIndex]),
                                state_Gindex_constraint_wrt_PositionBeforeEvent_variables);
                        }
                    }
                    this->dIndex_with_respect_to_PositionBeforeEvent.push_back(state_dIndex_with_respect_to_PositionBeforeEvent);
                    this->Gindex_constraint_wrt_PositionBeforeEvent_variables.push_back(state_Gindex_constraint_wrt_PositionBeforeEvent_variables);

                    //time variables
                    std::vector<size_t> state_dIndex_with_respect_to_PositionBeforeEvent_wrt_Time;
                    std::vector<size_t> state_Gindex_constraint_wrt_PositionBeforeEvent_time_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_PositionBeforeEvent_wrt_Time.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]),
                                state_Gindex_constraint_wrt_PositionBeforeEvent_time_variables);
                        }
                    }
                    this->dIndex_with_respect_to_PositionBeforeEvent_wrt_Time.push_back(state_dIndex_with_respect_to_PositionBeforeEvent_wrt_Time);
                    this->Gindex_constraint_wrt_PositionBeforeEvent_time_variables.push_back(state_Gindex_constraint_wrt_PositionBeforeEvent_time_variables);
                }

                // derivatives with respect to time variables that affect target body position
                // the target body is affected by ALL time variables (note that future phase flight time variables
                // have not been constructed yet, so the code below will correctly NOT include them)
                std::vector<size_t> timeVariables = this->myBoundaryEvent->get_Xindices_EventRightEpoch();
                for (size_t Xindex : timeVariables)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        Xindex,
                        this->Gindex_constraint_wrt_time_variables);
                }
                
            }//end calcbounds()

            void BoundaryTargetBodyDeticElevationConstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {
                //Step 1: get the state in ICRF
                math::Matrix<doubleType>& StateBeforeEventICRF = this->myBoundaryEvent->get_state_before_event();
                math::Matrix<doubleType> r_cb2sc_ICRF = StateBeforeEventICRF.getSubMatrix1D(0, 2);

                //Step 2: convert to constraint frame
                this->myUniverse->LocalFrame.construct_rotation_matrices(StateBeforeEventICRF(7), needG);
                math::Matrix<doubleType> R_from_ICRF_to_ConstraintFrame = this->myUniverse->LocalFrame.get_R(ReferenceFrame::ICRF, this->myReferenceFrame);
				this->r_cb2sc_constraint_frame = R_from_ICRF_to_ConstraintFrame * r_cb2sc_ICRF;

                //Step 3: compute deticElevation

                // compute the zenith angle: angle between spheroid surface normal and the target body of interest
				
                // We don't actually know what spheroid we're on yet
                // We need to determine that from the state vector and the flattening factor
                // This assumes that we are on some spheroid of unknown dimension, but with a flattening factor
                // defined by the universe central body
                double flattening = this->myUniverse->central_body.flattening_coefficient;
                double f2 = (1.0 - flattening) * (1.0 - flattening);
                doubleType rx = this->r_cb2sc_constraint_frame(0);
                doubleType ry = this->r_cb2sc_constraint_frame(1);
                doubleType rz = this->r_cb2sc_constraint_frame(2);
                doubleType spheroid_major_axis = sqrt(rx*rx + ry*ry + rz*rz / f2);
                doubleType spheroid_major_axis2 = spheroid_major_axis * spheroid_major_axis;

                // now compute the vector normal to the surface of the spheroid
                doubleType nx = 2.0 * rx / (spheroid_major_axis2);
                doubleType ny = 2.0 * ry / (spheroid_major_axis2);
                doubleType nz = 2.0 * rz / (spheroid_major_axis2 * f2);
                doubleType nnorm = sqrt(nx * nx + ny * ny + nz * nz);
                this->surface_normal_vec(0) = nx;
                this->surface_normal_vec(1) = ny;
                this->surface_normal_vec(2) = nz;

                // get the position vector of the target of interest w.r.t. the central body
                doubleType temp_body_state[12];
                this->target_body->locate_body(StateBeforeEventICRF(7),
                    temp_body_state,
                    needG,
                    *this->myOptions);

                for (size_t stateIndex : {0, 1, 2})
                {
                    this->r_cb2target_ICRF(stateIndex) = temp_body_state[stateIndex];
                    this->dr_cb2target_dt_ICRF(stateIndex) = temp_body_state[stateIndex + 6] _GETVALUE;
                }

                this->r_cb2target_constraint_frame = R_from_ICRF_to_ConstraintFrame * this->r_cb2target_ICRF;

                // get the position vector of the target of interest w.r.t. the spacecraft in ICRF
                this->r_sc2target_ICRF = this->r_cb2target_ICRF - r_cb2sc_ICRF;

                // get the position vector of the target of interest w.r.t. the spacecraft in the constraint frame
                this->r_sc2target_constraint_frame = R_from_ICRF_to_ConstraintFrame * this->r_sc2target_ICRF;

                // compute the angle between the surface normal vector and the spacecraft-target vector
                doubleType cos_zenith_angle = this->surface_normal_vec.dot(this->r_sc2target_constraint_frame) / (this->surface_normal_vec.norm() * this->r_sc2target_constraint_frame.norm());
                
                // trig. identity: cos(zenith) = sin(90 - zenith) = sin(elevation)
                this->sin_detic_elevation = cos_zenith_angle;
                this->detic_elevation = math::safe_asin(this->sin_detic_elevation);

                // populate NLP constraint vector
                F[Findex++] = this->sin_detic_elevation;

                // Step 4: compute partial derivatives if required
                if (needG)
                {   
                    // We are going to be adding components on top of components, so let's start by zeroing out all of our G entries
                    for (size_t stateIndex : {0, 1, 2})
                    {
                        for (size_t Gindex : this->Gindex_constraint_wrt_PositionBeforeEvent_variables[stateIndex])
                        {
                            G[Gindex] = 0.0;
                        }
                        for (size_t Gindex : this->Gindex_constraint_wrt_PositionBeforeEvent_time_variables[stateIndex])
                        {
                            G[Gindex] = 0.0;
                        }
                    }
                    for (size_t Gindex : this->Gindex_constraint_wrt_time_variables)
                    {
                        G[Gindex] = 0.0;
                    }

                    double x = rx _GETVALUE;
                    double y = ry _GETVALUE;
                    double z = rz _GETVALUE;
                    double x2 = x * x;
                    double y2 = y * y;
                    double z2 = z * z;
                    double xt = this->r_cb2target_constraint_frame(0) _GETVALUE;
                    double yt = this->r_cb2target_constraint_frame(1) _GETVALUE;
                    double zt = this->r_cb2target_constraint_frame(2) _GETVALUE;
                    double x_s2t = this->r_sc2target_constraint_frame(0) _GETVALUE;
                    double y_s2t = this->r_sc2target_constraint_frame(1) _GETVALUE;
                    double z_s2t = this->r_sc2target_constraint_frame(2) _GETVALUE;
                    double r_s2t = sqrt(x_s2t * x_s2t + y_s2t * y_s2t + z_s2t * z_s2t);
                    double fminus12 = (flattening - 1.0) * (flattening - 1.0);
                    double fminus14 = fminus12 * fminus12;
                    double a = spheroid_major_axis _GETVALUE;
                    double a2 = (spheroid_major_axis * spheroid_major_axis) _GETVALUE;
                    double a4 = a2 * a2;

                    // constraint partials w.r.t. spacecraft position (w.r.t. central body) in the constraint frame
                    double d_selev_dx = x * fminus12 * (z*(z - zt) + fminus12 * (x*(x - xt) + y * (y - yt))) / (a2*sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14)) *            (z2 + fminus14 * (x2 + y2))*r_s2t) - (2 * x - xt)                                                    / (a2*sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14))*r_s2t) - (-x + xt)*(z*(z - zt) + fminus12 * (x*(x - xt) + y * (y - yt))) / (a2*sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14))*fminus12*pow(r_s2t, 3.0));
                    double d_selev_dy = y * fminus12 * (z*(z - zt) + fminus12 * (x*(x - xt) + y * (y - yt))) / (a2*sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14)) *            (z2 + fminus14 * (x2 + y2))*r_s2t) - (2 * y - yt)                                                    / (a2*sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14))*r_s2t) - (-y + yt)*(z*(z - zt) + fminus12 * (x*(x - xt) + y * (y - yt))) / (a2*sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14))*fminus12*pow(r_s2t, 3.0));
                    double d_selev_dz = z            * (z*(z - zt) + fminus12 * (x*(x - xt) + y * (y - yt))) / (a2*sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14)) * fminus12 * (z2 + fminus14 * (x2 + y2))*r_s2t) - (-z + zt)*(z*(z - zt) + fminus12 * (x*(x - xt) + y * (y - yt))) / (a2*sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14))*fminus12*pow(r_s2t, 3.0)) - (2 * z - zt)                                 / (a2*sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14))*fminus12*r_s2t);

                    // all zero...the spheroid major axis does depend on the s/c position, but the partial of elevation w.r.t. the spheroid major axis is zero
                    // double d_selev_da
                    // double d_selev_da
                    // double d_selev_da

                    // constraint partials w.r.t. target body position (w.r.t. central body) in the constraint frame
                    double d_selev_dxt = x / (a2 * sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14))*r_s2t)          - (x - xt)*(z*(z - zt) + fminus12 * (x*(x - xt) + y * (y - yt))) / (a2*sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14))*fminus12*pow(r_s2t, 3.0));
                    double d_selev_dyt = y / (a2 * sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14))*r_s2t)          - (y - yt)*(z*(z - zt) + fminus12 * (x*(x - xt) + y * (y - yt))) / (a2*sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14))*fminus12*pow(r_s2t, 3.0));
                    double d_selev_dzt = z / (a2 * sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14))*fminus12*r_s2t) - (z - zt)*(z*(z - zt) + fminus12 * (x*(x - xt) + y * (y - yt))) / (a2*sqrt((z2 + fminus14 * (x2 + y2)) / (a4*fminus14))*fminus12*pow(r_s2t, 3.0));
                    

                    // derivatives of the spacecraft state w.r.t. non-time variables affecting position before boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent();
                    //position
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        double drxConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(0, stateIndex) _GETVALUE;
                        double dryConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(1, stateIndex) _GETVALUE;
                        double drzConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(2, stateIndex) _GETVALUE;

                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_PositionBeforeEvent[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_PositionBeforeEvent[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_PositionBeforeEvent_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (  d_selev_dx * drxConstraintFrame_dstateICRF
                                             + d_selev_dy * dryConstraintFrame_dstateICRF
                                             + d_selev_dz * drzConstraintFrame_dstateICRF)
                                             * std::get<2>(Derivatives_of_StateBeforeEvent[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }

                    // derivatives of the spacecraft position w.r.t. time variables affecting position before boundary event
                    // these are the implicit time partials of the spacecraft state
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        double drxConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(0, stateIndex) _GETVALUE;
                        double dryConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(1, stateIndex) _GETVALUE;
                        double drzConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(2, stateIndex) _GETVALUE;

                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_PositionBeforeEvent_wrt_Time[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_PositionBeforeEvent_wrt_Time[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_PositionBeforeEvent_time_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (  d_selev_dx * drxConstraintFrame_dstateICRF
                                             + d_selev_dy * dryConstraintFrame_dstateICRF
                                             + d_selev_dz * drzConstraintFrame_dstateICRF)
                                             * std::get<2>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }

                    // explicit time partials of the constraint
                    // the partials of this constraint w.r.t. epoch are the same for ALL previous and current flight time variables
                    math::Matrix<doubleType> dR_from_ICRF_to_ConstraintFrame_dt = this->myUniverse->LocalFrame.get_dRdt(ReferenceFrame::ICRF, this->myReferenceFrame);

                    // contributions due to target body motion
                    // let X be the position vector of the target body w.r.t. the central body
                    // dFdX_BCF * (dX_BCF/dX_ICRF * dX_ICRF/dt + d/dt(dX_BCF/dX_ICRF) * X_ICRF)
                    double time_derivative = (d_selev_dxt * (R_from_ICRF_to_ConstraintFrame(0, 0) * this->dr_cb2target_dt_ICRF(0)
                                                          + R_from_ICRF_to_ConstraintFrame(0, 1) * this->dr_cb2target_dt_ICRF(1)
                                                          + R_from_ICRF_to_ConstraintFrame(0, 2) * this->dr_cb2target_dt_ICRF(2)
                                                          + dR_from_ICRF_to_ConstraintFrame_dt(0, 0) * this->r_cb2target_ICRF(0)
                                                          + dR_from_ICRF_to_ConstraintFrame_dt(0, 1) * this->r_cb2target_ICRF(1)
                                                          + dR_from_ICRF_to_ConstraintFrame_dt(0, 2) * this->r_cb2target_ICRF(2))
                                           + d_selev_dyt * (R_from_ICRF_to_ConstraintFrame(1, 0) * this->dr_cb2target_dt_ICRF(0)
                                                          + R_from_ICRF_to_ConstraintFrame(1, 1) * this->dr_cb2target_dt_ICRF(1)
                                                          + R_from_ICRF_to_ConstraintFrame(1, 2) * this->dr_cb2target_dt_ICRF(2)
                                                          + dR_from_ICRF_to_ConstraintFrame_dt(1, 0) * this->r_cb2target_ICRF(0)
                                                          + dR_from_ICRF_to_ConstraintFrame_dt(1, 1) * this->r_cb2target_ICRF(1)
                                                          + dR_from_ICRF_to_ConstraintFrame_dt(1, 2) * this->r_cb2target_ICRF(2))
                                           + d_selev_dzt * (R_from_ICRF_to_ConstraintFrame(2, 0) * this->dr_cb2target_dt_ICRF(0)
                                                          + R_from_ICRF_to_ConstraintFrame(2, 1) * this->dr_cb2target_dt_ICRF(1)
                                                          + R_from_ICRF_to_ConstraintFrame(2, 2) * this->dr_cb2target_dt_ICRF(2)
                                                          + dR_from_ICRF_to_ConstraintFrame_dt(2, 0) * this->r_cb2target_ICRF(0)
                                                          + dR_from_ICRF_to_ConstraintFrame_dt(2, 1) * this->r_cb2target_ICRF(1)
                                                          + dR_from_ICRF_to_ConstraintFrame_dt(2, 2) * this->r_cb2target_ICRF(2))) _GETVALUE;
                    
                    // contributions due to spacecraft motion
                    // let X be the position vector of the spacecraft w.r.t. the central body
                    // dFdX_BCF * (dX_BCF/dX_ICRF * dX_ICRF/dt + d/dt(dX_BCF/dX_ICRF) * X_ICRF)
                    // dX_ICRF/dt = 0, the spacecraft position does not have any explicit time partials and implicit time partials are handled above
                    time_derivative += (d_selev_dx * (dR_from_ICRF_to_ConstraintFrame_dt(0, 0) * r_cb2sc_ICRF(0)
                                                    + dR_from_ICRF_to_ConstraintFrame_dt(0, 1) * r_cb2sc_ICRF(1)
                                                    + dR_from_ICRF_to_ConstraintFrame_dt(0, 2) * r_cb2sc_ICRF(2))
                                      + d_selev_dy * (dR_from_ICRF_to_ConstraintFrame_dt(1, 0) * r_cb2sc_ICRF(0)
                                                    + dR_from_ICRF_to_ConstraintFrame_dt(1, 1) * r_cb2sc_ICRF(1)
                                                    + dR_from_ICRF_to_ConstraintFrame_dt(1, 2) * r_cb2sc_ICRF(2))
                                      + d_selev_dz * (dR_from_ICRF_to_ConstraintFrame_dt(2, 0) * r_cb2sc_ICRF(0)
                                                    + dR_from_ICRF_to_ConstraintFrame_dt(2, 1) * r_cb2sc_ICRF(1)
                                                    + dR_from_ICRF_to_ConstraintFrame_dt(2, 2) * r_cb2sc_ICRF(2))) _GETVALUE;

                    for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_time_variables.size(); ++entryIndex)
                    {
                        size_t Gindex = this->Gindex_constraint_wrt_time_variables[entryIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                     * time_derivative;
                    }

                }//end derivatives
            }//end process_constraint()

            void BoundaryTargetBodyDeticElevationConstraint::output(std::ofstream& outputfile)
            {
                std::string framestring = this->myUniverse->central_body_name + " " + this->ReferenceFrameStrings[this->myReferenceFrame];

                outputfile << this->myBoundaryEvent->getName() << " detic elevation of " << this->target_body_name << " (deg, " << framestring << "): " << this->detic_elevation _GETVALUE * 180.0 / math::PI << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG
