// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2024 United States Government as represented by the
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
#include "BoundaryDeticElevationFromGroundStationConstraint.h"
#include "BodydeticConversions.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            BoundaryDeticElevationFromGroundStationConstraint::BoundaryDeticElevationFromGroundStationConstraint(const std::string& name,
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
                this->r_cb2sc_ICRF = math::Matrix<doubleType>(3, 1, 0.0);
                this->r_cb2gsbody_ICRF = math::Matrix<doubleType>(3, 1, 0.0);
                this->d_r_cb2gsbody_dt_ICRF = math::Matrix<double>(3, 1, 0.0);
                this->r_cb2gs_ICRF = math::Matrix<doubleType>(3, 1, 0.0);
                this->r_gsbody2gs_ICRF = math::Matrix<doubleType>(3, 1, 0.0);
                this->r_gsbody2gs_BCF = math::Matrix<doubleType>(3, 1, 0.0);
                this->r_gs2sc_ICRF = math::Matrix<doubleType>(3, 1, 0.0);

                this->surface_normal_vec = math::Matrix<doubleType> (3, 1, 0.0);

                this->myReferenceFrame = ReferenceFrame::TrueOfDate_BCF;
            }//end constructor

            void BoundaryDeticElevationFromGroundStationConstraint::calcbounds()
            {
                // parse the constraint definition
                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                // figure out what body the ground station is located on
                if (boost::algorithm::to_lower_copy(ConstraintDefinitionCell[3]) == "cb")
                {
                    this->ground_station_body = &this->myUniverse->central_body;
                    this->ground_station_body_name = this->ground_station_body->name;
                }
                else
                {
                    bool found_index = false;
                    for (auto const& body : this->myUniverse->bodies)
                    {
                        if (body.body_code == std::stoi(ConstraintDefinitionCell[3]))
                        {
                            found_index = true;
                        }
                    }
                    if (!found_index)
                    {
                        throw std::invalid_argument("The specified ground station body universe file index (" + ConstraintDefinitionCell[3] + ") for this detic elevation constraint is invalid. Constraint \"" + constraintDefinition + "\" in Journey " + std::to_string(journeyIndex));
                    }
                    else
                    {
                        int bodyIndex = std::stoi(ConstraintDefinitionCell[3]) - 1;
                        this->ground_station_body = &this->myUniverse->bodies[bodyIndex];
                        this->ground_station_body_name = this->ground_station_body->name;
                    }
                }

                // parse the lat./lon./alt.
                std::string latitude_string = ConstraintDefinitionCell[4];
                std::string longitude_string = ConstraintDefinitionCell[5];
                std::string altitude_string = ConstraintDefinitionCell[6];

                // this does not need to be a doubleType, but the conversion function requires that data type as an input
                math::Matrix<doubleType> LLA = math::Matrix<doubleType>(3, 1, 0.0);
                if (latitude_string.find("deg") < 1024)
                {
                    this->ground_station_latitude = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[4], "deg"));
                    LLA(0) = this->ground_station_latitude * math::deg2rad;
                }
                else
                {
                    throw std::invalid_argument("Please specify the ground station latitude units in degrees. Constraint \"" + constraintDefinition + "\" in Journey " + std::to_string(journeyIndex) + " Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
                if (longitude_string.find("deg") < 1024)
                {
                    this->ground_station_longitude = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[5], "deg"));
                    LLA(1) = this->ground_station_longitude * math::deg2rad;
                }
                else
                {
                    throw std::invalid_argument("Please specify the ground station longitude units in degrees. Constraint \"" + constraintDefinition + "\" in Journey " + std::to_string(journeyIndex) + " Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
                if (altitude_string.find("km") < 1024)
                {
                    this->ground_station_altitude = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[6], "km"));
                    LLA(2) = this->ground_station_altitude;
                }
                else
                {
                    throw std::invalid_argument("Please specify the ground station altitude units in km. Constraint \"" + constraintDefinition + "\" in Journey " + std::to_string(journeyIndex) + " Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }

                // compute the Cartesian position of the ground station in the body-fixed frame
                // the LLA2BCF_oblate function does not yet populate a partials matrix, but it requires one as an input 
                math::Matrix<doubleType> dBCFdLLA = math::Matrix<doubleType>(3, 3, 0.0);
                Astrodynamics::LLA2BCF_oblate(LLA,
                                              this->ground_station_body->radius,
                                              this->ground_station_body->flattening_coefficient,
                                              this->r_gsbody2gs_BCF,
                                              false,
                                              dBCFdLLA);

                // constraint bounds 
                this->Flowerbounds->push_back(sin(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[7], "deg")) * math::PI / 180.0));
                this->Fupperbounds->push_back(sin(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[8], "deg")) * math::PI / 180.0));
                this->Fdescriptions->push_back(prefix + this->ground_station_body_name + " ground station "
                                               + "lat: " + boost::erase_all_copy(ConstraintDefinitionCell[4], "deg") + " deg. "
                                               + "lon: " + boost::erase_all_copy(ConstraintDefinitionCell[5], "deg") + " deg. "
                                               + "alt: " + boost::erase_all_copy(ConstraintDefinitionCell[6], "km")  + " km "
                                               + "detic elevation constraint");

                // create sparsity pattern
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();

                // 
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

            void BoundaryDeticElevationFromGroundStationConstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {

                // Get the spacecraft state in ICRF
                math::Matrix<doubleType>& StateBeforeEventICRF = this->myBoundaryEvent->get_state_before_event();
                this->r_cb2sc_ICRF = StateBeforeEventICRF.getSubMatrix1D(0, 2);                

                // Get the position of the ground station body w.r.t. the central body
                // Check if the ground station is on the central body, or if it's on another body in the Universe
                if (this->ground_station_body_name == this->myUniverse->central_body_name)
                {
                    // the ground station body is the central body, so their relative position is zero
                    this->r_cb2gsbody_ICRF.assign_zeros();
                    this->d_r_cb2gsbody_dt_ICRF.assign_zeros();
                }
                else
                {
                    // get the position vector of the target of interest w.r.t. the central body
                    doubleType temp_body_state[12];
                    this->ground_station_body->locate_body(StateBeforeEventICRF(7),
                        temp_body_state,
                        needG,
                        *this->myOptions);

                    for (size_t stateIndex : {0, 1, 2})
                    {
                        this->r_cb2gsbody_ICRF(stateIndex) = temp_body_state[stateIndex];
                        this->d_r_cb2gsbody_dt_ICRF(stateIndex) = temp_body_state[stateIndex + 6] _GETVALUE;
                    }
                }

                // Compute the position of the ground station in ICRF
                this->ground_station_body->body_frame.construct_rotation_matrices(StateBeforeEventICRF(7), needG);
                //math::Matrix<doubleType> R_from_BCF_to_ICRF = this->ground_station_body->body_frame.get_R(ReferenceFrame::TrueOfDate_BCF, ReferenceFrame::ICRF);
                //math::Matrix<doubleType> R_from_ICRF_to_BCF = this->ground_station_body->body_frame.get_R(ReferenceFrame::ICRF, ReferenceFrame::TrueOfDate_BCF);
                math::Matrix<doubleType> R_from_BCF_to_ICRF = this->ground_station_body->body_frame.get_R_from_BCF_to_ICRF();
                math::Matrix<doubleType> R_from_ICRF_to_BCF = this->ground_station_body->body_frame.get_R_from_ICRF_to_BCF();
                this->r_gsbody2gs_ICRF = R_from_BCF_to_ICRF * this->r_gsbody2gs_BCF;

                // Compute the position of the spacecraft w.r.t. the ground station in ICRF
                this->r_gs2sc_ICRF = this->r_cb2sc_ICRF + this->r_cb2gsbody_ICRF - this->r_gsbody2gs_ICRF;
                
                //DEBUG
                math::Matrix<doubleType> r_cb2sc_BCF = R_from_ICRF_to_BCF * this->r_cb2sc_ICRF;

                // Compute the position of the spacecraft w.r.t. the ground station in TOD BCF
                this->r_gs2sc_BCF = R_from_ICRF_to_BCF * this->r_gs2sc_ICRF;
                doubleType r_gs2sc_BCF_norm = this->r_gs2sc_BCF.norm();

                // Compute the zenith vector at the location of the ground station 
                // First, compute the spheroid major axis size
                double flattening = this->ground_station_body->flattening_coefficient;
                double f2 = (1.0 - flattening) * (1.0 - flattening);
                doubleType rx = this->r_gsbody2gs_BCF(0);
                doubleType ry = this->r_gsbody2gs_BCF(1);
                doubleType rz = this->r_gsbody2gs_BCF(2);

                //TODO: this should just be the ground_station_body's equatorial radius
                doubleType spheroid_major_axis = sqrt(rx * rx + ry * ry + rz * rz / f2);
                doubleType spheroid_major_axis2 = spheroid_major_axis * spheroid_major_axis;

                // now compute the vector normal to the surface of the spheroid
                doubleType nx = 2.0 * rx / (spheroid_major_axis2);
                doubleType ny = 2.0 * ry / (spheroid_major_axis2);
                doubleType nz = 2.0 * rz / (spheroid_major_axis2 * f2);
                this->surface_normal_vec(0) = nx;
                this->surface_normal_vec(1) = ny;
                this->surface_normal_vec(2) = nz;
                doubleType surface_normal_vec_norm = this->surface_normal_vec.norm();

                // compute the angle between the surface normal vector and the ground station to spacecraft vector
                doubleType cos_zenith_angle = this->surface_normal_vec.dot(this->r_gs2sc_BCF) / (surface_normal_vec_norm * r_gs2sc_BCF_norm);

                // trig. identity: cos(zenith) = sin(90 - zenith) = sin(elevation)
                this->sin_detic_elevation = cos_zenith_angle;
                this->detic_elevation = math::safe_asin(this->sin_detic_elevation);

                // populate NLP constraint vector
                F[Findex++] = this->sin_detic_elevation;

                // Compute partial derivatives if required
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


                    // partials of sin(elevation) w.r.t. the position of the s/c w.r.t. the ground station in the constraint frame

                    // position vector components of the ground station w.r.t. the ground station body in the constraint frame
                    doubleType xg_BCF = this->r_gsbody2gs_BCF(0);
                    doubleType yg_BCF = this->r_gsbody2gs_BCF(1);
                    doubleType zg_BCF = this->r_gsbody2gs_BCF(2);

                    // position vector components of the spacecraft w.r.t. the ground station in the constraint frame
                    doubleType xg2sc_BCF = this->r_gs2sc_BCF(0);
                    doubleType yg2sc_BCF = this->r_gs2sc_BCF(1);
                    doubleType zg2sc_BCF = this->r_gs2sc_BCF(2);

                    double a = spheroid_major_axis _GETVALUE;
                    double a2 = a * a;

                    double d_sin_elev_dx_g2sc = 2.0 * (xg_BCF / (a2 * r_gs2sc_BCF_norm * surface_normal_vec_norm) - xg2sc_BCF * (((zg_BCF * zg2sc_BCF) / (a2 * (1.0 - flattening))) + ((xg_BCF * xg2sc_BCF) / (a2)) + ((yg_BCF * yg2sc_BCF) / (a2))) / (surface_normal_vec_norm * pow(r_gs2sc_BCF_norm, 3.0))) _GETVALUE;
                    double d_sin_elev_dy_g2sc = 2.0 * (yg_BCF / (a2 * r_gs2sc_BCF_norm * surface_normal_vec_norm) - yg2sc_BCF * (((zg_BCF * zg2sc_BCF) / (a2 * (1.0 - flattening))) + ((xg_BCF * xg2sc_BCF) / (a2)) + ((yg_BCF * yg2sc_BCF) / (a2))) / (surface_normal_vec_norm * pow(r_gs2sc_BCF_norm, 3.0))) _GETVALUE;
                    double d_sin_elev_dz_g2sc = 2.0 * (zg_BCF / (a2 * (1.0 - flattening) * r_gs2sc_BCF_norm * surface_normal_vec_norm) - zg2sc_BCF * (((zg_BCF * zg2sc_BCF) / (a2 * (1.0 - flattening))) + ((xg_BCF * xg2sc_BCF) / (a2)) + ((yg_BCF * yg2sc_BCF) / (a2))) / (surface_normal_vec_norm * pow(r_gs2sc_BCF_norm, 3.0))) _GETVALUE;
                    
                    // derivatives of the spacecraft state w.r.t. non-time variables affecting state before the boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent();
                    
                    // derivatives of the spacecraft state w.r.t. time variables affecting position before boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();
                    for (size_t stateIndex : {0, 1, 2})
                    {
                        // partials of the gs -> s/c vector in BCF w.r.t. the cb -> s/c vector in ICRF
                        double d_x_gs2sc_BCF_d_r_component_cb2sc_ICRF = R_from_ICRF_to_BCF(0, stateIndex) _GETVALUE;
                        double d_y_gs2sc_BCF_d_r_component_cb2sc_ICRF = R_from_ICRF_to_BCF(1, stateIndex) _GETVALUE;
                        double d_z_gs2sc_BCF_d_r_component_cb2sc_ICRF = R_from_ICRF_to_BCF(2, stateIndex) _GETVALUE;
                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_PositionBeforeEvent[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_PositionBeforeEvent[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_PositionBeforeEvent_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (d_sin_elev_dx_g2sc * d_x_gs2sc_BCF_d_r_component_cb2sc_ICRF
                                           + d_sin_elev_dy_g2sc * d_y_gs2sc_BCF_d_r_component_cb2sc_ICRF
                                           + d_sin_elev_dz_g2sc * d_z_gs2sc_BCF_d_r_component_cb2sc_ICRF) * std::get<2>(Derivatives_of_StateBeforeEvent[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex) * Gentry;
                        }

                        // these are the implicit time partials of the spacecraft state
                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_PositionBeforeEvent_wrt_Time[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_PositionBeforeEvent_wrt_Time[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_PositionBeforeEvent_time_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (d_sin_elev_dx_g2sc * d_x_gs2sc_BCF_d_r_component_cb2sc_ICRF
                                           + d_sin_elev_dy_g2sc * d_y_gs2sc_BCF_d_r_component_cb2sc_ICRF
                                           + d_sin_elev_dz_g2sc * d_z_gs2sc_BCF_d_r_component_cb2sc_ICRF) * std::get<2>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex) * Gentry;
                        }
                    }

                    // explicit time partials of the constraint
                    // the partials of this constraint w.r.t. epoch are the same for ALL previous and current flight time variables
                    math::Matrix<doubleType> dR_from_ICRF_to_BCF_dt = this->ground_station_body->body_frame.get_dRdt(ReferenceFrame::ICRF, ReferenceFrame::TrueOfDate_BCF);

                    // Let M be the rotation matrix from ICRF to BCF
                    // r_gs2sc_BCF = M * r_cb2sc_ICRF + M * r_cb2gbody_ICRF - r_gbody2gs_BCF
                    // d/dt(r_gs2sc_BCF) = dM/dt * r_cb2sc_ICRF 
                    //                   + M     * d/dt(r_cb2sc_ICRF) ---> handled in the for loop above
                    //                   + dM/dt * r_cb2gbody_ICRF
                    //                   + M * d/dt(r_cb2gbody_ICRF)  ---> explicit partials of the ground station body ephemeris
                    //                   - d/dt(r_gbody2gs_BCF)       ---> zero assuming the ground station does not move on the surface of the body

                    // contributions due to spacecraft motion
                    double time_derivative  = (d_sin_elev_dx_g2sc * (dR_from_ICRF_to_BCF_dt(0, 0) * this->r_cb2sc_ICRF(0)
                                                                   + dR_from_ICRF_to_BCF_dt(0, 1) * this->r_cb2sc_ICRF(1)
                                                                   + dR_from_ICRF_to_BCF_dt(0, 2) * this->r_cb2sc_ICRF(2))
                                             + d_sin_elev_dy_g2sc * (dR_from_ICRF_to_BCF_dt(1, 0) * this->r_cb2sc_ICRF(0)
                                                                   + dR_from_ICRF_to_BCF_dt(1, 1) * this->r_cb2sc_ICRF(1)
                                                                   + dR_from_ICRF_to_BCF_dt(1, 2) * this->r_cb2sc_ICRF(2))
                                             + d_sin_elev_dz_g2sc * (dR_from_ICRF_to_BCF_dt(2, 0) * this->r_cb2sc_ICRF(0)
                                                                   + dR_from_ICRF_to_BCF_dt(2, 1) * this->r_cb2sc_ICRF(1)
                                                                   + dR_from_ICRF_to_BCF_dt(2, 2) * this->r_cb2sc_ICRF(2))) _GETVALUE;

                    // contributions due to ground station body motion
                    // these will be zero if the ground station is located on the central body
                    time_derivative += (d_sin_elev_dx_g2sc * (R_from_ICRF_to_BCF(0, 0) * this->d_r_cb2gsbody_dt_ICRF(0)
                                                            + R_from_ICRF_to_BCF(0, 1) * this->d_r_cb2gsbody_dt_ICRF(1)
                                                            + R_from_ICRF_to_BCF(0, 2) * this->d_r_cb2gsbody_dt_ICRF(2)
                                                            + dR_from_ICRF_to_BCF_dt(0, 0) * this->r_cb2gsbody_ICRF(0)
                                                            + dR_from_ICRF_to_BCF_dt(0, 1) * this->r_cb2gsbody_ICRF(1)
                                                            + dR_from_ICRF_to_BCF_dt(0, 2) * this->r_cb2gsbody_ICRF(2))
                                      + d_sin_elev_dy_g2sc * (R_from_ICRF_to_BCF(1, 0) * this->d_r_cb2gsbody_dt_ICRF(0)
                                                            + R_from_ICRF_to_BCF(1, 1) * this->d_r_cb2gsbody_dt_ICRF(1)
                                                            + R_from_ICRF_to_BCF(1, 2) * this->d_r_cb2gsbody_dt_ICRF(2)
                                                            + dR_from_ICRF_to_BCF_dt(1, 0) * this->r_cb2gsbody_ICRF(0)
                                                            + dR_from_ICRF_to_BCF_dt(1, 1) * this->r_cb2gsbody_ICRF(1)
                                                            + dR_from_ICRF_to_BCF_dt(1, 2) * this->r_cb2gsbody_ICRF(2))
                                      + d_sin_elev_dz_g2sc * (R_from_ICRF_to_BCF(2, 0) * this->d_r_cb2gsbody_dt_ICRF(0)
                                                            + R_from_ICRF_to_BCF(2, 1) * this->d_r_cb2gsbody_dt_ICRF(1)
                                                            + R_from_ICRF_to_BCF(2, 2) * this->d_r_cb2gsbody_dt_ICRF(2)
                                                            + dR_from_ICRF_to_BCF_dt(2, 0) * this->r_cb2gsbody_ICRF(0)
                                                            + dR_from_ICRF_to_BCF_dt(2, 1) * this->r_cb2gsbody_ICRF(1)
                                                            + dR_from_ICRF_to_BCF_dt(2, 2) * this->r_cb2gsbody_ICRF(2))) _GETVALUE;

                    for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_time_variables.size(); ++entryIndex)
                    {
                        size_t Gindex = this->Gindex_constraint_wrt_time_variables[entryIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * time_derivative;
                    }
                }
            
            }//end process_constraint()

            void BoundaryDeticElevationFromGroundStationConstraint::output(std::ofstream& outputfile)
            {
                outputfile << this->myBoundaryEvent->getName() 
                           << " detic elevation of spacecraft w.r.t (" << this->ground_station_body_name << " station"
                           << " lat.: " << this->ground_station_latitude << " deg. "
                           << " lon.: " << this->ground_station_longitude << " deg. "
                           << " alt.: " << this->ground_station_altitude << " km): "
                           << this->detic_elevation _GETVALUE * 180.0 / math::PI << " deg." << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG
