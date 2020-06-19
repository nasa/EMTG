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

#include "BoundaryVelocitySphericalAzimuthConstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            BoundaryVelocitySphericalAzimuthConstraint::BoundaryVelocitySphericalAzimuthConstraint(const std::string& name,
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

            void BoundaryVelocitySphericalAzimuthConstraint::calcbounds()
            {
                // Step 1: parse the constraint definition
				// constraintDefinition is like p0_arrival_VelocitySphericalAzimuth_-10.9_-11.0_trueofdatebcf

                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                //Step 2: lower bound
                this->Flowerbounds->push_back(std::stod(ConstraintDefinitionCell[3]) * math::deg2rad);


                //Step 3: upper bound
                this->Fupperbounds->push_back(std::stod(ConstraintDefinitionCell[4]) * math::deg2rad);
                this->Fdescriptions->push_back(prefix + "velocity spherical Azimuth");

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

            void BoundaryVelocitySphericalAzimuthConstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {
                //Step 1: get the state in ICRF
                math::Matrix<doubleType>& StateBeforeEventICRF = this->myBoundaryEvent->get_state_before_event();

                //Step 2: convert to constraint frame
                this->myUniverse->LocalFrame.construct_rotation_matrices(StateBeforeEventICRF(7), needG);
                math::Matrix<doubleType> R_from_ICRF_to_ConstraintFrame = this->myUniverse->LocalFrame.get_R(ReferenceFrame::ICRF, this->myReferenceFrame);
				this->PositionBeforeEventConstraintFrame = R_from_ICRF_to_ConstraintFrame * StateBeforeEventICRF.getSubMatrix1D(0, 2);
				this->VelocityBeforeEventConstraintFrame = R_from_ICRF_to_ConstraintFrame * StateBeforeEventICRF.getSubMatrix1D(3, 5);


                //Step 3: compute velocity Azimuth
				//First, make sure that the probe is not at the pole because there is a singularity at the pole. 
				doubleType rx = this->PositionBeforeEventConstraintFrame(0);
				doubleType ry = this->PositionBeforeEventConstraintFrame(1);
				doubleType rz = this->PositionBeforeEventConstraintFrame(2);
				doubleType rxy = sqrt(rx * rx + ry * ry);
				if (rxy < math::SMALL)
				{
					throw std::runtime_error(this->prefix + "The probe is too close to the pole to calculate longitude and its derivatives.");
				}
				doubleType vx = this->VelocityBeforeEventConstraintFrame(0) ;
				doubleType vy = this->VelocityBeforeEventConstraintFrame(1) ;
				doubleType vz = this->VelocityBeforeEventConstraintFrame(2) ;

				doubleType r = sqrt(rx * rx + ry * ry + rz * rz);

				doubleType lambda = atan2(rz, rxy);
				doubleType longitude = atan2(ry, rx);

				doubleType clat = cos(lambda);
				doubleType slat = sin(lambda);
				doubleType clon = cos(longitude);
				doubleType slon = sin(longitude);

				doubleType vu = clat * clon * vx + clat * slon * vy + slat * vz;
				doubleType ve = -slon * vx + clon * vy;
				doubleType vn = -slat * clon * vx - slat * slon * vy + clat * vz;

				this->VelocitySphericalAzimuth = atan2(ve, vn);

                F[Findex++] = this->VelocitySphericalAzimuth;
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

                    //Step 4.2: derivatives of AZ w.r.t. CF states
					doubleType dr_drx = rx / r;
					doubleType dr_dry = ry / r;
					doubleType dr_drz = rz / r;

					doubleType dlambda_drx = -(rx * pow(rxy, -3.0) * rz) / (rz * rz / rxy / rxy + 1.0);
					doubleType dlambda_dry = -(ry * pow(rxy, -3.0) * rz) / (rz * rz / rxy / rxy + 1.0);
					doubleType dlambda_drz = 1.0 / (rxy * (rz * rz / rxy / rxy + 1.0));

					doubleType dlongitude_drx = -ry / rxy / rxy;
					doubleType dlongitude_dry = rx / rxy / rxy;

					doubleType dvu_drx = (-slat * clon * dlambda_drx - clat * slon * dlongitude_drx) * vx + (-slat * slon * dlambda_drx + clat * clon * dlongitude_drx) * vy + clat * dlambda_drx * vz;
					doubleType dvu_dry = (-slat * clon * dlambda_dry - clat * slon * dlongitude_dry) * vx + (-slat * slon * dlambda_dry + clat * clon * dlongitude_dry) * vy + clat * dlambda_dry * vz;
					doubleType dvu_drz = (-slat * clon * dlambda_drz)                                * vx + (-slat * slon * dlambda_drz)                                 * vy + clat * dlambda_drz * vz;
					doubleType dve_drx = -clon * dlongitude_drx * vx - slon * dlongitude_drx * vy;
					doubleType dve_dry = -clon * dlongitude_dry * vx - slon * dlongitude_dry * vy;
					doubleType dve_drz = 0.0;
					doubleType dvn_drx = (-clat * clon * dlambda_drx + slat * slon * dlongitude_drx) * vx + (-clat * slon * dlambda_drx - slat * clon * dlongitude_drx) * vy - slat * dlambda_drx * vz;
					doubleType dvn_dry = (-clat * clon * dlambda_dry + slat * slon * dlongitude_dry) * vx + (-clat * slon * dlambda_dry - slat * clon * dlongitude_dry) * vy - slat * dlambda_dry * vz;
					doubleType dvn_drz = (-clat * clon * dlambda_drz)                                * vx + (-clat * slon * dlambda_drz)                                * vy - slat * dlambda_drz * vz;

					doubleType dvu_dvx = clat * clon;
					doubleType dvu_dvy = clat * slon;
					doubleType dvu_dvz = slat;
					doubleType dve_dvx = -slon;
					doubleType dve_dvy = clon;
					doubleType dve_dvz = 0.0;
					doubleType dvn_dvx = -slat * clon;
					doubleType dvn_dvy = -slat * slon;
					doubleType dvn_dvz = clat;

					doubleType vne = sqrt(vn * vn + ve * ve);
					doubleType vne2 = vn / vne / vne;

					double dAz_drx = (vne2 * (dve_drx - dvn_drx * ve / vn)) _GETVALUE;
					double dAz_dry = (vne2 * (dve_dry - dvn_dry * ve / vn)) _GETVALUE;
					double dAz_drz = (vne2 * (dve_drz - dvn_drz * ve / vn)) _GETVALUE;
					double dAz_dvx = (vne2 * (dve_dvx - dvn_dvx * ve / vn)) _GETVALUE;
					double dAz_dvy = (vne2 * (dve_dvy - dvn_dvy * ve / vn)) _GETVALUE;
					double dAz_dvz = (vne2 * (dve_dvz - dvn_dvz * ve / vn)) _GETVALUE;

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

							double Gentry = ( dAz_drx * drxConstraintFrame_dstateICRF
											+ dAz_dry * dryConstraintFrame_dstateICRF
											+ dAz_drz * drzConstraintFrame_dstateICRF)
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

                            double Gentry = (  dAz_dvx * dvxConstraintFrame_dstateICRF
                                             + dAz_dvy * dvyConstraintFrame_dstateICRF
                                             + dAz_dvz * dvzConstraintFrame_dstateICRF)
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

                            double Gentry = (+ dAz_drx * drxConstraintFrame_dstateICRF
											 + dAz_dry * dryConstraintFrame_dstateICRF
											 + dAz_drz * drzConstraintFrame_dstateICRF
											 + dAz_dvx * dvxConstraintFrame_dstateICRF
                                             + dAz_dvy * dvyConstraintFrame_dstateICRF
                                             + dAz_dvz * dvzConstraintFrame_dstateICRF)
                                * std::get<2>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }


                    //Step 4.4: derivatives of longitude with respect to frame rotation due to change in epoch
                    math::Matrix<doubleType> dR_from_ICRF_to_ConstraintFrame_dt = this->myUniverse->LocalFrame.get_dRdt(ReferenceFrame::ICRF, this->myReferenceFrame);
					this->dPositionBeforeEventConstraintFrame_dt = dR_from_ICRF_to_ConstraintFrame_dt * StateBeforeEventICRF.getSubMatrix1D(0, 2);
					this->dVelocityBeforeEventConstraintFrame_dt = dR_from_ICRF_to_ConstraintFrame_dt * StateBeforeEventICRF.getSubMatrix1D(3, 5);

                    double timeDerivative = ( dAz_drx * this->dPositionBeforeEventConstraintFrame_dt(0)
											+ dAz_dry * this->dPositionBeforeEventConstraintFrame_dt(1)
											+ dAz_drz * this->dPositionBeforeEventConstraintFrame_dt(2)
											+ dAz_dvx * this->dVelocityBeforeEventConstraintFrame_dt(0)
											+ dAz_dvy * this->dVelocityBeforeEventConstraintFrame_dt(1)
											+ dAz_dvz * this->dVelocityBeforeEventConstraintFrame_dt(2)) _GETVALUE;

                    for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_time_variables.size(); ++entryIndex)
                    {
                        size_t Gindex = this->Gindex_constraint_wrt_time_variables[entryIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * timeDerivative;
                    }
                }//end derivatives
            }//end process_constraint()

            void BoundaryVelocitySphericalAzimuthConstraint::output(std::ofstream& outputfile)
            {
                std::string framestring = this->myUniverse->central_body_name + " " + this->ReferenceFrameStrings[this->myReferenceFrame];

                outputfile << this->myBoundaryEvent->getName() << " velocity spherical Azimuth (degrees, " << framestring << "): " << this->VelocitySphericalAzimuth _GETVALUE / math::deg2rad << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG