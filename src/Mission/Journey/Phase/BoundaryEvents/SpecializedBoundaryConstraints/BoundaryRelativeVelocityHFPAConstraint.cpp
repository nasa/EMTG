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

#include "BoundaryRelativeVelocityHFPAConstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            BoundaryRelativeVelocityHFPAConstraint::BoundaryRelativeVelocityHFPAConstraint(const std::string& name,
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

            void BoundaryRelativeVelocityHFPAConstraint::calcbounds()
            {
                // Step 1: parse the constraint definition
				// constraintDefinition is like p0_arrival_RelativeVHFPA_-10.9_-11.0_trueofdatebcf

                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                //Step 2: lower bound
                this->Flowerbounds->push_back(std::stod(ConstraintDefinitionCell[3]) * math::deg2rad);


                //Step 3: upper bound
                this->Fupperbounds->push_back(std::stod(ConstraintDefinitionCell[4]) * math::deg2rad);
                this->Fdescriptions->push_back(prefix + "relative velocity HFPA");

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

            void BoundaryRelativeVelocityHFPAConstraint::process_constraint(const std::vector<doubleType>& X,
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


                //Step 3: compute velocity HFPA
				//First, make sure that the probe is not at the pole because there is a singularity at the pole. 
				doubleType rx = this->PositionBeforeEventConstraintFrame(0);
				doubleType ry = this->PositionBeforeEventConstraintFrame(1);
				doubleType rz = this->PositionBeforeEventConstraintFrame(2);
				doubleType rxy = sqrt(rx * rx + ry * ry);
				if (rxy < math::SMALL)
				{
					throw std::runtime_error(this->prefix + "The probe is too close to the pole to make detic local coorinate for HFPA.");
				}
				doubleType R = this->myUniverse->central_body_radius;
				doubleType f = this->myUniverse->central_body.flattening_coefficient;
				doubleType f2 = (1.0 - f) * (1.0 - f);
				doubleType Wdot = this->myUniverse->LocalFrame.getWdot() / 86400.0;

				doubleType vx = this->VelocityBeforeEventConstraintFrame(0);
				doubleType vy = this->VelocityBeforeEventConstraintFrame(1);
				doubleType vz = this->VelocityBeforeEventConstraintFrame(2);
				doubleType vrx = vx + Wdot * ry;
				doubleType vry = vy - Wdot * rx;
				doubleType vrz = vz;

				doubleType r = sqrt(rx * rx + ry * ry + rz * rz);

				doubleType lambda = atan2(rz, rxy);
				doubleType longitude = atan2(ry, rx);
				doubleType tanlambda = tan(lambda);
				doubleType tanlambda2 = tanlambda * tanlambda;
				doubleType xa = (1.0 - f) * R / sqrt(tanlambda2 + f2);
				doubleType mua = atan2(tanlambda, f2);
				doubleType ra = xa / cos(lambda);
				doubleType l = r - ra;
				doubleType dellambda = mua - lambda;
				doubleType h = l * cos(dellambda);
				doubleType rhoa = R * f2 / pow(1.0 - (2.0 * f - f * f) * pow(sin(mua), 2.0), 1.5);
				doubleType delmu1 = l * sin(dellambda) / (rhoa + h);
				doubleType delmu = atan(delmu1);
				doubleType mu = mua - delmu;

				doubleType clat = cos(mu);
				doubleType slat = sin(mu);
				doubleType clon = cos(longitude);
				doubleType slon = sin(longitude);

				doubleType vu = clat * clon * vrx + clat * slon * vry + slat * vrz;
				doubleType ve = -slon * vrx + clon * vry;
				doubleType vn = -slat * clon * vrx - slat * slon * vry + clat * vrz;

				doubleType vne = sqrt(vn * vn + ve * ve);

				this->RelativeVelocityHFPA = atan2(vu, vne);

                F[Findex++] = this->RelativeVelocityHFPA;
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
					double dr_drx = (rx / r) _GETVALUE;
					double dr_dry = (ry / r) _GETVALUE;
					double dr_drz = (rz / r) _GETVALUE;

					double dlambda_drx = (-(rx * pow(rxy, -3.0) * rz) / (rz * rz / rxy / rxy + 1.0)) _GETVALUE;
					double dlambda_dry = (-(ry * pow(rxy, -3.0) * rz) / (rz * rz / rxy / rxy + 1.0)) _GETVALUE;
					double dlambda_drz = (1.0 / (rxy * (rz * rz / rxy / rxy + 1.0))) _GETVALUE;

					double dlongitude_drx = (-ry / rxy / rxy) _GETVALUE;
					double dlongitude_dry = (rx / rxy / rxy) _GETVALUE;

					double dxa_dlambda = (-(R * (1.0 - f) * pow(1.0 / cos(lambda), 2.0) * tanlambda) / pow(tanlambda2 + f2, 1.5)) _GETVALUE;
					double dxa_drx = dxa_dlambda * dlambda_drx;
					double dxa_dry = dxa_dlambda * dlambda_dry;
					double dxa_drz = dxa_dlambda * dlambda_drz;

					double dmua_dlambda = (pow(1.0 / cos(lambda), 2.0) / (f2 * (tanlambda2 / f2 / f2 + 1.0))) _GETVALUE;
					double dmua_drx = dmua_dlambda * dlambda_drx;
					double dmua_dry = dmua_dlambda * dlambda_dry;
					double dmua_drz = dmua_dlambda * dlambda_drz;

					double dra_drx = (dxa_drx / cos(lambda) + xa * tanlambda / cos(lambda) * dlambda_drx) _GETVALUE;
					double dra_dry = (dxa_dry / cos(lambda) + xa * tanlambda / cos(lambda) * dlambda_dry) _GETVALUE;
					double dra_drz = (dxa_drz / cos(lambda) + xa * tanlambda / cos(lambda) * dlambda_drz) _GETVALUE;

					double dl_drx = dr_drx - dra_drx;
					double dl_dry = dr_dry - dra_dry;
					double dl_drz = dr_drz - dra_drz;

					double ddellambda_drx = dmua_drx - dlambda_drx;
					double ddellambda_dry = dmua_dry - dlambda_dry;
					double ddellambda_drz = dmua_drz - dlambda_drz;

					double dh_drx = (dl_drx * cos(dellambda) - l * sin(dellambda) * ddellambda_drx) _GETVALUE;
					double dh_dry = (dl_dry * cos(dellambda) - l * sin(dellambda) * ddellambda_dry) _GETVALUE;
					double dh_drz = (dl_drz * cos(dellambda) - l * sin(dellambda) * ddellambda_drz) _GETVALUE;

					double drhoa_dmua = ((3.0 * R * f2 * (2.0 * f - f * f) * cos(mua) * sin(mua)) / pow(1.0 - (2.0 * f - f * f) * pow(sin(mua), 2.0), 2.5)) _GETVALUE;
					double drhoa_drx = drhoa_dmua * dmua_drx;
					double drhoa_dry = drhoa_dmua * dmua_dry;
					double drhoa_drz = drhoa_dmua * dmua_drz;

					double ddelmu1_drx = (l * cos(dellambda) / (rhoa + h) * ddellambda_drx + delmu1 * (dl_drx / l - (drhoa_drx + dh_drx) / (rhoa + h))) _GETVALUE;
					double ddelmu1_dry = (l * cos(dellambda) / (rhoa + h) * ddellambda_dry + delmu1 * (dl_dry / l - (drhoa_dry + dh_dry) / (rhoa + h))) _GETVALUE;
					double ddelmu1_drz = (l * cos(dellambda) / (rhoa + h) * ddellambda_drz + delmu1 * (dl_drz / l - (drhoa_drz + dh_drz) / (rhoa + h))) _GETVALUE;
					double ddelmu_drx = (1.0 / (delmu1 * delmu1 + 1) * ddelmu1_drx) _GETVALUE;
					double ddelmu_dry = (1.0 / (delmu1 * delmu1 + 1) * ddelmu1_dry) _GETVALUE;
					double ddelmu_drz = (1.0 / (delmu1 * delmu1 + 1) * ddelmu1_drz) _GETVALUE;

					double dmu_drx = dmua_drx - ddelmu_drx;
					double dmu_dry = dmua_dry - ddelmu_dry;
					double dmu_drz = dmua_drz - ddelmu_drz;

					double dvu_drx = ((-slat * clon * dmu_drx - clat * slon * dlongitude_drx) * vrx + (-slat * slon * dmu_drx + clat * clon * dlongitude_drx) * vry + clat * dmu_drx * vrz - clat * slon * Wdot) _GETVALUE;
					double dvu_dry = ((-slat * clon * dmu_dry - clat * slon * dlongitude_dry) * vrx + (-slat * slon * dmu_dry + clat * clon * dlongitude_dry) * vry + clat * dmu_dry * vrz + clat * clon * Wdot) _GETVALUE;
					double dvu_drz = ((-slat * clon * dmu_drz)                               * vrx + (-slat * slon * dmu_drz)                               * vry + clat * dmu_drz * vrz) _GETVALUE;
					double dve_drx = (-clon * dlongitude_drx * vrx - slon * dlongitude_drx * vry - clon * Wdot) _GETVALUE;
					double dve_dry = (-clon * dlongitude_dry * vrx - slon * dlongitude_dry * vry - slon * Wdot) _GETVALUE;
					double dve_drz = 0.0;
					double dvn_drx = ((-clat * clon * dmu_drx + slat * slon * dlongitude_drx) * vrx + (-clat * slon * dmu_drx - slat * clon * dlongitude_drx) * vry - slat * dmu_drx * vrz + slat * slon * Wdot) _GETVALUE;
					double dvn_dry = ((-clat * clon * dmu_dry + slat * slon * dlongitude_dry) * vrx + (-clat * slon * dmu_dry - slat * clon * dlongitude_dry) * vry - slat * dmu_dry * vrz - slat * clon * Wdot) _GETVALUE;
					double dvn_drz = ((-clat * clon * dmu_drz)                               * vrx + (-clat * slon * dmu_drz)                               * vry - slat * dmu_drz * vrz) _GETVALUE;

					double dvu_dvx = (clat * clon) _GETVALUE;
					double dvu_dvy = (clat * slon) _GETVALUE;
					double dvu_dvz = (slat)_GETVALUE;
					double dve_dvx = (-slon) _GETVALUE;
					double dve_dvy = (clon)_GETVALUE;
					double dve_dvz = 0.0;
					double dvn_dvx = (-slat * clon) _GETVALUE;
					double dvn_dvy = (-slat * slon) _GETVALUE;
					double dvn_dvz = (clat)_GETVALUE;

					double vne2 = (vu / vne / vne) _GETVALUE;
					double v2 = (vu * vu + vne * vne) _GETVALUE;

					double dHFPA_drx = (vne / v2 * (dvu_drx - vne2 * (vn * dvn_drx + ve * dve_drx))) _GETVALUE;
					double dHFPA_dry = (vne / v2 * (dvu_dry - vne2 * (vn * dvn_dry + ve * dve_dry))) _GETVALUE;
					double dHFPA_drz = (vne / v2 * (dvu_drz - vne2 * (vn * dvn_drz + ve * dve_drz))) _GETVALUE;
					double dHFPA_dvx = (vne / v2 * (dvu_dvx - vne2 * (vn * dvn_dvx + ve * dve_dvx))) _GETVALUE;
					double dHFPA_dvy = (vne / v2 * (dvu_dvy - vne2 * (vn * dvn_dvy + ve * dve_dvy))) _GETVALUE;
					double dHFPA_dvz = (vne / v2 * (dvu_dvz - vne2 * (vn * dvn_dvz + ve * dve_dvz))) _GETVALUE;

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

							double Gentry = ( dHFPA_drx * drxConstraintFrame_dstateICRF
											+ dHFPA_dry * dryConstraintFrame_dstateICRF
											+ dHFPA_drz * drzConstraintFrame_dstateICRF)
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

                            double Gentry = (  dHFPA_dvx * dvxConstraintFrame_dstateICRF
                                             + dHFPA_dvy * dvyConstraintFrame_dstateICRF
                                             + dHFPA_dvz * dvzConstraintFrame_dstateICRF)
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

                            double Gentry = (+ dHFPA_drx * drxConstraintFrame_dstateICRF
											 + dHFPA_dry * dryConstraintFrame_dstateICRF
											 + dHFPA_drz * drzConstraintFrame_dstateICRF
											 + dHFPA_dvx * dvxConstraintFrame_dstateICRF
                                             + dHFPA_dvy * dvyConstraintFrame_dstateICRF
                                             + dHFPA_dvz * dvzConstraintFrame_dstateICRF)
                                * std::get<2>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }


                    //Step 4.4: derivatives of longitude with respect to frame rotation due to change in epoch
                    math::Matrix<doubleType> dR_from_ICRF_to_ConstraintFrame_dt = this->myUniverse->LocalFrame.get_dRdt(ReferenceFrame::ICRF, this->myReferenceFrame);
					this->dPositionBeforeEventConstraintFrame_dt = dR_from_ICRF_to_ConstraintFrame_dt * StateBeforeEventICRF.getSubMatrix1D(0, 2);
					this->dVelocityBeforeEventConstraintFrame_dt = dR_from_ICRF_to_ConstraintFrame_dt * StateBeforeEventICRF.getSubMatrix1D(3, 5);

                    double timeDerivative = ( dHFPA_drx * this->dPositionBeforeEventConstraintFrame_dt(0)
											+ dHFPA_dry * this->dPositionBeforeEventConstraintFrame_dt(1)
											+ dHFPA_drz * this->dPositionBeforeEventConstraintFrame_dt(2)
											+ dHFPA_dvx * this->dVelocityBeforeEventConstraintFrame_dt(0)
											+ dHFPA_dvy * this->dVelocityBeforeEventConstraintFrame_dt(1)
											+ dHFPA_dvz * this->dVelocityBeforeEventConstraintFrame_dt(2)) _GETVALUE;

                    for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_time_variables.size(); ++entryIndex)
                    {
                        size_t Gindex = this->Gindex_constraint_wrt_time_variables[entryIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * timeDerivative;
                    }
                }//end derivatives
            }//end process_constraint()

            void BoundaryRelativeVelocityHFPAConstraint::output(std::ofstream& outputfile)
            {
                std::string framestring = this->myUniverse->central_body_name + " " + this->ReferenceFrameStrings[this->myReferenceFrame];

                outputfile << this->myBoundaryEvent->getName() << " relative velocity HFPA (degrees, " << framestring << "): " << this->RelativeVelocityHFPA _GETVALUE / math::deg2rad << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG