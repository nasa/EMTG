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

#include "BoundaryDeticAltitudeConstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            BoundaryDeticAltitudeConstraint::BoundaryDeticAltitudeConstraint(const std::string& name,
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

            void BoundaryDeticAltitudeConstraint::calcbounds()
            {
                // Step 1: parse the constraint definition
				// constraintDefinition is like p0_arrival_DeticAltitude_200.0_200.1_J2000BCF

                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                //Step 2: lower bound
                this->Flowerbounds->push_back(std::stod(ConstraintDefinitionCell[3]) );


                //Step 3: upper bound
                this->Fupperbounds->push_back(std::stod(ConstraintDefinitionCell[4]) );
                this->Fdescriptions->push_back(prefix + "detic Altitude");

                //Step 5: sparsity pattern
                //derivatives with respect to anything influencing the boundary event's right-hand velocity vector
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

                //derivatives with respect to time variables that affect current epoch
                std::vector<size_t> timeVariables = this->myBoundaryEvent->get_Xindices_EventRightEpoch();
                for (size_t Xindex : timeVariables)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        Xindex,
                        this->Gindex_constraint_wrt_time_variables);
                }
            }//end calcbounds()

            void BoundaryDeticAltitudeConstraint::process_constraint(const std::vector<doubleType>& X,
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

                //Step 3: compute deticAltitude magnitude
				//First, make sure that the probe is not at the pole because there is a singularity at the pole. 
				doubleType rx = this->PositionBeforeEventConstraintFrame(0);
				doubleType ry = this->PositionBeforeEventConstraintFrame(1);
				doubleType rz = this->PositionBeforeEventConstraintFrame(2);
				doubleType rxy = sqrt(rx * rx + ry * ry);
				if (rxy < math::SMALL)
				{
					throw std::runtime_error(this->prefix + "The probe is too close to the pole to calculate detic altitude and its derivatives.");
				}

				doubleType R = this->myUniverse->central_body_radius;
				doubleType f = this->myUniverse->central_body.flattening_coefficient;
				doubleType f2 = (1.0 - f) * (1.0 - f) ;
				doubleType e2 = 1 - f2;

				doubleType r = this->PositionBeforeEventConstraintFrame.norm();

				doubleType lambda = atan2(rz, rxy);
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

				doubleType N = R / sqrt(1.0 - e2 * slat * slat);

				this->DeticAltitude = rxy * clat + (rz + e2 * N * slat) * slat - N;

                F[Findex++] = this->DeticAltitude;
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
					doubleType dr_drx = rx / r;
					doubleType dr_dry = ry / r;
					doubleType dr_drz = rz / r;

					doubleType dlambda_drx = -(rx * pow(rxy, -3.0) * rz) / (rz * rz / rxy / rxy + 1.0);
					doubleType dlambda_dry = -(ry * pow(rxy, -3.0) * rz) / (rz * rz / rxy / rxy + 1.0);
					doubleType dlambda_drz = 1.0 / (rxy * (rz * rz / rxy / rxy + 1.0));

					doubleType dxa_dlambda = -(R * (1.0 - f) * pow(1.0 / cos(lambda), 2.0) * tanlambda) / pow(tanlambda2 + f2, 1.5);
					doubleType dxa_drx = dxa_dlambda * dlambda_drx;
					doubleType dxa_dry = dxa_dlambda * dlambda_dry;
					doubleType dxa_drz = dxa_dlambda * dlambda_drz;

					doubleType dmua_dlambda = pow(1.0 / cos(lambda), 2.0) / (f2 * (tanlambda2 / f2 / f2 + 1.0));
					doubleType dmua_drx = dmua_dlambda * dlambda_drx;
					doubleType dmua_dry = dmua_dlambda * dlambda_dry;
					doubleType dmua_drz = dmua_dlambda * dlambda_drz;

					doubleType dra_drx = dxa_drx / cos(lambda) + xa * tanlambda / cos(lambda) * dlambda_drx;
					doubleType dra_dry = dxa_dry / cos(lambda) + xa * tanlambda / cos(lambda) * dlambda_dry;
					doubleType dra_drz = dxa_drz / cos(lambda) + xa * tanlambda / cos(lambda) * dlambda_drz;

					doubleType dl_drx = dr_drx - dra_drx;
					doubleType dl_dry = dr_dry - dra_dry;
					doubleType dl_drz = dr_drz - dra_drz;

					doubleType ddellambda_drx = dmua_drx - dlambda_drx;
					doubleType ddellambda_dry = dmua_dry - dlambda_dry;
					doubleType ddellambda_drz = dmua_drz - dlambda_drz;

					doubleType dh_drx = dl_drx * cos(dellambda) - l * sin(dellambda) * ddellambda_drx;
					doubleType dh_dry = dl_dry * cos(dellambda) - l * sin(dellambda) * ddellambda_dry;
					doubleType dh_drz = dl_drz * cos(dellambda) - l * sin(dellambda) * ddellambda_drz;

					doubleType drhoa_dmua = (3.0 * R * f2 * (2.0 * f - f * f) * cos(mua) * sin(mua)) / pow(1.0 - (2.0 * f - f * f) * pow(sin(mua), 2.0), 2.5);
					doubleType drhoa_drx = drhoa_dmua * dmua_drx;
					doubleType drhoa_dry = drhoa_dmua * dmua_dry;
					doubleType drhoa_drz = drhoa_dmua * dmua_drz;

					doubleType ddelmu1_drx = l * cos(dellambda) / (rhoa + h) * ddellambda_drx + delmu1 * (dl_drx / l - (drhoa_drx + dh_drx) / (rhoa + h));
					doubleType ddelmu1_dry = l * cos(dellambda) / (rhoa + h) * ddellambda_dry + delmu1 * (dl_dry / l - (drhoa_dry + dh_dry) / (rhoa + h));
					doubleType ddelmu1_drz = l * cos(dellambda) / (rhoa + h) * ddellambda_drz + delmu1 * (dl_drz / l - (drhoa_drz + dh_drz) / (rhoa + h));
					doubleType ddelmu_drx = 1.0 / (delmu1 * delmu1 + 1) * ddelmu1_drx;
					doubleType ddelmu_dry = 1.0 / (delmu1 * delmu1 + 1) * ddelmu1_dry;
					doubleType ddelmu_drz = 1.0 / (delmu1 * delmu1 + 1) * ddelmu1_drz;

					doubleType dmu_drx = (dmua_drx - ddelmu_drx) ;
					doubleType dmu_dry = (dmua_dry - ddelmu_dry) ;
					doubleType dmu_drz = (dmua_drz - ddelmu_drz) ;

					doubleType dN_drx = R * e2 * slat / pow(1.0 - e2 * slat * slat, 1.5) * clat * dmu_drx;
					doubleType dN_dry = R * e2 * slat / pow(1.0 - e2 * slat * slat, 1.5) * clat * dmu_dry;
					doubleType dN_drz = R * e2 * slat / pow(1.0 - e2 * slat * slat, 1.5) * clat * dmu_drz;

					double ddetalt_drx = (rx / rxy * clat - rxy * slat * dmu_drx + (e2 * dN_drx * slat + e2 * N * clat * dmu_drx) * slat + (rz + e2 * N * slat) * clat * dmu_drx - dN_drx) _GETVALUE;
					double ddetalt_dry = (ry / rxy * clat - rxy * slat * dmu_dry + (e2 * dN_dry * slat + e2 * N * clat * dmu_dry) * slat + (rz + e2 * N * slat) * clat * dmu_dry - dN_dry) _GETVALUE;
					double ddetalt_drz = (-rxy * slat * dmu_drz + (1.0 + e2 * dN_drz * slat + e2 * N * clat * dmu_drz) * slat + (rz + e2 * N * slat) * clat * dmu_drz - dN_drz) _GETVALUE;

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

							double Gentry = ( ddetalt_drx * drxConstraintFrame_dstateICRF
											+ ddetalt_dry * dryConstraintFrame_dstateICRF
											+ ddetalt_drz * drzConstraintFrame_dstateICRF)
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

                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_StateBeforeEvent_wrt_Time[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_StateBeforeEvent_wrt_Time[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_StateBeforeEvent_time_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (+ ddetalt_drx * drxConstraintFrame_dstateICRF
											 + ddetalt_dry * dryConstraintFrame_dstateICRF
											 + ddetalt_drz * drzConstraintFrame_dstateICRF)
                                * std::get<2>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }


                    //Step 4.4: derivatives of altitude with respect to frame rotation due to change in epoch
                    math::Matrix<doubleType> dR_from_ICRF_to_ConstraintFrame_dt = this->myUniverse->LocalFrame.get_dRdt(ReferenceFrame::ICRF, this->myReferenceFrame);
					this->dPositionBeforeEventConstraintFrame_dt = dR_from_ICRF_to_ConstraintFrame_dt * StateBeforeEventICRF.getSubMatrix1D(0, 2);

                    double timeDerivative = ( ddetalt_drx * this->dPositionBeforeEventConstraintFrame_dt(0)
											+ ddetalt_dry * this->dPositionBeforeEventConstraintFrame_dt(1)
											+ ddetalt_drz * this->dPositionBeforeEventConstraintFrame_dt(2)) _GETVALUE;

                    for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_time_variables.size(); ++entryIndex)
                    {
                        size_t Gindex = this->Gindex_constraint_wrt_time_variables[entryIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * timeDerivative;
                    }
                }//end derivatives
            }//end process_constraint()

            void BoundaryDeticAltitudeConstraint::output(std::ofstream& outputfile)
            {
                std::string framestring = this->myUniverse->central_body_name + " " + this->ReferenceFrameStrings[this->myReferenceFrame];

                outputfile << this->myBoundaryEvent->getName() << " detic Altitude (km, " << framestring << "): " << this->DeticAltitude _GETVALUE << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG