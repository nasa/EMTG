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

#include "BoundaryOrbitalEnergyConstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            BoundaryOrbitalEnergyConstraint::BoundaryOrbitalEnergyConstraint(const std::string& name,
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
                this->rvec_cb2sc = math::Matrix<doubleType>(3, 1, 0.0);
                this->vvec_cb2sc = math::Matrix<doubleType>(3, 1, 0.0);
                this->d_rvec_cb2sc_dt = math::Matrix<double>(3, 1, 0.0);
                this->d_vvec_cb2sc_dt = math::Matrix<double>(3, 1, 0.0);

                this->myReferenceFrame = ReferenceFrame::ICRF;
            }

            void BoundaryOrbitalEnergyConstraint::calcbounds()
            {
                // parse the constraint definition
				// constraintDefinition is like j1p0_departure_orbitalenergy_18.1km2s2_22.5km2s2
                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);


                //create the constraint
                //figure out the lower and upper bounds
                if (ConstraintDefinitionCell[3].find("km2s2") < 1024)
                {
                    this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[3], "km2s2")) * ((this->myUniverse->TU * this->myUniverse->TU) / (this->myUniverse->LU * this->myUniverse->LU)));
                }
                else
				{
					return throw std::logic_error("Please specify km2s2 for orbital energy constraint bounds: " + this->Fdescriptions->back() + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				}
                

                if (ConstraintDefinitionCell[4].find("km2s2") < 1024)
                {
                    this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[4], "km2s2")) * ((this->myUniverse->TU * this->myUniverse->TU) / (this->myUniverse->LU * this->myUniverse->LU)));
                }
                else
				{
					return throw std::logic_error("Please specify km2s2 for orbital energy constraint bounds: " + this->Fdescriptions->back() + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				}
                
                this->Fdescriptions->push_back(prefix + "orbital energy constraint w.r.t. " + this->myUniverse->central_body_name);

                // sparsity pattern
				// derivatives with respect to anything influencing the boundary event's right-hand six-state vector
				std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent();
				std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();

				for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
				{
					// POSITION
					// non-time variables
					std::vector<size_t> state_dIndex_energy_position_wrt_StateBeforeEvent;
					std::vector<size_t> state_Gindex_energy_position_wrt_StateBeforeEvent_variables;
					for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforeEvent.size(); ++dIndex)
					{
						if (std::get<1>(Derivatives_of_StateBeforeEvent[dIndex]) == stateIndex)
						{
							state_dIndex_energy_position_wrt_StateBeforeEvent.push_back(dIndex);

							this->create_sparsity_entry(this->Fdescriptions->size() - 1,
								std::get<0>(Derivatives_of_StateBeforeEvent[dIndex]),
								state_Gindex_energy_position_wrt_StateBeforeEvent_variables);
						}
					}
					this->dIndex_energy_position_wrt_StateBeforeEvent.push_back(state_dIndex_energy_position_wrt_StateBeforeEvent);
					this->Gindex_energy_position_wrt_StateBeforeEvent_variables.push_back(state_Gindex_energy_position_wrt_StateBeforeEvent_variables);

					// time variables
					std::vector<size_t> state_dIndex_energy_position_wrt_StateBeforeEvent_wrt_Time;
					std::vector<size_t> state_Gindex_energy_position_wrt_StateBeforeEvent_time_variables;
					for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
					{
						if (std::get<1>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]) == stateIndex)
						{
							state_dIndex_energy_position_wrt_StateBeforeEvent_wrt_Time.push_back(dIndex);

							this->create_sparsity_entry(this->Fdescriptions->size() - 1,
								std::get<0>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]),
								state_Gindex_energy_position_wrt_StateBeforeEvent_time_variables);
						}
					}
					this->dIndex_energy_position_wrt_StateBeforeEvent_wrt_Time.push_back(state_dIndex_energy_position_wrt_StateBeforeEvent_wrt_Time);
					this->Gindex_energy_position_wrt_StateBeforeEvent_time_variables.push_back(state_Gindex_energy_position_wrt_StateBeforeEvent_time_variables);

					// VELOCITY
					// non-time variables
					std::vector<size_t> state_dIndex_energy_velocity_with_respect_to_StateBeforeEvent;
					std::vector<size_t> state_Gindex_energy_velocity_wrt_StateBeforeEvent_variables;
					for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforeEvent.size(); ++dIndex)
					{
						if (std::get<1>(Derivatives_of_StateBeforeEvent[dIndex]) == stateIndex + 3)
						{
							state_dIndex_energy_velocity_with_respect_to_StateBeforeEvent.push_back(dIndex);

							this->create_sparsity_entry(this->Fdescriptions->size() - 1,
								std::get<0>(Derivatives_of_StateBeforeEvent[dIndex]),
								state_Gindex_energy_velocity_wrt_StateBeforeEvent_variables);
						}
					}
					this->dIndex_energy_velocity_wrt_StateBeforeEvent.push_back(state_dIndex_energy_velocity_with_respect_to_StateBeforeEvent);
					this->Gindex_energy_velocity_wrt_StateBeforeEvent_variables.push_back(state_Gindex_energy_velocity_wrt_StateBeforeEvent_variables);

					// time variables
					std::vector<size_t> state_dIndex_energy_velocity_wrt_StateBeforeEvent_wrt_Time;
					std::vector<size_t> state_Gindex_energy_velocity_wrt_StateBeforeEvent_time_variables;
					for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
					{
						if (std::get<1>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]) == stateIndex + 3)
						{
							state_dIndex_energy_velocity_wrt_StateBeforeEvent_wrt_Time.push_back(dIndex);

							this->create_sparsity_entry(this->Fdescriptions->size() - 1,
								std::get<0>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]),
								state_Gindex_energy_velocity_wrt_StateBeforeEvent_time_variables);
						}
					}
					this->dIndex_energy_velocity_wrt_StateBeforeEvent_wrt_Time.push_back(state_dIndex_energy_velocity_wrt_StateBeforeEvent_wrt_Time);
					this->Gindex_energy_velocity_wrt_StateBeforeEvent_time_variables.push_back(state_Gindex_energy_velocity_wrt_StateBeforeEvent_time_variables);
				}
            }//end calcbounds()

            void BoundaryOrbitalEnergyConstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {
                // get the state relative to the central body
                // grab the state before event to get the current central body relative state
                // state after event will grab the current central body's central body (body that the central body is orbiting) relative state
                math::Matrix<doubleType>& StateRelativeToCentralBody = this->myBoundaryEvent->get_state_before_event();
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->rvec_cb2sc(stateIndex) = StateRelativeToCentralBody(stateIndex);
                    this->vvec_cb2sc(stateIndex) = StateRelativeToCentralBody(stateIndex + 3);
                }

                doubleType r = this->rvec_cb2sc.norm();
                doubleType v = this->vvec_cb2sc.norm();

                // compute the specific orbital energy
                this->orbital_energy = 0.5 * v * v - this->myUniverse->mu / r;
                F[Findex++] = this->orbital_energy * ((this->myUniverse->TU * this->myUniverse->TU) / (this->myUniverse->LU * this->myUniverse->LU));

                // derivatives
                if (needG)
                {
                    // zero-out all relevant Jacobian entries
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
					{
						// POSITION
						// non-time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_energy_position_wrt_StateBeforeEvent_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_energy_position_wrt_StateBeforeEvent_variables[stateIndex][entryIndex]] = 0.0;
						}

						// time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_energy_position_wrt_StateBeforeEvent_time_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_energy_position_wrt_StateBeforeEvent_time_variables[stateIndex][entryIndex]] = 0.0;
						}

						// VELOCITY
						// non-time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_energy_velocity_wrt_StateBeforeEvent_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_energy_velocity_wrt_StateBeforeEvent_variables[stateIndex][entryIndex]] = 0.0;
						}

						// time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_energy_velocity_wrt_StateBeforeEvent_time_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_energy_velocity_wrt_StateBeforeEvent_time_variables[stateIndex][entryIndex]] = 0.0;
						}
					}
                    
                    double denergy_dr = (this->myUniverse->mu / (r * r)) _GETVALUE;
                    double denergy_dv = v _GETVALUE;
                    double dr_drvec[3] = {(this->rvec_cb2sc(0) / r) _GETVALUE, (this->rvec_cb2sc(1) / r) _GETVALUE, (this->rvec_cb2sc(2) / r) _GETVALUE};
                    double dv_dvvec[3] = {(this->vvec_cb2sc(0) / v) _GETVALUE, (this->vvec_cb2sc(1) / v) _GETVALUE, (this->vvec_cb2sc(2) / v) _GETVALUE};
                    
                    // first for non-time decision variables
					std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent();

					for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
					{
						// POSITION
						for (size_t entryIndex = 0; entryIndex < this->dIndex_energy_position_wrt_StateBeforeEvent[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_energy_position_wrt_StateBeforeEvent[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_energy_position_wrt_StateBeforeEvent_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);
							
							double Gentry = (denergy_dr * dr_drvec[stateIndex]) * std::get<2>(Derivatives_of_StateBeforeEvent[dIndex]);

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry * ((this->myUniverse->TU * this->myUniverse->TU) / (this->myUniverse->LU * this->myUniverse->LU));
						}

						// VELOCITY
						for (size_t entryIndex = 0; entryIndex < this->dIndex_energy_velocity_wrt_StateBeforeEvent[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_energy_velocity_wrt_StateBeforeEvent[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_energy_velocity_wrt_StateBeforeEvent_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);

							double Gentry = (denergy_dv * dv_dvvec[stateIndex]) * std::get<2>(Derivatives_of_StateBeforeEvent[dIndex]);

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry * ((this->myUniverse->TU * this->myUniverse->TU) / (this->myUniverse->LU * this->myUniverse->LU));
						}
					}

					// now for the time decision variables
					std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforeEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();

					for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
					{
						// POSITION
						for (size_t entryIndex = 0; entryIndex < this->dIndex_energy_position_wrt_StateBeforeEvent_wrt_Time[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_energy_position_wrt_StateBeforeEvent_wrt_Time[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_energy_position_wrt_StateBeforeEvent_time_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);

							double Gentry = (denergy_dr * dr_drvec[stateIndex]) * std::get<2>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry / this->myUniverse->LU;
						}

						// VELOCITY
						for (size_t entryIndex = 0; entryIndex < this->dIndex_energy_velocity_wrt_StateBeforeEvent_wrt_Time[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_energy_velocity_wrt_StateBeforeEvent_wrt_Time[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_energy_velocity_wrt_StateBeforeEvent_time_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);

							double Gentry = (denergy_dv * dv_dvvec[stateIndex]) * std::get<2>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry / this->myUniverse->LU;
						}
					}
                    
                    
                }//end derivatives
            }//end process_constraint()

            void BoundaryOrbitalEnergyConstraint::output(std::ofstream& outputfile)
            {
                std::string bodystring = this->myUniverse->central_body_name;

                outputfile << this->myBoundaryEvent->getName() << " orbital energy w.r.t. " << bodystring << " [km^2/s^2]: " << (this->orbital_energy) _GETVALUE << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG