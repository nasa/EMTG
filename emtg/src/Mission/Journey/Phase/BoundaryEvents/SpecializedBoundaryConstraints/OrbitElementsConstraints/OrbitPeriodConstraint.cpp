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

#include "OrbitPeriodConstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {
			OrbitPeriodConstraint::OrbitPeriodConstraint(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                BoundaryEventBase* myBoundaryEvent,
                const std::string& constraintDefinition) :
				OrbitElementConstraintBase::OrbitElementConstraintBase(name,
                    journeyIndex,
                    phaseIndex,
                    stageIndex,
                    Universe,
                    mySpacecraft,
                    myOptions,
                    myBoundaryEvent,
                    constraintDefinition)
            {
				// indicate to the boundary event that it should compute the orbital elements of the spacecraft
				this->myBoundaryEvent->setComputeOrbitElements(true);
                
                this->myReferenceFrame = ReferenceFrame::ICRF;

                this->myBoundaryEvent->add_orbit_element_reference_frame(this->myReferenceFrame);
            }

			void OrbitPeriodConstraint::calcbounds()
			{
				// Step 1: parse the constraint definition
				// constraintDefinition is like j0p3_arrival_orbit_period_lbkm_ubkm
				std::vector<std::string> ConstraintDefinitionCell;
				std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
				boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

				// Step 2: create the constraint
				// figure out the lower and upper bounds - depends on units
				this->Fdescriptions->push_back(prefix + " " + ConstraintDefinitionCell[1] + " " + ConstraintDefinitionCell[2]);

				std::string lower_bound = ConstraintDefinitionCell[3];
				if (ConstraintDefinitionCell[3].find("sec") < 1024)
				{
					this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(lower_bound, "sec")) / this->myUniverse->TU);
				}
				else if (lower_bound.find("min") < 1024)
				{
					this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(lower_bound, "min")) * 60.0 / this->myUniverse->TU);
				}
				else if (lower_bound.find("hr") < 1024)
				{
					this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(lower_bound, "hr")) * 3600.0 / this->myUniverse->TU);
				}
				else if (lower_bound.find("d") < 1024)
				{
					this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(lower_bound, "d")) * 86400.0 / this->myUniverse->TU);
				}
				else if (lower_bound.find("yr") < 1024)
				{
					this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(lower_bound, "yr")) * 365.25 * 86400.0 / this->myUniverse->TU);
				}
				else
				{
					return throw std::logic_error("Please specify either sec, hr, d or yr for orbit period constraint bounds: " + this->Fdescriptions->back() + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				}

				std::string upper_bound = ConstraintDefinitionCell[4];
				if (ConstraintDefinitionCell[4].find("sec") < 1024)
				{
					this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(upper_bound, "sec")) / this->myUniverse->TU);
				}
				else if (upper_bound.find("min") < 1024)
				{
					this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(upper_bound, "min")) * 60.0 / this->myUniverse->TU);
				}
				else if (upper_bound.find("hr") < 1024)
				{
					this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(upper_bound, "hr")) * 3600.0 / this->myUniverse->TU);
				}
				else if (upper_bound.find("d") < 1024)
				{
					this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(upper_bound, "d")) * 86400.0 / this->myUniverse->TU);
				}
				else if (upper_bound.find("yr") < 1024)
				{
					this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(upper_bound, "yr")) * 365.25 * 86400.0 / this->myUniverse->TU);
				}
				else
				{
					return throw std::logic_error("Please specify either sec, hr, d or yr for orbit period constraint bounds: " + this->Fdescriptions->back() + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				}


				// Step 3: sparsity pattern
				// derivatives with respect to anything influencing the boundary event's right-hand six-state vector
				std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();
				std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

				for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
				{
					// POSITION
					// non-time variables
					std::vector<size_t> state_dIndex_orbit_period_position_wrt_StateAfterEvent;
					std::vector<size_t> state_Gindex_orbit_period_position_wrt_StateAfterEvent_variables;
					for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent.size(); ++dIndex)
					{
						if (std::get<1>(Derivatives_of_StateAfterEvent[dIndex]) == stateIndex)
						{
							state_dIndex_orbit_period_position_wrt_StateAfterEvent.push_back(dIndex);

							this->create_sparsity_entry(this->Fdescriptions->size() - 1,
								std::get<0>(Derivatives_of_StateAfterEvent[dIndex]),
								state_Gindex_orbit_period_position_wrt_StateAfterEvent_variables);
						}
					}
					this->dIndex_orbit_period_position_wrt_StateAfterEvent.push_back(state_dIndex_orbit_period_position_wrt_StateAfterEvent);
					this->Gindex_orbit_period_position_wrt_StateAfterEvent_variables.push_back(state_Gindex_orbit_period_position_wrt_StateAfterEvent_variables);

					// time variables
					std::vector<size_t> state_dIndex_orbit_period_position_wrt_StateAfterEvent_wrt_Time;
					std::vector<size_t> state_Gindex_orbit_period_position_wrt_StateAfterEvent_time_variables;
					for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
					{
						if (std::get<1>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]) == stateIndex)
						{
							state_dIndex_orbit_period_position_wrt_StateAfterEvent_wrt_Time.push_back(dIndex);

							this->create_sparsity_entry(this->Fdescriptions->size() - 1,
								std::get<0>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]),
								state_Gindex_orbit_period_position_wrt_StateAfterEvent_time_variables);
						}
					}
					this->dIndex_orbit_period_position_wrt_StateAfterEvent_wrt_Time.push_back(state_dIndex_orbit_period_position_wrt_StateAfterEvent_wrt_Time);
					this->Gindex_orbit_period_position_wrt_StateAfterEvent_time_variables.push_back(state_Gindex_orbit_period_position_wrt_StateAfterEvent_time_variables);

					// VELOCITY
					// non-time variables
					std::vector<size_t> state_dIndex_orbit_period_velocity_with_respect_to_StateAfterEvent;
					std::vector<size_t> state_Gindex_orbit_period_velocity_wrt_StateAfterEvent_variables;
					for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent.size(); ++dIndex)
					{
						if (std::get<1>(Derivatives_of_StateAfterEvent[dIndex]) == stateIndex + 3)
						{
							state_dIndex_orbit_period_velocity_with_respect_to_StateAfterEvent.push_back(dIndex);

							this->create_sparsity_entry(this->Fdescriptions->size() - 1,
								std::get<0>(Derivatives_of_StateAfterEvent[dIndex]),
								state_Gindex_orbit_period_velocity_wrt_StateAfterEvent_variables);
						}
					}
					this->dIndex_orbit_period_velocity_wrt_StateAfterEvent.push_back(state_dIndex_orbit_period_velocity_with_respect_to_StateAfterEvent);
					this->Gindex_orbit_period_velocity_wrt_StateAfterEvent_variables.push_back(state_Gindex_orbit_period_velocity_wrt_StateAfterEvent_variables);

					// time variables
					std::vector<size_t> state_dIndex_orbit_period_velocity_wrt_StateAfterEvent_wrt_Time;
					std::vector<size_t> state_Gindex_orbit_period_velocity_wrt_StateAfterEvent_time_variables;
					for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
					{
						if (std::get<1>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]) == stateIndex + 3)
						{
							state_dIndex_orbit_period_velocity_wrt_StateAfterEvent_wrt_Time.push_back(dIndex);

							this->create_sparsity_entry(this->Fdescriptions->size() - 1,
								std::get<0>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]),
								state_Gindex_orbit_period_velocity_wrt_StateAfterEvent_time_variables);
						}
					}
					this->dIndex_orbit_period_velocity_wrt_StateAfterEvent_wrt_Time.push_back(state_dIndex_orbit_period_velocity_wrt_StateAfterEvent_wrt_Time);
					this->Gindex_orbit_period_velocity_wrt_StateAfterEvent_time_variables.push_back(state_Gindex_orbit_period_velocity_wrt_StateAfterEvent_time_variables);
				}
			} // end calcbounds()

            void OrbitPeriodConstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {
                // Step 1: get the orbital elements from the boundary event
				this->boundary_orbit_elements = this->myBoundaryEvent->get_orbit_elements_after_event(ReferenceFrame::ICRF);
				this->semimajor_axis = this->boundary_orbit_elements(0);
				doubleType sma_squared = this->semimajor_axis * this->semimajor_axis;
				doubleType sma_cubed = sma_squared * this->semimajor_axis;

				// compute the period of the orbit
				this->orbit_period = math::TwoPI * sqrt(sma_cubed / this->myUniverse->mu);

                // Step 2: apply the constraint
                F[Findex++] = this->orbit_period / this->myUniverse->TU;

                // Step 3: derivatives
                if (needG)
                {

					//Step 4.1: for safety, start by zeroing out all relevant G entries
					for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
					{
						// POSITION
						// non-time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_orbit_period_position_wrt_StateAfterEvent_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_orbit_period_position_wrt_StateAfterEvent_variables[stateIndex][entryIndex]] = 0.0;
						}

						// time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_orbit_period_position_wrt_StateAfterEvent_time_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_orbit_period_position_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex]] = 0.0;
						}

						// VELOCITY
						// non-time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_orbit_period_velocity_wrt_StateAfterEvent_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_orbit_period_velocity_wrt_StateAfterEvent_variables[stateIndex][entryIndex]] = 0.0;
						}

						// time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_orbit_period_velocity_wrt_StateAfterEvent_time_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_orbit_period_velocity_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex]] = 0.0;
						}
					}

					// Compute the actual Jacobian entries
					// we'll need the partials of the orbital elements w.r.t. the boundary Cartesian state
					this->boundary_orbit_elements_Jacobian = this->myBoundaryEvent->get_orbit_element_Jacobian_after_event(ReferenceFrame::ICRF);

					// partial of period w.r.t. SMA
					double dperiod_dSMA = (3.0 * math::PI * sma_squared / (this->myUniverse->mu * sqrt(sma_cubed / this->myUniverse->mu))) _GETVALUE;

					// first for non-time decision variables
					std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();

					for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
					{
						// POSITION
						for (size_t entryIndex = 0; entryIndex < this->dIndex_orbit_period_position_wrt_StateAfterEvent[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_orbit_period_position_wrt_StateAfterEvent[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_orbit_period_position_wrt_StateAfterEvent_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);
							
							double Gentry = dperiod_dSMA * this->boundary_orbit_elements_Jacobian(0, stateIndex) _GETVALUE * std::get<2>(Derivatives_of_StateAfterEvent[dIndex]);

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry / this->myUniverse->TU;
						}

						// VELOCITY
						for (size_t entryIndex = 0; entryIndex < this->dIndex_orbit_period_velocity_wrt_StateAfterEvent[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_orbit_period_velocity_wrt_StateAfterEvent[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_orbit_period_velocity_wrt_StateAfterEvent_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);

							double Gentry = dperiod_dSMA * this->boundary_orbit_elements_Jacobian(0, stateIndex + 3) _GETVALUE * std::get<2>(Derivatives_of_StateAfterEvent[dIndex]);

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry / this->myUniverse->TU;
						}
					}

					// now for the time decision variables
					std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

					for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
					{
						// POSITION
						for (size_t entryIndex = 0; entryIndex < this->dIndex_orbit_period_position_wrt_StateAfterEvent_wrt_Time[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_orbit_period_position_wrt_StateAfterEvent_wrt_Time[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_orbit_period_position_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);

							double Gentry = dperiod_dSMA * this->boundary_orbit_elements_Jacobian(0, stateIndex) _GETVALUE * std::get<2>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]);

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry / this->myUniverse->TU;
						}

						// VELOCITY
						for (size_t entryIndex = 0; entryIndex < this->dIndex_orbit_period_velocity_wrt_StateAfterEvent_wrt_Time[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_orbit_period_velocity_wrt_StateAfterEvent_wrt_Time[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_orbit_period_velocity_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);

							double Gentry = dperiod_dSMA * this->boundary_orbit_elements_Jacobian(0, stateIndex + 3) _GETVALUE * std::get<2>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]);

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry / this->myUniverse->TU;
						}
					}
                }
            } // end process_constraint()

            void OrbitPeriodConstraint::output(std::ofstream& outputfile)
            {
                std::string framestring = this->myUniverse->central_body_name + " " + this->ReferenceFrameStrings[this->myReferenceFrame];

                outputfile << this->myBoundaryEvent->getName() << " orbit period (days, " << framestring << "): " << this->orbit_period _GETVALUE / 86400.0 << std::endl;
            } // end output()

        } // end namespace SpecializedConstraints
    } // end namespace BoundaryEvents
} //end namespace EMTG