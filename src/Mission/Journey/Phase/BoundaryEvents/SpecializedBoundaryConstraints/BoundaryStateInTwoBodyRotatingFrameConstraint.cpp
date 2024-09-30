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

#include "BoundaryStateInTwoBodyRotatingFrameConstraint.h"
#include "StateInTwoBodyRotatingFrame.h"

//#define DEBUG_BOUNDARY_STATE_IN_TWO_BODY_ROTATING_FRAME

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

			BoundaryStateInTwoBodyRotatingFrameConstraint::BoundaryStateInTwoBodyRotatingFrameConstraint(const std::string& name,
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
                this->myReferenceFrame = ReferenceFrame::ICRF;
            }

			void BoundaryStateInTwoBodyRotatingFrameConstraint::calcbounds()
			{
				// parse the constraint definition
				// constraintDefinition is like j1p0_departure+_stateintwobodyrotatingframe_x_cb_2_body2_900km_1000km
				std::vector<std::string> ConstraintDefinitionCell;
				std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
				boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

				SetMyBoundaryEventBody(ConstraintDefinitionCell[1], 0); // the 0 indicates before/after event is set via the ConstraintDefinition string

				//// arrival or departure?
				//if (ConstraintDefinitionCell[1].find("departure") < 1024)
				//{
				//	this->departureOrArrival = 0;
				//}
				//else if (ConstraintDefinitionCell[1].find("arrival") < 1024)
				//{
				//	this->departureOrArrival = 1;
				//}
				//else
				//{
				//	throw std::logic_error(this->name + " invalid boundary departure/arrival input for boundary constraint on StateInTwoBodyRotatingFrame."
				//		+ " Valid values are: departure, arrival"
				//		+ " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				//}

				// before/after event?
				//if (ConstraintDefinitionCell[1].find("-") < 1024)
				//{
				//	this->getStateAndDerivativesIndex = -1; // before the event
				//}
				//else if (ConstraintDefinitionCell[1].find("+") < 1024)
				//{
				//	this->getStateAndDerivativesIndex = 1; // after the event
				//}
				//else
				//{
				//	throw std::logic_error(this->name + " invalid boundary before/after input for boundary constraint on StateInTwoBodyRotatingFrame."
				//		+ " Valid values are: -, +"
				//		+ " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				//}

				// which state element is being constrained?
				if (ConstraintDefinitionCell[3] == "x")
				{
					this->stateConstrainedIndex = 0;
				}
				else if (ConstraintDefinitionCell[3] == "y")
				{
					this->stateConstrainedIndex = 1;
				}
				else if (ConstraintDefinitionCell[3] == "z")
				{
					this->stateConstrainedIndex = 2;
				}
				else if (ConstraintDefinitionCell[3] == "vx")
				{
					this->stateConstrainedIndex = 3;
				}
				else if (ConstraintDefinitionCell[3] == "vy")
				{
					this->stateConstrainedIndex = 4;
				}
				else if (ConstraintDefinitionCell[3] == "vz")
				{
					this->stateConstrainedIndex = 5;
				}
				else
				{
					throw std::logic_error(this->name + " invalid state element input for boundary constraint on StateInTwoBodyRotatingFrame."
						+ " Valid values are: x, y, z, vx, vy, vz."
						+ " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				}

				// body 1 can't be the same as body 2
				if (ConstraintDefinitionCell[4] == ConstraintDefinitionCell[5])
				{
					throw std::logic_error(this->name + " Body1 and Body2 cannot be the same for boundary constraint on StateInTwoBodyRotatingFrame."
						+ " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				}

				// body 1
				if (ConstraintDefinitionCell[4] == "cb")
				{
					this->body1IsCB = true;
					this->body1Name = this->myUniverse->central_body_name;
				}
				else
				{
					this->body1IsCB = false;

					int bodyIndex = std::stoi(ConstraintDefinitionCell[4]) - 1;
					this->BoundaryStateInTwoBodyRotatingFrameConstraint::myBody1 = &this->myUniverse->bodies[bodyIndex];
					this->body1Name = this->BoundaryStateInTwoBodyRotatingFrameConstraint::myBody1->name;
				}

				// body 2
				if (ConstraintDefinitionCell[5] == "cb")
				{
					this->body2IsCB = true;
					this->body2Name = this->myUniverse->central_body_name;
				}
				else
				{
					this->body2IsCB = false;

					int bodyIndex = std::stoi(ConstraintDefinitionCell[5]) - 1;
					this->BoundaryStateInTwoBodyRotatingFrameConstraint::myBody2 = &this->myUniverse->bodies[bodyIndex];
					this->body2Name = this->BoundaryStateInTwoBodyRotatingFrameConstraint::myBody2->name;
				}

				// origin of the reference frame
				// currently, only valid option is body 2
				if (ConstraintDefinitionCell[6] == "body2")
				{
					this->originOfConstraintFrame = 2;
					this->nameOfOriginOfConstraintFrame = this->body2Name;
				}
				else
				{
					throw std::logic_error(this->name + " invalid constraint frame origin input for boundary constraint on StateInTwoBodyRotatingFrame."
						+ " Valid values are: body2"
						+ " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				}


				// create the constraint
				//figure out the lower and upper bounds - depends on units

				size_t lbIndex = 7;
				size_t ubIndex = 8;

				// distance units
				if (this->stateConstrainedIndex < 3)
				{
					if (ConstraintDefinitionCell[lbIndex].find("km") < 1024)
					{
						this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[lbIndex], "km")) / this->myUniverse->LU);
					}
					else if (ConstraintDefinitionCell[lbIndex].find("au") < 1024)
					{
						this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[lbIndex], "au")) * this->myOptions->AU / this->myUniverse->LU);
					}
					else if (ConstraintDefinitionCell[lbIndex].find("lu") < 1024)
					{
						this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[lbIndex], "lu")) * this->myUniverse->LU / this->myUniverse->LU);
					}
					else
					{
						return throw std::logic_error("Please specify either km, nmi, AU or LU for position boundary state in two-body rotating frame constraint bounds: " + this->Fdescriptions->back() + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
					}

					if (ConstraintDefinitionCell[ubIndex].find("km") < 1024)
					{
						this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[ubIndex], "km")) / this->myUniverse->LU);
					}
					else if (ConstraintDefinitionCell[ubIndex].find("au") < 1024)
					{
						this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[ubIndex], "au")) * this->myOptions->AU / this->myUniverse->LU);
					}
					else if (ConstraintDefinitionCell[ubIndex].find("lu") < 1024)
					{
						this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[ubIndex], "lu")) * this->myUniverse->LU / this->myUniverse->LU);
					}
					else
					{
						return throw std::logic_error("Please specify either km, nmi, AU or LU for position boundary state in two-body rotating frame constraint bounds: " + this->Fdescriptions->back() + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
					}
				}
				else
				{
					if (ConstraintDefinitionCell[lbIndex].find("km/s") < 1024)
					{
						this->Flowerbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[lbIndex], "km/s")) / (this->myUniverse->LU / this->myUniverse->TU));
					}
					else
					{
						return throw std::logic_error("Please specify km/s for velocity boundary state in two-body rotating frame constraint bounds: " + this->Fdescriptions->back() + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
					}

					if (ConstraintDefinitionCell[ubIndex].find("km/s") < 1024)
					{
						this->Fupperbounds->push_back(std::stod(boost::erase_all_copy(ConstraintDefinitionCell[ubIndex], "km/s")) / (this->myUniverse->LU / this->myUniverse->TU));
					}
					else
					{
						return throw std::logic_error("Please specify km/s for velocity boundary state in two-body rotating frame constraint bounds: " + this->Fdescriptions->back() + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
					}
				}

				// description
				std::string plusOrMinus;
				switch (this->getStateAndDerivativesIndex)
				{
				case (-1):
				{
					plusOrMinus = "-";
					break;
				}
				default:
				{
					plusOrMinus = "+";
					break;
				}
				}
				this->Fdescriptions->push_back(prefix + "constraint on " + ConstraintDefinitionCell[3] + plusOrMinus +
					+ " in two-body rotating reference frame based on rotation of " + this->body2Name
					+ " about " + this->body1Name + "; centered at " + this->nameOfOriginOfConstraintFrame);

                // sparsity pattern
				// if it is a constraint on an element of the position vector, 
				// then there are derivatives w/r/t/ position and time.
				// if it is a constraint on an element of the velocity vector, 
				// then there are derivatives w/r/t/ position, velocity, and time.

				std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAroundEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeOrAfterEvent(this->getStateAndDerivativesIndex);
				std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAroundEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeOrAfterEvent_wrt_Time(this->getStateAndDerivativesIndex);

                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
					// position: always a dependency
                    //non-time variables
					std::vector<size_t> state_dIndex_constraint_position_wrt_StateAroundEvent;
					std::vector<size_t> state_Gindex_constraint_position_wrt_StateAroundEvent_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAroundEvent.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAroundEvent[dIndex]) == stateIndex)
                        {
                            state_dIndex_constraint_position_wrt_StateAroundEvent.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAroundEvent[dIndex]),
                                state_Gindex_constraint_position_wrt_StateAroundEvent_variables);
                        }
                    }
					this->dIndex_constraint_position_wrt_StateAroundEvent.push_back(state_dIndex_constraint_position_wrt_StateAroundEvent);
					this->Gindex_constraint_position_wrt_StateAroundEvent_variables.push_back(state_Gindex_constraint_position_wrt_StateAroundEvent_variables);

                    //time variables
                    std::vector<size_t> state_dIndex_constraint_position_wrt_StateAroundEvent_wrt_Time;
                    std::vector<size_t> state_Gindex_constraint_position_wrt_StateAroundEvent_time_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAroundEvent_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAroundEvent_wrt_Time[dIndex]) == stateIndex)
                        {
                            state_dIndex_constraint_position_wrt_StateAroundEvent_wrt_Time.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAroundEvent_wrt_Time[dIndex]),
                                state_Gindex_constraint_position_wrt_StateAroundEvent_time_variables);
                        }
                    }
                    this->dIndex_constraint_position_wrt_StateAroundEvent_wrt_Time.push_back(state_dIndex_constraint_position_wrt_StateAroundEvent_wrt_Time);
                    this->Gindex_constraint_position_wrt_StateAroundEvent_time_variables.push_back(state_Gindex_constraint_position_wrt_StateAroundEvent_time_variables);

					// velocity: only a dependency if stateConstrainedIndex > 2
					if (this->stateConstrainedIndex > 2)
					{
						// non-time variables
						std::vector<size_t> state_dIndex_constraint_velocity_wrt_StateAroundEvent;
						std::vector<size_t> state_Gindex_constraint_velocity_wrt_StateAroundEvent_variables;
						for (size_t dIndex = 0; dIndex < Derivatives_of_StateAroundEvent.size(); ++dIndex)
						{
							if (std::get<1>(Derivatives_of_StateAroundEvent[dIndex]) == stateIndex + 3)
							{
								state_dIndex_constraint_velocity_wrt_StateAroundEvent.push_back(dIndex);

								this->create_sparsity_entry(this->Fdescriptions->size() - 1,
									std::get<0>(Derivatives_of_StateAroundEvent[dIndex]),
									state_Gindex_constraint_velocity_wrt_StateAroundEvent_variables);
							}
						}
						this->dIndex_constraint_velocity_wrt_StateAroundEvent.push_back(state_dIndex_constraint_velocity_wrt_StateAroundEvent);
						this->Gindex_constraint_velocity_wrt_StateAroundEvent_variables.push_back(state_Gindex_constraint_velocity_wrt_StateAroundEvent_variables);

						// time variables
						std::vector<size_t> state_dIndex_constraint_velocity_wrt_StateAroundEvent_wrt_Time;
						std::vector<size_t> state_Gindex_constraint_velocity_wrt_StateAroundEvent_time_variables;
						for (size_t dIndex = 0; dIndex < Derivatives_of_StateAroundEvent_wrt_Time.size(); ++dIndex)
						{
							if (std::get<1>(Derivatives_of_StateAroundEvent_wrt_Time[dIndex]) == stateIndex + 3)
							{
								state_dIndex_constraint_velocity_wrt_StateAroundEvent_wrt_Time.push_back(dIndex);

								this->create_sparsity_entry(this->Fdescriptions->size() - 1,
									std::get<0>(Derivatives_of_StateAroundEvent_wrt_Time[dIndex]),
									state_Gindex_constraint_velocity_wrt_StateAroundEvent_time_variables);
							}
						}
						this->dIndex_constraint_velocity_wrt_StateAroundEvent_wrt_Time.push_back(state_dIndex_constraint_velocity_wrt_StateAroundEvent_wrt_Time);
						this->Gindex_constraint_velocity_wrt_StateAroundEvent_time_variables.push_back(state_Gindex_constraint_velocity_wrt_StateAroundEvent_time_variables);
					}
                }

                // derivatives with respect to time variables that affect body positions
				// at least one of the 2 bodies is not the central body, so this dependency always exists
                std::vector<size_t> timeVariables = this->myBoundaryEvent->get_Xindices_EventRightEpoch();
                for (size_t Xindex : timeVariables)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        Xindex,
                        this->Gindex_constraint_wrt_time_variables);
                }
            }//end calcbounds()

            void BoundaryStateInTwoBodyRotatingFrameConstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {
				// get the inputs we need to call the math routine

				// state of s/c w/r/t/ central body
				math::Matrix<doubleType> stateSCWrtCBInertial(6, 1), augmentedStateSCWrtCBInertial(8, 1); // augmented also has mass and time
				augmentedStateSCWrtCBInertial = this->myBoundaryEvent->get_state_before_or_after_event(this->getStateAndDerivativesIndex);
				stateSCWrtCBInertial = augmentedStateSCWrtCBInertial.getSubMatrix1D(0, 5);
				doubleType time = augmentedStateSCWrtCBInertial(7);

#ifdef DEBUG_BOUNDARY_STATE_IN_TWO_BODY_ROTATING_FRAME
				size_t gsadIndex = 100;
				time.setDerivative(gsadIndex, 1.0);
				gsadIndex = gsadIndex + 1;
				stateSCWrtCBInertial(0).setDerivative(gsadIndex, 1.0);
				gsadIndex = gsadIndex + 1;
				stateSCWrtCBInertial(1).setDerivative(gsadIndex, 1.0);
				gsadIndex = gsadIndex + 1;
				stateSCWrtCBInertial(2).setDerivative(gsadIndex, 1.0);
				gsadIndex = gsadIndex + 1;
				stateSCWrtCBInertial(3).setDerivative(gsadIndex, 1.0);
				gsadIndex = gsadIndex + 1;
				stateSCWrtCBInertial(4).setDerivative(gsadIndex, 1.0);
				gsadIndex = gsadIndex + 1;
				stateSCWrtCBInertial(5).setDerivative(gsadIndex, 1.0);
				gsadIndex = gsadIndex + 1;
#endif

				// state of b1 w/r/t/ central body
				math::Matrix<doubleType> stateB1WrtCBInertial(6, 1), dStateB1WrtCBInertialDt(6, 1);
				if (this->body1IsCB)
				{
					for (size_t i = 0; i < 6; ++i)
					{
						stateB1WrtCBInertial(i) = 0.0;
						dStateB1WrtCBInertialDt(i) = 0.0;
					}
				}
				else
				{
					doubleType temp_body_state[12];
					this->BoundaryStateInTwoBodyRotatingFrameConstraint::myBody1->locate_body(time,
						temp_body_state,
						needG,
						*this->myOptions);
					for (size_t i = 0; i < 6; ++i)
					{
						stateB1WrtCBInertial(i) = temp_body_state[i];
						dStateB1WrtCBInertialDt(i) = temp_body_state[i + 6];
					}
				}

				// state of b2 w/r/t/ central body
				math::Matrix<doubleType> stateB2WrtCBInertial(6, 1), dStateB2WrtCBInertialDt(6, 1);
				if (this->body2IsCB)
				{
					for (size_t i = 0; i < 6; ++i)
					{
						stateB2WrtCBInertial(i) = 0.0;
						dStateB2WrtCBInertialDt(i) = 0.0;
					}
				}
				else
				{
					doubleType temp_body_state[12];
					this->BoundaryStateInTwoBodyRotatingFrameConstraint::myBody2->locate_body(time,
						temp_body_state,
						needG,
						*this->myOptions);
					for (size_t i = 0; i < 6; ++i)
					{
						stateB2WrtCBInertial(i) = temp_body_state[i];
						dStateB2WrtCBInertialDt(i) = temp_body_state[i + 6];
					}
				}

				// call the routine to get the state in the rotating frame
				math::Matrix<doubleType> stateSCWrtB2Rot(6, 1), dStateSCWrtB2RotDStateInertial(6, 6), dStateSCWrtB2RotDt(6, 1);
				Astrodynamics::TwoBodyRotatingFrame::CalculateStateInRotatingFrameWrtB2(stateB1WrtCBInertial, dStateB1WrtCBInertialDt,
					stateB2WrtCBInertial, dStateB2WrtCBInertialDt,
					stateSCWrtCBInertial,
					stateSCWrtB2Rot,
					needG, dStateSCWrtB2RotDStateInertial, dStateSCWrtB2RotDt);

#ifdef DEBUG_BOUNDARY_STATE_IN_TWO_BODY_ROTATING_FRAME
				// check state partials
				for (size_t i = 0; i < 6; ++i)
				{
					for (size_t j = 0; j < 6; ++j)
					{
						std::cout << "F(" << i << ", " << j << ") AD = " << stateSCWrtB2Rot(i).getDerivative(j + 1) << "\n";
					}
				}

				for (size_t i = 0; i < 6; ++i)
				{
					for (size_t j = 0; j < 6; ++j)
					{
						std::cout << "F(" << i << ", " << j << ") Ana. = " << dStateSCWrtB2RotDStateInertial(i, j) _GETVALUE << "\n";
					}
				}

				for (size_t i = 0; i < 6; ++i)
				{
					for (size_t j = 0; j < 6; ++j)
					{
						std::cout << "F(" << i << ", " << j << ") AD - Ana. = " << stateSCWrtB2Rot(i).getDerivative(j + 101) - dStateSCWrtB2RotDStateInertial(i, j) _GETVALUE << "\n";
					}
				}
#endif

				this->constrainedStateValue = stateSCWrtB2Rot(this->stateConstrainedIndex);

				// scale factor
				doubleType scaleFactor;
				if (this->stateConstrainedIndex < 3) // position constraint
				{
					scaleFactor = this->myUniverse->LU;
				}
				else // velocity constraint
				{
					scaleFactor = this->myUniverse->LU / this->myUniverse->TU;
				}

				// scale and add to constraint vector
				F[Findex++] = this->constrainedStateValue / scaleFactor;

				// add derivatives to overall derivatives array
				if (needG)
				{
					// start by zeroing out all relevant G entries
					for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
					{
						// POSITION
						// non-time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_position_wrt_StateAroundEvent_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_constraint_position_wrt_StateAroundEvent_variables[stateIndex][entryIndex]] = 0.0;
						}

						// time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_position_wrt_StateAroundEvent_time_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_constraint_position_wrt_StateAroundEvent_time_variables[stateIndex][entryIndex]] = 0.0;
						}

						// VELOCITY
						if (this->stateConstrainedIndex > 2)
						{
							// non-time
							for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_velocity_wrt_StateAroundEvent_variables[stateIndex].size(); ++entryIndex)
							{
								G[this->Gindex_constraint_velocity_wrt_StateAroundEvent_variables[stateIndex][entryIndex]] = 0.0;
							}

							// time
							for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_velocity_wrt_StateAroundEvent_time_variables[stateIndex].size(); ++entryIndex)
							{
								G[this->Gindex_constraint_velocity_wrt_StateAroundEvent_time_variables[stateIndex][entryIndex]] = 0.0;
							}
						}
					}

					//explicit time
					for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_time_variables.size(); ++entryIndex)
					{
						G[this->Gindex_constraint_wrt_time_variables[entryIndex]] = 0.0;
					}

					// derivatives with respect to non-time variables affecting state around boundary event
					std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAroundEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeOrAfterEvent(this->getStateAndDerivativesIndex);

					for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
					{
						// position
						for (size_t entryIndex = 0; entryIndex < this->dIndex_constraint_position_wrt_StateAroundEvent[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_constraint_position_wrt_StateAroundEvent[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_constraint_position_wrt_StateAroundEvent_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);

							double Gentry = (dStateSCWrtB2RotDStateInertial(this->stateConstrainedIndex, stateIndex) * std::get<2>(Derivatives_of_StateAroundEvent[dIndex])) _GETVALUE;

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry
								/ scaleFactor _GETVALUE;
						}

						// velocity
						if (this->stateConstrainedIndex > 2)
						{
							for (size_t entryIndex = 0; entryIndex < this->dIndex_constraint_velocity_wrt_StateAroundEvent[stateIndex].size(); ++entryIndex)
							{
								size_t dIndex = this->dIndex_constraint_velocity_wrt_StateAroundEvent[stateIndex][entryIndex];
								size_t Gindex = this->Gindex_constraint_velocity_wrt_StateAroundEvent_variables[stateIndex][entryIndex];
								size_t Xindex = this->jGvar->operator[](Gindex);

								double Gentry = (dStateSCWrtB2RotDStateInertial(this->stateConstrainedIndex, stateIndex + 3) * std::get<2>(Derivatives_of_StateAroundEvent[dIndex])) _GETVALUE;

								G[Gindex] += this->X_scale_factors->operator[](Xindex)
									* Gentry
									/ scaleFactor _GETVALUE;
							}
						}
					}

					// derivatives with respect to time variables affecting state around boundary event
					std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAroundEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeOrAfterEvent_wrt_Time(this->getStateAndDerivativesIndex);

					for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
					{
						// POSITION
						for (size_t entryIndex = 0; entryIndex < this->dIndex_constraint_position_wrt_StateAroundEvent_wrt_Time[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_constraint_position_wrt_StateAroundEvent_wrt_Time[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_constraint_position_wrt_StateAroundEvent_time_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);

							double Gentry = (dStateSCWrtB2RotDStateInertial(this->stateConstrainedIndex, stateIndex) * std::get<2>(Derivatives_of_StateAroundEvent_wrt_Time[dIndex])) _GETVALUE;

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry / scaleFactor _GETVALUE;
						}

						// VELOCITY
						if (this->stateConstrainedIndex > 2)
						{
							for (size_t entryIndex = 0; entryIndex < this->dIndex_constraint_velocity_wrt_StateAroundEvent_wrt_Time[stateIndex].size(); ++entryIndex)
							{
								size_t dIndex = this->dIndex_constraint_velocity_wrt_StateAroundEvent_wrt_Time[stateIndex][entryIndex];
								size_t Gindex = this->Gindex_constraint_velocity_wrt_StateAroundEvent_time_variables[stateIndex][entryIndex];
								size_t Xindex = this->jGvar->operator[](Gindex);

								double Gentry = (dStateSCWrtB2RotDStateInertial(this->stateConstrainedIndex, stateIndex + 3) * std::get<2>(Derivatives_of_StateAroundEvent_wrt_Time[dIndex])) _GETVALUE;

								G[Gindex] += this->X_scale_factors->operator[](Xindex)
									* Gentry / scaleFactor _GETVALUE;
							}
						}
					}

					// explicit dependencies on time from the states of body1 and body2
					// derivatives are in: dStateSCWrtB2RotDStateInertial, dStateSCWrtB2RotDt
					for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_time_variables.size(); ++entryIndex)
					{
						size_t Gindex = this->Gindex_constraint_wrt_time_variables[entryIndex];
						size_t Xindex = this->jGvar->operator[](Gindex);
						double Gentry = dStateSCWrtB2RotDt(this->stateConstrainedIndex) _GETVALUE;
						G[Gindex] += this->X_scale_factors->operator[](Xindex)
							* Gentry
							/ scaleFactor _GETVALUE;
					}
				}
            }//end process_constraint()

            void BoundaryStateInTwoBodyRotatingFrameConstraint::output(std::ofstream& outputfile)
            {
                //std::string framestring = this->myUniverse->central_body_name + " " + this->ReferenceFrameStrings[this->myReferenceFrame];
				std::vector<std::string> ConstraintDefinitionCell;
				std::string plusOrMinus;
				switch (this->getStateAndDerivativesIndex)
				{
				case (-1):
				{
					plusOrMinus = "-";
					break;
				}
				default:
				{
					plusOrMinus = "+";
					break;
				}
				}

				std::string units;
				if (this->stateConstrainedIndex > 2)
				{
					units = "km/s";
				}
				else
				{
					units = "km";
				}


				std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
				boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);
				outputfile << this->myBoundaryEvent->getName() << " " << ConstraintDefinitionCell[3] << plusOrMinus << " in frame rotating with " << this->body2Name << " about " << this->body1Name << "; centered at " << this->nameOfOriginOfConstraintFrame << " (" << units << "): " << this->constrainedStateValue _GETVALUE << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG