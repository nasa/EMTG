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

#include "BoundaryVelocityDeclinationAnybodyConstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            BoundaryVelocityDeclinationAnybodyConstraint::BoundaryVelocityDeclinationAnybodyConstraint(const std::string& name,
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
				this->VelocityAroundEventConstraintFrame = math::Matrix<doubleType>(3, 1, 0.0);
				this->dVelocityAroundEventConstraintFrame_dt = math::Matrix<doubleType>(3, 1, 0.0);
				this->V_cbEvent_wrt_cbJourney = math::Matrix<doubleType>(3, 1, 0.0);
				this->dV_cbEvent_wrt_cbJourney_dt = math::Matrix<doubleType>(3, 1, 0.0);
				this->V_cbConstraint_wrt_cbJourney = math::Matrix<doubleType>(3, 1, 0.0);
				this->dV_cbConstraint_wrt_cbJourney_dt = math::Matrix<doubleType>(3, 1, 0.0);
                this->myReferenceFrame = constraint_reference_frame;
            }//end constructor

            void BoundaryVelocityDeclinationAnybodyConstraint::calcbounds()
            {
                // Step 1: parse the constraint definition
				// constraintDefinition is like p0_departure+_velocitydeclinationanybody_bodyID_84.0deg_85.0deg_j2000bcf
                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

				// arrival or departure?
				bool classTypeComboOK = false;
				if (ConstraintDefinitionCell[1].find("departure") < 1024)
				{
					//this->departureOrArrival = 0; // moved setting to SetMyBoundaryEventBody

					// check that the departure class/type combination is supported
					if (this->myJourneyOptions->departure_type == EMTG::DepartureType::LAUNCH_OR_DIRECT_INSERTION)
					{
						if (this->myJourneyOptions->departure_class == EMTG::BoundaryClass::EphemerisPegged ||
							this->myJourneyOptions->departure_class == EMTG::BoundaryClass::FreePoint ||
							this->myJourneyOptions->departure_class == EMTG::BoundaryClass::Periapse)
						{
							classTypeComboOK = true;
						}
					}
					else if (this->myJourneyOptions->departure_type == EMTG::DepartureType::FREE_DIRECT_DEPARTURE)
					{
						if (this->myJourneyOptions->departure_class == EMTG::BoundaryClass::FreePoint)
						{
							classTypeComboOK = true;
						}
					}
					else if (this->myJourneyOptions->departure_type == EMTG::DepartureType::FLYBY)
					{
						if (this->myJourneyOptions->departure_class == EMTG::BoundaryClass::EphemerisPegged)
						{
							classTypeComboOK = true;
						}
					}

					if (!classTypeComboOK)
					{
						throw std::logic_error(this->name + " invalid departure class/type combination for VelocityDeclinationAnybody."
							+ " Departure type is " + EMTG::DepartureTypeStrings[this->myJourneyOptions->departure_type] + " and departure class is " + EMTG::BoundaryClassStrings[this->myJourneyOptions->departure_class] + "."
							+ " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
					}
				}
				else if (ConstraintDefinitionCell[1].find("arrival") < 1024)
				{
					//this->departureOrArrival = 1; // moved setting to SetMyBoundaryEventBody

					// check that the arrival class/type combination is supported
					if (this->myJourneyOptions->arrival_type == EMTG::ArrivalType::INSERTION_INTO_PARKING)
					{
						if (this->myJourneyOptions->arrival_class == EMTG::BoundaryClass::EphemerisPegged)
						{
							classTypeComboOK = true;
						}
					}
					else if (this->myJourneyOptions->arrival_type == EMTG::ArrivalType::CHEM_RENDEZVOUS)
					{
						if (this->myJourneyOptions->arrival_class == EMTG::BoundaryClass::FreePoint)
						{
							classTypeComboOK = true;
						}
					}
					else if (this->myJourneyOptions->arrival_type == EMTG::ArrivalType::INTERCEPT)
					{
						if (this->myJourneyOptions->arrival_class == EMTG::BoundaryClass::EphemerisPegged ||
							this->myJourneyOptions->arrival_class == EMTG::BoundaryClass::FreePoint ||
							this->myJourneyOptions->arrival_class == EMTG::BoundaryClass::EphemerisReferenced ||
							this->myJourneyOptions->arrival_class == EMTG::BoundaryClass::Periapse)
						{
							classTypeComboOK = true;
						}
					}
					else if (this->myJourneyOptions->arrival_type == EMTG::ArrivalType::LT_RENDEZVOUS)
					{
						if (this->myJourneyOptions->arrival_class == EMTG::BoundaryClass::FreePoint)
						{
							classTypeComboOK = true;
						}
					}

					if (!classTypeComboOK)
					{
						throw std::logic_error(this->name + " invalid arrival class/type combination for VelocityDeclinationAnybody."
							+ " Arrival type is " + EMTG::ArrivalTypeStrings[this->myJourneyOptions->arrival_type] + " and arrival class is " + EMTG::BoundaryClassStrings[this->myJourneyOptions->arrival_class] + "."
							+ " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
					}
				}
				else
				{
					throw std::logic_error(this->name + " invalid boundary departure/arrival input for boundary constraint on VelocityDeclinationAnybody."
						+ " Valid values are: departure, arrival"
						+ " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				}

				// before/after event?
				// now handled in SetMyBoundaryEventBody
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
				//	throw std::logic_error(this->name + " invalid boundary before/after input for boundary constraint on VelocityDeclinationAnybody."
				//		+ " Valid values are: -, +"
				//		+ " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				//}

				// when we grab states for the s/c at the boundary, what body will they be with respect to?
				// it will be either the CB of the journey or the boundary event body

				// it is the CB of the journey UNLESS we want the state after the event
				// at an ephemeris-referenced intercept with bounded v infinity
				SetMyBoundaryEventBody(ConstraintDefinitionCell[1], 0); // the 0 indicates before/after event is set via the ConstraintDefinition string

				// SetMyBoundaryEventBody sets getStateAndDerivativesIndex,
				// so we can validate it now
				if (this->getStateAndDerivativesIndex != -1 && this->getStateAndDerivativesIndex != 1)
				{
					throw std::logic_error(this->name + " invalid boundary before/after input for boundary constraint on VelocityDeclinationAnybody."
						+ " Valid values are: -, +"
						+ " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				}

				//if (ConstraintDefinitionCell[1].find("arrival") < 1024
				//	&& this->getStateAndDerivativesIndex == 1
				//	&& this->myJourneyOptions->arrival_class == EMTG::BoundaryClass::EphemerisReferenced
				//	&& this->myJourneyOptions->arrival_type == EMTG::ArrivalType::INTERCEPT)
				//{
				//	this->myBoundaryEventBody = this->myBoundaryEvent->getBody();
				//}
				//else
				//{
				//	this->myBoundaryEventBody = nullptr;
				//}

				// body w/r/t/ which velocity and velocity declination is to be calculated
				if (ConstraintDefinitionCell[3] == "cb")
				{
					this->bodyIsCBOfJourney = true;
					this->bodyName = this->myUniverse->central_body_name;
				}
				else
				{
					this->bodyIsCBOfJourney = false;

					int bodyIndex = std::stoi(ConstraintDefinitionCell[3]) - 1;
					this->BoundaryVelocityDeclinationAnybodyConstraint::myBody = &this->myUniverse->bodies[bodyIndex];
					this->bodyName = this->BoundaryVelocityDeclinationAnybodyConstraint::myBody->name;
				}

                // lower bound
				std::string lower_bound = ConstraintDefinitionCell[4];
				if (lower_bound.find("deg") < 1024)
				{
					double lower_bound_float = std::stod(boost::erase_all_copy(lower_bound, "deg"));
					this->Flowerbounds->push_back(lower_bound_float * math::deg2rad);
				}
				else if (lower_bound.find("rad") < 1024)
				{
					double lower_bound_float = std::stod(boost::erase_all_copy(lower_bound, "rad"));
					this->Flowerbounds->push_back(lower_bound_float);
				}
				else
				{
					return throw std::logic_error("Please specify either deg or rad for VelocityDeclinationAnybody constraint lower bound. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				}

				// upper bound
				std::string upper_bound = ConstraintDefinitionCell[5];
				if (upper_bound.find("deg") < 1024)
				{
					double upper_bound_float = std::stod(boost::erase_all_copy(upper_bound, "deg"));
					this->Fupperbounds->push_back(upper_bound_float * math::deg2rad);
				}
				else if (upper_bound.find("rad") < 1024)
				{
					double upper_bound_float = std::stod(boost::erase_all_copy(upper_bound, "rad"));
					this->Fupperbounds->push_back(upper_bound_float);
				}
				else
				{
					return throw std::logic_error("Please specify either deg or rad for VelocityDeclinationAnybody constraint upper bound. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
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
				this->Fdescriptions->push_back(prefix + "constraint on velocity declination"  + plusOrMinus +
					+ " with respect to " + this->bodyName + " in frame " + ConstraintDefinitionCell[6]);


                //Step 5: sparsity pattern
                //derivatives with respect to anything influencing the boundary event's right-hand velocity vector
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAroundEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeOrAfterEvent(this->getStateAndDerivativesIndex);
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAroundEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeOrAfterEvent_wrt_Time(this->getStateAndDerivativesIndex);

                for (size_t stateIndex = 3; stateIndex < 6; ++stateIndex)
                {
                    //non-time variables
                    std::vector<size_t> state_dIndex_with_respect_to_StateAroundEvent;
                    std::vector<size_t> state_Gindex_constraint_wrt_StateAroundEvent_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAroundEvent.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAroundEvent[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_StateAroundEvent.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAroundEvent[dIndex]),
                                state_Gindex_constraint_wrt_StateAroundEvent_variables);
                        }
                    }
                    this->dIndex_with_respect_to_StateAroundEvent.push_back(state_dIndex_with_respect_to_StateAroundEvent);
                    this->Gindex_constraint_wrt_StateAroundEvent_variables.push_back(state_Gindex_constraint_wrt_StateAroundEvent_variables);

                    //time variables
                    std::vector<size_t> state_dIndex_with_respect_to_StateAroundEvent_wrt_Time;
                    std::vector<size_t> state_Gindex_constraint_wrt_StateAroundEvent_time_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAroundEvent_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAroundEvent_wrt_Time[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_StateAroundEvent_wrt_Time.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAroundEvent_wrt_Time[dIndex]),
                                state_Gindex_constraint_wrt_StateAroundEvent_time_variables);
                        }
                    }
                    this->dIndex_with_respect_to_StateAroundEvent_wrt_Time.push_back(state_dIndex_with_respect_to_StateAroundEvent_wrt_Time);
                    this->Gindex_constraint_wrt_StateAroundEvent_time_variables.push_back(state_Gindex_constraint_wrt_StateAroundEvent_time_variables);
                }

                //derivatives with respect to time variables that affect current epoch
				// can we get rid of these if it is a central body constraint?
                std::vector<size_t> timeVariables = this->myBoundaryEvent->get_Xindices_EventRightEpoch();
                for (size_t Xindex : timeVariables)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        Xindex,
                        this->Gindex_constraint_wrt_time_variables);
                }
            }//end calcbounds()

            void BoundaryVelocityDeclinationAnybodyConstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {
				// get the state relative to the central body of boundary event
				math::Matrix<doubleType>& SpacecraftStateRelativeToCentralBodyOfEvent = this->myBoundaryEvent->get_state_before_or_after_event(this->getStateAndDerivativesIndex);
				math::Matrix<doubleType>& SpacecraftStateRelativeToCentralBodyOfEventBefore = this->myBoundaryEvent->get_state_before_or_after_event(-1);
				math::Matrix<doubleType>& SpacecraftStateRelativeToCentralBodyOfEventAfter = this->myBoundaryEvent->get_state_before_or_after_event(1);
				math::Matrix<doubleType> SpacecraftVelocityRelativeToCentralBodyOfJourney(3, 1, 0.0);

				// if boundary event body is not necessarily the same as the central body of the journey ...
				if (this->myBoundaryEventBody)
				{
					doubleType temp_body_state[12];

					// get state of boundary event central body w/r/t/ journey central body
					this->myBoundaryEventBody->locate_body(SpacecraftStateRelativeToCentralBodyOfEvent(7),
						temp_body_state,
						needG,
						*this->myOptions);

					// compute the velocity of the event body w/r/t/ the journey central body
					for (size_t stateIndex : {0, 1, 2})
					{
						this->V_cbEvent_wrt_cbJourney(stateIndex) = temp_body_state[stateIndex + 3];
						this->dV_cbEvent_wrt_cbJourney_dt(stateIndex) = temp_body_state[stateIndex + 9];
					}
				}
				else
				{
					for (size_t stateIndex : {0, 1, 2})
					{
						this->V_cbEvent_wrt_cbJourney(stateIndex) = 0.0;
						this->dV_cbEvent_wrt_cbJourney_dt(stateIndex) = 0.0;
					}
				}

				// also need velocity of constraint body w/r/t/ central body of journey unless
				// the constraint body is the central body of the journey
				if (this->bodyIsCBOfJourney)
				{
					for (size_t stateIndex : {0, 1, 2})
					{
						this->V_cbConstraint_wrt_cbJourney(stateIndex) = 0.0;
						this->dV_cbConstraint_wrt_cbJourney_dt(stateIndex) = 0.0;
					}
				}
				else
				{
					doubleType temp_body_state[12];

					// get state of boundary event central body w/r/t/ journey central body
					this->BoundaryVelocityDeclinationAnybodyConstraint::myBody->locate_body(SpacecraftStateRelativeToCentralBodyOfEvent(7),
						temp_body_state,
						needG,
						*this->myOptions);

					// compute the velocity of the event body w/r/t/ the journey central body
					for (size_t stateIndex : {0, 1, 2})
					{
						this->V_cbConstraint_wrt_cbJourney(stateIndex) = temp_body_state[stateIndex + 3];
						this->dV_cbConstraint_wrt_cbJourney_dt(stateIndex) = temp_body_state[stateIndex + 9];
					}
				}

				// get state of s/c w/r/t/ constraint central body:
				// s = spacecraft
				// c = constraint body
				// j = journey body
				// e = boundary event body
				// / = "with respect to"
				// v_s/c = v_s/e + v_e/j - v_c/j
				// only actually updating velocity because that's all we use
				math::Matrix<doubleType> SpacecraftVelocityRelativeToConstraintBody(3, 1, 0.0);
				for (size_t stateIndex : {0, 1, 2})
				{
					SpacecraftVelocityRelativeToConstraintBody(stateIndex) = 
						SpacecraftStateRelativeToCentralBodyOfEvent(stateIndex + 3) + 
						this->V_cbEvent_wrt_cbJourney(stateIndex) -
						this->V_cbConstraint_wrt_cbJourney(stateIndex);
				}

				// we now have velocity w/r/t/ the constraint body in ICRF
				// now need to convert to reference frame of choice
				math::Matrix<doubleType> R_from_ICRF_to_ConstraintFrame(3, 3, 0.0);
				if (this->bodyIsCBOfJourney)
				{
					// if the constraint body is the central body of the journey, its rotations are stored in universe's LocalFrame
					this->myUniverse->LocalFrame.construct_rotation_matrices(SpacecraftStateRelativeToCentralBodyOfEvent(7), needG);
					R_from_ICRF_to_ConstraintFrame = this->myUniverse->LocalFrame.get_R(ReferenceFrame::ICRF, this->myReferenceFrame);
				}
				else
				{
					// if the constraint body is not the central body of the journey, its rotations are stored in
					// myBody's body_frame
					this->myBody->body_frame.construct_rotation_matrices(SpacecraftStateRelativeToCentralBodyOfEvent(7), needG);
					R_from_ICRF_to_ConstraintFrame = this->myBody->body_frame.get_R(ReferenceFrame::ICRF, this->myReferenceFrame);
				}
                
				// rotate into constraint frame
                this->VelocityAroundEventConstraintFrame = R_from_ICRF_to_ConstraintFrame * SpacecraftVelocityRelativeToConstraintBody;


                // compute declination
                doubleType vxy_ConstraintFrame = sqrt(this->VelocityAroundEventConstraintFrame(0) * this->VelocityAroundEventConstraintFrame(0) + this->VelocityAroundEventConstraintFrame(1) * this->VelocityAroundEventConstraintFrame(1));
                this->Declination = atan2(this->VelocityAroundEventConstraintFrame(2), vxy_ConstraintFrame);

				// add to constraint vector
                F[Findex++] = this->Declination;

                //Step 4: derivatives
                if (needG)
                {
                    //Step 4.1: we are going to be adding components on top of components, so let's start by zeroing out all of our G entries
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_StateAroundEvent_variables[stateIndex].size(); ++entryIndex)
                        {
                            G[this->Gindex_constraint_wrt_StateAroundEvent_variables[stateIndex][entryIndex]] = 0.0;
                        }
                        for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_StateAroundEvent_time_variables[stateIndex].size(); ++entryIndex)
                        {
                            G[this->Gindex_constraint_wrt_StateAroundEvent_time_variables[stateIndex][entryIndex]] = 0.0;
                        }
                    }
                    for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_time_variables.size(); ++entryIndex)
                    {
                        G[this->Gindex_constraint_wrt_time_variables[entryIndex]] = 0.0;
                    }

                    // derivatives of declination w.r.t. states in constraint frame
                    double dDeclination_dvzConstraintFrame = (vxy_ConstraintFrame / (this->VelocityAroundEventConstraintFrame(2) * this->VelocityAroundEventConstraintFrame(2) + vxy_ConstraintFrame * vxy_ConstraintFrame))_GETVALUE;
                    double dDeclination_dvxy_ConstraintFrame = (-this->VelocityAroundEventConstraintFrame(2) / (this->VelocityAroundEventConstraintFrame(2) * this->VelocityAroundEventConstraintFrame(2) + vxy_ConstraintFrame * vxy_ConstraintFrame))_GETVALUE;
                    double dDeclination_dvxConstraintFrame = dDeclination_dvxy_ConstraintFrame * (this->VelocityAroundEventConstraintFrame(0) / vxy_ConstraintFrame) _GETVALUE;
                    double dDeclination_dvyConstraintFrame = dDeclination_dvxy_ConstraintFrame * (this->VelocityAroundEventConstraintFrame(1) / vxy_ConstraintFrame) _GETVALUE;

                    // derivatives with respect to non-time variables affecting state around boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAroundEvent = this->myBoundaryEvent->get_Derivatives_of_StateBeforeOrAfterEvent(this->getStateAndDerivativesIndex);

                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        double dvxConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(0, stateIndex) _GETVALUE;
                        double dvyConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(1, stateIndex) _GETVALUE;
                        double dvzConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(2, stateIndex) _GETVALUE;

                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_StateAroundEvent[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_StateAroundEvent[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_StateAroundEvent_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (dDeclination_dvxConstraintFrame * dvxConstraintFrame_dstateICRF + dDeclination_dvyConstraintFrame * dvyConstraintFrame_dstateICRF + dDeclination_dvzConstraintFrame * dvzConstraintFrame_dstateICRF)
                                * std::get<2>(Derivatives_of_StateAroundEvent[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }

                    // derivatives with respect to time variables affecting state around boundary event
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAroundEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateBeforeOrAfterEvent_wrt_Time(this->getStateAndDerivativesIndex);

                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        double dvxConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(0, stateIndex) _GETVALUE;
                        double dvyConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(1, stateIndex) _GETVALUE;
                        double dvzConstraintFrame_dstateICRF = R_from_ICRF_to_ConstraintFrame(2, stateIndex) _GETVALUE;

                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_StateAroundEvent_wrt_Time[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_StateAroundEvent_wrt_Time[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_StateAroundEvent_time_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = (dDeclination_dvxConstraintFrame * dvxConstraintFrame_dstateICRF + dDeclination_dvyConstraintFrame * dvyConstraintFrame_dstateICRF + dDeclination_dvzConstraintFrame * dvzConstraintFrame_dstateICRF)
                                * std::get<2>(Derivatives_of_StateAroundEvent_wrt_Time[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;
                        }
                    }


                    // derivatives of declination with respect to frame rotation due to change in epoch
					math::Matrix<doubleType> dR_from_ICRF_to_ConstraintFrame_dt(3, 3, 0.0);
					if (this->bodyIsCBOfJourney)
					{
						// if the constraint body is the central body of the journey, its rotations are stored in universe's LocalFrame
						dR_from_ICRF_to_ConstraintFrame_dt = this->myUniverse->LocalFrame.get_dRdt(ReferenceFrame::ICRF, this->myReferenceFrame);
					}
					else
					{
						// if the constraint body is not the central body of the journey, its rotations are stored in
						// myBody's body_frame
						dR_from_ICRF_to_ConstraintFrame_dt = this->myBody->body_frame.get_dRdt(ReferenceFrame::ICRF, this->myReferenceFrame);
					}
                    this->dVelocityAroundEventConstraintFrame_dt = dR_from_ICRF_to_ConstraintFrame_dt * SpacecraftVelocityRelativeToConstraintBody + 
						R_from_ICRF_to_ConstraintFrame * dV_cbEvent_wrt_cbJourney_dt - 
						R_from_ICRF_to_ConstraintFrame * dV_cbConstraint_wrt_cbJourney_dt;

                    double timeDerivative = (dDeclination_dvxConstraintFrame * this->dVelocityAroundEventConstraintFrame_dt(0)
                        + dDeclination_dvyConstraintFrame * this->dVelocityAroundEventConstraintFrame_dt(1)
                        + dDeclination_dvzConstraintFrame * this->dVelocityAroundEventConstraintFrame_dt(2)) _GETVALUE;

                    for (size_t entryIndex = 0; entryIndex < this->Gindex_constraint_wrt_time_variables.size(); ++entryIndex)
                    {
                        size_t Gindex = this->Gindex_constraint_wrt_time_variables[entryIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * timeDerivative;
                    }
                }//end derivatives
            }//end process_constraint()

            void BoundaryVelocityDeclinationAnybodyConstraint::output(std::ofstream& outputfile)
            {
				std::vector<std::string> ConstraintDefinitionCell;
				std::string plusOrMinus;
				std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
				boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

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

				outputfile << this->myBoundaryEvent->getName() << " " << "velocity declination" << plusOrMinus << " relative to " << this->bodyName <<"'s " << this->ReferenceFrameStrings[this->myReferenceFrame] << " frame (degrees): " << this->Declination _GETVALUE / math::deg2rad << std::endl;

                //std::string framestring = this->myUniverse->central_body_name + " " + this->ReferenceFrameStrings[this->myReferenceFrame];

                //outputfile << this->myBoundaryEvent->getName() << " velocity declination [anybody] (degrees " << framestring << "): " << this->Declination _GETVALUE / math::deg2rad  << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG