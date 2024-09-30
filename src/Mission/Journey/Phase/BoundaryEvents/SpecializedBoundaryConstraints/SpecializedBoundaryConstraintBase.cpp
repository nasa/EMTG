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

#include "SpecializedBoundaryConstraintBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            SpecializedBoundaryConstraintBase::SpecializedBoundaryConstraintBase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                BoundaryEventBase* myBoundaryEvent,
                const std::string& constraintDefinition ) :
                SpecializedBoundaryConstraintBase()
            {

                this->name = name;
                this->journeyIndex = journeyIndex;
                this->phaseIndex = phaseIndex;
                this->stageIndex = stageIndex;
                this->myUniverse = Universe;
                this->mySpacecraft = mySpacecraft;
                this->myOptions = myOptions;
                this->myJourneyOptions = &this->myOptions->Journeys[this->journeyIndex];
                this->myBoundaryEvent = myBoundaryEvent;
                this->constraintDefinition = constraintDefinition;
                this->myReferenceFrame = ReferenceFrame::ICRF;
                this->ReferenceFrameStrings = std::vector<std::string> ({ "ICRF", "J2000_BCI", "J2000_BCF", "TrueOfDateBCI", "TrueOfDate_BCF", "PrincipleAxes", "Topocentric", "Polar" });

				this->state_cbEvent_wrt_cbJourney = math::Matrix<doubleType>(6, 1, 0.0);
				this->d_state_cbEvent_wrt_cbJourney_dt = math::Matrix<double>(6, 1, 0.0);
				this->SpacecraftStateRelativeToCentralBodyOfJourney = math::Matrix<doubleType>(8, 1, 0.0);
            }

            void SpecializedBoundaryConstraintBase::setup_calcbounds(
                std::vector<double>* Xupperbounds,
                std::vector<double>* Xlowerbounds,
                std::vector<double>* X_scale_factors,
                std::vector<double>* Fupperbounds,
                std::vector<double>* Flowerbounds,
				std::vector<double>* F_scale_factors,
                std::vector<std::string>* Xdescriptions,
                std::vector<std::string>* Fdescriptions,
                std::vector<size_t>* iGfun,
                std::vector<size_t>* jGvar,
                std::vector<std::string>* Gdescriptions,
                std::vector<size_t>* iAfun,
                std::vector<size_t>* jAvar,
                std::vector<std::string>* Adescriptions,
                std::vector<double>* A)
            {
                this->prefix = this->name + ": ";

                this->sparsey_thing::setup_calcbounds(Xupperbounds,
                    Xlowerbounds,
                    X_scale_factors,
                    Fupperbounds,
                    Flowerbounds,
					F_scale_factors,
                    Xdescriptions,
                    Fdescriptions,
                    iGfun,
                    jGvar,
                    Gdescriptions,
                    iAfun,
                    jAvar,
                    Adescriptions,
                    A);
            }

			/**
			* Depending on the departure/arrival class/type, the state returned by BoundaryEvent's get_state_before_event() or get_state_after_event()
			* may or may not be what the user expects. The purpose of this method is do all the checks and end up with states that are
			* in ICRF w/r/t/ the central body of the journey.
			* This routine currently only handles state; does not yet handle derivatives.
			*/
			void SpecializedBoundaryConstraintBase::GetStateRelativeToCentralBodyOfJourney(const bool& needG)
			{
				math::Matrix<doubleType> SpacecraftStateRelativeToCentralBodyOfEvent;

				// EphemerisReferencedArrivalInterior needs to be handled if we are grabbing the state AFTER the event
				if (this->myBoundaryEvent->getIsEphemerisReferencedArrivalInterior()
					&& this->getStateAndDerivativesIndex == 1)
				{
					SpacecraftStateRelativeToCentralBodyOfEvent = this->myBoundaryEvent->get_state_after_event_central_body_unchanged();
				}
				else
				{
					SpacecraftStateRelativeToCentralBodyOfEvent = this->myBoundaryEvent->get_state_before_or_after_event(this->getStateAndDerivativesIndex);
				}

				// handle case where boundary event body is not the same as the central body of the journey
				// this will be the case if this->myBoundaryEventBody is not NULL, as set by SetMyBoundaryEventBody()
				if (this->myBoundaryEventBody)
				{
					doubleType temp_body_state[12];

					// get state of boundary event central body w/r/t/ journey central body
					this->myBoundaryEventBody->locate_body(SpacecraftStateRelativeToCentralBodyOfEvent(7),
						temp_body_state,
						needG,
						*this->myOptions);

					// compute the position and velocity of the event body w/r/t/ the journey central body
					for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
					{
						this->state_cbEvent_wrt_cbJourney(stateIndex) = temp_body_state[stateIndex];
						this->d_state_cbEvent_wrt_cbJourney_dt(stateIndex) = temp_body_state[stateIndex + 6] _GETVALUE;
					}
				}
				else
				{
					for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
					{
						this->state_cbEvent_wrt_cbJourney(stateIndex) = 0.0;
						this->d_state_cbEvent_wrt_cbJourney_dt(stateIndex) = 0.0;
					}
				}

				for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
				{
					this->SpacecraftStateRelativeToCentralBodyOfJourney(stateIndex) = SpacecraftStateRelativeToCentralBodyOfEvent(stateIndex) + this->state_cbEvent_wrt_cbJourney(stateIndex);
				}

				// time and mass do not change even if the central body changes
				for (size_t stateIndex : {6, 7})
				{
					this->SpacecraftStateRelativeToCentralBodyOfJourney(stateIndex) = SpacecraftStateRelativeToCentralBodyOfEvent(stateIndex);
				}
			}

			/**
			* Method for determining what the boundary event body of the boundary constraint is.
			* Also sets the class variables departureOrArrival and getStateAndDerivativesIndex.
			@param arrivalOrDepartureString The string of the boundary event definition that contains
			the text "arrival" or "departure"
			@param hardCodedBeforeOrAfterEvent Integer.
			0: getStateAndDerivativesIndex is taken from arrivalOrDepartureString
			-1: getStateAndDerivativesIndex is set to before the event
			1: getStateAndDerivativesIndex is set to after the event
			*/
			void SpecializedBoundaryConstraintBase::SetMyBoundaryEventBody(std::string arrivalOrDepartureString, int hardCodedBeforeOrAfterEvent)
			{
				/*
				reason: at a boundary, the boundary state can be referenced to a body that is not the
				central body of the journey.example : a journey whose central body is the sun, but which arrives
				at an ephemeris - referenced point on the soi of a planet.in that case, if the user grabs the state
				before the boundary, emtg gives it back sun - centered.however, if the user grabs the state after the
				boundary, emtg gives it back planet - centered.it will be the opposite for a departure.currently,
				the process_constraint method grabs states after the boundary.so, for a departure, the central body
				of the constraint is the central body of the journey. for a arrival, the central body of the constraint
				is the central body of the boundary event.
				*/
				// when we grab states for the s/c at the boundary, what body will they be with respect to?
				// it will be either the CB of the journey or the boundary event body
				// also, we have to keep track of if we have an EphemerisReferencedInterceptInterior arrival
				// if that is the case, we need to use a different state getter later on

				// it is the CB of the journey UNLESS we want the state after the event
				// at an ephemeris-referenced intercept with bounded v infinity
				if (hardCodedBeforeOrAfterEvent == -1)
				{
					this->getStateAndDerivativesIndex = -1; // before the event
				}
				else if (hardCodedBeforeOrAfterEvent == 1)
				{
					this->getStateAndDerivativesIndex = 1; // after the event
				}
				else
				{
					if (arrivalOrDepartureString.find("-") < 1024)
					{
						this->getStateAndDerivativesIndex = -1; // before the event
					}
					else if (arrivalOrDepartureString.find("+") < 1024)
					{
						this->getStateAndDerivativesIndex = 1; // after the event
					}
					else
					{
						// if not hard-coded, wwe NEED to find + or -
						throw std::logic_error(this->name + " invalid boundary before/after input for boundary constraint "
							+ this->constraintDefinition + "."
							+ " Valid values are: -, +"
							+ " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
					}
				}

				if (arrivalOrDepartureString.find("arrival") < 1024)
				{
					this->departureOrArrival = 1;
				}
				else if (arrivalOrDepartureString.find("departure") < 1024)
				{
					this->departureOrArrival = 0;
				}
				else
				{
					throw std::logic_error(this->name + " invalid boundary departure/arrival input for boundary constraint "
						+ this->constraintDefinition + "."
						+ " Valid values are: departure, arrival"
						+ " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
				}


				if (arrivalOrDepartureString.find("arrival") < 1024
					&& this->getStateAndDerivativesIndex == 1
					&& this->myJourneyOptions->arrival_class == EMTG::BoundaryClass::EphemerisReferenced
					&& this->myJourneyOptions->arrival_type == EMTG::ArrivalType::INTERCEPT)
				{
					this->myBoundaryEventBody = this->myBoundaryEvent->getBody();
				}
				else
				{
					// nothing for departure because departure state is always w/r/t/ central body
					// of journey because we are grabbing the state after the boundary event.
					// the boundary event body is just the CB
					this->myBoundaryEventBody = nullptr;
				}

				
			}
        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG