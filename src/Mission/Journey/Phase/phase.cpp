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

//EMTGv9 phase class

#include "phase.h"

#include "EphemerisPeggedDeparture/EphemerisPeggedLaunchDirectInsertion.h"
#include "EphemerisPeggedDeparture/EphemerisPeggedFreeDirectDeparture.h"
#include "EphemerisPeggedDeparture/EphemerisPeggedZeroTurnFlyby.h"
#include "EphemerisPeggedDeparture/EphemerisPeggedUnpoweredFlyby.h"
#include "EphemerisPeggedDeparture/EphemerisPeggedPoweredFlyby.h"
#include "EphemerisPeggedDeparture/EphemerisPeggedSpiralDeparture.h"

#include "EphemerisPeggedArrival/EphemerisPeggedLTRendezvous.h"
#include "EphemerisPeggedArrival/EphemerisPeggedFlybyIn.h"
#include "EphemerisPeggedArrival/EphemerisPeggedIntercept.h"
#include "EphemerisPeggedArrival/EphemerisPeggedChemRendezvous.h"
#include "EphemerisPeggedArrival/EphemerisPeggedOrbitInsertion.h"
#include "EphemerisPeggedArrival/EphemerisPeggedArrivalWithVinfinity.h"
#include "EphemerisPeggedArrival/EphemerisPeggedMomentumTransfer.h"
#include "EphemerisPeggedArrival/EphemerisPeggedSpiralArrival.h"

#include "FreePointDeparture/FreePointDirectInsertion.h"
#include "FreePointDeparture/FreePointFreeDirectDeparture.h"

#include "FreePointArrival/FreePointLTRendezvous.h"
#include "FreePointArrival/FreePointChemRendezvous.h"
#include "FreePointArrival/FreePointIntercept.h"

#include "EphemerisReferencedDeparture/Exterior/EphemerisReferencedFreeDirectDepartureExterior.h"

#include "EphemerisReferencedDeparture/Interior/EphemerisReferencedFreeDirectDepartureInterior.h"

#include "EphemerisReferencedArrival/Exterior/EphemerisReferencedLTRendezvousExterior.h"
#include "EphemerisReferencedArrival/Exterior/EphemerisReferencedInterceptExterior.h"

#include "EphemerisReferencedArrival/Interior/EphemerisReferencedLTRendezvousInterior.h"
#include "EphemerisReferencedArrival/Interior/EphemerisReferencedInterceptInterior.h"

#include "PeriapseDeparture/PeriapseLaunch.h"

#include "PeriapseArrival/PeriapseFlybyIn.h"

#include "EMTG_solver_utilities.h"
#include "EMTG_math.h"
#include "mjd_to_mdyhms.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <exception>

namespace EMTG
{
    namespace Phases
    {
        phase::phase() :
            name("dummy_phase"),
            stateVectorNames(std::vector<std::string>({ "x", "y", "z", "xdot", "ydot", "zdot", "mass", "epoch" })),
            journeyIndex(0),
            phaseIndex(0),
            stageIndex(0),
            PhaseTotalDeterministicDeltav(0.0),
            virtual_electric_propellant_used(0.0),
            virtual_chemical_fuel_used(0.0),
            virtual_chemical_oxidizer_used(0.0),
            electric_propellant_used(0.0),
            chemical_fuel_used(0.0),
            chemical_oxidizer_used(0.0),
            hasElectricManeuver(false),
            hasMonopropManeuver(false),
            hasBipropManeuver(false),
            isFirstPhaseInMission(false),
            isFirstPhaseInJourney(false),
            isLastPhaseInJourney(false),
            StageBeforePhase(false),
            hasInitialTCM(false),
            initial_TCM_magnitude(0.0),
            hasInitialCoast(false),
            InitialCoastDuration(0.0),
            hasTerminalCoast(false),
            TerminalCoastDuration(0.0),
            numMatchConstraints(7),
            PhaseDutyCycle(1.0),
            myArrivalEvent(nullptr),
            myDepartureEvent(nullptr)
        {
        }//end default constructor   

        phase::phase(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            phase* previousPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            HardwareModels::LaunchVehicle* myLaunchVehicle,
            missionoptions* myOptions,
            const size_t& numStatesToPropagate,
            const size_t& numMatchConstraints) :
            phase()
        {

            this->name = name;
            this->journeyIndex = journeyIndex;
            this->phaseIndex = phaseIndex;
            this->mySpacecraft = mySpacecraft;
            this->myLaunchVehicle = myLaunchVehicle;
            this->myJourneyOptions = &myOptions->Journeys[this->journeyIndex];

            if (this->myJourneyOptions->override_PropagatorType)
            {
                this->myPropagatorType = this->myJourneyOptions->propagatorType;
            }
            else
            {
                this->myPropagatorType = myOptions->propagatorType;
            }

            //initialize the writey thing
            this->writey_thing::initialize(myOptions, Universe);

            this->EphemerisOutputResolution = 2.0 * math::PI * sqrt(this->myUniverse->central_body.radius * this->myUniverse->central_body.radius * this->myUniverse->central_body.radius / this->myUniverse->central_body.mu) / 100.0;

            //continuity constraint setup
            this->numStatesToPropagate = numStatesToPropagate;
            this->numMatchConstraints = numMatchConstraints;
            this->continuity_constraint_scale_factors.resize(this->numMatchConstraints, 1);
            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                this->continuity_constraint_scale_factors(stateIndex) = this->myUniverse->continuity_constraint_scale_factors(stateIndex);

            //initialize the derivative truth table to all true for now - the individual phase types can change this
            this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates.resize(this->numMatchConstraints, std::vector<size_t>(this->numStatesToPropagate, true));
            this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates.resize(this->numMatchConstraints, std::vector<size_t>(this->numStatesToPropagate, true));
            this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime.resize(this->numMatchConstraints, true);

            //state holders
            this->state_at_beginning_of_phase.resize(this->numStatesToPropagate, 1, 0.0);
            this->state_at_end_of_phase.resize(this->numStatesToPropagate, 1, 0.0);
            this->state_after_initial_TCM.resize(this->numStatesToPropagate, 1, 0.0);

            //some flags
            if (this->phaseIndex == 0)
            {
                this->isFirstPhaseInJourney = true;

                if (this->journeyIndex == 0)
                    this->isFirstPhaseInMission = true;
                else if (this->myOptions->Journeys[this->journeyIndex - 1].stage_after_arrival)
                    this->StageBeforePhase = false;
            }
            
            if (this->phaseIndex == this->myJourneyOptions->number_of_phases - 1)
                this->isLastPhaseInJourney = true;

            if (this->stageIndex >= this->mySpacecraft->getNumberOfStages())
            {
                throw std::invalid_argument(this->name + " has invalid stage index " + std::to_string(this->stageIndex) + ". The spacecraft has run out of stages. Check your journey staging options. If you want to debug, place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }


            //create departure event
            this->PreviousPhaseArrivalEvent = this->isFirstPhaseInMission ? 
                NULL
                :
                previousPhase->getArrivalEvent();

            //throw an error if the previous phase ended in an ephemeris-referenced whathaveyou and this phase doesn't start with free point free direct departure
            if (this->journeyIndex > 0)
            {
                //this is not the first journey, then we know that there was a previous journey and the following is safe 
                if (this->myOptions->Journeys[this->journeyIndex - 1].arrival_class == BoundaryClass::EphemerisReferenced)
                {
                    if (!(this->myJourneyOptions->departure_class == BoundaryClass::FreePoint && this->myJourneyOptions->departure_type == DepartureType::FREE_DIRECT_DEPARTURE))
                    {
                        throw std::invalid_argument("Journey " + std::to_string(this->journeyIndex) + " needs to start with a FreePointFreeDirectDeparture because Journey " + std::to_string(this->journeyIndex - 1) + " ends with an EphemerisReferencedArrival.");
                    }
                }
            }

            if (this->isFirstPhaseInJourney)
            {
                if (this->myJourneyOptions->departure_class == BoundaryClass::EphemerisPegged)
                {

                    if (this->myJourneyOptions->sequence[phaseIndex] < 1)
                    {
                        throw std::invalid_argument("j" + std::to_string(this->journeyIndex) + "p" + std::to_string(this->phaseIndex) + ": you have selected an ephemeris pegged departure event that occurs at either the central body or the universe's sphere of influence. Neither of these is allowed. Ephemeris pegged boundary events must occur at a body in the universe.");
                    }

                    switch (this->myJourneyOptions->departure_type)
                    {
                        case LAUNCH_OR_DIRECT_INSERTION:
                            this->myDepartureEvent = new BoundaryEvents::EphemerisPeggedLaunchDirectInsertion(name + "EphemerisPeggedLaunchDirectInsertion",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myLaunchVehicle,
                                this->myOptions,
                                PreviousPhaseArrivalEvent);
                            break;
                        case DEPART_PARKING_ORBIT:
                            throw std::invalid_argument("EphemerisPeggedDepartParkingOrbit not implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                            break;
                            /* case DEPART_PARKING_ORBIT:
                                    this->myDepartureEvent = new BoundaryEvents::OrbitDeparture2D(name + "OrbitDeparture2D",
                                        this->journeyIndex,
                                        this->phaseIndex,
                                        stageIndex,
                                        this->myUniverse,
                                        this->mySpacecraft,
                                        this->myOptions,
                                        PreviousPhaseArrivalEvent);
                                    break;
                                    */
                        case FREE_DIRECT_DEPARTURE:
                            this->myDepartureEvent = new BoundaryEvents::EphemerisPeggedFreeDirectDeparture(name + "EphemerisPeggedFreeDirectDeparture",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions,
                                PreviousPhaseArrivalEvent);
                            break;
                        case FLYBY:
                            if (this->myJourneyOptions->enable_periapse_burns)
                            {
								this->myDepartureEvent = new BoundaryEvents::EphemerisPeggedPoweredFlyby(name + "EphemerisPeggedPoweredFlyby",
									this->journeyIndex,
									this->phaseIndex,
									stageIndex,
									this->myUniverse,
									this->mySpacecraft,
									this->myOptions,
									dynamic_cast<BoundaryEvents::EphemerisPeggedArrivalWithVinfinity*>(PreviousPhaseArrivalEvent));
                            }
                            else
                            {
                                this->myDepartureEvent = new BoundaryEvents::EphemerisPeggedUnpoweredFlyby(name + "EphemerisPeggedUnpoweredFlyby",
                                    this->journeyIndex,
                                    this->phaseIndex,
                                    stageIndex,
                                    this->myUniverse,
                                    this->mySpacecraft,
                                    this->myOptions,
									dynamic_cast<BoundaryEvents::EphemerisPeggedArrivalWithVinfinity*>(PreviousPhaseArrivalEvent));
                            }
                            break;
                        case FLYBY_FIXED_VINF:
                            throw std::invalid_argument("EphemerisPeggedFlybyFixedVinf not implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                            break;
                        case SPIRAL_ESCAPE:
                        {

                            this->myDepartureEvent = new BoundaryEvents::EphemerisPeggedSpiralDeparture(name + "EphemerisPeggedSpiralDeparture",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions,
                                PreviousPhaseArrivalEvent);

                            break;
                        }
                        case ZERO_TURN_FLYBY:
							this->myDepartureEvent = new BoundaryEvents::EphemerisPeggedZeroTurnFlyby(name + "EphemerisPeggedZeroTurnFlyby",
								this->journeyIndex,
								this->phaseIndex,
								stageIndex,
								this->myUniverse,
								this->mySpacecraft,
								this->myOptions,
								dynamic_cast<BoundaryEvents::EphemerisPeggedArrivalWithVinfinity*>(PreviousPhaseArrivalEvent));
                            break;
                        default:
                            throw std::invalid_argument("Undefined ephemeris pegged departure type '" + DepartureTypeStrings[this->myJourneyOptions->departure_type] + "', j" + std::to_string(this->journeyIndex) + "p" + std::to_string(this->phaseIndex) + ".Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }//end switch/case on boundary types
                }//end ephemeris pegged departures
                else if (this->myJourneyOptions->departure_class == BoundaryClass::FreePoint)
                {
                    switch (this->myJourneyOptions->departure_type)
                    {

                        case LAUNCH_OR_DIRECT_INSERTION:
                            this->myDepartureEvent = new BoundaryEvents::FreePointDirectInsertion(name + "FreePointDirectInsertion",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions,
                                PreviousPhaseArrivalEvent);
                            break;
                        case FREE_DIRECT_DEPARTURE:
                            this->myDepartureEvent = new BoundaryEvents::FreePointFreeDirectDeparture(name + "FreePointFreeDirectDeparture",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions,
                                PreviousPhaseArrivalEvent);
                            break;
                        default:
                            throw std::invalid_argument("Undefined free point departure type '" + DepartureTypeStrings[this->myJourneyOptions->departure_type] + "', j" + std::to_string(this->journeyIndex) + "p" + std::to_string(this->phaseIndex) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }
                }
                else if (this->myJourneyOptions->departure_class == BoundaryClass::EphemerisReferenced)
                {
                    if (this->myJourneyOptions->destination_list[0] == -1) //point on the edge of the universe
                    {
                        switch (this->myJourneyOptions->departure_type)
                        {
                        case FREE_DIRECT_DEPARTURE:
                            this->myDepartureEvent = new BoundaryEvents::EphemerisReferencedFreeDirectDepartureInterior(name + "EphemerisPeggedFreeDirectDepartureInterior",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions,
                                PreviousPhaseArrivalEvent);
                            break;
                        default:
                            throw std::invalid_argument("Undefined ephemeris-referenced departure type '" + DepartureTypeStrings[this->myJourneyOptions->departure_type] + "', j" + std::to_string(this->journeyIndex) + "p" + std::to_string(this->phaseIndex) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                        }
                    }
                    else //central body or body in the universe
                    {
                        switch (this->myJourneyOptions->departure_type)
                        {
                        case FREE_DIRECT_DEPARTURE:
                            this->myDepartureEvent = new BoundaryEvents::EphemerisReferencedFreeDirectDepartureExterior(name + "EphemerisReferencedFreeDirectDepartureExterior",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions,
                                PreviousPhaseArrivalEvent);
                            break;
                        default:
                            throw std::invalid_argument("Undefined ephemeris-referenced departure type '" + DepartureTypeStrings[this->myJourneyOptions->departure_type] + "', j" + std::to_string(this->journeyIndex) + "p" + std::to_string(this->phaseIndex) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                        }
                    }
                }//end ephemeris_referenced
                else if (this->myJourneyOptions->departure_class == BoundaryClass::Periapse)
                {
                switch (this->myJourneyOptions->departure_type)
                {
                    case LAUNCH_OR_DIRECT_INSERTION:
                        this->myDepartureEvent = new BoundaryEvents::PeriapseLaunch(name + "PeriapseLaunch",
                            this->journeyIndex,
                            this->phaseIndex,
                            stageIndex,
                            this->myUniverse,
                            this->mySpacecraft,
                            this->myLaunchVehicle,
                            this->myOptions,
                            PreviousPhaseArrivalEvent);
                        break;
                    default:
                        std::invalid_argument("Undefined periapse departure type '" + DepartureTypeStrings[this->myJourneyOptions->departure_type] + "', j" + std::to_string(this->journeyIndex) + "p" + std::to_string(this->phaseIndex));
                }//end switch
                }//end periapse arrivals
            }//end first phase in journey departures
            else
            {
                if (this->myJourneyOptions->enable_periapse_burns)
                {
                    this->myDepartureEvent = new BoundaryEvents::EphemerisPeggedPoweredFlyby(name + "EphemerisPeggedPoweredFlyby",
                        this->journeyIndex,
                        this->phaseIndex,
                        stageIndex,
                        this->myUniverse,
                        this->mySpacecraft,
                        this->myOptions,
						dynamic_cast<BoundaryEvents::EphemerisPeggedArrivalWithVinfinity*>(PreviousPhaseArrivalEvent));
                }
                else
                {
                    this->myDepartureEvent = new BoundaryEvents::EphemerisPeggedUnpoweredFlyby(name + "EphemerisPeggedUnpoweredFlyby",
                        this->journeyIndex,
                        this->phaseIndex,
                        stageIndex,
                        this->myUniverse,
                        this->mySpacecraft,
                        this->myOptions,
						dynamic_cast<BoundaryEvents::EphemerisPeggedArrivalWithVinfinity*>(PreviousPhaseArrivalEvent));
                }
            }

            //staging
            if (this->isFirstPhaseInJourney
                && this->myJourneyOptions->stage_after_departure)
                ++stageIndex;

            this->stageIndex = stageIndex;

            if (this->stageIndex >= this->mySpacecraft->getNumberOfStages())
            {
                throw std::invalid_argument(this->name + " has invalid stage index " + std::to_string(this->stageIndex) + ". The spacecraft has run out of stages. Check your journey staging options. If you want to debug, place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            if (this->isLastPhaseInJourney
                && this->myJourneyOptions->stage_before_arrival)
                ++stageIndex;

            //I don't need another stageIndex check here because there is one in ArrivalEvent

            //create arrival event
            if (!this->isLastPhaseInJourney)
            {
                this->myArrivalEvent = new BoundaryEvents::EphemerisPeggedFlybyIn(name + "EphemerisPeggedFlybyIn",
                    this->journeyIndex,
                    this->phaseIndex,
                    stageIndex,
                    this->myUniverse,
                    this->mySpacecraft,
                    this->myOptions);
            }
            else
            {
                if (this->myJourneyOptions->arrival_class == BoundaryClass::FreePoint)
                {
                    switch (this->myJourneyOptions->arrival_type)
                    {
                        case CHEM_RENDEZVOUS:
                            this->myArrivalEvent = new BoundaryEvents::FreePointChemRendezvous(name + "FreePointChemRendezvous",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions);
                            break;
                        case INTERCEPT:
                            this->myArrivalEvent = new BoundaryEvents::FreePointIntercept(name + "FreePointIntercept",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions);
                            break;
                        case LT_RENDEZVOUS:
                            this->myArrivalEvent = new BoundaryEvents::FreePointLTRendezvous(name + "FreePointLTRendezvous",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions);
                            break;
                        default:
                            throw std::invalid_argument("Undefined free point arrival type '" + ArrivalTypeStrings[this->myJourneyOptions->arrival_type] + "', j" + std::to_string(this->journeyIndex) + "p" + std::to_string(this->phaseIndex) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }
                }//end free point arrivals
                else if (this->myJourneyOptions->arrival_class == BoundaryClass::EphemerisPegged)
                {
                    if (this->myJourneyOptions->sequence[phaseIndex + 1] < 1)
                    {
                        throw std::invalid_argument("j" + std::to_string(this->journeyIndex) + "p" + std::to_string(this->phaseIndex) + ": you have selected an ephemeris pegged arrival event that occurs at either the central body or the universe's sphere of influence. Neither of these is allowed. Ephemeris pegged boundary events must occur at a body in the universe.");
                    }

                    switch (this->myJourneyOptions->arrival_type)
                    {
                        case INSERTION_INTO_PARKING:
                            this->myArrivalEvent = new BoundaryEvents::EphemerisPeggedOrbitInsertion(name + "EphemerisPeggedOrbitInsertion",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions);
                            break;
                        case CHEM_RENDEZVOUS:
                            this->myArrivalEvent = new BoundaryEvents::EphemerisPeggedChemRendezvous(name + "EphemerisPeggedChemRendezvous",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions);
                            break;
                        case INTERCEPT:
                            this->myArrivalEvent = new BoundaryEvents::EphemerisPeggedIntercept(name + "EphemerisPeggedIntercept",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions);
                            break;
						case MOMENTUM_EXCHANGE:
							this->myArrivalEvent = new BoundaryEvents::EphemerisPeggedMomentumTransfer(name + "EphemerisPeggedMomentumTransfer",
								this->journeyIndex,
								this->phaseIndex,
								stageIndex,
								this->myUniverse,
								this->mySpacecraft,
								this->myOptions);
							break;
                        case LT_RENDEZVOUS:
                            this->myArrivalEvent = new BoundaryEvents::EphemerisPeggedLTRendezvous(name + "EphemerisPeggedLTRendezvous",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions);
                            break;
                        case CHEM_MATCH_VINF:
                            throw std::invalid_argument("EphemerisPeggedChemMatchVinf not implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                            break;
                        case LT_MATCH_VINF:
                            throw std::invalid_argument("EphemerisPeggedLTMatchVinf not implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                            break;
                        case SPIRAL_CAPTURE:
                        {
                            this->myArrivalEvent = new BoundaryEvents::EphemerisPeggedSpiralArrival(name + "EphemerisPeggedSpiralArrival",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions);

                            break;
                        }
                        default:
                            throw std::invalid_argument("Undefined ephemeris pegged arrival type '" + ArrivalTypeStrings[this->myJourneyOptions->arrival_type] + "', j" + std::to_string(this->journeyIndex) + "p" + std::to_string(this->phaseIndex) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }
                }//end ephemeris pegged arrivals
                else if (this->myJourneyOptions->arrival_class == BoundaryClass::EphemerisReferenced)
                {
                    if (this->myJourneyOptions->destination_list[1] == -1) //this is an INTERIOR boundary, i.e. we are going to a point on the edge of the universe
                    { 
                        switch (this->myJourneyOptions->arrival_type)
                        {
                            case CHEM_RENDEZVOUS:
                                throw std::invalid_argument("EphemerisReferencedChemRendezvousInterior not implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                                break;
                                /*
                                case CHEM_RENDEZVOUS:
                                this->myArrivalEvent = new BoundaryEvents::EphemerisReferencedChemRendezvousInterior(name + "EphemerisReferencedChemRendezvousInterior",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions);
                                break;*/
                            case INTERCEPT:
                                this->myArrivalEvent = new BoundaryEvents::EphemerisReferencedInterceptInterior(name + "EphemerisReferencedInterceptInterior",
                                    this->journeyIndex,
                                    this->phaseIndex,
                                    stageIndex,
                                    this->myUniverse,
                                    this->mySpacecraft,
                                    this->myOptions);
                                break;
                            case LT_RENDEZVOUS:
                                this->myArrivalEvent = new BoundaryEvents::EphemerisReferencedLTRendezvousInterior(name + "EphemerisReferencedLTRendezvousInterior",
                                    this->journeyIndex,
                                    this->phaseIndex,
                                    stageIndex,
                                    this->myUniverse,
                                    this->mySpacecraft,
                                    this->myOptions);
                                break;
                            case CHEM_MATCH_VINF:
                                throw std::invalid_argument("EphemerisReferencedChemMatchVinfInterior not implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                                break;
                            case LT_MATCH_VINF:
                                throw std::invalid_argument("EphemerisReferencedLTMatchVinfInterior not implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                                break;
                            default:
                                throw std::invalid_argument("Undefined ephemeris referenced interior arrival type '" + ArrivalTypeStrings[this->myJourneyOptions->arrival_type] + "', j" + std::to_string(this->journeyIndex) + "p" + std::to_string(this->phaseIndex) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                        }//end switch
                    }
                    else //exterior arrival
                    {
                        switch (this->myJourneyOptions->arrival_type)
                        {
                            case CHEM_RENDEZVOUS:
                                throw std::invalid_argument("EphemerisReferencedChemRendezvousExterior not implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                                break;
                            /*
                            case CHEM_RENDEZVOUS:
                                this->myArrivalEvent = new BoundaryEvents::EphemerisReferencedChemRendezvousExterior(name + "EphemerisReferencedChemRendezvousExterior",
                                    this->journeyIndex,
                                    this->phaseIndex,
                                    stageIndex,
                                    this->myUniverse,
                                    this->mySpacecraft,
                                    this->myOptions);
                                break;*/
                            case INTERCEPT:
                                this->myArrivalEvent = new BoundaryEvents::EphemerisReferencedInterceptExterior(name + "EphemerisReferencedInterceptExterior",
                                    this->journeyIndex,
                                    this->phaseIndex,
                                    stageIndex,
                                    this->myUniverse,
                                    this->mySpacecraft,
                                    this->myOptions);
                                break;
                            case LT_RENDEZVOUS:
                                this->myArrivalEvent = new BoundaryEvents::EphemerisReferencedLTRendezvousExterior(name + "EphemerisReferencedLTRendezvousExterior",
                                    this->journeyIndex,
                                    this->phaseIndex,
                                    stageIndex,
                                    this->myUniverse,
                                    this->mySpacecraft,
                                    this->myOptions);
                                break;
                            case CHEM_MATCH_VINF:
                                throw std::invalid_argument("EphemerisReferencedChemMatchVinfExterior not implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                                break;
                            case LT_MATCH_VINF:
                                throw std::invalid_argument("EphemerisReferencedLTMatchVinfExterior not implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                                break;
                            default:
                                throw std::invalid_argument("Undefined ephemeris referenced exterior arrival type '" + ArrivalTypeStrings[this->myJourneyOptions->arrival_type] + "', j" + std::to_string(this->journeyIndex) + "p" + std::to_string(this->phaseIndex) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                        }//end switch
                    }//end exterior ephemeris-referenced arrivals
                }//end ephemeris referenced arrivals
                else if (this->myJourneyOptions->arrival_class == BoundaryClass::Periapse)
                {
                    switch (this->myJourneyOptions->arrival_type)
                    {
                        case INTERCEPT:
                            this->myArrivalEvent = new BoundaryEvents::PeriapseFlybyIn(name + "PeriapseFlybyIn",
                                this->journeyIndex,
                                this->phaseIndex,
                                stageIndex,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions);
                            break;
                        default:
                            throw std::invalid_argument("Undefined periapse arrival type '" + ArrivalTypeStrings[this->myJourneyOptions->arrival_type] + "', j" + std::to_string(this->journeyIndex) + "p" + std::to_string(this->phaseIndex) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }//end switch
                }//end periapse arrivals
                
            }//end arrivals

            //initial TCM stuff
            if (this->isFirstPhaseInJourney && this->myOptions->TCM_post_launch > 0.0)
            {
                this->hasInitialTCM = true;
                this->initial_TCM_magnitude = this->myOptions->TCM_post_launch;
                this->hasMonopropManeuver = true;
            }

            //forced coast stuff

            if (this->isFirstPhaseInJourney && this->myJourneyOptions->forced_initial_coast > math::SMALL)
            {
                this->hasInitialCoast = true;
                this->InitialCoastDuration = this->myJourneyOptions->forced_initial_coast;
            }
            else if (this->isFirstPhaseInMission && this->myOptions->forced_post_launch_coast > math::SMALL)
            {
                this->hasInitialCoast = true;
                this->InitialCoastDuration = this->myOptions->forced_post_launch_coast;
            }
            else if ((!this->isFirstPhaseInJourney
                || this->myJourneyOptions->departure_type == DepartureType::FLYBY 
                || this->myJourneyOptions->departure_type == DepartureType::FLYBY_FIXED_VINF 
                || this->myJourneyOptions->departure_type == DepartureType::ZERO_TURN_FLYBY)
                && this->myOptions->forced_post_flyby_coast > math::SMALL)
            {
                this->hasInitialCoast = true;
                this->InitialCoastDuration = this->myOptions->forced_post_flyby_coast;
            }
            if (!this->isLastPhaseInJourney && this->myOptions->forced_pre_flyby_coast > math::SMALL)
            {
                this->hasTerminalCoast = true;
                this->TerminalCoastDuration = this->myOptions->forced_pre_flyby_coast;
            }
            else if (this->isLastPhaseInJourney
                && (   this->myJourneyOptions->arrival_type == ArrivalType::INTERCEPT
                    || this->myJourneyOptions->arrival_type == ArrivalType::LT_MATCH_VINF
                    || this->myJourneyOptions->arrival_type == ArrivalType::CHEM_RENDEZVOUS)
                && this->myJourneyOptions->forced_terminal_coast > math::SMALL)
            {
                this->hasTerminalCoast = true;
                this->TerminalCoastDuration = this->myJourneyOptions->forced_terminal_coast;
            }

            //electric thruster stuff
            this->PhaseDutyCycle = (this->myJourneyOptions->override_duty_cycle ?
                this->myJourneyOptions->duty_cycle
                : this->myOptions->engine_duty_cycle);
        }//end main constructor

        //destructor
        phase::~phase()
        {
            delete this->myDepartureEvent;
            delete this->myArrivalEvent;
        }//end destructor

        void phase::setup_calcbounds(std::vector<double>* Xupperbounds,
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

            this->myDepartureEvent->setup_calcbounds(Xupperbounds,
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

            this->myArrivalEvent->setup_calcbounds(Xupperbounds,
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

        }//end setup_calcbounds()

        //******************************************calcbounds methods
        void phase::calcbounds_phase_flight_time()
        {
            double SMA1, SMA2, ECC1, ECC2, pseudoa1, pseudoa2, T1, T2;

            if (this->myDepartureEvent->getLeftBoundaryIsABody())
            {
                SMA1 = this->myDepartureEvent->getBody()->SMA;
                ECC1 = this->myDepartureEvent->getBody()->ECC;

				//if the departure body is the central body, use the radius of the body as SMA and set ECC = 0
				if (this->myDepartureEvent->getBody()->spice_ID == this->myUniverse->central_body_SPICE_ID)
				{
					SMA2 = this->myDepartureEvent->getBody()->radius;
					ECC2 = 0.0;
				}
            }
            else if (this->journeyIndex > 0 || this->phaseIndex > 0)//free point, but passed in from the previous phase
            {
                SMA1 = this->myUniverse->r_SOI;
                ECC1 = 1.0e-8;
            }
            else //begin at fixed or free point in space
            {
                math::Matrix<doubleType> temp_coordinates(6, 1);
                math::Matrix<doubleType> temp_elements(6, 1);
                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                {
                    if (this->myJourneyOptions->departure_elements_vary_flag[stateIndex])
                        temp_coordinates(stateIndex) = this->myJourneyOptions->departure_elements_bounds[2 * stateIndex + 1];
                    else
                        temp_coordinates(stateIndex) = this->myJourneyOptions->departure_elements[stateIndex];
                }

                Astrodynamics::inertial2COE(temp_coordinates, this->myUniverse->mu, temp_elements);
                SMA1 = temp_elements(0) _GETVALUE;
                ECC1 = temp_elements(1) _GETVALUE;
            }

            if (this->myArrivalEvent->getLeftBoundaryIsABody())
            {
                SMA2 = this->myArrivalEvent->getBody()->SMA;
                ECC2 = this->myArrivalEvent->getBody()->ECC;

				//if the arrival body is the central body, use the radius of the body as SMA and set ECC = 0
				if (this->myArrivalEvent->getBody()->spice_ID == this->myUniverse->central_body_SPICE_ID)
				{
					SMA2 = this->myArrivalEvent->getBody()->radius;
					ECC2 = 0.0;
				}
			}
            else //end at fixed or free point in space
            {
                math::Matrix<doubleType> temp_coordinates(6, 1);
                math::Matrix<doubleType> temp_elements(6, 1);
                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                {
                    if (this->myJourneyOptions->arrival_elements_vary_flag[stateIndex])
                        temp_coordinates(stateIndex) = this->myJourneyOptions->arrival_elements_bounds[2 * stateIndex + 1];
                    else
                        temp_coordinates(stateIndex) = this->myJourneyOptions->arrival_elements[stateIndex];
                }

                Astrodynamics::inertial2COE(temp_coordinates, this->myUniverse->mu, temp_elements);
                SMA2 = temp_elements(0) _GETVALUE;
                ECC2 = temp_elements(1) _GETVALUE;
            }

            if (ECC1 < 1.0)
                pseudoa1 = SMA1 * (1 + ECC1);
            else
                pseudoa1 = this->myUniverse->r_SOI / 5.0;

            if (ECC2 < 1.0)
                pseudoa2 = SMA2 * (1 + ECC2);
            else
                pseudoa2 = this->myUniverse->r_SOI / 5.0;

            T1 = 2 * math::PI*sqrt(pseudoa1*pseudoa1*pseudoa1 / this->myUniverse->mu);// pseudo-period of body 1 in days
            T2 = 2 * math::PI*sqrt(pseudoa2*pseudoa2*pseudoa2 / this->myUniverse->mu);// pseudo-period of body 2 in days

            double forced_coast_this_phase = 1.0; //always have a buffer day
            if (this->hasInitialCoast)
                forced_coast_this_phase += this->InitialCoastDuration;
            if (this->hasTerminalCoast)
                forced_coast_this_phase += this->TerminalCoastDuration;

            if (this->myDepartureEvent->getLeftBoundaryIsABody() && this->myArrivalEvent->getLeftBoundaryIsABody()
                && (this->myDepartureEvent->getBody() == this->myArrivalEvent->getBody()))
            {
                double lowerbound_temp = T1 * 0.5;
                Xlowerbounds->push_back(lowerbound_temp > forced_coast_this_phase ? lowerbound_temp : forced_coast_this_phase);
                Xupperbounds->push_back(T1 * 100.0);
            }
            else if (!this->myDepartureEvent->getLeftBoundaryIsABody() && !this->myArrivalEvent->getLeftBoundaryIsABody()) //for transfers between two free or fixed orbits
            {
                double lowerbound_temp = 1.0;
                Xlowerbounds->push_back(lowerbound_temp > forced_coast_this_phase ? lowerbound_temp : forced_coast_this_phase);
                Xupperbounds->push_back(fmax(T1, T2) * 20.0);
            }
            else if (fabs(pseudoa1 - pseudoa2) < 0.1 * pseudoa1) //if the bodies have very similar semi-major axes
            {
                //lower bound is the same for all non-resonant phases
                double lowerbound_temp = 0.001 * fmin(T1, T2);


                //upper-bound is the larger of the two periods
                lowerbound_temp = lowerbound_temp > forced_coast_this_phase ? lowerbound_temp : forced_coast_this_phase;
                Xlowerbounds->push_back(lowerbound_temp > 10.0 *  this->myUniverse->TU ? 10.0 *  this->myUniverse->TU : lowerbound_temp);
                Xupperbounds->push_back(10.0 * fmax(T1, T2));
            }
            else
            {
                //lower bound is the same for all non-resonant phases
                double lowerbound_temp = 0.01 * fmin(T1, T2);

                lowerbound_temp = lowerbound_temp > forced_coast_this_phase ? lowerbound_temp : forced_coast_this_phase;
                Xlowerbounds->push_back(lowerbound_temp > 10.0 *  this->myUniverse->TU ? 10.0 *  this->myUniverse->TU : lowerbound_temp);

                if (fmax(pseudoa1, pseudoa2) / this->myUniverse->LU < 2.0) //outermost body is an inner body with a < 2 LU
                    Xupperbounds->push_back(10.0 * fmax(T1, T2) < 45.0 * this->myUniverse->TU ? 45.0 * this->myUniverse->TU : 10.0 * fmax(T1, T2));

                else //outermost body is an outer body
                    Xupperbounds->push_back(10.0 * fmax(T1, T2) < 45.0 * this->myUniverse->TU ? 45.0 * this->myUniverse->TU : 10.0 * fmax(T1, T2));
            }

            //make sure the upper bound on the phase flight time does not exceed the upper bound on the mission flight time
            if (this->myOptions->global_timebounded && Xupperbounds->back() > this->myOptions->total_flight_time_bounds[1])
                Xupperbounds->back() = this->myOptions->total_flight_time_bounds[1];

            //make sure the upper bound on the phase flight time does not exceed the upper bound on the journey flight time or journey aggregate flight time
            if ((this->myJourneyOptions->timebounded == 1 || this->myJourneyOptions->timebounded == 3) && Xupperbounds->back() > this->myJourneyOptions->flight_time_bounds[1])
                Xupperbounds->back() = this->myJourneyOptions->flight_time_bounds[1];

            Xdescriptions->push_back(prefix + "phase flight time");
            X_scale_factors->push_back(this->myUniverse->TU);
            this->Xindex_PhaseFlightTime = Xdescriptions->size() - 1;

            //compute the synodic period of the boundary points, for use in the MBH synodic period perturbation
            //these are "true" periods, not the pseudo-periods used for computing the bounds
            T1 = 2 * math::PI*sqrt(fabs(SMA1*SMA1*SMA1) / this->myUniverse->mu);
            T2 = 2 * math::PI*sqrt(fabs(SMA2*SMA2*SMA2) / this->myUniverse->mu);

            if ( fabs(T1 - T2) > 1.0 )//this is in seconds so it's pretty restrictive
                this->synodic_period = 1.0 / (fabs(1.0 / T1 - 1.0 / T2));
            else
                this->synodic_period = T1; //effectively allows you to hop resonances

        }//end calcbounds_phase_flight_time()

        void phase::calcbounds_phase_left_boundary()
        {
            if (this->isFirstPhaseInMission)
            {
                //there are no pre-existing time variables, so pass in an empty vector
                this->myDepartureEvent->calcbounds(std::vector<size_t> {});
            }
            else
            {
                //there are pre-existing time variables in the previous phase's arrival event
                this->myDepartureEvent->calcbounds(this->PreviousPhaseArrivalEvent->get_Xindices_EventRightEpoch());
            }

            this->X_index_of_first_decision_variable_in_this_phase = this->Xdescriptions->size();
            this->F_index_of_first_constraint_in_this_phase = this->Fdescriptions->size();

            //staging constraint, if applicable
            if (this->StageBeforePhase)
            {
                //set active stage
                this->mySpacecraft->setActiveStage(this->stageIndex);

                //create the staging mass constraint
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterDeparture =
                    this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();

                this->Flowerbounds->push_back(this->mySpacecraft->getCurrentDryMass() / this->myJourneyOptions->maximum_mass);
                this->Fupperbounds->push_back(1.0);
                this->Fdescriptions->push_back(prefix + "staging mass constraint");

                //has derivative with respect to current left boundary mass ONLY
                size_t Xindex_left_mass;
                for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterDeparture.size(); ++dIndex)
                {
                    if (std::get<1>(Derivatives_of_StateAfterDeparture[dIndex]) == 6)
                    {
                        Xindex_left_mass = std::get<0>(Derivatives_of_StateAfterDeparture[dIndex]);
                        break;
                    }
                }
                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    Xindex_left_mass,
                    Gindex_StageMass_LeftBoundaryMass);
            }

            //this is a good opportunity to get the X index of the launch epoch
            for (size_t Xindex = 0; Xindex < Xdescriptions->size(); ++Xindex)
            {
                if (Xdescriptions->at(Xindex).find("epoch") < 1024)
                {
                    this->Xindex_LaunchDate = Xindex;
                    break;
                }
            }
        }//end calcbounds_phase_left_boundary()

        void phase::calcbounds_phase_right_boundary()
        {
            std::vector<size_t> Departure_Xindices_Epoch = this->myDepartureEvent->get_Xindices_EventRightEpoch();

            Departure_Xindices_Epoch.push_back(this->Xindex_PhaseFlightTime);

            this->myArrivalEvent->calcbounds(Departure_Xindices_Epoch);
        }//end calcbounds_phase_right_boundary()
        
        //******************************************process methods
        
        //function to process phase flight time
        void phase::process_phase_flight_time(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            this->PhaseFlightTime = X[Xindex++];

            this->LaunchDate = X[this->Xindex_LaunchDate];
        }//end process_phase_flight_time()

        //function to process phase boundaries
        void
            phase::process_phase_left_boundary(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //reset the phase
            this->reset_phase(X, Xindex, F, Findex, G, needG);

            //call the events themselves
            this->myDepartureEvent->reset_ETM();
            this->myDepartureEvent->process_event(X, Xindex, F, Findex, G, needG);

            //extract the states
            math::Matrix<doubleType>& DepartureState = this->myDepartureEvent->get_state_after_event();
            for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
                this->state_at_beginning_of_phase(stateIndex) = DepartureState(stateIndex);
            this->state_after_initial_TCM = this->state_at_beginning_of_phase;

            //process staging if appropriate
            if (this->StageBeforePhase)
            {
                F[Findex++] = this->state_at_beginning_of_phase(6) * this->myUniverse->continuity_constraint_scale_factors(6);

                size_t Gindex = this->Gindex_StageMass_LeftBoundaryMass;
                size_t Xindex = this->jGvar->operator[](Gindex);
                
                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                    * 1.0
                    * this->myUniverse->continuity_constraint_scale_factors(6);
            }

            //process the initial TCM
            if (this->hasInitialTCM)
            {
                this->mySpacecraft->computeChemicalPropulsionPerformance(this->initial_TCM_magnitude, this->state_at_beginning_of_phase(6), true, PropulsionSystemChoice::Monoprop); //monoprop

                this->chemical_fuel_used += this->mySpacecraft->getChemFuelConsumedThisManeuver();

                this->state_after_initial_TCM(6) -= this->chemical_fuel_used;

                this->TCMTM(6, 6) = (this->state_after_initial_TCM(6) / this->state_at_beginning_of_phase(6) )_GETVALUE;
            }
        }//end process_phase_left_boundary()

        void
            phase::process_phase_right_boundary(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //call the events themselves
            this->myArrivalEvent->reset_ETM();
            this->myArrivalEvent->process_event(X, Xindex, F, Findex, G, needG);

            //extract the states
            math::Matrix<doubleType>& ArrivalState = this->myArrivalEvent->get_state_before_event();
            for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
                this->state_at_end_of_phase(stateIndex) = ArrivalState(stateIndex);
        }//end process_phase_right_boundary()

        void phase::reset_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //get the value of the launch date
            this->LaunchDate = X[this->Xindex_LaunchDate];

            //reset staging
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //reset the propellant tanks
            this->electric_propellant_used = 0.0;
            this->chemical_fuel_used = 0.0;
            this->chemical_oxidizer_used = 0.0;
            this->PhaseTotalDeterministicDeltav = 0.0;
            this->PhaseTotalStatisticalDeltav = 0.0;
        }//end reset_phase()

        //***************************************************common output segments
        void phase::output_initial_TCM(std::ofstream& outputfile,
            size_t& eventcount)
        {
            //where are we?
            math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
            if (this->myUniverse->central_body_SPICE_ID == 10)
            {
                R_sc_Sun = this->state_after_initial_TCM.getSubMatrix1D(0, 2);
            }
            else
            {
                //where is the central body relative to the sun?
                doubleType central_body_state_and_derivatives[12];
                this->myUniverse->locate_central_body(this->state_after_initial_TCM(7),
                    central_body_state_and_derivatives,
                    *this->myOptions,
                    false);

                math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                }

                R_sc_Sun = this->state_after_initial_TCM.getSubMatrix1D(0, 2) + R_CB_Sun;
            }

            //what are the available and active power?
            doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
            this->mySpacecraft->computePowerState(r_sc_sun_AU, this->state_after_initial_TCM(7));

            this->output_power = this->mySpacecraft->getAvailablePower();
            this->output_active_power = this->mySpacecraft->getEPActivePower();

            if (this->hasInitialTCM)
            {
                phase::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "TCM",//event_type
                    "deep-space",//event_location
                    0.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    this->state_after_initial_TCM,//state
                    math::Matrix<doubleType>(3, 1, 0.0),//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    this->initial_TCM_magnitude,//dVmag
                    0.0,//Thrust
                    this->mySpacecraft->getMonopropIsp(),//Isp
                    this->output_power,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    this->output_active_power,
                    "none");//active_power)
            }
        }//end output_initial_TCM()
    }//end namespace Phases
}//end namespace EMTG