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

#include "arrival.h"
#include "SpecializedBoundaryConstraintFactory.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        //specialized constructor
        ArrivalEvent::ArrivalEvent(const std::string& name,
                                const size_t& journeyIndex,
                                const size_t& phaseIndex,
                                size_t& stageIndex,
                                Astrodynamics::universe* Universe,
                                HardwareModels::Spacecraft* mySpacecraft,
                                missionoptions* myOptions)
        {
            this->initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);
        }

        void ArrivalEvent::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            //specialized classes turn on TCMs as appropriate - right now that's only the ones that arrive with a nonzero v-infinity
            this->hasTCM = false;

            if (this->phaseIndex == this->myOptions->Journeys[this->journeyIndex].number_of_phases - 1)
            {
                this->isLastEventInJourney = true;

                if (this->journeyIndex == this->myOptions->number_of_journeys - 1)
                    this->isLastEventInMission = true;
                else
                    this->isLastEventInMission = false;
            }
            else
            {
                this->isLastEventInJourney = false;
                this->isLastEventInMission = false;
            }

            //staging
            if (this->isLastEventInJourney
                && this->myOptions->Journeys[this->journeyIndex].stage_after_arrival)
                ++stageIndex;

            if (this->stageIndex >= this->mySpacecraft->getNumberOfStages())
            {
                throw std::invalid_argument(this->name + " has invalid stage index " + std::to_string(this->stageIndex) + ". The spacecraft has run out of stages. Check your journey staging options. If you want to debug, place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            //boundary constraints
            this->construct_journey_boundary_constraints();

            //journey-end maneuver propellant
            this->journey_end_propellant_used = 0.0;

            //final mass increment
            this->final_mass_increment = 0.0;
        }//end initialize()

        void ArrivalEvent::construct_journey_boundary_constraints()
        {
            std::vector<std::string> constraintsToAdd = this->myJourneyOptions->BoundaryConstraintDefinitions;
            this->construct_boundary_constraints(constraintsToAdd);
        }

        void ArrivalEvent::construct_boundary_constraints(std::vector<std::string> givenConstraints)
        {
            //first construct this event's tag
            std::string Tag = "p" + std::to_string(this->phaseIndex)
                + "_arrival";

            std::string TagEnd = "pEnd_arrival";

            //clear the current constraint vector
            this->mySpecializedConstraints.clear();

            //are we adding constraints from journeyOptions or from an input vector?
            std::vector<std::string>* constraintsToAdd = &givenConstraints;

            //now, loop over constraints to see if they are relevant
            for (std::string& constraint : *constraintsToAdd)
            {
                if (constraint.find("#") != 0) //don't create a constraint if it is commented out
                {
                    if (constraint.find(Tag) < 1024 || (constraint.find(TagEnd) < 1024 && this->isLastEventInJourney))
                    {
                        if (constraint.find("monoprop") < 1024)
                        {
                            this->ChemicalManeuverType = PropulsionSystemChoice::Monoprop;
                        }
                        else if (constraint.find("biprop") < 1024)
                        {
                            this->ChemicalManeuverType = PropulsionSystemChoice::Biprop;
                        }
                        else
                        {
                            this->mySpecializedConstraints.push_back(BoundaryEvents::SpecializedConstraints::create_boundary_event_constraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                (BoundaryEventBase*)this,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions,
                                constraint,
                                "arrival"));
                        }
                    }
                }//end constraint comment check
            }//end loop over constraints
        }//end construct_boundary_constraints


        //******************************************calcbounds methods
        void ArrivalEvent::calcbounds_event_right_side()
        {
            if (this->isLastEventInJourney)
            {
                //stage mass
                if (this->myOptions->Journeys[this->journeyIndex].stage_after_arrival || this->isLastEventInMission)
                {
                    size_t Xindex_mass_after_event;
                    for (size_t Xindex = this->X_index_of_first_decision_variable_in_this_event; Xindex < this->Xdescriptions->size(); ++Xindex)
                    {
                        if (this->Xdescriptions->at(Xindex).find("event left state mass") < 1024)
                        {
                            Xindex_mass_after_event = Xindex;
                            break;
                        }
                    }
                    this->mySpacecraft->setXindex_StageFinalMass(this->stageIndex, Xindex_mass_after_event);
                    this->mySpacecraft->setXscale_StageFinalMass(this->stageIndex, this->X_scale_factors->at(Xindex_mass_after_event));
                }
            }
        }//end calcbounds_event_right_side()

         //******************************************process methods
        void ArrivalEvent::process_post_arrival_mass_increment()
        {
            if (this->isLastEventInJourney)
            {
                this->final_mass_increment = this->myJourneyOptions->fixed_ending_mass_increment;
                this->state_after_event(6) += this->final_mass_increment;
            }
        }//end process_post_arrival_mass_increment()

        void ArrivalEvent::process_post_arrival_deltav()
        {
            //note that this delta-v happens AFTER any staging that might occur at the end of the arrival

            //Step 1: perform the delta-v
            this->mySpacecraft->computeChemicalPropulsionPerformance(this->myJourneyOptions->journey_end_deltav,
                this->state_after_event(6),
                false,
                this->myJourneyOptions->journey_end_propulsion_system);

            if (this->myJourneyOptions->journey_end_propulsion_system == PropulsionSystemChoice::Monoprop)
            {
                this->journey_end_propellant_used = this->mySpacecraft->getChemFuelConsumedThisManeuver();
                this->state_after_event(6) -= this->journey_end_propellant_used;

                //I really should make this come out of the monoprop tank, too
                //will need to add to virtual chemical tank?

                this->ETM(6, 6) *= (this->state_after_event(6) / (this->state_after_event(6) + this->journey_end_propellant_used))_GETVALUE;
            }
            else// if (this->myJourneyOptions->journey_end_propulsion_system == PropulsionSystemChoice::Biprop)
            {
                throw std::invalid_argument("At the current time, journey_end_deltav has to be done with a monoprop system. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            //Step 2: derivatives - just need to modify mass w.r.t. mass
        }//end process_post_arrival_deltav()
        
        //******************************************output methods
        void ArrivalEvent::output_mass_increment(std::ofstream& outputfile)
        {
            outputfile << "Journey final mass increment: " << this->myJourneyOptions->fixed_ending_mass_increment << " kg" << std::endl;
        }//end output_mass_increment()

        void ArrivalEvent::output_post_arrival_maneuver(std::ofstream& outputfile)
        {
            std::vector<std::string> ManeuverTypeNames({ "Monoprop", "Biprop", "Electric" });
            outputfile << "Journey post-arrival delta-v (" << ManeuverTypeNames[this->myJourneyOptions->journey_end_propulsion_system] << "): " << this->myJourneyOptions->journey_end_deltav << " km/s" << std::endl;
            outputfile << "Journey post-arrival propellant consumed: " << this->journey_end_propellant_used _GETVALUE << " kg" << std::endl;
            outputfile << std::endl;
        }//end output_mass_increment()
    }//end namespace BoundaryEvents
}//end namespace EMTG