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

#include "PeriapseArrival.h"
#include "bplane.h"
#include "StateRepresentationFactory.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        PeriapseArrival::PeriapseArrival(const std::string& name,
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
        }//end constructor

        void PeriapseArrival::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {

            //set this boundary's state representation
            this->myStateRepresentationEnum = myOptions->PeriapseBoundaryStateRepresentation;

            //periapse arrival cannot use IncomingBplane
            if (this->myStateRepresentationEnum == StateRepresentation::OutgoingBplane)
            {
                std::cout << "In Journey " << this->journeyIndex << "'s arrival event, the state representation is set to "
                    << StateRepresentationStrings[this->myStateRepresentationEnum == StateRepresentation::OutgoingBplane]
                    << ". PeriapseArrival automatically switches this to "
                    << StateRepresentationStrings[this->myStateRepresentationEnum == StateRepresentation::IncomingBplane] << std::endl;
                this->myStateRepresentationEnum = StateRepresentation::IncomingBplane;
            }

            this->myStateRepresentation = Astrodynamics::CreateStateRepresentation(this->myStateRepresentationEnum, Universe->mu);

            this->PeriapseBoundary::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);

            this->ArrivalEvent::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);

            //set periapse distance bounds

            if (this->myJourneyOptions->PeriapseArrival_override_altitude)
            {
                this->periapseDistanceBounds[0] = this->myUniverse->central_body.radius + this->myJourneyOptions->PeriapseArrival_altitude_bounds[0];
                this->periapseDistanceBounds[1] = this->myUniverse->central_body.radius + this->myJourneyOptions->PeriapseArrival_altitude_bounds[1];
            }
            else
            {
                this->periapseDistanceBounds[0] = this->myUniverse->central_body.radius + this->myUniverse->central_body.minimum_safe_flyby_altitude;
                if (this->myUniverse->central_body.mass < 1.0e+25)
                    this->periapseDistanceBounds[1] = 10.0 * this->myUniverse->central_body.radius;
                else
                    this->periapseDistanceBounds[1] = 300.0 * this->myUniverse->central_body.radius;
            }
        }//end initialize()
        
        //******************************************calcbounds methods
        void PeriapseArrival::calcbounds_event_left_side(const std::vector<double>& RadiusBounds,
                                                         const std::vector<double>& VelocityMagnitudeBounds,
                                                         std::vector<size_t> timeVariables)
        {            
            this->PeriapseBoundary::calcbounds_event_left_side(RadiusBounds,
                VelocityMagnitudeBounds,
                timeVariables);
        }//end calcbounds_event_left_side()

        void PeriapseArrival::calcbounds_event_right_side()
        {
            //base classes
            this->PeriapseBoundary::calcbounds_event_right_side();

            this->ArrivalEvent::calcbounds_event_right_side();
        }//end calcbounds_event_right_side

        void PeriapseArrival::calcbounds_specialized_constraints()
        {
            this->PeriapseBoundary::calcbounds_specialized_constraints();
        }

        //**************************************process functions
        void PeriapseArrival::process_event_left_side(const std::vector<doubleType>& X,
                                    size_t& Xindex,
                                    std::vector<doubleType>& F,
                                    size_t& Findex,
                                    std::vector<double>& G,
                                    const bool& needG)
        {

            //epoch
            this->BoundaryEventBase::process_left_epoch(X, Xindex, F, Findex, G, needG);

            //base ephemeris pegged boundary
            this->PeriapseBoundary::process_event_left_side(X, Xindex, F, Findex, G, needG);

            //save the raw ephemeris state, i.e. before any maneuvering happens
            //we are cheating here, taking advantage of the fact that we "know" ephemeris pegged arrivals have no time width
            this->state_after_event_raw = this->state_before_event;
        }//end process_event_left_side()
              

        void PeriapseArrival::process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //Step 1: copy the state
            this->state_after_event = this->state_before_event;
            this->state_after_event(6) -= this->chemical_fuel_used + this->chemical_oxidizer_used + this->electric_propellant_used;

            //Step 2: copy the derivative entries, but recognize that we don't want to kill off the derivative entries that are special to StateAfterEvent and not part of StateBefore Event
            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent.size(); ++dIndex)
                this->Derivatives_of_StateAfterEvent[dIndex] = this->Derivatives_of_StateBeforeEvent[dIndex];

            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
                this->Derivatives_of_StateAfterEvent_wrt_Time[dIndex] = this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex];

            this->EventRightEpoch = this->state_after_event(7);

            //Step 3: perform arrival mass increment
            this->process_post_arrival_mass_increment();

            //Step 4: perform post-arrival delta-v
            this->process_post_arrival_deltav();
            
            //Step 5: adjust the mass derivative as necessary
            std::get<2>(this->Derivatives_of_StateAfterEvent[this->dIndex_mass_wrt_encodedMass]) *= this->ETM(6, 6);
        }//end process_event_right_side()

        void PeriapseArrival::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            //The base-class version of this method assumes that there is NOT an arrival maneuver,
            //and therefore does not write any lines in the maneuver spec
            //sometimes we use PeriapseArrival just to place a periapse distance constraint, i.e. about the sun. We don't want to generate maneuvers and targets for that
            //we might do this in other places, too, but for now let's just lock it out if we're at the sun

            if (this->myUniverse->central_body.spice_ID != 10
                && haveManeuverNeedTarget)
            {
                //reset the flag that tells us we need a target spec line
                haveManeuverNeedTarget = false;

                //Step 1: initialize a target spec object
                Astrodynamics::bplane myBplane(this->myUniverse->central_body.mu);
                myBplane.define_bplane(this->state_before_event);

                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->state_after_event(7),
                    this->state_after_event,
                    myBplane.getBdotR(),
                    myBplane.getBdotT());

                //Step 2: write target spec object
                myTargetSpecLine.write(target_spec_file);
            }
        }//end output_maneuver_and_target_spec()
    }//end namespace BoundaryEvents
}//end namespace EMTG