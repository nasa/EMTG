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

#include "FreePointArrival.h"
#include "StateRepresentationFactory.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        FreePointArrival::FreePointArrival(const std::string& name,
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

        void FreePointArrival::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            //we have to do this up front, before delegating
            this->AllowStateToPropagate = myOptions->Journeys[journeyIndex].AllowJourneyFreePointArrivalToPropagate;

            this->FreePointBoundary::initialize(name,
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

            this->myStateRepresentationEnum = this->myJourneyOptions->arrival_elements_state_representation;

            //create the state representation
            this->myStateRepresentation = Astrodynamics::CreateStateRepresentation(this->myStateRepresentationEnum, this->myUniverse->mu);

            this->myEncodedReferenceFrame = this->myJourneyOptions->arrival_elements_frame;

            this->ReferenceEpoch = this->myJourneyOptions->arrival_elements_reference_epoch;

            //are we using an object-referenced frame? If so, let's make a body
            if (this->myEncodedReferenceFrame == ReferenceFrame::ObjectReferenced)
            {
                if (this->myOptions->Journeys[journeyIndex].destination_list[1] < 1)
                {
                    throw std::invalid_argument(this->name + " reference body must be set to a body in the universe. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
                else
                    this->myBody = &myUniverse->bodies[this->myOptions->Journeys[journeyIndex].destination_list[1] - 1];
            }
        }//end initialize()
        
        //******************************************calcbounds methods
        void FreePointArrival::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            this->X_index_of_first_decision_variable_in_this_event = this->Xdescriptions->size();

            std::vector< std::tuple<double, double> > StateBounds;
            //Step 1: state bounds
            //Step 1.1: position and velocity
            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
            {
                if (this->myJourneyOptions->arrival_elements_vary_flag[stateIndex])
                    StateBounds.push_back({ this->myJourneyOptions->arrival_elements_bounds[2 * stateIndex], this->myJourneyOptions->arrival_elements_bounds[2 * stateIndex + 1] });
                else
                    StateBounds.push_back({ this->myJourneyOptions->arrival_elements[stateIndex] - 1.0e-13, this->myJourneyOptions->arrival_elements[stateIndex] + 1.0e-13 });
            }

            //Step 1.2: mass
            StateBounds.push_back({ 1.0e-13, this->myJourneyOptions->maximum_mass });

            //Step 2: base free point boundary
            FreePointBoundary::calcbounds_event_left_side(StateBounds, timeVariables);
        }//end calcbounds_event_left_side()

        void FreePointArrival::calcbounds_event_right_side()
        {
            //base classes
            FreePointBoundary::calcbounds_event_right_side();

            ArrivalEvent::calcbounds_event_right_side();
        }//end calcbounds_event_right_side

        void FreePointArrival::calcbounds_specialized_constraints()
        {
            this->ArrivalEvent::calcbounds_specialized_constraints();
        }

        //**************************************process functions
        void FreePointArrival::
            process_event_left_side(const std::vector<doubleType>& X,
                                    size_t& Xindex,
                                    std::vector<doubleType>& F,
                                    size_t& Findex,
                                    std::vector<double>& G,
                                    const bool& needG)
        {
            //base ephemeris pegged boundary
            FreePointBoundary::process_event_left_side(X, Xindex, F, Findex, G, needG);

            //save the raw ephemeris state, i.e. before any maneuvering happens
            //we are cheating here, taking advantage of the fact that we "know" ephemeris pegged arrivals have no time width
            this->state_after_event_raw = this->state_before_event;
        }//end process_event_left_side()

        

        void FreePointArrival::
            process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //Step 1: base class
            this->FreePointBoundary::process_event_right_side(X, Xindex, F, Findex, G, needG);

            //note, all of the time derivatives need to have a sign flip since this is an arrival event
            /*
			if (this->AllowStateToPropagate && this->myPropagatorType == IntegratedPropagator)
            {
                for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                {
                    for (size_t Xepoch_index = 0; Xepoch_index < this->Xindices_EventLeftEpoch.size(); ++Xepoch_index)
                    {
                        std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[this->dIndex_StateBeforeEvent_wrt_Time[stateIndex][Xepoch_index]]) *= -1.0;
                    }
                }
            }
			*/

            //Step 2: decrement propellant
            this->state_after_event(6) -= this->chemical_fuel_used + this->chemical_oxidizer_used + this->electric_propellant_used;

            //Step 3: perform arrival mass increment
            this->process_post_arrival_mass_increment();

            //Step 4: perform post-arrival delta-v
            this->process_post_arrival_deltav();
            
            //Step 5: adjust the mass derivative as necessary
            std::get<2>(this->Derivatives_of_StateAfterEvent[this->dIndex_mass_wrt_encodedMass]) *= this->ETM(6, 6);
        }//end process_event_right_side()

        void FreePointArrival::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            //The base-class version of this method assumes that there is NOT an arrival maneuver,
            //and therefore does not write any lines in the maneuver spec
            if (haveManeuverNeedTarget && this->myJourneyOptions->FreePointArrival_print_target_spec)
            {

                //reset the flag that tells us we need a target spec line
                haveManeuverNeedTarget = false;

                //Step 1: initialize a target spec object
                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->state_after_event(7),
                    this->state_after_event);

                //Step 2: write target spec object
                myTargetSpecLine.write(target_spec_file);
            }
        }//end output_maneuver_and_target_spec()
    }//end namespace BoundaryEvents
}//end namespace EMTG