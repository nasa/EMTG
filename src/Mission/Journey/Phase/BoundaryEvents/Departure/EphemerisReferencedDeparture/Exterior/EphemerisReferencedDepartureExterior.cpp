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

#include "EphemerisReferencedDepartureExterior.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisReferencedDepartureExterior::EphemerisReferencedDepartureExterior(const std::string& name,
                                           const size_t& journeyIndex,
                                           const size_t& phaseIndex,
                                           size_t& stageIndex,
                                           Astrodynamics::universe* Universe,
                                           HardwareModels::Spacecraft* mySpacecraft,
                                           missionoptions* myOptions,
                                           ArrivalEvent* PreviousPhaseArrivalEvent)
        {
            this->EphemerisReferencedDepartureExterior::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions,
                PreviousPhaseArrivalEvent);

        }//end constructor


        void EphemerisReferencedDepartureExterior::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            ArrivalEvent* PreviousPhaseArrivalEvent)
        {
            this->EphemerisReferencedDeparture::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions,
                PreviousPhaseArrivalEvent);

            if (this->myJourneyOptions->sequence[phaseIndex] == 0) //launching from the central body
            {
                this->myBody = &myUniverse->central_body;
                this->semi_axis_a = this->myJourneyOptions->departure_ellipsoid_axes[0];
                this->semi_axis_b = this->myJourneyOptions->departure_ellipsoid_axes[1];
                this->semi_axis_c = this->myJourneyOptions->departure_ellipsoid_axes[2];
            }
            else //body in the universe
            {
                this->myBody = &myUniverse->bodies[this->myJourneyOptions->sequence[phaseIndex] - 1];
                this->semi_axis_a = this->myJourneyOptions->departure_ellipsoid_axes[0];
                this->semi_axis_b = this->myJourneyOptions->departure_ellipsoid_axes[1];
                this->semi_axis_c = this->myJourneyOptions->departure_ellipsoid_axes[2];

                if (this->semi_axis_a == 0.0 && this->semi_axis_b == this->semi_axis_a && this->semi_axis_b == this->semi_axis_c)
                {
                    this->semi_axis_a = this->myBody->r_SOI;
                    this->semi_axis_b = this->myBody->r_SOI;
                    this->semi_axis_c = this->myBody->r_SOI;

                    std::cout << "User inputed 0.0 for journey departure ellipsoid axes in journey " + std::to_string(this->journeyIndex) + ". This instructs EMTG to use the body's sphere of influence of " + std::to_string(this->myBody->r_SOI) + " km." << std::endl;
                }
            }

            this->LeftBoundaryIsABody = true;
        }//end initialize()
        
        //******************************************calcbounds methods
        void EphemerisReferencedDepartureExterior::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            //Step 1: base EphemerisReferencedDeparture::calcbounds_event_left_side()
            this->EphemerisReferencedDeparture::calcbounds_event_left_side(timeVariables);

            //Step 2: additional derivative with respect to time
            for (size_t varIndex = 0; varIndex < this->Xindices_EventLeftEpoch.size(); ++varIndex)
            {
                size_t Xindex = this->Xindices_EventLeftEpoch[varIndex];

                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                {
                    //all state variables except mass in an EphemerisPeggedBoundary event have a derivative with respect to epoch
                    //we'll put in a dummy derivative of 0.0 for now, and later, when the event is processed, we'll do it right
                    this->Derivatives_of_StateBeforeEvent_wrt_Time.push_back({ Xindex, stateIndex, 0.0 });
                }
            }
        }//end calcbounds_event_left_side()

        void EphemerisReferencedDepartureExterior::calcbounds_event_right_side()
        {
            //in an EphemerisReferencedDepartureExterior, the right-hand side is identical to the left-hand side, minus the state of the body
            //in other words, it's the interface state directly.
            //one can, after calling this, override the derivatives of mass if one wants to

            this->Derivatives_of_StateAfterEvent = this->Derivatives_of_StateBeforeEvent;
            this->Derivatives_of_StateAfterEvent_wrt_Time = this->Derivatives_of_StateBeforeEvent_wrt_Time;

            //Xindices_EventRightEpoch
            this->Xindices_EventRightEpoch = this->Xindices_EventLeftEpoch;
        }//end calcbounds_event_right_side

        //**************************************process functions
        void EphemerisReferencedDepartureExterior::process_event_left_side(const std::vector<doubleType>& X,
                                                                         size_t& Xindex,
                                                                         std::vector<doubleType>& F,
                                                                         size_t& Findex,
                                                                         std::vector<double>& G,
                                                                         const bool& needG)
        {
            //Step 0: assume that the velocity entries have already been populated by whatever function calls this one

            //Step 1: base EphemerisReferencedDeparture
            this->EphemerisReferencedDeparture::process_event_left_side(X, Xindex, F, Findex, G, needG);

            //Step 2: add the state vector of the body
            doubleType body_state[12];

            this->myBody->locate_body(this->state_before_event(7),
                body_state,
                needG,
                *this->myOptions);

            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
            {
                this->state_before_event(stateIndex) += body_state[stateIndex];
                this->boundary_state(stateIndex) = body_state[stateIndex]; //this isn't strictly right, but it will achieve the objective of getting PyEMTG to plot the body orbit
            }

            if (needG)
            {
                for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
                {
                    size_t stateIndex = std::get<1>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                    if (stateIndex < 6)
                        std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]) = body_state[6 + stateIndex]_GETVALUE;
                }
            }

        }//end process_event_left_side()        

        void EphemerisReferencedDepartureExterior::process_event_right_side(const std::vector<doubleType>& X,
                                                                          size_t& Xindex,
                                                                          std::vector<doubleType>& F,
                                                                          size_t& Findex,
                                                                          std::vector<double>& G,
                                                                          const bool& needG)
        {
            this->state_after_event = this->state_before_event;
            this->state_after_event(6) -= this->chemical_fuel_used + this->chemical_oxidizer_used + this->electric_propellant_used;

            //base class
            this->EphemerisReferencedDeparture::process_event_right_side(X, Xindex, F, Findex, G, needG);
        }//end process_event_right_side()
    }//end namespace BoundaryEvents
}//end namespace EMTG