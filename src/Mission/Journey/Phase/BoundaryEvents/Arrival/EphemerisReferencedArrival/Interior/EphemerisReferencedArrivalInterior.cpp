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

#include "EphemerisReferencedArrivalInterior.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisReferencedArrivalInterior::EphemerisReferencedArrivalInterior(const std::string& name,
                                           const size_t& journeyIndex,
                                           const size_t& phaseIndex,
                                           size_t& stageIndex,
                                           Astrodynamics::universe* Universe,
                                           HardwareModels::Spacecraft* mySpacecraft,
                                           missionoptions* myOptions)
        {
            this->EphemerisReferencedArrivalInterior::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);
        }//end constructor

        void EphemerisReferencedArrivalInterior::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            this->EphemerisReferencedArrival::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);

            //body needs to be the central body
            this->LeftBoundaryIsABody = false;

            this->myBody = &myUniverse->central_body;
            this->semi_axis_a = this->myJourneyOptions->arrival_ellipsoid_axes[0];
            this->semi_axis_b = this->myJourneyOptions->arrival_ellipsoid_axes[1];
            this->semi_axis_c = this->myJourneyOptions->arrival_ellipsoid_axes[2];

            if (this->semi_axis_a == 0.0 && this->semi_axis_b == this->semi_axis_a && this->semi_axis_b == this->semi_axis_c)
            {
                this->semi_axis_a = this->myBody->r_SOI;
                this->semi_axis_b = this->myBody->r_SOI;
                this->semi_axis_c = this->myBody->r_SOI;

                std::cout << "User inputed 0.0 for journey arrival ellipsoid axes in journey " + std::to_string(this->journeyIndex) + ". This instructs EMTG to use the body's sphere of influence of " + std::to_string(this->myBody->r_SOI) + " km." << std::endl;
            }
        }//end initialize

        //******************************************calcbounds methods
        void EphemerisReferencedArrivalInterior::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            //Step 1: base EphemerisReferencedArrival::calcbounds_event_left_side()
            this->EphemerisReferencedArrival::calcbounds_event_left_side(timeVariables);
        }//end calcbounds_event_left_side()

        void EphemerisReferencedArrivalInterior::calcbounds_event_right_side()
        {
            //in an EphemerisReferencedArrivalInterior, the right-hand side is identical to the left-hand side, plus the state of the central body
            //one can, after calling this, override the derivatives of mass if one wants to

            this->Derivatives_of_StateAfterEvent = this->Derivatives_of_state_on_interface_cartesian;
            this->Derivatives_of_StateAfterEvent_wrt_Time = this->Derivatives_of_state_on_interface_cartesian_wrt_Time;

            //Step 2: additional derivative with respect to time
            for (size_t varIndex = 0; varIndex < this->Xindices_EventLeftEpoch.size(); ++varIndex)
            {
                size_t Xindex = this->Xindices_EventLeftEpoch[varIndex];

                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                {
                    this->Derivatives_of_StateAfterEvent_wrt_Time.push_back({ Xindex, stateIndex, 1.0 });
                }
            }

            //Step 3: Xindices_EventRightEpoch
            this->Xindices_EventRightEpoch = this->Xindices_EventLeftEpoch;
        }//end calcbounds_event_right_side

        //**************************************process functions
        void EphemerisReferencedArrivalInterior::process_event_left_side(const std::vector<doubleType>& X,
                                                                         size_t& Xindex,
                                                                         std::vector<doubleType>& F,
                                                                         size_t& Findex,
                                                                         std::vector<double>& G,
                                                                         const bool& needG)
        {
            //Step 0: assume that the velocity entries have already been populated by whatever function calls this one

            //Step 1: base EphemerisReferencedArrival
            this->EphemerisReferencedArrival::process_event_left_side(X, Xindex, F, Findex, G, needG);

        }//end process_event_left_side()        

        void EphemerisReferencedArrivalInterior::process_event_right_side(const std::vector<doubleType>& X,
                                                                          size_t& Xindex,
                                                                          std::vector<doubleType>& F,
                                                                          size_t& Findex,
                                                                          std::vector<double>& G,
                                                                          const bool& needG)
        {
            this->state_after_event = this->state_on_interface_cartesian;
            this->state_after_event(6) -= this->chemical_fuel_used + this->chemical_oxidizer_used + this->electric_propellant_used;

            //Step 1: base class
            this->EphemerisReferencedArrival::process_event_right_side(X, Xindex, F, Findex, G, needG);

            //Step 2: compute the state offset vector
            //Step 2.1: central body with respect to the Sun
            doubleType offset_state[12];

            this->myUniverse->locate_central_body(this->state_after_event(7),
                offset_state,
                *this->myOptions,
                needG);

            //Step 2.2: next universe's central body with respect to the sun
            if (!this->isLastEventInMission)
            {
                doubleType next_universe_central_body_state[12];
                
                this->myUniverse->get_nextUniverse()->locate_central_body(this->state_after_event(7),
                    next_universe_central_body_state,
                    *this->myOptions,
                    needG);

                for (size_t stateIndex = 0; stateIndex < 12; ++stateIndex)
                    offset_state[stateIndex] -= next_universe_central_body_state[stateIndex];
            }

            //Step 3: offset the state
            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
            {
                this->state_after_event(stateIndex) += offset_state[stateIndex];
            }

            if (needG)
            {
                for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
                {
                    size_t stateIndex = std::get<1>(this->Derivatives_of_StateAfterEvent_wrt_Time[dIndex]);

                    if (stateIndex < 6)
                        std::get<2>(this->Derivatives_of_StateAfterEvent_wrt_Time[dIndex]) = offset_state[6 + stateIndex]_GETVALUE;
                }
            }
        }//end process_event_right_side()
    }//end namespace BoundaryEvents
}//end namespace EMTG