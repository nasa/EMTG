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

#include "EphemerisReferencedDepartureInterior.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisReferencedDepartureInterior::EphemerisReferencedDepartureInterior(const std::string& name,
                                           const size_t& journeyIndex,
                                           const size_t& phaseIndex,
                                           size_t& stageIndex,
                                           Astrodynamics::universe* Universe,
                                           HardwareModels::Spacecraft* mySpacecraft,
                                           missionoptions* myOptions,
                                           ArrivalEvent* PreviousPhaseArrivalEvent)
        {
            this->EphemerisReferencedDepartureInterior::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions,
                PreviousPhaseArrivalEvent);
        }//end constructor

        void EphemerisReferencedDepartureInterior::initialize(const std::string& name,
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

            //body needs to be the central body
            this->LeftBoundaryIsABody = false;

            this->myBody = &myUniverse->central_body;
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
        }//end initialize

        //******************************************calcbounds methods
        void EphemerisReferencedDepartureInterior::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            //Step 1: base EphemerisReferencedDeparture::calcbounds_event_left_side()
            this->EphemerisReferencedDeparture::calcbounds_event_left_side(timeVariables);
        }//end calcbounds_event_left_side()

        void EphemerisReferencedDepartureInterior::calcbounds_event_right_side()
        {
            //in an EphemerisReferencedDepartureInterior, the right-hand side is identical to the left-hand side, plus the state of the central body
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
        void EphemerisReferencedDepartureInterior::process_event_left_side(const std::vector<doubleType>& X,
                                                                         size_t& Xindex,
                                                                         std::vector<doubleType>& F,
                                                                         size_t& Findex,
                                                                         std::vector<double>& G,
                                                                         const bool& needG)
        {
            //Step 0: assume that the velocity entries have already been populated by whatever function calls this one

            //Step 1: base EphemerisReferencedDeparture
            this->EphemerisReferencedDeparture::process_event_left_side(X, Xindex, F, Findex, G, needG);

        }//end process_event_left_side()        

        void EphemerisReferencedDepartureInterior::process_event_right_side(const std::vector<doubleType>& X,
                                                                          size_t& Xindex,
                                                                          std::vector<doubleType>& F,
                                                                          size_t& Findex,
                                                                          std::vector<double>& G,
                                                                          const bool& needG)
        {
            this->state_after_event = this->state_on_interface_cartesian;
            this->state_after_event(6) -= this->chemical_fuel_used + this->chemical_oxidizer_used + this->electric_propellant_used;

            //base class
            this->EphemerisReferencedDeparture::process_event_right_side(X, Xindex, F, Findex, G, needG);
        }//end process_event_right_side()
    }//end namespace BoundaryEvents
}//end namespace EMTG