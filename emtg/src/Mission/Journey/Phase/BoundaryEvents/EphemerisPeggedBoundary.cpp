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

#include "EphemerisPeggedBoundary.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedBoundary::EphemerisPeggedBoundary(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            EphemerisPeggedBoundary()
        {
            this->initialize(name,
                             journeyIndex,
                             phaseIndex,
                             stageIndex,
                             Universe,
                             mySpacecraft,
                             myOptions);
        }

        void EphemerisPeggedBoundary::initialize(const std::string& name,
                                const size_t& journeyIndex,
                                const size_t& phaseIndex,
                                size_t& stageIndex,
                                Astrodynamics::universe* Universe,
                                HardwareModels::Spacecraft* mySpacecraft,
                                missionoptions* myOptions)
        {
            //base class
            this->BoundaryEventBase::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);

            this->LeftBoundaryIsABody = true;
        }//end initialize()

        //calcbounds methods

        void EphemerisPeggedBoundary::calcbounds_event_left_side(const std::vector<double>& MassBounds, const std::vector<double>& EpochBounds, std::vector<size_t> timeVariables)
        {
            //Step 1: mass variable
            this->Xlowerbounds->push_back(MassBounds[0]);
            this->Xupperbounds->push_back(MassBounds[1]);
            this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(6));
            this->Xdescriptions->push_back(prefix + "event left state mass");
            this->Xindex_mass = this->Xdescriptions->size() - 1;
            this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(Xindex_mass, 6, 1.0));
            this->dIndex_mass_wrt_encodedMass = this->Derivatives_of_StateBeforeEvent.size() - 1;


            //Step 2: left epoch
            if (this->isFirstEventInMission)
            {
                this->Xlowerbounds->push_back(EpochBounds[0]);
                this->Xupperbounds->push_back(EpochBounds[1]);
                this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(7));
                this->Xdescriptions->push_back(prefix + "event left state epoch");
                timeVariables.insert(timeVariables.begin(), this->Xdescriptions->size() - 1);
            }

            this->calculate_dependencies_left_epoch(timeVariables);

            //all state variables except mass in an EphemerisPeggedBoundary event have a derivative with respect to epoch
            //we'll put in a dummy derivative of 0.0 for now, and later, when the event is processed, we'll do it right
            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
            {
                for (size_t Xepoch_index = 0; Xepoch_index < this->Xindices_EventLeftEpoch.size(); ++Xepoch_index)
                {
                    this->Derivatives_of_StateBeforeEvent_wrt_Time.push_back(std::make_tuple(this->Xindices_EventLeftEpoch[Xepoch_index], stateIndex, 0.0));
                }
            }
        }//end calcbounds_event_left_side

        void EphemerisPeggedBoundary::calcbounds_event_right_side()
        {
            this->Derivatives_of_StateAfterEvent = this->Derivatives_of_StateBeforeEvent;
            this->Derivatives_of_StateAfterEvent_wrt_Time = this->Derivatives_of_StateBeforeEvent_wrt_Time;

            this->Xindices_EventRightEpoch = this->Xindices_EventLeftEpoch;
        }//end calcbounds_event_right_side

        //process methods
        
        void EphemerisPeggedBoundary::process_event_left_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: extract mass variable
            this->state_before_event(6) = X[Xindex++];
            this->boundary_state(6) = this->state_before_event(6);

            //Step 2: left epoch
            this->process_left_epoch(X, Xindex, F, Findex, G, needG);
            this->boundary_state(7) = this->state_before_event(7);

            //Step 3: ephemeris lookup
            doubleType body_state[12];
            this->myBody->locate_body(this->EventLeftEpoch,
                body_state,
                needG,
                *this->myOptions);

            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
            {
                this->state_before_event(stateIndex) = body_state[stateIndex];
                this->boundary_state(stateIndex) = body_state[stateIndex];
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
        }//end process_event_right_side

        
        void EphemerisPeggedBoundary::process_event_right_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->state_after_event.shallow_copy(this->state_before_event);

            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent.size(); ++dIndex)
                this->Derivatives_of_StateAfterEvent[dIndex] = this->Derivatives_of_StateBeforeEvent[dIndex];

            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
                this->Derivatives_of_StateAfterEvent_wrt_Time[dIndex] = this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex];

            this->EventRightEpoch = this->EventLeftEpoch;
        }//end process_event_right_side
    }//close namespace events
}//close namespace EMTG