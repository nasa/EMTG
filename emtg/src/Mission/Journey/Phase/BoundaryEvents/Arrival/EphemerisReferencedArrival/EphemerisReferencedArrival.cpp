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

#include "EphemerisReferencedArrival.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisReferencedArrival::EphemerisReferencedArrival(const std::string& name,
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
        
        

        void EphemerisReferencedArrival::initialize(const std::string& name,
                                const size_t& journeyIndex,
                                const size_t& phaseIndex,
                                size_t& stageIndex,
                                Astrodynamics::universe* Universe,
                                HardwareModels::Spacecraft* mySpacecraft,
                                missionoptions* myOptions)
        {
            this->EphemerisReferencedBoundary::initialize(name,
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
        }//end initialize()
        
        //******************************************calcbounds methods
        void EphemerisReferencedArrival::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            //in an EphemerisReferencedArrival, the left side of the event is the interface state. Additional time derivatives may be necessary depending on whether the event is Interior or Exterior
            
            this->Derivatives_of_StateBeforeEvent = this->Derivatives_of_state_on_interface_cartesian;

            this->Derivatives_of_StateBeforeEvent_wrt_Time = this->Derivatives_of_state_on_interface_cartesian_wrt_Time;
        }//end calcbounds_event_left_side()

        void EphemerisReferencedArrival::calcbounds_specialized_constraints()
        {
            this->ArrivalEvent::calcbounds_specialized_constraints();
        }

        //**************************************process functions
        void EphemerisReferencedArrival::process_event_left_side(const std::vector<doubleType>& X,
                                                                 size_t& Xindex,
                                                                 std::vector<doubleType>& F,
                                                                 size_t& Findex,
                                                                 std::vector<double>& G,
                                                                 const bool& needG)
        {
            this->state_before_event = this->state_on_interface_cartesian;

            //copy the derivative entries, but recognize that we don't want to kill off the derivative entries that are special to StateBeforeEvent and not part of state_on_interface_cartesian
            for (size_t dIndex = 0; dIndex < this->Derivatives_of_state_on_interface_cartesian.size(); ++dIndex)
                this->Derivatives_of_StateBeforeEvent[dIndex] = this->Derivatives_of_state_on_interface_cartesian[dIndex];

            for (size_t dIndex = 0; dIndex < this->Derivatives_of_state_on_interface_cartesian_wrt_Time.size(); ++dIndex)
                this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex] = this->Derivatives_of_state_on_interface_cartesian_wrt_Time[dIndex];

            //save the raw ephemeris state, i.e. before any maneuvering happens
            //we are cheating here, taking advantage of the fact that we "know" ephemeris pegged arrivals have no time width
            this->state_after_event_raw = this->state_before_event;
        }//end process_event_left_side()

        void EphemerisReferencedArrival::process_event_right_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: copy the derivative entries, but recognize that we don't want to kill off the derivative entries that are special to StateAfterEvent and not part of state_on_interface_cartesian
            for (size_t dIndex = 0; dIndex < this->Derivatives_of_state_on_interface_cartesian.size(); ++dIndex)
                this->Derivatives_of_StateAfterEvent[dIndex] = this->Derivatives_of_state_on_interface_cartesian[dIndex];

            for (size_t dIndex = 0; dIndex < this->Derivatives_of_state_on_interface_cartesian_wrt_Time.size(); ++dIndex)
                this->Derivatives_of_StateAfterEvent_wrt_Time[dIndex] = this->Derivatives_of_state_on_interface_cartesian_wrt_Time[dIndex];

            this->EventRightEpoch = this->state_after_event(7);

            //Step 2: perform arrival mass increment
            this->process_post_arrival_mass_increment();

            //Step 3: perform post-arrival delta-v
            this->process_post_arrival_deltav();

            //Step 4: adjust the mass derivative as necessary
            std::get<2>(this->Derivatives_of_StateAfterEvent[this->dIndex_mass_wrt_encodedMass]) *= this->ETM(6, 6);
        }//end process_event_right_side()
    }//end namespace BoundaryEvents
}//end namespace EMTG