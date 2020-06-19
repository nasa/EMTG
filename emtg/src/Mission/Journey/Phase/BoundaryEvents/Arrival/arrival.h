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

#pragma once

#include "BoundaryEventBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class ArrivalEvent : virtual public BoundaryEventBase
        {
        public:
            //default constructor
            ArrivalEvent() : BoundaryEventBase::BoundaryEventBase() {};

            //specialized constructor
            ArrivalEvent(const std::string& name,
                         const size_t& journeyIndex,
                         const size_t& phaseIndex,
                         size_t& stageIndex,
                         Astrodynamics::universe* Universe,
                         HardwareModels::Spacecraft* mySpacecraft,
                         missionoptions* myOptions);

            void initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            void construct_journey_boundary_constraints();
            virtual void construct_boundary_constraints(std::vector<std::string> givenConstraints = {});

            virtual ~ArrivalEvent() {};

            virtual void output_mass_increment(std::ofstream& outputfile);
            virtual void output_post_arrival_maneuver(std::ofstream& outputfile);

            virtual void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget) = 0;

            inline doubleType getFinalMassIncrement() const { return this->final_mass_increment; }

            inline math::Matrix<doubleType>& get_state_after_event_raw() { return this->state_after_event_raw; }

        protected:

            virtual void calcbounds_event_right_side() = 0;

            virtual void process_post_arrival_mass_increment();

            virtual void process_post_arrival_deltav();
            
            //fields
            bool isLastEventInJourney;
            bool isLastEventInMission;

            doubleType journey_end_propellant_used;
            doubleType final_mass_increment;

            math::Matrix<doubleType> state_after_event_raw;
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG