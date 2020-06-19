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

//MGAnDSMs maneuver cone angle exclusion constraint
//3-27-2018

#pragma once

#include "doubleType.h"
#include "MGAnDSMs_maneuver_constraint.h"

#include "MGAnDSMs_subphase.h"

namespace EMTG
{
    namespace Phases
    {
        //forward declare parent class
        class MGAnDSMs_subphase;

        class MGAnDSMs_maneuver_cone_angle_exclusion_constraint : public MGAnDSMs_maneuver_constraint
        {
        public:
            //constructor
            MGAnDSMs_maneuver_cone_angle_exclusion_constraint() {};
            MGAnDSMs_maneuver_cone_angle_exclusion_constraint(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& subphaseIndex,
                const size_t& stageIndex,
                MGAnDSMs_subphase* mySubPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                const std::string& ConstraintDefinition);

            //clone
            virtual MGAnDSMs_maneuver_cone_angle_exclusion_constraint* clone() const { return new MGAnDSMs_maneuver_cone_angle_exclusion_constraint(*this); }

            virtual ~MGAnDSMs_maneuver_cone_angle_exclusion_constraint() {};

            //calcbounds goes in the specialized phase
            virtual void calcbounds();

            //process goes in the specialized phase
            virtual void process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

        protected:
            std::vector<size_t> Gindex_wrt_ManeuverComponents;
            std::vector<size_t> Gindex_wrt_time;
            size_t reference_body_index; //-2 for central body
            double minimum_cone_angle;
        };
    }//close namespace Phases
}//close namespace EMTG