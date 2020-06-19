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

//MGAnDSMs maneuver epoch constraint
//9-22-2017

#pragma once

#include "doubleType.h"

#include "MGAnDSMs_maneuver_constraint.h"

namespace EMTG
{
    namespace Phases
    {
        //forward declare parent class
        class MGAnDSMs_subphase;

        class MGAnDSMs_maneuver_epoch_constraint : public MGAnDSMs_maneuver_constraint
        {
        public:
            //constructor
            MGAnDSMs_maneuver_epoch_constraint() {};
            MGAnDSMs_maneuver_epoch_constraint(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& subphaseIndex,
                const size_t& stageIndex,
                MGAnDSMs_subphase* mySubPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                const std::string& ConstraintDefinition);

            virtual ~MGAnDSMs_maneuver_epoch_constraint() {};

            std::vector<size_t> Gindex_wrt_BurnIndices;
            std::vector<size_t> Gindex_wrt_TimeVariables;
            size_t Xindex_PhaseFlightTime;
        };
    }//close namespace Phases
}//close namespace EMTG