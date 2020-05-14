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

//boundary constraint factory for ephemeris pegged boundary events
//11-1-2017

#pragma once

#include "doubleType.h"

#include "SpecializedBoundaryConstraintBase.h"

#include "BoundaryEventBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class BoundaryEventBase;

        namespace SpecializedConstraints
        {
            //forward declare stuff
            class SpecializedBoundaryConstraintBase;

            SpecializedBoundaryConstraintBase* create_boundary_event_constraint(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                BoundaryEventBase* myBoundaryEvent,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                const std::string& ConstraintDefinition,
                const std::string& EventDefinition);
        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG