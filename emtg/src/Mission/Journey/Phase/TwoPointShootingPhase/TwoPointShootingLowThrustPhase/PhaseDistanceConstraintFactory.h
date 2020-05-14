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

//phase distance constraint factory
//for TwoPointShootingLowThrustPhase
//Jacob Englander 1/18/2018

#include <tuple>
#include <string>

#include "missionoptions.h"
#include "universe.h"

namespace EMTG
{
    namespace Phases
    {
        std::tuple<int, double, double> CreatePhaseDistanceConstraint(const std::string& ConstraintDefinition,
            const missionoptions* myOptions,
            const Astrodynamics::universe* myUniverse);
    }//end namespace Phases
}//end namespace EMTG