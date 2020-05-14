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

#include "PhaseDistanceConstraintFactory.h"

#include <vector>

#include "boost/algorithm/string.hpp"

namespace EMTG
{
    namespace Phases
    {
        std::tuple<int, double, double> CreatePhaseDistanceConstraint(const std::string& ConstraintDefinition,
            const missionoptions* myOptions,
            const Astrodynamics::universe* myUniverse)
        {

            std::vector<std::string> ConstraintDefinitionCell;
            boost::split(ConstraintDefinitionCell, ConstraintDefinition, boost::is_any_of("_"), boost::token_compress_on);

            int bodyIndex;
            if (boost::to_lower_copy(ConstraintDefinitionCell[1]) == "cb")
                bodyIndex = -2;
            else
                bodyIndex = std::stoi(ConstraintDefinitionCell[1]);

            double lowerBound, upperBound;

            if (ConstraintDefinitionCell[2].find("km") < 1024)
            {
                lowerBound = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[2], "km"));
            }
            else if (ConstraintDefinitionCell[2].find("au") < 1024)
            {
                lowerBound = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[2], "au")) * myOptions->AU;
            }
            else if (ConstraintDefinitionCell[2].find("lu") < 1024)
            {
                lowerBound = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[2], "lu")) * myUniverse->LU;
            }
            else
            {
                lowerBound = std::stod(ConstraintDefinitionCell[2]);
            }

            if (ConstraintDefinitionCell[3].find("km") < 1024)
            {
                upperBound = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[3], "km"));
            }
            else if (ConstraintDefinitionCell[3].find("au") < 1024)
            {
                upperBound = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[3], "au")) * myOptions->AU;
            }
            else if (ConstraintDefinitionCell[3].find("lu") < 1024)
            {
                upperBound = std::stod(boost::erase_all_copy(ConstraintDefinitionCell[3], "lu")) * myUniverse->LU;
            }
            else
            {
                upperBound = std::stod(ConstraintDefinitionCell[3]);
            }

            return std::tuple<int, double, double>({ bodyIndex, lowerBound, upperBound });
        }
    }//end namespace Phases
}//end namespace EMTG