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

#include "MinimizeDeltavObjective.h"

namespace EMTG
{
    namespace ObjectiveFunctions
{

        MinimizeDeltavObjective::MinimizeDeltavObjective(Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            Mission* myMission) :
            ObjectiveFunctionBase(Universe,
                mySpacecraft,
                myOptions,
                myMission)
        {
            this->name = "MinimizeDeltavObjectiveFunction";
        }//end constructor

        void MinimizeDeltavObjective::calcbounds()
        {
            //oddly enough, no bounds - the phases took care of it
        }//end calcbounds()

        void MinimizeDeltavObjective::process(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->ObjectiveValue = 0.0;
            for (size_t journeyIndex = 0; journeyIndex < this->myOptions->number_of_journeys; ++journeyIndex)
                this->ObjectiveValue += this->myMission->getJourney(journeyIndex)->getDeterministicDeltav();

            F[0] = this->ObjectiveValue;
        }//end process()

        void MinimizeDeltavObjective::output(std::ofstream& outputfile)
        {
            outputfile << "Objective function is \"minimize delta-v\"" << std::endl;
            outputfile << "J = " << this->ObjectiveValue _GETVALUE << std::endl;
        }
    }//end namespace ObjectiveFunctions
}//end namespace EMTG