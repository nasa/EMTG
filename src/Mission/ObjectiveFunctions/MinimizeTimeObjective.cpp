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

#include "MinimizeTimeObjective.h"

namespace EMTG
{
    namespace ObjectiveFunctions
{

        MinimizeTimeObjective::MinimizeTimeObjective(Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            Mission* myMission) :
            ObjectiveFunctionBase(Universe,
                mySpacecraft,
                myOptions,
                myMission)
        {
            this->name = "MaximizeMassObjectiveFunction";
        }//end constructor

        void MinimizeTimeObjective::calcbounds()
        {
            //this objective function has derivatives with respect to all time variables
            std::vector<size_t> timeVariables = this->myMission->getJourney(this->myOptions->number_of_journeys - 1)->getLastPhase()->getArrivalEvent()->get_Xindices_EventRightEpoch();

            for (size_t Xindex : timeVariables)
            {
                if (!(Xdescriptions->operator[](Xindex).find("epoch") < 1024))//we want to exclude launch epoch, which is the only thing called "epoch"
                {
                    this->iAfun->push_back(0);
                    this->jAvar->push_back(Xindex);
                    this->Adescriptions->push_back("Derivative of objective function F[0] with respect to X["
                        + std::to_string(this->jAvar->back()) + "]: " + this->Xdescriptions->operator[](Xindex));
                    this->A->push_back(this->X_scale_factors->operator[](Xindex) / this->myUniverse->LU);
                    this->time_variables_Xindices.push_back(Xindex);
                }
            }
        }//end calcbounds()

        void MinimizeTimeObjective::process(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            doubleType TotalFlightTime = 0.0;
            for (size_t timeIndex = 0; timeIndex < this->time_variables_Xindices.size(); ++timeIndex)
            {
                TotalFlightTime += X[this->time_variables_Xindices[timeIndex]];
            }

            this->ObjectiveValue = TotalFlightTime / this->myUniverse->LU;
            F[0] = ObjectiveValue;
        }//end process()

        void MinimizeTimeObjective::output(std::ofstream& outputfile)
        {
            outputfile << "Objective function is \"minimize flight time\"" << std::endl;
            outputfile << "J = " << this->ObjectiveValue _GETVALUE << std::endl;
        }
    }//end namespace ObjectiveFunctions
}//end namespace EMTG