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

#include "ArriveAsLateAsPossibleObjective.h"

namespace EMTG
{
    namespace ObjectiveFunctions
{

        ArriveAsLateAsPossibleObjective::ArriveAsLateAsPossibleObjective(Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            Mission* myMission) :
            ObjectiveFunctionBase(Universe,
                mySpacecraft,
                myOptions,
                myMission)
        {
            this->name = "ArriveAsEarlyAsPossibleObjective";

            size_t number_of_phases_in_objective_journey = this->myMission->getJourney(this->myOptions->objective_journey)->getNumberOfPhases();

            this->myArrivalEvent = this->myMission->getJourney(this->myOptions->objective_journey)->getPhase(number_of_phases_in_objective_journey - 1)->getArrivalEvent();

            this->ObjectiveScale = 1.0e-3;
        }//end constructor

        void ArriveAsLateAsPossibleObjective::calcbounds()
        {
            //this objective function has derivatives with respect to all time variables that affect the arrival epoch
            this->Xindices_ArrivalEpoch = this->myArrivalEvent->get_Xindices_EventLeftEpoch();

            for (size_t Xindex : this->Xindices_ArrivalEpoch)
            {
                this->create_sparsity_entry(0,
                    Xindex,
                    this->Gindex_derivatives_of_objective_function);
            }
        }//end calcbounds()

        void ArriveAsLateAsPossibleObjective::process(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->ObjectiveValue = 0.0;

            for (size_t Xindex : this->Xindices_ArrivalEpoch)
            {
                this->ObjectiveValue += -this->ObjectiveScale * X[Xindex] * this->myUniverse->continuity_constraint_scale_factors(7);
            }
            F[0] = this->ObjectiveValue;

            if (needG)
            {
                for (size_t Gindex : this->Gindex_derivatives_of_objective_function)
                {
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = -this->ObjectiveScale * this->X_scale_factors->operator[](Xindex)
                        * 1.0
                        * this->myUniverse->continuity_constraint_scale_factors(7);
                }
            }
        }//end process()

        void ArriveAsLateAsPossibleObjective::output(std::ofstream& outputfile)
        {

            outputfile << "Objective function is \"arrive as late as possible\"" << std::endl;
            outputfile << "J = " << this->ObjectiveValue _GETVALUE << std::endl;
        }
    }//end namespace ObjectiveFunctions
}//end namespace EMTG