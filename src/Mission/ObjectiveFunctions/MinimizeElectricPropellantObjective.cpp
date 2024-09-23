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

#include "MinimizeElectricPropellantObjective.h"

namespace EMTG
{
    namespace ObjectiveFunctions
{

        MinimizeElectricPropellantobjective::MinimizeElectricPropellantobjective(Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            Mission* myMission) :
            ObjectiveFunctionBase(Universe,
                mySpacecraft,
                myOptions,
                myMission)
        {
            this->name = "MinimizeElectricPropellantObjectiveFunction";

            size_t number_of_phases_in_final_journey = this->myMission->getJourney(this->myOptions->number_of_journeys - 1)->getNumberOfPhases();

            this->myArrivalEvent = this->myMission->getJourney(this->myOptions->number_of_journeys - 1)->getPhase(number_of_phases_in_final_journey - 1)->getArrivalEvent();

            this->ObjectiveScale = 1.0;
                
        }//end constructor

        void MinimizeElectricPropellantobjective::calcbounds()
        {
            //electric propellant
            std::vector<size_t> Xindices_tank = this->mySpacecraft->getGlobalElectricPropellantTank_Xindices();
            this->nElectricPropellantTanks = Xindices_tank.size();

            for (size_t virtualTankIndex = 0; virtualTankIndex < this->nElectricPropellantTanks; ++virtualTankIndex)
            {
                size_t Xindex = Xindices_tank[virtualTankIndex];

                this->create_sparsity_entry(0,
                    Xindex,
                    this->Gindex_dVirtualElectricPropellant);
            }
        }//end calcbounds()

        void MinimizeElectricPropellantobjective::process(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->mySpacecraft->computePropellantState(X, this->myOptions->electric_propellant_margin, this->myOptions->chemical_propellant_margin);


            this->ObjectiveValue = this->ObjectiveScale
                                    * (this->mySpacecraft->getGlobalElectricPropellantUsed() + this->mySpacecraft->getElectricPropellantMargin())
                                    * this->myUniverse->continuity_constraint_scale_factors(6);

            F[0] = this->ObjectiveValue;

            if (needG)
            {
                //chemical fuel - opposite sign from mass derivatives
                for (size_t dTankIndex = 0; dTankIndex < this->nElectricPropellantTanks; ++dTankIndex)
                {
                    size_t Gindex = this->Gindex_dVirtualElectricPropellant[dTankIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->ObjectiveScale * this->X_scale_factors->operator[](Xindex)
                        * (1.0 + this->myOptions->electric_propellant_margin)
                        * this->myUniverse->continuity_constraint_scale_factors(6);
                }
            }
        }//end process()

        void MinimizeElectricPropellantobjective::output(std::ofstream& outputfile)
        {
            outputfile << "Objective function is \"minimize electric propellant\"" << std::endl;
            outputfile << "J = " << this->ObjectiveValue _GETVALUE << std::endl;
        }
    }//end namespace ObjectiveFunctions
}//end namespace EMTG