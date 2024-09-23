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

#include "MinimizeInitialImpulseObjective.h"

namespace EMTG
{
    namespace ObjectiveFunctions
{

        MinimizeInitialImpulseObjective::MinimizeInitialImpulseObjective(Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            Mission* myMission) :
            ObjectiveFunctionBase(Universe,
                mySpacecraft,
                myOptions,
                myMission)
        {
            this->name = "MinimizeInitialImpulseObjectiveFunction";

            this->myDepartureEvent = this->myMission->getJourney(0)->getPhase(0)->getDepartureEvent();

            this->ObjectiveScale = 1.0;

            //does the first departure event have an impulse? If not, we need to throw an error
            if (!this->myDepartureEvent->get_hasManeuver())
            {
                throw std::invalid_argument(this->name + " only works if the first departure event in your mission has a maneuver and therefore has an initial impulse. Right now the only boundary types that do this are 'EphemerisPeggedLaunchDirectInsertion,' 'FreePointDirectInsertion,' and 'PeriapseLaunch.' You have apparently not selected one of those.");
            }
        }//end constructor

        void MinimizeInitialImpulseObjective::calcbounds()
        {
            //there is only one derivative, with respect to the departure event's initial impulse magnitude
            this->Xindex_initial_impulse = this->myDepartureEvent->getXindex_of_initial_impulse();

            this->create_sparsity_entry(0,
                this->Xindex_initial_impulse,
                this->Gindex_derivatives_of_objective_function);
        }//end calcbounds()

        void MinimizeInitialImpulseObjective::process(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->ObjectiveValue = this->ObjectiveScale * X[this->Xindex_initial_impulse];
            F[0] = this->ObjectiveValue;

            if (needG)
            {
                size_t Gindex = this->Gindex_derivatives_of_objective_function.front();

                G[Gindex] = this->ObjectiveScale;
            }
        }//end process()

        void MinimizeInitialImpulseObjective::output(std::ofstream& outputfile)
        {
            outputfile << "Objective function is \"minimize initial impulse\"" << std::endl;
            outputfile << "J = " << this->ObjectiveValue _GETVALUE << std::endl;
        }
    }//end namespace ObjectiveFunctions
}//end namespace EMTG