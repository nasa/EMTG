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

#include "MinimizeWaypointTrackingErrorObjective.h"

namespace EMTG
{
    namespace ObjectiveFunctions
{

        MinimizeWaypointTrackingErrorObjective::MinimizeWaypointTrackingErrorObjective(Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            Mission* myMission) :
            ObjectiveFunctionBase(Universe,
                mySpacecraft,
                myOptions,
                myMission)
        {
            this->name = "MinimizeWaypointTrackingErrorObjective";

            this->ObjectiveScale = 1.0;

            //this->myOptions->myWaypointTracker.fit_splines();                
        }//end constructor

        void MinimizeWaypointTrackingErrorObjective::calcbounds()
        {
            //this objective function has derivatives with respect to all virtual waypoint error variables
            for (size_t Xindex = 0; Xindex < this->Xdescriptions->size(); ++Xindex)
            {
                if (this->Xdescriptions->operator[](Xindex).find("virtual mahalanobis distance") < 1024)
                {
                    this->Xindex_virtual_mahalanobis_distance.push_back(Xindex);

                    this->create_sparsity_entry(0,
                        Xindex,
                        this->Gindex_virtual_mahalanobis_distance);
                }
            }
        }//end calcbounds()

        void MinimizeWaypointTrackingErrorObjective::process(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->ObjectiveValue = 0.0;
            
            for (size_t errorIndex : this->Xindex_virtual_mahalanobis_distance)
                this->ObjectiveValue += this->ObjectiveScale * X[errorIndex];

            F[0] = this->ObjectiveValue;

            if (needG)
            {
                for (size_t Gindex : this->Gindex_virtual_mahalanobis_distance)
                {
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * this->ObjectiveScale;
                }
            }
        }//end process()

        void MinimizeWaypointTrackingErrorObjective::output(std::ofstream& outputfile)
        {
            outputfile << "Objective function is \"minimize waypoint tracking error\"" << std::endl;
            outputfile << "J = " << this->ObjectiveValue _GETVALUE << std::endl;
        }
    }//end namespace ObjectiveFunctions
}//end namespace EMTG