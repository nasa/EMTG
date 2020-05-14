
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

#include "EphemerisPeggedArrivalWithVinfinity.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedArrivalWithVinfinity::EphemerisPeggedArrivalWithVinfinity(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            EphemerisPeggedArrival::EphemerisPeggedArrival(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions)
        {
            this->Vinfinity_in.resize(3, 1, 0.0);

            if (this->isLastEventInJourney)
            {
                if (this->myOptions->Journeys[this->journeyIndex].journey_end_TCM > math::SMALL)
                {
                    this->hasTCM = true;
                    this->TCM_magnitude = this->myOptions->Journeys[this->journeyIndex].journey_end_TCM;
                }
            }
            else if (this->myOptions->TCM_pre_flyby > math::SMALL)
            {
                this->hasTCM = true;
                this->TCM_magnitude = this->myOptions->TCM_pre_flyby;
            }
            else
                this->hasTCM = false;
        }

        //******************************************calcbounds methods
        void EphemerisPeggedArrivalWithVinfinity::
            calcbounds_event_main(const std::vector< std::tuple<double, double> >& vinfBounds)
        {
			this->Vinfinity_upperbound = std::get<1>(vinfBounds[0]);

            //encode v-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                Xlowerbounds->push_back(std::get<0>(vinfBounds[Vindex]));
                Xupperbounds->push_back(std::get<1>(vinfBounds[Vindex]));
                Xdescriptions->push_back(prefix + "V_infinity_" + this->stateNames[Vindex]);
                X_scale_factors->push_back(this->myUniverse->LU / this->myUniverse->TU);
                this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple<size_t, size_t, double>(Xdescriptions->size() - 1, Vindex + 3, 1.0));
                this->dIndex_VbeforeEvent_dVinfinity_in.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
            }
        }//end calcbounds_event_main()

        //******************************************process methods
        void EphemerisPeggedArrivalWithVinfinity::process_event_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //extract v-infinity and add to state_before_event
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->Vinfinity_in(Vindex) = X[Xindex++];
                this->state_before_event(3 + Vindex) += this->Vinfinity_in(Vindex);
            }
        }//end process_event_left_side()
    }//end namespace BoundaryEvents
}//end namespace EMTG