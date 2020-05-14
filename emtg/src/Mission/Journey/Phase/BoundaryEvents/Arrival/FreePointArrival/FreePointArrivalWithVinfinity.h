
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

#pragma once

#include "FreePointArrival.h"

//small class that serves as an abstract base for arrivals that need a v-infinity vector

namespace EMTG
{
    namespace BoundaryEvents
    {
        class FreePointArrivalWithVinfinity : public FreePointArrival
        {
        public:
            //specialized constructor
            FreePointArrivalWithVinfinity(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

        protected:
            void calcbounds_event_main() {}; //stub because otherwise it's pure abstract

            virtual void calcbounds_event_main(const std::vector< std::tuple<double, double> >& vinfBounds);

            virtual void calcbounds_virtual_propellant_constraints() = 0;

            virtual void process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_virtual_propellant_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            //fields
            math::Matrix<doubleType> Vinfinity_in;
            std::vector<size_t> dIndex_VbeforeEvent_dVinfinity_in;

            
            //TCM propellant
            size_t Gindices_dVirtualChemicalFuel_TCM_dLeftMass;

            std::vector<size_t> Gindices_incoming_velocity_magnitude_wrt_Vinfinity_in;
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG