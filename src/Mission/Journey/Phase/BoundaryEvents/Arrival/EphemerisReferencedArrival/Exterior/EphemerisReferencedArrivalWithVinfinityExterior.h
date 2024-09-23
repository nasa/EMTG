
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

#include "EphemerisReferencedArrivalExterior.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class EphemerisReferencedArrivalWithVinfinityExterior : virtual public EphemerisReferencedArrivalExterior
        {
        public:
            //specialized constructor
            EphemerisReferencedArrivalWithVinfinityExterior(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

        protected:
            void calcbounds_event_interface_state(const std::vector<double>& vinf_max,
                const std::vector<double>& RAbounds,
                const std::vector<double>& DECbounds,
                std::vector<double>& MassBounds,
                const std::vector<double>& EpochBounds,
                std::vector<size_t> timeVariables);

            void process_event_interface_state(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void calcbounds_event_left_side(std::vector<size_t> timeVariables);

            virtual void calcbounds_event_main() {};//doesn't need to do anything - just a stub

            void calcbounds_virtual_propellant_constraints() {};//doesn't need to do anything - just a stub

            void calcbounds_event_right_side(); //derivatives of the right-side state with respect to event main decision variables

            void calcbounds_virtual_deltav_constraint() {};//doesn't need to do anything - just a stub

            void process_virtual_propellant_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) {};//doesn't need to do anything - just a stub

            void process_virtual_deltav_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) {};//doesn't need to do anything - just a stub

            virtual void process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) {};//doesn't need to do anything - just a stub

            //fields
            std::vector<size_t> dIndex_interface_velocity_wrt_v;//vx, vy, vz
            std::vector<size_t> dIndex_interface_velocity_wrt_vRA;//vx, vy, vz
            std::vector<size_t> dIndex_interface_velocity_wrt_vDEC;//vx, vy, vz
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG