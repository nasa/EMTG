
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

#include "EphemerisReferencedDepartureInterior.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class EphemerisReferencedFreeDirectDepartureInterior : virtual public EphemerisReferencedDepartureInterior
        {
        public:
            //specialized constructor
            EphemerisReferencedFreeDirectDepartureInterior(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                ArrivalEvent* PreviousPhaseArrivalEvent);

            //output
            void output(std::ofstream& outputfile,
                const double& launchdate,
                size_t& eventcount);

            //process
            void process_event(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //calcbounds
            void calcbounds(std::vector<size_t> timeVariables);

        private:
            void calcbounds_event_interface_state(const std::vector<double>& RAbounds,
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

            void calcbounds_event_main() {};//doesn't need to do anything - just a stub

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

            void process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) {};//doesn't need to do anything - just a stub

            void process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //no fields, because this event is completely trivial
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG