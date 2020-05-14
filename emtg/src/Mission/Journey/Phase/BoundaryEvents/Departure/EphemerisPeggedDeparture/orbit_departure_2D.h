
// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2014 - 2017 United States Government as represented by the
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

#include "departure.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class OrbitDeparture2D : public DepartureEvent
        {
        public:
            //specialized constructor
            OrbitDeparture2D(const std::string& name,
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
            void calcbounds();

        private:
            void calcbounds_event_main();

            void calcbounds_event_right_side(); //derivatives of the right-side linkage constraints with respect to event main decision variables

            void calcbounds_virtual_propellant_constraints();

            void calcbounds_deltav_contribution();

            void process_virtual_propellant_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_deltav_contribution(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);


            //fields
            doubleType C3;
            doubleType RLA;
            doubleType DLA;
            doubleType dvDeparture;

            //derivative entries
            double ddvDeparture_dVinfinityOut;
            double dm_dVinfinityOut;

            //index holders
            std::vector<size_t> Gindices_DLA_must_be_compatible_with_INC;

            //derivatives of state after event linkage constraint
            std::vector<size_t> Gindices_dStateAfterEventLinkage_dVinfinity;
            std::vector<size_t> Gindices_dStateAfterEventLinkage_dRLA;
            std::vector<size_t> Gindices_dStateAfterEventLinkage_dDLA;

            //delta-v
            size_t Gindex_dVirtualEventTotalDeltav_dVinfinity;

            //propellant
            size_t Gindices_dVirtualChemicalFuel_dLeftMass;
            size_t Gindices_dVirtualChemicalFuel_dVinfinity;
            size_t Gindices_dVirtualChemicalOxidizer_dLeftMass;
            size_t Gindices_dVirtualChemicalOxidizer_dVinfinity;

            double dChemicalFuel_dVinfinity;
            double dChemicalOxidizer_dVinfinity;
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG