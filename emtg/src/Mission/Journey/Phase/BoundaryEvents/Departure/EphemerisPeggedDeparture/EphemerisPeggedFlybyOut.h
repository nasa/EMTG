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

#include "EphemerisPeggedDeparture.h"
#include "Arrival/EphemerisPeggedArrival/EphemerisPeggedArrivalWithVinfinity.h"
#include "bplane.h"

//small class that serves as an abstract base for flyby departures
//encodes a v-infinity vector
//does NOT include the flyby constraints

namespace EMTG
{
    namespace BoundaryEvents
    {
        class EphemerisPeggedFlybyOut : public EphemerisPeggedDeparture
        {
        public:
            //specialized constructor
            EphemerisPeggedFlybyOut(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                EphemerisPeggedArrivalWithVinfinity* PreviousPhaseArrivalEvent);

        protected:
            void calcbounds_event_left_side(std::vector<size_t> timeVariables);

            virtual void calcbounds_event_main(const std::vector< std::tuple<double, double> >& vinfBounds);

            void calcbounds_event_main() {};//this is a stub because in this event and its children, we have to pass in an arguement

            void process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            math::Matrix<doubleType> calculate_flyby_periapse_state();
            virtual void output_periapse_state(size_t& flybyIndex, std::ofstream& outputfile);

            //fields
            math::Matrix<doubleType> Vinfinity_out;
            std::vector<size_t> Xindices_Vinfinity_out;

            std::vector<size_t> dIndex_PreviousEventPosVelEpoch_PreviousTimeVariables;
            size_t dIndex_PreviousEventMass_PreviousEventMass;
            std::vector<size_t> dIndex_CurrentEventPosVelEpoch_PreviousTimeVariables;
            size_t dIndex_CurrentEventMass_PreviousEventMass;

            std::vector<size_t> dIndex_VbeforeEvent_dVinfinity_out;

            math::Matrix<doubleType> Vinfinity_in;
            std::vector<size_t> Xindices_Vinfinity_in;

            doubleType FlybyTurnAngle;
            doubleType FlybyAltitude;
            doubleType BdotR, BdotT;

            Astrodynamics::bplane myBplane;

			EphemerisPeggedArrivalWithVinfinity* PreviousPhaseArrivalEvent;//pointer to previous phase arrival event
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG