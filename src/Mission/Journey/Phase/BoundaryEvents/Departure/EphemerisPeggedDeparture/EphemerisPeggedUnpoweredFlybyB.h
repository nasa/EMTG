
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

#include "EphemerisPeggedFlybyout.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class EphemerisPeggedUnpoweredFlyby : public EphemerisPeggedFlybyOut
        {
        public:
            //specialized constructor
            EphemerisPeggedUnpoweredFlyby(const std::string name,
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
            
            math::Matrix<doubleType> get_periapse_state();

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
            void calcbounds_event_main(const std::vector< std::tuple<double, double> >& vinfBounds);
            
            void calcbounds_virtual_propellant_constraints() {};//stub
            
            void calcbounds_virtual_deltav_constraint() {};//stub

            void process_virtual_deltav_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) {};//stub

            void process_virtual_propellant_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) {};//stub

            void process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //fields
            doubleType FlybyPeriapseDistance;

            std::vector< std::vector<size_t> > Gindices_Vinfinity_magnitude_match; //first incoming, then outgoing
            std::vector< std::vector<size_t> > Gindices_turn_angle; //first incoming, then outgoing


            std::vector< std::vector<size_t> > Gindices_TurnAngleConstraint_with_respect_to_Vinfinity; //first incoming, then outgoing
            size_t Gindices_TurnAngleConstraint_with_respect_to_FlybyPeriapseDistance;
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG