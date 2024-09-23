
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

#include "EphemerisPeggedFlybyOut.h"

//class for zero-turn flybys

namespace EMTG
{
    namespace BoundaryEvents
    {
        class EphemerisPeggedZeroTurnFlyby : public EphemerisPeggedFlybyOut
        {
        public:
            //specialized constructor
            EphemerisPeggedZeroTurnFlyby(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
				EphemerisPeggedArrivalWithVinfinity* PreviousPhaseArrivalEvent);

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
            void calcbounds(std::vector<size_t> timeVariables);

            //output
            void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

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
            std::vector< std::vector<size_t> > Gindices_Vinfinity_match; //first incoming, then outgoing
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG