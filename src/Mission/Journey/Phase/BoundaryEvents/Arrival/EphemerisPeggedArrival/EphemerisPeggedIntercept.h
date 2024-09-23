
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

#include "EphemerisPeggedFlybyIn.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class EphemerisPeggedIntercept : public EphemerisPeggedFlybyIn
        {
        public:
            //specialized constructor
            EphemerisPeggedIntercept(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //output
            virtual void output(std::ofstream& outputfile,
                const double& launchdate,
                size_t& eventcount);

            //process
            virtual void process_event(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //calcbounds
            void calcbounds(std::vector<size_t> timeVariables);
            
            //output
            void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

        protected:
            virtual void calcbounds_event_main();

            virtual void process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            math::Matrix<doubleType> get_periapse_state(); //this is private because we only use it to compute b-plane targets

            //fields
            std::vector<size_t> Gindices_incoming_velocity_magnitude_wrt_Vinfinity_in;
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG