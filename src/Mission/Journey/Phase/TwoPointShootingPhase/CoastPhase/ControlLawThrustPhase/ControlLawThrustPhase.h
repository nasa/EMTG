// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2024 United States Government as represented by the
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

//EMTGv9 ControlLawThrustPhase
//Jacob Englander 9-5-2020

#pragma once

#include "CoastPhase.h"

namespace EMTG
{
    namespace Phases
    {
        class ControlLawThrustPhase : public CoastPhase
        {
        public:
            //constructor
            ControlLawThrustPhase() : CoastPhase::CoastPhase() {};
            ControlLawThrustPhase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions);

            virtual ~ControlLawThrustPhase() {};

            virtual void initialize();

            //clone
            virtual ControlLawThrustPhase* clone() const { return new ControlLawThrustPhase(*this); }

            //output
            void output(std::ofstream& outputfile,
                size_t& eventcount);

            virtual void output_STMs() { this->CoastPhase::output_STMs(); };

            virtual void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

            void output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file);

            void process_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

        protected:
            //calcbounds
            void calcbounds_match_point_constraints();

            void calcbounds_virtual_propellant_tanks();

            //process

            void process_forward_half_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_backward_half_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_match_point_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_virtual_propellant_tanks(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);
        };//close class ControlLawThrustPhase
    }//close namespace Phases
}//close namespace EMTG