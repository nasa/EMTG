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

//EMTGv9 SundmanCoastPhase
//Jacob Englander 4-5-2018

#pragma once

#include "CoastPhase.h"
#include "IntegratedPropagator.h"
#include "IntegrationScheme.h"
#include "SpacecraftAccelerationModel.h"
#include "SundmanSpacecraftEOM.h"

namespace EMTG
{
    namespace Phases
    {
        class SundmanCoastPhase : public CoastPhase
        {
        public:
            //constructor
            SundmanCoastPhase() : CoastPhase::CoastPhase() {};
            SundmanCoastPhase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions);

            virtual ~SundmanCoastPhase() {};

            virtual void initialize();

            //clone
            virtual SundmanCoastPhase* clone() const { return new SundmanCoastPhase(*this); }

            //output
            void output(std::ofstream& outputfile,
                size_t& eventcount);

            virtual void output_STMs() { this->CoastPhase::output_STMs(); };

            virtual void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

            void output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file);

        protected:
            //calcbounds
            void calcbounds_phase_flight_time();

            void calcbounds_match_point_constraints();

            //process
            void process_phase_flight_time(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

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

            //fields
            doubleType SundmanIndependentVariable;
            size_t Xindex_SundmanIndependentVariable;
            std::vector<bool> TruthTable_MatchConstraints_Derivative_wrt_SundmanIndependentVariable;
            std::vector<size_t> G_indices_match_point_constraints_wrt_SundmanIndependentVariable;
            double dStepSize_dPropagationVariable;
		private:
			Astrodynamics::SundmanSpacecraftEOM myEOM;
        };
    }//close namespace Phases
}//close namespace EMTG