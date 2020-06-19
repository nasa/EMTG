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

//forward subphase for EMTGv9 MGAnDSMs
//Jacob Englander 8-24-2017

#pragma once

#include "MGAnDSMs_subphase.h"

namespace EMTG
{
    namespace Phases
    {
        class Forward_MGAnDSMs_subphase : public MGAnDSMs_subphase
        {
        public:
            Forward_MGAnDSMs_subphase() : MGAnDSMs_subphase::MGAnDSMs_subphase() {};
            Forward_MGAnDSMs_subphase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& subphaseIndex,
                const size_t& stageIndex,
                MGAnDSMs_phase* myPhase,
                MGAnDSMs_subphase* previousSubPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //clone
            virtual Forward_MGAnDSMs_subphase* clone() const { return new Forward_MGAnDSMs_subphase(*this); }

            void output(std::ofstream& outputfile,
                size_t& eventcount);
            
            void output_ephemeris(std::ofstream& outputfile, std::ofstream & acceleration_model_file);

            void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

            void calcbounds(std::vector<size_t>& timeVariables);

            void process_subphase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void configure_propagator(); 
            
            inline math::Matrix<double> get_dMassAfterTCM_dDSMcomponents() const { return this->dMassAfterTCM_dDSMcomponents; }
            inline bool getBordersBoundary() const { return this->isFirstSubPhase; };

        protected:
            void calcbounds_subphase_time();

            void process_subphase_time(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //fields
            bool isFirstSubPhase;
            math::Matrix<double> dMassAfterTCM_dDSMcomponents;

            double dFuelConsumedDSM_dMassAtDSM;
            double dOxidizerConsumedDSM_dMassAtDSM;
            double dMassAfterDSM_dMassAtDSM; 
            double dFuelConsumedDSM_ddeltav;
            double dOxidizerConsumedDSM_ddeltav;
            double dMassAfterDSM_ddeltav;
            double dFuelConsumedTCM_dMassAtTCM;
            double dOxidizerConsumedTCM_dMassAtTCM;
            double dMassAfterTCM_dMassAtTCM;
            double dFuelConsumedTCM_ddeltav;
            double dOxidizerConsumedTCM_ddeltav;
            double dMassAfterTCM_ddeltav;

        };//end class Forward_MGAnDSMs_subphase
    }//end namespace Phases
}//end namespace EMTG