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


//PSBI "last step" for EMTGv9
//Jacob Englander 2/8/2019

#include "PSBIlaststep.h"

namespace EMTG
{
    namespace Phases
    {
        PSBIlaststep::PSBIlaststep() {};

        PSBIlaststep::PSBIlaststep(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stepIndex,
            const size_t& stageIndex,
            ParallelShootingStep* previousStep,
            ParallelShootingPhase* myPhase,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            this->initialize(name,
                journeyIndex,
                phaseIndex,
                stepIndex,
                stageIndex,
                previousStep,
                myPhase,
                myUniverse,
                mySpacecraft,
                myOptions);
        }

        void PSBIlaststep::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stepIndex,
            const size_t& stageIndex,
            ParallelShootingStep* previousStep,
            ParallelShootingPhase* myPhase,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            this->PSBIstep::initialize(name,
                journeyIndex,
                phaseIndex,
                stepIndex,
                stageIndex,
                previousStep,
                myPhase,
                myUniverse,
                mySpacecraft,
                myOptions);
        }//end initialize

         //master calcbounds
        void PSBIlaststep::calcbounds_step()
        {
            this->calcbounds_step_left_state();

            this->calcbounds_step_control();

            this->calcbounds_step_left_match_point_constraints();

            this->calcbounds_step_main();

            this->calcbounds_distance_constraints();

            this->calcbounds_maneuver_constraints();

            this->calcbounds_step_right_match_point_constraints();
        }//end calcbounds_step

         //master process
        void PSBIlaststep::process_step(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_step_left_state(X, Xindex, F, Findex, G, needG);

            this->process_step_control(X, Xindex, F, Findex, G, needG);

            this->process_step_left_match_point_constraints(X, Xindex, F, Findex, G, needG);

            this->process_step_main(X, Xindex, F, Findex, G, needG);

            this->process_distance_constraints(X, Xindex, F, Findex, G, needG);

            this->process_maneuver_constraints(X, Xindex, F, Findex, G, needG);

            if (needG)
                this->process_derivative_tuples(X, Xindex, F, Findex, G, needG);

            this->process_step_right_match_point_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_step
    }//end namespace Phases
}//end namespace EMTG