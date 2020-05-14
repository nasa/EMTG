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

//base parallel shooting "first step" for EMTGv9
//Jacob Englander 2-21-2018

#pragma once

#include "ParallelShootingStep.h"

namespace EMTG
{
    namespace Phases
    {
        class ParallelShootingFirstStep : virtual public ParallelShootingStep
        {
        public:
            //constructor
            ParallelShootingFirstStep();
            ParallelShootingFirstStep(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                const size_t& stageIndex,
                ParallelShootingPhase* myPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //destructor
            virtual ~ParallelShootingFirstStep() {};
            
        protected:
            //calcbounds
            void calcbounds_step_left_match_point_constraints();

            //process
            void process_step_left_match_point_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //fields
            std::vector<size_t> ListOfVariablesAffectingPreviousStepRightState;
            std::vector< std::vector< std::tuple<size_t, size_t> > > DerivativesOfPreviousStepRightStateByVariable; //[varIndex][entryIndex]{stateIndex, dIndex}
            std::vector< std::vector<size_t> > Gindices_LeftMatchPoint_wrt_PreviousStepRightState; //[varIndex][entryIndex]
            std::vector<size_t> ListOfTimeVariablesAffectingPreviousStepRightState;
            std::vector< std::vector< std::tuple<size_t, size_t> > > DerivativesOfPreviousStepRightStateByTimeVariable; //[varIndex][entryIndex]{stateIndex, dIndex}
            std::vector< std::vector<size_t> > Gindices_LeftMatchPoint_wrt_PreviousStepRightStateTime; //[varIndex][entryIndex]
        };
    }//close namespace Phases
}//close namespace EMTG