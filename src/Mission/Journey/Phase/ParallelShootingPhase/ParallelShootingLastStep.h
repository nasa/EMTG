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

//base parallel shooting "last step" for EMTGv9
//Jacob Englander 2-21-2018

#pragma once

#include "ParallelShootingStep.h"

namespace EMTG
{
    namespace Phases
    {
        class ParallelShootingLastStep : virtual public ParallelShootingStep
        {
        public:
            //constructor
            ParallelShootingLastStep();
            ParallelShootingLastStep(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& subphaseIndex,
                const size_t& stageIndex,
                ParallelShootingStep* previousStep,
                ParallelShootingPhase* myPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //destructor
            virtual ~ParallelShootingLastStep() {};

        protected:
            //calcbounds
            void calcbounds_step_right_match_point_constraints();

            //process
            void process_step_right_match_point_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //right match point
            std::vector<size_t> Findices_right_match_point_constraints;


            std::vector<size_t> ListOfVariablesAffectingPhaseTerminalCoastLeftState;
            std::vector< std::vector< std::tuple<size_t, size_t> > > DerivativesOfPhaseTerminalCoastLeftStateByVariable; //[varIndex][entryIndex]{stateIndex, dIndex}
            std::vector< std::vector<size_t> > Gindices_RightMatchPoint_wrt_PhaseTerminalCoastLeftState; //[varIndex][entryIndex]
            std::vector<size_t> ListOfTimeVariablesAffectingPhaseTerminalCoastLeftState;
            std::vector< std::vector< std::tuple<size_t, size_t> > > DerivativesOfPhaseTerminalCoastLeftStateByTimeVariable; //[varIndex][entryIndex]{stateIndex, dIndex}
            std::vector< std::vector<size_t> > Gindices_RightMatchPoint_wrt_PhaseTerminalCoastLeftStateTime; //[varIndex][entryIndex]

            std::vector< std::vector<size_t> > Gindices_StepRightMatchPoint_wrt_CurrentStepRightStateVariables;//stateIndex, varIndex
            std::vector< std::vector<size_t> > Gindices_StepRightMatchPoint_wrt_CurrentStepRightStateTimeVariables;//stateIndex, varIndex          

            std::vector< std::vector< std::vector<size_t> > > Gindices_StepRightMatchPoint_wrt_CurrentStepControlVariables;//[subStepIndex][controlIndex][entryIndex]
        };//end class ParallelShootingLastStep
    }//close namespace Phases
}//close namespace EMTG