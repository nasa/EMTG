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

//PSBI step factory
//Jacob Englander 2/8/2019

#include "PSBIstep_factory.h"

#include "PSBIstep.h"
#include "PSBIfirststep.h"
#include "PSBIlaststep.h"
#include "PSBIOneStepToRuleThemAll.h"


namespace EMTG
{
    namespace Phases
    {
        ParallelShootingStep* create_PSBI_step(const std::string& name,
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
            if (myPhase->get_num_steps() == 1)
            {
                return new PSBIOneStepToRuleThemAll(name,
                    journeyIndex,
                    phaseIndex,
                    stepIndex,
                    stageIndex,
                    myPhase,
                    myUniverse,
                    mySpacecraft,
                    myOptions);
            }
            else if (stepIndex == 0)
            {
                return new PSBIfirststep(name,
                    journeyIndex,
                    phaseIndex,
                    stepIndex,
                    stageIndex,
                    myPhase,
                    myUniverse,
                    mySpacecraft,
                    myOptions);
            }
            else if (stepIndex == myPhase->get_num_steps() - 1)
            {
                return new PSBIlaststep(name,
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
            else
            {
                return new PSBIstep(name,
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
            
        }
    }//end namespace Phases
}//end namespace EMTG