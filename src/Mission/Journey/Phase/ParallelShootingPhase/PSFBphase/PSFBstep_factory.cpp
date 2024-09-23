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

//PSFB step factory
//Jacob Englander 3/19/2018

#include "PSFBstep_factory.h"

#include "PSFBstep.h"
#include "PSFBfirststep.h"
#include "PSFBlaststep.h"
#include "PSFBOneStepToRuleThemAll.h"

#include "PSFB_HifiDuty_step.h"
#include "PSFB_HifiDuty_firststep.h"
#include "PSFB_HifiDuty_laststep.h"
#include "PSFB_HifiDuty_OneStepToRuleThemAll.h"


namespace EMTG
{
    namespace Phases
    {
        ParallelShootingStep* create_PSFB_step(const std::string& name,
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
            if (myOptions->duty_cycle_type == DutyCycleType::Averaged)
            {
                if (myPhase->get_num_steps() == 1)
                {
                    return new PSFBOneStepToRuleThemAll(name,
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
                    return new PSFBfirststep(name,
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
                    return new PSFBlaststep(name,
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
                    return new PSFBstep(name,
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
            }//end averaged duty cycle step types
            else if (myOptions->duty_cycle_type == DutyCycleType::Realistic)
            {
                if (myPhase->get_num_steps() == 1)
                {
                    return new PSFB_HifiDuty_OneStepToRuleThemAll(name,
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
                    return new PSFB_HifiDuty_firststep(name,
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
                    return new PSFB_HifiDuty_laststep(name,
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
                    return new PSFB_HifiDuty_step(name,
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
            else
            {
                throw std::invalid_argument("PSFBstep_factory: Unrecognized duty cycle type in journey " + std::to_string(journeyIndex) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
        }
    }//end namespace Phases
}//end namespace EMTG