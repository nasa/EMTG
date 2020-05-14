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

//PSFB "last step" for EMTGv9
//Jacob Englander 2-26-2018

#pragma once

#include "PSFBstep.h"
#include "ParallelShootingLastStep.h"

namespace EMTG
{
    namespace Phases
    {
        class PSFBlaststep : virtual public ParallelShootingLastStep, virtual public PSFBstep
        {
        public:
            //constructor
            PSFBlaststep();
            PSFBlaststep(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& subphaseIndex,
                const size_t& stageIndex,
                ParallelShootingStep* previousStep,
                ParallelShootingPhase* myPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            virtual void initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                const size_t& stageIndex,
                ParallelShootingStep* previousStep,
                ParallelShootingPhase* myPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //clone
            virtual PSFBlaststep* clone() const { return new PSFBlaststep(*this); }

            //destructor
            virtual ~PSFBlaststep() {};

            virtual void calcbounds_step();

            virtual void process_step(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);
        };//end class ParallelShootingLastStep
    }//close namespace Phases
}//close namespace EMTG