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

//PSFB "first step" for EMTGv9
//High-fidelity duty cycle
//Jacob Englander 3-19-2018

#pragma once

#include "PSFBfirststep.h"
#include "PSFB_HifiDuty_step.h"

namespace EMTG
{
    namespace Phases
    {
        class PSFB_HifiDuty_firststep : virtual public ParallelShootingFirstStep, virtual public PSFB_HifiDuty_step
        {
        public:
            //constructor
            PSFB_HifiDuty_firststep();
            PSFB_HifiDuty_firststep(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                const size_t& stageIndex,
                ParallelShootingPhase* myPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            virtual void initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                const size_t& stageIndex,
                ParallelShootingPhase* myPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //clone
            virtual PSFB_HifiDuty_firststep* clone() const { return new PSFB_HifiDuty_firststep(*this); }

            //destructor
            virtual ~PSFB_HifiDuty_firststep() {};

            virtual void calcbounds_step();

            virtual void process_step(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);
        };
    }//close namespace Phases
}//close namespace EMTG