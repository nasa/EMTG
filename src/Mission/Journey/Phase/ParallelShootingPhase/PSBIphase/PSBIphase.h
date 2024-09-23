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

//EMTGv9 parallel shooting with bounded impulse (PSBI) phase
//Jacob Englander 2/8/2019

#pragma once

#include "ParallelShootingPhase.h"

#include "IntegratedPropagator.h"
#include "IntegrationScheme.h"
#include "SpacecraftAccelerationModel.h"
#include "TimeDomainSpacecraftEOM.h"

#include "PSBIstep.h"

#include "boost/ptr_container/ptr_vector.hpp"

namespace EMTG
{
    namespace Phases
    {
        class PSBIstep;

        class PSBIphase : public ParallelShootingPhase
        {
        public:
            //constructor
            PSBIphase() {};
            PSBIphase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions);

            virtual ~PSBIphase();

            //clone
            virtual PSBIphase* clone() const { return new PSBIphase(*this); }

        protected:
            //we have to overload the initial and terminal coast propagations because the Kepler propagator does not propagate mass
            virtual void process_phase_initial_coast(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_phase_terminal_coast(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);
        };
    }//close namespace Phases
}//close namespace EMTG