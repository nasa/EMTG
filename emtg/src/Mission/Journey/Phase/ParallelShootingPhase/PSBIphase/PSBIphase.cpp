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

#include "PSBIphase.h"

#include "PropagatorFactory.h"
#include "PSBIstep_factory.h"

namespace EMTG
{
    namespace Phases
    {
        PSBIphase::PSBIphase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions) :
            ParallelShootingPhase::ParallelShootingPhase(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                previousPhase,
                Universe,
                mySpacecraft,
                myLaunchVehicle,
                myOptions)
        {

            //integrators for initial and terminal coasts
            if (this->hasInitialCoast)
            {
                this->InitialCoastPropagatorObject = CreatePropagator(this->myOptions,
                    this->myUniverse,
                    6,
                    this->state_after_initial_TCM,
                    this->state_after_initial_coast,
                    this->STM_initial_coast,
                    this->InitialCoast_dStatedIndependentVariable,
                    &this->ForcedCoast_dStepSize_dPropagationVariable);
            }

            if (this->hasTerminalCoast)
            {
                this->TerminalCoastPropagatorObject = CreatePropagator(this->myOptions,
                    this->myUniverse,
                    6,
                    this->state_at_end_of_phase,
                    this->state_before_terminal_coast,
                    this->STM_terminal_coast,
                    this->TerminalCoast_dStatedIndependentVariable,
                    &this->ForcedCoast_dStepSize_dPropagationVariable);
            }

            //steps
            ParallelShootingStep* previousStep = NULL;
            for (size_t stepIndex = 0; stepIndex < this->num_steps; ++stepIndex)
            {
                this->mySteps.push_back(create_PSBI_step(this->name + "_Step" + std::to_string(stepIndex),
                    this->journeyIndex,
                    this->phaseIndex,
                    stepIndex,
                    this->stageIndex,
                    previousStep,
                    this,
                    this->myUniverse,
                    this->mySpacecraft,
                    this->myOptions));

                previousStep = &this->mySteps.back();
            }
        }//end constructor

        PSBIphase::~PSBIphase()
        {
            //forced coast STMs and such
            if (this->hasInitialCoast)
            {
                delete this->InitialCoastPropagatorObject;
            }

            if (this->hasTerminalCoast)
            {
                delete this->TerminalCoastPropagatorObject;
            }
        }//end destructor

        void PSBIphase::process_phase_initial_coast(const std::vector<doubleType>& X,
                                                    size_t& Xindex,
                                                    std::vector<doubleType>& F,
                                                    size_t& Findex,
                                                    std::vector<double>& G,
                                                    const bool& needG)
        {
            this->ParallelShootingPhase::process_phase_initial_coast(X, Xindex, F, Findex, G, needG);

            if (this->hasInitialCoast)
            {
                //copy over mass, epoch, chemical fuel, and electric propellant
                for (size_t stateIndex : {6, 8, 9})
                {
                    this->state_after_initial_coast(stateIndex) = this->state_after_initial_TCM(stateIndex);
                    this->STM_Augmented_initial_coast(stateIndex, stateIndex) = 1.0;
                }
            }
        }//end process_phase_initial_coast

        void PSBIphase::process_phase_terminal_coast(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->ParallelShootingPhase::process_phase_terminal_coast(X, Xindex, F, Findex, G, needG);

            //copy over mass, epoch, chemical fuel, and electric propellant
            if (this->hasTerminalCoast)
            {
                for (size_t stateIndex : {6, 8, 9})
                {
                    this->state_before_terminal_coast(stateIndex) = this->state_at_end_of_phase(stateIndex);
                    this->STM_Augmented_terminal_coast(stateIndex, stateIndex) = 1.0;
                }
            }
        }//end process_phase_terminal_coast
    }//close namespace Phases
}//close namespace EMTG