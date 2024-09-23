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

//backward MGAnDSMs maneuver epoch constraint with respect to next event
//9-26-2017


#include "MGAnDSMs_backward_maneuver_epoch_constraint_wrt_next_event.h"

#include "EMTG_solver_utilities.h"

#include "MGAnDSMs_subphase.h"

#include "boost/algorithm/string/split.hpp"

#include <exception>

namespace EMTG
{
    namespace Phases
    {
        MGAnDSMs_backward_maneuver_epoch_constraint_wrt_next_event::MGAnDSMs_backward_maneuver_epoch_constraint_wrt_next_event(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& subphaseIndex,
            const size_t& stageIndex,
            MGAnDSMs_subphase* mySubPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            const std::string& ConstraintDefinition) :
            MGAnDSMs_maneuver_epoch_constraint::MGAnDSMs_maneuver_epoch_constraint(name,
                journeyIndex,
                phaseIndex,
                subphaseIndex,
                stageIndex,
                mySubPhase,
                Universe,
                mySpacecraft,
                myOptions,
                ConstraintDefinition)
        {
        }//end constructor

        //************************************calcbounds methods
        void MGAnDSMs_backward_maneuver_epoch_constraint_wrt_next_event::calcbounds()
        {
            //Step 0: can we even create this constraint?
            if (this->mySubPhase->getBordersBoundary())
            {
                throw std::invalid_argument("You cannot constrain the epoch of an MGAnDSMs phase's final backward subphase maneuver with respect to the next maneuver. (1) There is no next maneuver, you're at the end of the phase. (2) There is no maneuver in the final backward subphase anyway, by definition.\n" + this->ConstraintDefinition + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            //Step 1: parse the constraint definition
            std::vector<std::string> ConstraintDefinitionCell;
            boost::split(ConstraintDefinitionCell, ConstraintDefinition, boost::is_any_of("_"), boost::token_compress_on);

            //Step 2: create the constraint
            this->lowerBound = std::stod(ConstraintDefinitionCell[3]) / 100.0;
            this->upperBound = std::stod(ConstraintDefinitionCell[4]) / 100.0;
            this->Flowerbounds->push_back(this->lowerBound - this->upperBound);
            this->Fupperbounds->push_back(0.0);
            this->Fdescriptions->push_back(prefix + "maneuver epoch relative to next event constraint");

            //sparsity pattern - derivatives with respect to phase flight time
            this->Xindex_PhaseFlightTime = this->mySubPhase->get_timeVariables().back();

            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                this->Xindex_PhaseFlightTime,
                this->Gindex_wrt_TimeVariables);

            //derivatives with respect to burn index - THE NEXT SUBPHASE's current burn index because that's the measurement between now and the next event
            std::vector<size_t>& Xindex_burnIndex = this->mySubPhase->getXindex_burnIndex();
            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                Xindex_burnIndex[this->subphaseIndex - 1],
                this->Gindex_wrt_BurnIndices);
        }//end calcbounds()

         //******************************************process methods
        void MGAnDSMs_backward_maneuver_epoch_constraint_wrt_next_event::process_constraint(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: evaluate the constraint
            doubleType SubPhaseTime = this->mySubPhase->get_previousSubphase()->getSubPhaseTime();

            F[Findex++] = SubPhaseTime / 86400.0 / 100.0 - this->upperBound;

            //Step 2: derivatives
            if (needG)
            {
                size_t Gindex_BurnIndex = this->Gindex_wrt_BurnIndices.back();
                size_t Xindex_BurnIndex = this->jGvar->operator[](Gindex_BurnIndex);
                size_t Gindex_PhaseFlightTime = this->Gindex_wrt_TimeVariables.back();
                size_t Xindex_PhaseFlightTime = this->jGvar->operator[](Gindex_PhaseFlightTime);

                doubleType PhaseFlightTime = X[Xindex_PhaseFlightTime];
                doubleType burnIndex = X[Xindex_BurnIndex];

                //with respect to burn index                        
                G[Gindex_BurnIndex] = this->X_scale_factors->operator[](Xindex_BurnIndex)
                    * PhaseFlightTime _GETVALUE
                    / 86400.0 / 100.0;

                //phase flight time
                G[Gindex_PhaseFlightTime] = this->X_scale_factors->operator[](Xindex_PhaseFlightTime)
                    * burnIndex _GETVALUE
                    / 86400.0 / 100.0;
            }//end derivatives
        }//end process_constraint()
    }//close namespace Phases
}//close namespace EMTG