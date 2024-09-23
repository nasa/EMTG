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

//forward MGAnDSMs maneuver epoch constraint with respect to left boundary
//9-26-2017

#include "MGAnDSMs_forward_maneuver_epoch_constraint_wrt_left_boundary.h"

#include "EMTG_solver_utilities.h"

#include "MGAnDSMs_subphase.h"

#include "boost/algorithm/string/split.hpp"

namespace EMTG
{
    namespace Phases
    {
        MGAnDSMs_forward_maneuver_epoch_constraint_wrt_left_boundary::MGAnDSMs_forward_maneuver_epoch_constraint_wrt_left_boundary(const std::string& name,
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
        void MGAnDSMs_forward_maneuver_epoch_constraint_wrt_left_boundary::calcbounds()
        {
            //Step 1: parse the constraint definition
            std::vector<std::string> ConstraintDefinitionCell;
            boost::split(ConstraintDefinitionCell, ConstraintDefinition, boost::is_any_of("_"), boost::token_compress_on);

            //Step 2: create the constraint
            this->lowerBound = std::stod(ConstraintDefinitionCell[3]) / 100.0;
            this->upperBound = std::stod(ConstraintDefinitionCell[4]) / 100.0;
            this->Flowerbounds->push_back(this->lowerBound - this->upperBound);
            this->Fupperbounds->push_back(0.0);
            this->Fdescriptions->push_back(prefix + "maneuver epoch relative to left boundary constraint");

            //sparsity pattern - derivatives with respect to phase flight time
            this->Xindex_PhaseFlightTime = this->mySubPhase->get_timeVariables().back();

            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                this->Xindex_PhaseFlightTime,
                this->Gindex_wrt_TimeVariables);

            //derivatives with respect to burn index - all burn indices between now and the left boundary
            std::vector<size_t>& Xindex_burnIndex = this->mySubPhase->getXindex_burnIndex();
            for (size_t entry = 0; entry < Xindex_burnIndex.size(); ++entry)
            {
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    Xindex_burnIndex[entry],
                    this->Gindex_wrt_BurnIndices);
            }            
        }//end calcbounds()

         //******************************************process methods
        void MGAnDSMs_forward_maneuver_epoch_constraint_wrt_left_boundary::process_constraint(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: apply the constraint
            doubleType PhaseFlightTime = X[this->Xindex_PhaseFlightTime];
            std::vector<size_t>& Xindex_burnIndex = this->mySubPhase->getXindex_burnIndex();
            doubleType SumOfBurnIndices = 0.0;

            for (size_t BurnIndexIndex = 0; BurnIndexIndex < Xindex_burnIndex.size(); ++BurnIndexIndex)
                SumOfBurnIndices += X[Xindex_burnIndex[BurnIndexIndex]];

            doubleType TimeSinceLeftBoundary = SumOfBurnIndices * PhaseFlightTime;

            F[Findex++] = TimeSinceLeftBoundary / 86400.0 / 100.0 - this->upperBound;

            //Step 2: derivatives
            if (needG)
            {
                //derivative with respect to phase flight time
                size_t Gindex = this->Gindex_wrt_TimeVariables.front();
                size_t Xindex = this->Xindex_PhaseFlightTime;

                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                    * SumOfBurnIndices _GETVALUE
                    / 86400.0 / 100.0;

                //derivative with respect to burn indices
                for (size_t BurnIndexIndex = 0; BurnIndexIndex < Xindex_burnIndex.size(); ++BurnIndexIndex)
                {
                    size_t Gindex = this->Gindex_wrt_BurnIndices[BurnIndexIndex];
                    size_t Xindex = Xindex_burnIndex[BurnIndexIndex];

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * PhaseFlightTime _GETVALUE
                        / 86400.0 / 100.0;
                }
            }
        }//end process_constraint()
    }//close namespace Phases
}//close namespace EMTG