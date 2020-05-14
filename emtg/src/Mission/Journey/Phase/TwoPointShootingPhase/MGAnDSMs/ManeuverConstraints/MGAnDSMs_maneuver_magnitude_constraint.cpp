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

//MGAnDSMs maneuver epoch constraint
//9-22-2017

#include "MGAnDSMs_maneuver_magnitude_constraint.h"

#include "EMTG_solver_utilities.h"

#include "boost/algorithm/string/split.hpp"

#include "MGAnDSMs_subphase.h"

namespace EMTG
{
    namespace Phases
    {
        MGAnDSMs_maneuver_magnitude_constraint::MGAnDSMs_maneuver_magnitude_constraint(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& subphaseIndex,
            const size_t& stageIndex,
            MGAnDSMs_subphase* mySubPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            const std::string& ConstraintDefinition) :
            MGAnDSMs_maneuver_magnitude_constraint()
        {
            this->journeyIndex = journeyIndex;
            this->phaseIndex = phaseIndex;
            this->subphaseIndex = subphaseIndex;
            this->stageIndex = stageIndex;
            this->Universe = Universe;
            this->mySpacecraft = mySpacecraft;
            this->myOptions = myOptions;
            this->myJourneyOptions = &(this->myOptions->Journeys[this->journeyIndex]);
            this->mySubPhase = mySubPhase;
            this->ConstraintDefinition = ConstraintDefinition;
            this->name = this->mySubPhase->getName();
        }//end constructor

         //************************************calcbounds methods
        void MGAnDSMs_maneuver_magnitude_constraint::calcbounds()
        {
            //Step 1: parse the constraint definition
            std::vector<std::string> ConstraintDefinitionCell;
            boost::split(ConstraintDefinitionCell, ConstraintDefinition, boost::is_any_of("_"), boost::token_compress_on);

            this->lowerBound = std::stod(ConstraintDefinitionCell[2]);
            this->upperBound = std::stod(ConstraintDefinitionCell[3]);

            if (lowerBound < 0.0)
                lowerBound = 0.0;

            if (upperBound < 1.0e-8)
                upperBound = 1.0e-8;

            //Step 2: create the constraint
            this->Flowerbounds->push_back(this->lowerBound - this->upperBound);
            this->Fupperbounds->push_back(0.0);
            this->Fdescriptions->push_back(prefix + "maneuver magnitude constraint");

            //sparsity pattern - derivatives with DSM components
            std::vector<size_t>& Xindex_DSM_components = this->mySubPhase->getXindex_DSM_components();
            
            for (size_t entry = 0; entry < Xindex_DSM_components.size(); ++entry)
            {
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    Xindex_DSM_components[entry],
                    this->Gindex_wrt_ManeuverComponents);
            }
        }//end calcbounds()

         //******************************************process methods
        void MGAnDSMs_maneuver_magnitude_constraint::process_constraint(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: evaluate the constraint
            math::Matrix<doubleType>& myDSM = this->mySubPhase->getDSM();
            doubleType DSM_magnitude = this->mySubPhase->getDSMmagnitude();

            F[Findex++] = DSM_magnitude - this->upperBound;

            //Step 2: derivatives
            if (needG)
            {
                
                for (size_t DSMcomponentIndex = 0; DSMcomponentIndex < this->Gindex_wrt_ManeuverComponents.size(); ++DSMcomponentIndex)
                {
                    size_t Gindex = this->Gindex_wrt_ManeuverComponents[DSMcomponentIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * (myDSM(DSMcomponentIndex) / DSM_magnitude)_GETVALUE;
                }
            }//end derivatives
        }//end process_constraint()
    }//close namespace Phases
}//close namespace EMTG