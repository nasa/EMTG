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

//MGAnDSMs maneuver cone angle exclusion constraint
//3-27-2018

#include "MGAnDSMs_maneuver_cone_angle_exclusion_constraint.h"

#include "EMTG_solver_utilities.h"

#include "boost/algorithm/string/split.hpp"

namespace EMTG
{
    namespace Phases
    {
        MGAnDSMs_maneuver_cone_angle_exclusion_constraint::MGAnDSMs_maneuver_cone_angle_exclusion_constraint(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& subphaseIndex,
            const size_t& stageIndex,
            MGAnDSMs_subphase* mySubPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            const std::string& ConstraintDefinition) :
            MGAnDSMs_maneuver_cone_angle_exclusion_constraint()
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

            //parse the constraint
            std::vector<std::string> ConstraintDefinitionCell;
            boost::split(ConstraintDefinitionCell, ConstraintDefinition, boost::is_any_of("_"), boost::token_compress_on);

            std::string reference_body_name;
            if (boost::to_lower_copy(ConstraintDefinitionCell[2]) == "cb")
            {
                this->reference_body_index = -2;
                reference_body_name = this->Universe->central_body_name;
            }                
            else
            {
                this->reference_body_index = std::stoi(ConstraintDefinitionCell[2]);
                reference_body_name = this->Universe->bodies[this->reference_body_index].name;
            }

            this->minimum_cone_angle = std::stod(ConstraintDefinitionCell[3]);

            if (this->minimum_cone_angle < 0.0)
                this->minimum_cone_angle = 0.0;

            this->name = this->mySubPhase->getName() + " maneuver cone angle exclusion constraint relative to " + reference_body_name;
        }//end constructor

         //************************************calcbounds methods
        void MGAnDSMs_maneuver_cone_angle_exclusion_constraint::calcbounds()
        {
            //Step 1: create the constraint
            this->Flowerbounds->push_back(this->minimum_cone_angle * math::deg2rad);
            this->Fupperbounds->push_back(math::PI);
            this->Fdescriptions->push_back(prefix + "cone angle ");

            //Step 2: sparsity pattern
            //Step 2.1: derivatives with respect to DSM components
            std::vector<size_t>& Xindex_DSM_components = this->mySubPhase->getXindex_DSM_components();
            
            for (size_t entry = 0; entry < Xindex_DSM_components.size(); ++entry)
            {
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    Xindex_DSM_components[entry],
                    this->Gindex_wrt_ManeuverComponents);
            }

            //Steo 2.2: derivatives with respect to non-time variables that affect the spacecraft state

            //Steo 2.3: derivatives with respect to time variables that affect the spacecraft state

        }//end calcbounds()

         //******************************************process methods
        void MGAnDSMs_maneuver_cone_angle_exclusion_constraint::process_constraint(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            
        }//end process_constraint()
    }//close namespace Phases
}//close namespace EMTG