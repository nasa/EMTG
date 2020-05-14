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

//MGAnDSMs maneuver constraint factory
//9-26-2017

#include "MGAnDSMs_maneuver_constraint_factory.h"

#include "MGAnDSMs_subphase.h"

#include "MGAnDSMs_maneuver_magnitude_constraint.h"
#include "MGAnDSMs_maneuver_cone_angle_exclusion_constraint.h"

#include "MGAnDSMs_forward_maneuver_epoch_constraint_wrt_previous_event.h"
#include "MGAnDSMs_forward_maneuver_epoch_constraint_wrt_next_event.h"
#include "MGAnDSMs_forward_maneuver_epoch_constraint_wrt_left_boundary.h"
#include "MGAnDSMs_forward_maneuver_epoch_constraint_wrt_right_boundary.h"
#include "MGAnDSMs_forward_maneuver_epoch_constraint_absolute.h"

#include "MGAnDSMs_backward_maneuver_epoch_constraint_wrt_previous_event.h"
#include "MGAnDSMs_backward_maneuver_epoch_constraint_wrt_next_event.h"
#include "MGAnDSMs_backward_maneuver_epoch_constraint_wrt_left_boundary.h"
#include "MGAnDSMs_backward_maneuver_epoch_constraint_wrt_right_boundary.h"
#include "MGAnDSMs_backward_maneuver_epoch_constraint_absolute.h"

namespace EMTG
{
    namespace Phases
    {
        MGAnDSMs_maneuver_constraint* create_MGAnDSMs_maneuver_constraint(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& subphaseIndex,
            const size_t& stageIndex,
            MGAnDSMs_subphase* mySubPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            const std::string& ConstraintDefinition,
            const std::string& HalfPhaseDefinition)
        {
            std::vector<std::string> ConstraintDefinitionCell;
            boost::split(ConstraintDefinitionCell, ConstraintDefinition, boost::is_any_of("_"), boost::token_compress_on);

            if (boost::to_lower_copy(ConstraintDefinitionCell[1]).find("epoch") < 1024)
            {
                if (HalfPhaseDefinition == "Forward")
                {
                    if (ConstraintDefinitionCell[2].find("prev") < 1024)
                    {
                        return new MGAnDSMs_forward_maneuver_epoch_constraint_wrt_previous_event(name,
                            journeyIndex,
                            phaseIndex,
                            subphaseIndex,
                            stageIndex,
                            mySubPhase,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("next") < 1024)
                    {
                        return new MGAnDSMs_forward_maneuver_epoch_constraint_wrt_next_event(name,
                            journeyIndex,
                            phaseIndex,
                            subphaseIndex,
                            stageIndex,
                            mySubPhase,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("lboundary") < 1024)
                    {
                        return new MGAnDSMs_forward_maneuver_epoch_constraint_wrt_left_boundary(name,
                            journeyIndex,
                            phaseIndex,
                            subphaseIndex,
                            stageIndex,
                            mySubPhase,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("rboundary") < 1024)
                    {
                        return new MGAnDSMs_forward_maneuver_epoch_constraint_wrt_right_boundary(name,
                            journeyIndex,
                            phaseIndex,
                            subphaseIndex,
                            stageIndex,
                            mySubPhase,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("abs") < 1024)
                    {
                        return new MGAnDSMs_forward_maneuver_epoch_constraint_absolute(name,
                            journeyIndex,
                            phaseIndex,
                            subphaseIndex,
                            stageIndex,
                            mySubPhase,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            ConstraintDefinition);
                    }
                    else
                    {
                        std::cout << "Invalid maneuver constraint, " << ConstraintDefinition << std::endl;
                        return NULL;
                    }
                }//end forward
                else //backward half-phase
                {
                    if (ConstraintDefinitionCell[2].find("prev") < 1024)
                    {
                        return new MGAnDSMs_backward_maneuver_epoch_constraint_wrt_previous_event(name,
                            journeyIndex,
                            phaseIndex,
                            subphaseIndex,
                            stageIndex,
                            mySubPhase,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("next") < 1024)
                    {
                        return new MGAnDSMs_backward_maneuver_epoch_constraint_wrt_next_event(name,
                            journeyIndex,
                            phaseIndex,
                            subphaseIndex,
                            stageIndex,
                            mySubPhase,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("lboundary") < 1024)
                    {
                        return new MGAnDSMs_backward_maneuver_epoch_constraint_wrt_left_boundary(name,
                            journeyIndex,
                            phaseIndex,
                            subphaseIndex,
                            stageIndex,
                            mySubPhase,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("rboundary") < 1024)
                    {
                        return new MGAnDSMs_backward_maneuver_epoch_constraint_wrt_right_boundary(name,
                            journeyIndex,
                            phaseIndex,
                            subphaseIndex,
                            stageIndex,
                            mySubPhase,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("abs") < 1024)
                    {
                        return new MGAnDSMs_backward_maneuver_epoch_constraint_absolute(name,
                            journeyIndex,
                            phaseIndex,
                            subphaseIndex,
                            stageIndex,
                            mySubPhase,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            ConstraintDefinition);
                    }
                    else
                    {
                        std::cout << "Invalid maneuver constraint, " << ConstraintDefinition << std::endl;
                        return NULL;
                    }
                }//end backward maneuver epoch constraints
            }//end maneuver epoch constraints
            else if (boost::to_lower_copy(ConstraintDefinitionCell[1]).find("magnitude") < 1024)
            {
                return new MGAnDSMs_maneuver_magnitude_constraint(name,
                    journeyIndex,
                    phaseIndex,
                    subphaseIndex,
                    stageIndex,
                    mySubPhase,
                    Universe,
                    mySpacecraft,
                    myOptions,
                    ConstraintDefinition);
            }
            else if (boost::to_lower_copy(ConstraintDefinitionCell[1]).find("coneangleexclusion") < 1024)
            {
                return new MGAnDSMs_maneuver_cone_angle_exclusion_constraint(name,
                    journeyIndex,
                    phaseIndex,
                    subphaseIndex,
                    stageIndex,
                    mySubPhase,
                    Universe,
                    mySpacecraft,
                    myOptions,
                    ConstraintDefinition);
            }
            else
            {
                std::cout << "Invalid maneuver constraint, " << ConstraintDefinition << std::endl;
                return NULL;
            }
        }
    }//end namespace Phases
}//end namespace EMTG