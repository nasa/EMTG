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

//objective function factory
//Jacob Englander 4/1/2018

#include "ObjectiveFunctionFactory.h"

#include "MinimizeDeltavObjective.h"
#include "MaximizeMassObjective.h"
#include "MaximizeDryMassObjective.h"
#include "MaximizeLog10DryMassObjective.h"
#include "MaximizeLogeDryMassObjective.h"
#include "MinimizeTimeObjective.h"
#include "MaximizeInitialMassObjective.h"
#include "MaximizeLogeMassObjective.h"
#include "MaximizeLog10MassObjective.h"
#include "ArriveAsEarlyAsPossibleObjective.h"
#include "ArriveAsLateAsPossibleObjective.h"
#include "DepartAsEarlyAsPossibleObjective.h"
#include "DepartAsLateAsPossibleObjective.h"
#include "MinimizeChemicalFuelObjective.h"
#include "MinimizeChemicalOxidizerObjective.h"
#include "MinimizeElectricPropellantObjective.h"
#include "MinimizeTotalPropellantObjective.h"
#include "MinimizeWaypointTrackingErrorObjective.h"
#include "MinimizeInitialImpulseObjective.h"
#include "MaximizeDistanceFromCentralBodyObjective.h"

namespace EMTG
{

    namespace ObjectiveFunctions
    {
        ObjectiveFunctionBase* create_objective_function(Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            Mission* myMission)
        {
            switch (myOptions->objective_type)
            {
            case ObjectiveFunctionType::MINIMIZE_DELTAV:
            {
                return new ObjectiveFunctions::MinimizeDeltavObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::MAXIMIZE_MASS:
            {
                return new ObjectiveFunctions::MaximizeMassObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::MAXIMIZE_INITIAL_MASS:
            {
                return new ObjectiveFunctions::MaximizeInitialMassObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::MINIMIZE_TIME:
            {
                return new ObjectiveFunctions::MinimizeTimeObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
                break;
            }
            case ObjectiveFunctionType::MAXIMIZE_LOGE_MASS:
            {
                return new ObjectiveFunctions::MaximizeLogeMassObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::MAXIMIZE_LOG10_MASS:
            {
                return new ObjectiveFunctions::MaximizeLog10MassObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::ARRIVE_AS_EARLY_AS_POSSIBLE:
            {
                return new ObjectiveFunctions::ArriveAsEarlyAsPossibleObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::ARRIVE_AS_LATE_AS_POSSIBLE:
            {
                return new ObjectiveFunctions::ArriveAsLateAsPossibleObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::DEPART_AS_EARLY_AS_POSSIBLE:
            {
                return new ObjectiveFunctions::DepartAsEarlyAsPossibleObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::DEPART_AS_LATE_AS_POSSIBLE:
            {
                return new ObjectiveFunctions::DepartAsLateAsPossibleObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::MAXIMIZE_DRY_MASS:
            {
                return new ObjectiveFunctions::MaximizeDryMassObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::MAXIMIZE_LOG10_DRY_MASS:
            {
                return new ObjectiveFunctions::MaximizeLog10DryMassObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::MAXIMIZE_LOGE_DRY_MASS:
            {
                return new ObjectiveFunctions::MaximizeLogeDryMassObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::MINIMIZE_CHEMICAL_FUEL:
            {
                return new ObjectiveFunctions::MinimizeChemicalFuelObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::MINIMIZE_CHEMICAL_OXIDIZER:
            {
                return new ObjectiveFunctions::MinimizeChemicalOxidizerObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::MINIMIZE_ELECTRIC_PROPELLANT:
            {
                return new ObjectiveFunctions::MinimizeElectricPropellantobjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            case ObjectiveFunctionType::MINIMIZE_TOTAL_PROPELLANT:
            {
                return new ObjectiveFunctions::MinimizeTotalPropellantObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
            /*case ObjectiveFunctionType::MINIMIZE_WAYPOINT_TRACKING_ERROR:
            {
                return new ObjectiveFunctions::MinimizeWaypointTrackingErrorObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }*/
            case ObjectiveFunctionType::MINIMIZE_INITIAL_IMPULSE:
            {
                return new ObjectiveFunctions::MinimizeInitialImpulseObjective(Universe,
                    mySpacecraft,
                    myOptions,
                    myMission);
            }
			case ObjectiveFunctionType::MAXIMIZE_FINAL_PERIAPSIS:
			{
				return new ObjectiveFunctions::MaximizeDistanceFromCentralBodyObjective(Universe,
					mySpacecraft,
					myOptions,
					myMission);
			}
            default:
            {
                throw std::invalid_argument("Undefined objective function type! Ask Jacob to implement your objective function. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            }
        }//end create_objective_function()
    }//end namespace ObjectiveFunctions
}//end namespace EMTG