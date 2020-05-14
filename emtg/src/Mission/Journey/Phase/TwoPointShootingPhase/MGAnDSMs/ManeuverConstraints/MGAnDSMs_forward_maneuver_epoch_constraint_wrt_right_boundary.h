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

//forward MGAnDSMs maneuver epoch constraint with respect to right boundary
//9-26-2017

#pragma once

#include "doubleType.h"

#include "MGAnDSMs_maneuver_epoch_constraint.h"

namespace EMTG
{
    namespace Phases
    {
        //forward declare parent class
        class MGAnDSMs_subphase;

        class MGAnDSMs_forward_maneuver_epoch_constraint_wrt_right_boundary : public MGAnDSMs_maneuver_epoch_constraint
        {
        public:
            //constructor
            MGAnDSMs_forward_maneuver_epoch_constraint_wrt_right_boundary() {};
            MGAnDSMs_forward_maneuver_epoch_constraint_wrt_right_boundary(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& subphaseIndex,
                const size_t& stageIndex,
                MGAnDSMs_subphase* mySubPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                const std::string& ConstraintDefinition);

            //clone
            virtual MGAnDSMs_forward_maneuver_epoch_constraint_wrt_right_boundary* clone() const { return new MGAnDSMs_forward_maneuver_epoch_constraint_wrt_right_boundary(*this); }

            //calcbounds goes in the specialized phase
            void calcbounds();

            //process goes in the specialized phase
            void process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);
        };
    }//close namespace Phases
}//close namespace EMTG