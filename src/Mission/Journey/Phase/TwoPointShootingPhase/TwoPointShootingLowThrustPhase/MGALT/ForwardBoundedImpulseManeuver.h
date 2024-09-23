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

//Jacob Englander 6/25/2017

#pragma once

#include "BoundedImpulseManeuver.h"




namespace EMTG 
{ 
    namespace Phases
    {
        class ForwardBoundedImpulseManeuver : public BoundedImpulseManeuver
        {
        public:
            //constructor
            ForwardBoundedImpulseManeuver() : BoundedImpulseManeuver() {};
            ForwardBoundedImpulseManeuver(const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //methods
            void process_maneuver(const bool& needMTM);

            //get
            inline math::Matrix<double> get_dMassAfterManeuver_dThrottleComponents() const { return this->dMassAfterManeuver_dThrottleComponents; }

        protected:
            //fields
            math::Matrix<double> dMassAfterManeuver_dThrottleComponents;
        };
    } //close namespace Astrodynamics
} //close namespace EMTG