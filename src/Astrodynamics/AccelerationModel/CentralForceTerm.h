// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2014 - 2017 United States Government as represented by the
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

//central force term

#pragma once

#include "AccelerationModelTerm.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        class CentralForceTerm : public AccelerationModelTerm
        {
        public:
            //constructor
            CentralForceTerm(missionoptions* myOptions,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& STMrows,
                const size_t& STMcolumns);

            //methods
            void computeAccelerationTerm(const math::Matrix<doubleType>& spacecraft_state_relative_to_central_body,
                const math::Matrix <double> & dspacecraft_state_relative_to_central_bodydTOF,
                const doubleType& launch_epoch,
                const math::Matrix<doubleType>& control,
                const doubleType & epoch_step_left,
                std::vector <double> & depoch_left_segmentdTOF,
                const double & c,
                const doubleType & h,
                const double & dhdTOF,
                const bool& generate_derivatives);
        };
    }
}