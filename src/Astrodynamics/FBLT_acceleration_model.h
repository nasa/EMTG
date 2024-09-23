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

//force model for EMTGv9
//Donald Ellison August 1st 2015

#ifndef FBLT_ACCELERATION_MODEL
#define FBLT_ACCELERATION_MODEL

#include "doubleType.h"
#include "universe.h"
#include "missionoptions.h"
#include "Spacecraft.h"

namespace EMTG { 
    namespace Astrodynamics {

        int FBLT_acceleration_model(const EMTG::missionoptions & options,
            Astrodynamics::universe & Universe,
            HardwareModels::Spacecraft * mySpacecraft,
            std::vector <doubleType> & spacecraft_state_relative_to_central_body,
            EMTG::math::Matrix <double> & dspacecraft_state_relative_to_central_bodydTOF,
            const doubleType & epoch_step_left,
            std::vector <double> & depoch_left_segmentdTOF,
            const double & c,
            const doubleType & h,
            const double & dhdTOF,
            const doubleType & launch_epoch,
            const std::vector <doubleType> & control,
            const doubleType & u_command,
            const double& DutyCycle,
            EMTG::math::Matrix <double> & dfdTOF,
            doubleType & max_thrust,
            doubleType & max_mass_flow_rate,
            doubleType & Isp,
            doubleType & power,
            doubleType & active_power,
            int & number_of_active_engines,
            EMTG::math::Matrix<doubleType> & fx,
            std::vector <doubleType> & acceleration_vector,
            const bool & generate_derivatives,
            const int& j,
            const int& p,
            const int& step);

    }
} //close namespace

#endif