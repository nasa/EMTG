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

//differential equations of motion for EMTGv9
//Donald Ellison August 1st 2015

#include <vector>

#include "missionoptions.h"
#include "EMTG_Matrix.h"
#include "doubleType.h"
#include "Spacecraft.h"
#include "universe.h"

#ifndef FBLTEOM
#define FBLTEOM

namespace EMTG {
    namespace Astrodynamics {
        namespace EOM {
            //cartesian coordinate equations of motion for a spacecraft with a thrust term

            void EOM_inertial_continuous_thrust(std::vector <doubleType> & spacecraft_state,
                EMTG::math::Matrix <double> & dxdTOF,
                const doubleType & epoch_step_left,
                std::vector <double> & depoch_step_leftdTOF,
                const double & c,
                doubleType & h,
                const double & dhdTOF,
                const doubleType & launch_epoch,
                const std::vector <doubleType> & u,
                const doubleType& u_command,
                const double& DutyCycle,
                std::vector <doubleType> & f, // EOM gradient vector
                EMTG::math::Matrix <double> & dfdTOF,
                doubleType & thrust,
                doubleType & mdot,
                doubleType & Isp,
                doubleType & power,
                doubleType & active_power,
                int & number_of_active_engines,
                int & STMrows,
                int & STMcolumns,
                const missionoptions& options,
                EMTG::Astrodynamics::universe& Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                const bool& GenerateDerivatives,
                const int& j,
                const int& p,
                const int& step,
                const double& mu,
                const double& LU,
                const double& TU,
                const double& MU);

        } //end EOM namespace
    } //end Astrodynamics namespace
} //end EMTG namespace

#endif