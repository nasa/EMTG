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

//event testbed class

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <vector>

#include "EMTG_Matrix.h"

#include "mathTestbed.h"

void mathTestbed(EMTG::missionoptions& options,
    std::vector< EMTG::Astrodynamics::universe > TheUniverse,
    EMTG::HardwareModels::Spacecraft& mySpacecraft,
    EMTG::HardwareModels::LaunchVehicle& myLaunchVehicle,
    std::mt19937 RNG,
    std::uniform_real_distribution<> UniformDouble)
{
    size_t n = 14;
    EMTG::math::Matrix<double> A(n, n, 0.0);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
            A(i, j) = UniformDouble(RNG);
    }

    A.print_to_file("tests/A.txt");

    EMTG::math::Matrix<double> Ainv = A.inverse();

    Ainv.print_to_file("tests/Ainv.txt");
}//end main