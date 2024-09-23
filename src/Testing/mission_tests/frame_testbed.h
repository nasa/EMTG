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

//event testbed

#pragma once

#include "missionoptions.h"
#include "universe.h"
#include "frame.h"

#include <random>
#include <sstream>


void frameTestbed(EMTG::missionoptions& options,
    std::vector< EMTG::Astrodynamics::universe > TheUniverse,
    std::mt19937 RNG,
    std::uniform_real_distribution<> UniformDouble);

void testRotation(EMTG::ReferenceFrame InputFrame,
                  EMTG::ReferenceFrame OutputFrame,
                  std::ofstream& outputfile,
                  const double& alpha0,
                  const double& alphadot,
                  const double& delta0,
                  const double& deltadot,
                  const double& W0,
                  const double& Wdot,
                  const double& theta1_0,
                  const double& theta1dot_0,
                  const double& theta2_0,
                  const double& theta2dot_0,
                  const double& theta3_0,
                  const double& theta3dot_0,
                  const double& semi_axis_a,
                  const double& semi_axis_b,
                  const double& semi_axis_c,
                  const doubleType& ETepoch,
                  EMTG::math::Matrix<doubleType> inputVector,
                  EMTG::math::Matrix<doubleType> referenceVector);