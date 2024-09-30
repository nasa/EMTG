// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2024 United States Government as represented by the
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

//header file containing EMTG macro definitions
//Noble Hatten 08/02/2021

#pragma once

/*
STM_TEST_TURN_OFF_SCALING:
Define this if executing the acceleration_model_tests STM_test executable.
Otherwise, comment out this definition.
What it does: Changes how state derivatives are calculated in universe.cpp
(around line 425) and body.cpp (around line 273).
STM_test does not want derivatives scaled, but the rest of EMTG does.
Turned off by default.
*/
//#define STM_TEST_TURN_OFF_SCALING