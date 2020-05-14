
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

//file to define doubleType
//depending on a preprocessor switch, this redefines doubleType
//no more templates!

#pragma once

#ifdef AD_INSTRUMENTATION

//needed for GSAD A4
//#include "GSAD_A5.h"
#include "GSAD_2B.h"

#define doubleType GSAD::adouble
#define _GETVALUE .getValue()

#define _SETVALUE(x,y) .setValue(x,y)

#else

#define doubleType double
#define _GETVALUE 
#define _SETVALUE(x,y)

#endif