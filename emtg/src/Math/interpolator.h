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

//EMTG interpolation class header file
//Jacob Englander 6/14/2013

#pragma once

#include <vector>
#include <utility>

#include "EMTG_math.h"

namespace EMTG
{ 
    namespace math
    {

        class interpolator
        {
        public:
            //default constructor
            interpolator(void);

            //constructor for a known data table
            interpolator (std::vector < std::pair<double, double> > InputTable);

            //destructor
            virtual ~interpolator(void);

            //fields
            std::vector< std::pair<double, double> > DataTable;

            //methods
            double interpolate(double x);
        };

    }
}//close namespace EMTG::math