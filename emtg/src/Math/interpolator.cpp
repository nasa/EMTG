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

#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>

#include "interpolator.h"
#include "EMTG_math.h"

namespace EMTG { namespace math {

    //default constructor
    interpolator::interpolator(void)
    {
    }

    //constructor for a known data table
    interpolator::interpolator(std::vector < std::pair<double, double> > InputTable)
    {
        DataTable = InputTable;
    }

    //destructor
    interpolator::~interpolator(void)
    {
    }

    //function to interpolate the data table
    double interpolator::interpolate(double x)
    {
        // Check if x is out of bound
        if ( x > DataTable.back().first )
            return -LARGE;
        if ( x < DataTable.front().first )
            return LARGE;

        std::vector<std::pair<double, double> >::iterator it, it2;

        it = lower_bound(DataTable.begin(), DataTable.end(), std::make_pair(x, -LARGE));

        // Corner case

        if (it == DataTable.begin()) 
            return it->second;

        it2 = it;
        --it2;

        return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
    }

}} //close namespace EMTG::math