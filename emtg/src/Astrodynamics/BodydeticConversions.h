// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2014 - 2016 United States Government as represented by the
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

// header file for bodydetic conversions
// Noble Hatten 11/19/2019

#ifndef BODYDETIC_CONVERSIONS_H
#define BODYDETIC_CONVERSIONS_H

#include "EMTG_math.h"
#include "EMTG_Matrix.h"
#include "doubleType.h"

#include<vector>

namespace EMTG
{
	namespace Astrodynamics
	{
        void LLA2BCF_oblate(const math::Matrix<doubleType>& LLA,
            const doubleType& Re,
            const doubleType& f,
            math::Matrix<doubleType>& rbcf,
            const bool& generateDerivatives,
            math::Matrix<doubleType>& dBCFdLLA);
		
        void BCF2LLA_oblate(const math::Matrix<doubleType>& rbcf,
            const doubleType& Re,
            const doubleType& f,
            math::Matrix<doubleType>& LLA,
            const bool& generateDerivatives,
            math::Matrix<doubleType>& dLLAdBCF);

	}//end namespace Astrodynamics
}//end namespace EMTG

#endif // BODYDETIC_CONVERSIONS_H