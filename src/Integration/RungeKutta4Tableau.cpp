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

#include "RungeKutta4Tableau.h"

namespace EMTG {
	namespace Integration
	{
		// standard constructor
		RungeKutta4Tableau::RungeKutta4Tableau()
		{
			this->setHasVariableStepCoefficients(false); // fixed-step method
			this->num_stages = 4;
			this->resizeArrays(); // set sizes of coefficient arrays

			// populate A
			this->A.assign_zeros(); // mostly zeros
			this->A(1, 0) = 1.0 / 2.0;
			this->A(2, 1) = 1.0 / 2.0;
			this->A(3, 2) = 1.0;

			// populate bUpper
			this->bUpper(0) = 1.0 / 6.0;
			this->bUpper(1) = 1.0 / 3.0;
			this->bUpper(2) = 1.0 / 3.0;
			this->bUpper(3) = 1.0 / 6.0;

			// populate c
			this->c(0, 0) = 0.0;
			this->c(1, 0) = 1.0 / 2.0;
			this->c(2, 0) = 1.0 / 2.0;
			this->c(3, 0) = 1.0;
			
		}
	} // end integration namespace
} // end EMTG namespace