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

#include "RungeKuttaTableau.h"

namespace EMTG {
	namespace Integration
	{
		// standard constructor
		RungeKuttaTableau::RungeKuttaTableau()
		{
		}

		// destructor
		RungeKuttaTableau::~RungeKuttaTableau()
		{
		}

		// methods
		void RungeKuttaTableau::resizeArrays()
		{
			this->A.resize(this->getNumStages(), this->getNumStages());
			this->bUpper.resize(this->getNumStages(), 1);
			this->bLower.resize(this->getNumStages(), 1);
			this->c.resize(this->getNumStages(), 1);
		}

		void RungeKuttaTableau::writeTableauToFile(const std::string & fileName)
		{
			std::ofstream f;
			f.open(fileName);

            size_t num_stages = this->getNumStages();

			// A 
			for (size_t i = 0; i < num_stages; ++i)
			{
				for (size_t j = 0; j < num_stages; ++j)
				{
					f << std::fixed << std::setprecision(16) << "A(" << i << ", " << j << ") = " << this->getA()(i, j) << "\n";
				}
			}

			// bUpper
			for (size_t i = 0; i < num_stages; ++i)
			{
				f << "bUpper(" << i << ") = " << this->getBUpper()(i, 0) << "\n";
			}

			// bLower
			if (this->getHasVariableStepCoefficients())
			{
				for (size_t i = 0; i < num_stages; ++i)
				{
					f << "bLower(" << i << ") = " << this->getBLower()(i, 0) << "\n";
				}
			}

			// c 
			for (size_t i = 0; i < num_stages; ++i)
			{
				f << "c(" << i << ") = " << this->getC()(i, 0) << "\n";
			}

			f.close();
		}
	} // end integration namespace
} // end EMTG namespace