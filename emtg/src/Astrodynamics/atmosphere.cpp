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

//source file for EMTG atmosphere class
#include "atmosphere.h"

namespace EMTG
{
	namespace Astrodynamics
	{
		atmosphere::atmosphere()
		{
		}

		//constructor to load a data file
		atmosphere::atmosphere(const size_t& j,
			const std::string & atmosphere_file,
			const missionoptions& options,
			const doubleType& alpha)
		{
		}

		//destructor
		atmosphere::~atmosphere() {}

		/* 
		getter for density
		@return density atmospheric density in kg/m^3
		*/
		doubleType atmosphere::getDensity(void)
		{
			return this->density;
		}

		// set alpha parameter (km)
		void atmosphere::setAlpha(const doubleType& alpha)
		{
			this->alpha = alpha;
		}

		// get alpha parameter (km)
		doubleType atmosphere::getAlpha(void)
		{
			return this->alpha;
		}
	}//close namespace Astrodynamics
}//close namespace EMTG