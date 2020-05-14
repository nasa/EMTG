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

//header file for EMTG ExponentialAtmosphere class

#ifndef EXPONENTIAL_ATMOSPHERE_H
#define EXPONENTIAL_ATMOSPHERE_H

#include "atmosphere.h"
#include "doubleType.h"

#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

#include "missionoptions.h"
#include "frame.h"

#include "EMTG_math.h"
#include "orbit_element_conversions.h"

#include "SpiceUsr.h"

#ifdef SPLINE_EPHEM
#include "SplineEphem_universe.h"
#endif


namespace EMTG
{
	namespace Astrodynamics
	{
		class ExponentialAtmosphere : public atmosphere
		{
		public:
			//constructor
			ExponentialAtmosphere() {}; //default constructor

			ExponentialAtmosphere(const size_t & j,
				const std::string & atmosphere_file,
				const missionoptions & options,
				const doubleType& alpha = 0.5);

			//destructor
			virtual ~ExponentialAtmosphere();


			//**************************************
			//methods
			void load_atmosphere_data(const size_t& j, const std::string & atmospherefile, const missionoptions& options);
			void computeScaleHeights(void); // method to convert densities at altitudes to scale heights
			void computeDensity(const doubleType& h, const bool & generate_derivatives); // exponential atmosphere
			size_t computeAltitudeRegion(const doubleType& h); // method to calculate what altitude region the spacecraft is in
			void computeWeighting(const doubleType xbar[2], const doubleType& h,
				const doubleType scaleHeights[2],
				doubleType& scaleHeightPrime,
				doubleType& DscaleHeightPrimeDh); // compute weighting parameters necessary for continuity
			doubleType getDdensityDh(void); // get density derivative
			//**************************************
			//members
			
		private:
			//**************************************
			//methods
			
			
		};//end class ExponentialAtmosphere

	}//close namespace Astrodynamics
}//close namespace EMTG

#endif // EXPONENTIAL_ATMOSPHERE_H