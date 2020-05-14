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

//header file for EMTG atmosphere class

#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include "doubleType.h"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "missionoptions.h"
#include "EMTG_math.h"

#include "boost/ptr_container/ptr_vector.hpp"
#include "boost/algorithm/string.hpp"

#ifdef AD_INSTRUMENTATION
#include <set>
#endif

namespace EMTG
{
	namespace Astrodynamics
	{
		class atmosphere
		{
		public:
			//constructors
			atmosphere(); //default constructor
			atmosphere(const size_t& j,
				const std::string & atmospherefile,
				const missionoptions& options,
				const doubleType& alpha);

			//destructor
			virtual ~atmosphere();

			//**************************************
			//methods

			// get/set alpha parameter (km)
			void setAlpha(const doubleType& alpha);
			doubleType getAlpha(void);

			// get density
			doubleType getDensity(void);

			// get density derivative for exponential atmosphere
			virtual doubleType getDdensityDh(void) = 0;

			// function to get density
			virtual void computeDensity(const doubleType& h, const bool & generate_derivatives) = 0; // exponential atmosphere

			//**************************************
			//fields

			//the following fields are read in
			EMTG::AtmosphereModel model_type; // possibilities: "Exponential"

			//the following fields are computed internal to the class

		protected:
			//function to load new data into the atmosphere
			virtual void load_atmosphere_data(const size_t& j, const std::string & atmospherefile, const missionoptions& options) = 0;

			// method to convert densities at altitudes to scale heights
			virtual void computeScaleHeights(void) = 0;

			// method to calculate what altitude region the spacecraft is in
			virtual size_t computeAltitudeRegion(const doubleType& h) = 0;

			virtual void computeWeighting(const doubleType xbar[2], const doubleType& h,
				const doubleType scaleHeights[2],
				doubleType& scaleHeightPrime,
				doubleType& DscaleHeightPrimeDh) = 0; // compute weighting parameters necessary for continuity

			//**************************************
			//fields
			doubleType alpha; // km; parameter that controls size of smoothing region for exponential atmosphere
			doubleType density; // density in kg/m^3
			doubleType DdensityDh; // d [density] / d [altitude] in [kg/m^3] / km

			// used for exponential atmosphere
			std::vector<doubleType> altitudes;
			std::vector<doubleType> densities;
			std::vector<doubleType> scaleHeights;
			size_t nDataPoints;
			size_t altitudeRegion;

		};//end class atmosphere

	}//close namespace Astrodynamics
}//close namespace EMTG

#endif // ATMOSPHERE_H
