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

//source file for EMTG ExponentialAtmosphere class

#include "atmosphere.h"
#include "ExponentialAtmosphere.h"


namespace EMTG
{
	namespace Astrodynamics
	{

		// constructors
		//constructor to load a data file
		ExponentialAtmosphere::ExponentialAtmosphere(const size_t& j,
			const std::string & atmosphere_file,
			const missionoptions& options,
			const doubleType& alpha)
		{
			this->load_atmosphere_data(j, atmosphere_file, options);

			this->computeScaleHeights(); // convert densities at altitudes to scale heights
			this->setAlpha(alpha); // default value of alpha is 0.5
		}

		// destructor
		ExponentialAtmosphere::~ExponentialAtmosphere() {}

		// methods

		// load data from file
		void ExponentialAtmosphere::load_atmosphere_data(const size_t& j, const std::string & atmospherefile, const missionoptions& options)
		{
            std::string file_to_open = options.universe_folder + "/atmosphere_files/" + atmospherefile;
			std::ifstream inputfile(file_to_open.c_str());
			int linenumber = 0;
			std::string choice;
			std::string peek;
			char dump_buffer[1024];

			if (!inputfile.is_open())
			{
				throw std::invalid_argument("Failure to read " + file_to_open + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
			}

			while (!inputfile.eof())
			{
				peek = inputfile.peek();
				if (peek == "#" || peek == "\r" || peek == "\n") {
					//comment or blank line, do not parse
					inputfile.getline(dump_buffer, 1024);
					++linenumber;
				}
				else
				{
					inputfile >> choice;

					if (choice == "ModelType")
					{
                        std::string temp_type;
                        inputfile >> temp_type;
                            
                        if (temp_type == "Exponential")
                        {
                            this->model_type = EMTG::AtmosphereModel::Exponential;
                        }
                        else
                        {
                            throw std::invalid_argument("Atmosphere model type [" + temp_type + "] in " + atmospherefile + " is incompatible with AtmosphericDensityModelKey in .emtgopt file" + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                        }
					}
					else if (choice == "begin_data_list")
					{
						double tempAltitude;
						double tempDensity;
						std::string tempname;

						do
						{
							tempname.clear();

							//if not yet reached end of list:
							inputfile >> tempname;
							if (tempname == "end_data_list")
								break;
							tempAltitude = std::stod(tempname);
							inputfile >> tempDensity;

							// add to arrays: .push_back
							this->altitudes.push_back(tempAltitude);
							this->densities.push_back(tempDensity);
						} while (true);
						choice.clear();
					}

				}
			}

			this->nDataPoints = this->altitudes.size();
		}

		// method to convert densities at altitudes to scale heights
		void ExponentialAtmosphere::computeScaleHeights(void)
		{
			for (size_t i = 0; i < this->nDataPoints - 1; ++i)
			{
				this->scaleHeights.push_back((this->altitudes[i]-this->altitudes[i+1]) / log(this->densities[i+1] / this->densities[i]));
			}
			this->scaleHeights.push_back(this->scaleHeights.back()); // assume second-to-highest scale height is same as highest scale height

		}

		/* function to get density and derivative of density for atmosphere that depends only on altitude
		derivative is so simple, just always compute it regardless of value of generate_derivatives
		@param h altitude (km)
		*/
		void ExponentialAtmosphere::computeDensity(const doubleType& h, const bool & generate_derivatives)
		{
			// get altitude region
			this->altitudeRegion = this->computeAltitudeRegion(h);
			doubleType hi, scaleHeighti, densityi, hip1, scaleHeightip1, densityip1;
			doubleType scaleHeightim1, expTerm;
			doubleType scaleHeightPrimei, DscaleHeightPrimeDhi;
			doubleType xbar[2], scaleHeights[2];

			// get necessary coefficients
			if (this->altitudeRegion < this->nDataPoints - 1)
			{
				// not above nominal altitude max
				hi = this->altitudes[this->altitudeRegion]; //self.hVec[i]
				scaleHeighti = this->scaleHeights[this->altitudeRegion]; //self.scaleHeights[i]
				densityi = this->densities[this->altitudeRegion]; //self.densityVec[i]
				hip1 = this->altitudes[this->altitudeRegion + 1]; // self.hVec[i + 1]
				scaleHeightip1 = this->scaleHeights[this->altitudeRegion + 1]; // self.scaleHeights[i + 1]
				densityip1 = this->densities[this->altitudeRegion + 1]; //self.densityVec[i + 1]
			}
			else
			{
				// above nominal altitude max
				hi = this->altitudes[this->nDataPoints-2]; //self.hVec[-2]
				scaleHeighti = this->scaleHeights[this->nDataPoints - 2];  //self.scaleHeights[-2]
				densityi = this->densities[this->nDataPoints - 2]; //self.densityVec[-2]
				hip1 = this->altitudes[this->nDataPoints - 1]; //self.hVec[-1]
				scaleHeightip1 = this->scaleHeights[this->nDataPoints - 1]; //self.scaleHeights[-1]
				densityip1 = this->densities[this->nDataPoints - 1]; //self.densityVec[-1]
			}

			if (h < (hi + this->alpha) && this->altitudeRegion > 0)
			{
				// near lower boundary of layer and not in lowermost layer
				scaleHeightim1 = this->scaleHeights[this->altitudeRegion - 1]; //self.scaleHeights[i - 1]

				// boundaries of weighting region
				xbar[0] = hi - this->alpha;
				xbar[1] = hi + this->alpha;

				// get weighted scale heights
				scaleHeights[0] = scaleHeightim1;
				scaleHeights[1] = scaleHeighti;
				this->computeWeighting(xbar, h, scaleHeights, scaleHeightPrimei, DscaleHeightPrimeDhi);

				// get density
				expTerm = exp((hi - h) / scaleHeightPrimei);
				this->density = densityi * expTerm;

				// derivative
				this->DdensityDh = (-this->density / scaleHeightPrimei) * (1. + ((hi - h) / scaleHeightPrimei) * DscaleHeightPrimeDhi); // d[rho(h)] / d[h]
			}
			else if (h > (hip1 - this->alpha) && this->altitudeRegion < this->nDataPoints - 2)
			{
				// near upper boundary of layer and not in or above uppermost layer(might need to be - 2 instead of - 1)
				scaleHeightip1 = this->scaleHeights[this->altitudeRegion + 1];
				xbar[0] = hip1 - this->alpha;
				xbar[1] = hip1 + this->alpha;
				
				// get weighted scale heights
				scaleHeights[0] = scaleHeighti;
				scaleHeights[1] = scaleHeightip1;
				this->computeWeighting(xbar, h, scaleHeights, scaleHeightPrimei, DscaleHeightPrimeDhi);

				// get density
				expTerm = exp((hip1 - h) / scaleHeightPrimei);
				this->density = densityip1 * expTerm;

				// derivative
				this->DdensityDh = (-this->density / scaleHeightPrimei) * (1. + ((hip1 - h) / scaleHeightPrimei) * DscaleHeightPrimeDhi); // d[rho(h)] / d[h]
			}
			else
			{
				// not near boundary or near lowermost or uppermost boundary or above uppermost boundary, do not need to do weighting
				scaleHeightPrimei = scaleHeighti;
				expTerm = exp((hi - h) / scaleHeightPrimei);
				density = densityi * expTerm;

				// derivative
				this->DdensityDh = -this->density / scaleHeightPrimei; // d[rho(h)] / d[h]
			}
		}

		/*
		method to calculate what altitude region the spacecraft is in
		@param h altitude (km)
		*/
		size_t ExponentialAtmosphere::computeAltitudeRegion(const doubleType& h)
		{
			if (h < this->altitudes[0])
			{
				// altitude technically too low; use workaround where we use scale height for lowest region
				return 0;
			}
			else if (h >= this->altitudes[this->nDataPoints - 1])
			{
				// altitude technically too high; use workaround where we use scale height for highest region
				return this->nDataPoints - 1;
			}
			else
			{
				// simple loop through every altitude region
				for (size_t i = 0; i < this->nDataPoints - 1; ++i)
				{
					if (h >= this->altitudes[i] && h < this->altitudes[i+1])
						return i;
				}

                // this just gets rid of a "not all paths return a value" compiler warning and potential crash. This line should never be reached.
                return 0;
			}
		}

		/*
		compute weighting parameters necessary for continuity
		@param xbar 2-element vector: lower and upper boundaries of weighting region (altitudes, in km)
		@param h altitude (km)
		@param scaleHeights 2-element vector: scale heights for lower and upper regions
		@param scaleHeightPrime continuous scale height
		@param DscaleHeightPrimeDh d [scaleHeightPrime] / d [altitude]
		always: him1 < h < hip1
		if h is near him1, then scaleHeights[0] = scaleHeightsim1 and scaleHeights[1] = scaleHeightsi
		if h is near hip1, then scaleHeights[0] = scaleHeightsi and scaleHeights[1] = scaleHeightsip1
		*/
		void ExponentialAtmosphere::computeWeighting(const doubleType xbar[2], const doubleType& h,
			const doubleType scaleHeights[2],
			doubleType& scaleHeightPrime,
			doubleType& DscaleHeightPrimeDh)
		{
			doubleType xbarDiff = xbar[1] - xbar[0];
			doubleType x = (h - xbar[0]) / xbarDiff; // h scaled within weighting region
			doubleType x2 = x * x;
			doubleType x3 = x * x2;

			// third - order continuity:
			doubleType x4 = x * x3;
			doubleType t1 = -20. * x3 + 70. * x2 - 84. * x + 35.;
			doubleType w1 = x4 * t1;
			doubleType dw1 = 4. * x3 * t1 + x4 * (-60. * x2 + 140. * x - 84.);
			doubleType scaleHeightDiff = scaleHeights[1] - scaleHeights[0];
			scaleHeightPrime = scaleHeights[0] + w1 * scaleHeightDiff;
			DscaleHeightPrimeDh = (scaleHeightDiff / xbarDiff) * dw1;
		}

		/*
		getter for derivative of density with respect to altitude
		@return DdensityDh derivative of density with respect to altitude ([kg/m^3] / km)
		*/
		doubleType ExponentialAtmosphere::getDdensityDh(void)
		{
			return this->DdensityDh;
		}
	}//close namespace Astrodynamics
}//close namespace EMTG