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

//central force term

#ifndef SOLAR_RADIATION_PRESSURE_TERM_H
#define SOLAR_RADIATION_PRESSURE_TERM_H

#include "SpacecraftAccelerationModel.h"
#include "SpacecraftAccelerationModelTerm.h"
#include "doubleType.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        class SolarRadiationPressureTerm : public SpacecraftAccelerationModelTerm
        {
        public:
            //constructor
            SolarRadiationPressureTerm(SpacecraftAccelerationModel * acceleration_model_in);

            virtual ~SolarRadiationPressureTerm();

            //methods
            virtual void computeAccelerationTerm() override;
            virtual void computeAccelerationTerm(const bool & generate_derivatives) override;
            void populateInstrumentationFile(std::ofstream & acceleration_model_file) override;
            inline SolarRadiationPressureTerm* clone() const override { return new SolarRadiationPressureTerm(*this); }

            //fields
            double Cr; // coefficient of reflectivity [0,2]
            double sc_area; // spacecraft area exposed to sunlight km^2
            double K; // percentage of Sun seen by spacecraft [0,1]
            double solar_flux;  // W/m^2 = kg/s^3 --> equal to 3.846e26 Watts / (4 * pi * r_km^2)
            double speed_of_light_vac;  // speed of light in a vacuum km/s
            double SRP_coeff;
            doubleType SRP_force_norm;

        };
    } // end Astrodynamics namespace
} // end EMTG namespace

#endif