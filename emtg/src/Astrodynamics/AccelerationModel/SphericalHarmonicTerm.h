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

#ifndef SPHERICAL_HARMONIC_TERM_H
#define SPHERICAL_HARMONIC_TERM_H

#include "SpacecraftAccelerationModelTerm.h"
#include "body.h"
#include "SpacecraftAccelerationModel.h"


namespace EMTG
{
    namespace Astrodynamics
    {
        class CentralBodyGravityTerm;
        class SphericalHarmonicTerm : public SpacecraftAccelerationModelTerm
        {
        public:
            // constructor
            SphericalHarmonicTerm(SpacecraftAccelerationModel * acceleration_model_in, body * my_body_in, CentralBodyGravityTerm * parent_gravity_term_in, const size_t & degree_in, const size_t & order_in);

            virtual ~SphericalHarmonicTerm();

            // methods
            virtual void computeAccelerationTerm() override;
            virtual void computeAccelerationTerm(const bool & generate_derivatives) override;
            virtual void populateInstrumentationFile(std::ofstream & acceleration_model_file) override;
            inline SphericalHarmonicTerm* clone() const override { return new SphericalHarmonicTerm(*this); }


        protected:
            // fields
            body * my_body;
            CentralBodyGravityTerm * parent_gravity_term;
            size_t degree;
            size_t order;

            math::Matrix<doubleType> r_body2sc_BCF;
            doubleType r_body2sc_BCF_norm;
            math::Matrix<doubleType> spherical_harmonic_term_acceleration_BCF;

        };
    } // end Astrodynamics namespace
} // end EMTG namespace
#endif // SPHERICAL_HARMONIC_TERM_H