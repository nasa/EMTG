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

// gravity term

#ifndef CENTRAL_BODY_GRAVITY_TERM_H
#define CENTRAL_BODY_GRAVITY_TERM_H

#include "boost/ptr_container/ptr_vector.hpp"

#include "SpacecraftAccelerationModel.h"
#include "SpacecraftAccelerationModelTerm.h"
#include "body.h"
#include "CentralBody.h"
#include "GravityTerm.h"
#include "SphericalHarmonicTerm.h"
#include "doubleType.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        class CentralBodyGravityTerm : public GravityTerm
        {
            friend class SphericalHarmonicTerm;

        public:
            // constructor
            CentralBodyGravityTerm(SpacecraftAccelerationModel * acceleration_model_in, CentralBody * central_body_in);
            ~CentralBodyGravityTerm();

            // methods
            inline CentralBodyGravityTerm* clone() const override { return new CentralBodyGravityTerm(*this); }
            virtual void computeFrameDragAcceleration() override;
            virtual void computeAccelerationTerm() override;
            virtual void computeAccelerationTerm(const bool & generate_derivatives) override;
            virtual void computePointMassGravityTimeDerivatives() override;
            virtual void populateInstrumentationFile(std::ofstream & acceleration_model_file) override;
            
        protected:
            // fields
            boost::ptr_vector< SphericalHarmonicTerm > gravitational_harmonic_terms;

        };
    }
}

#endif // CENTRAL_BODY_GRAVITY_TERM_H