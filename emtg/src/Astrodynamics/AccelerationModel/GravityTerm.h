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

#ifndef GRAVITY_TERM_H
#define GRAVITY_TERM_H

#include "SpacecraftAccelerationModel.h"
#include "SpacecraftAccelerationModelTerm.h"
#include "body.h"
#include "doubleType.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        class GravityTerm : public SpacecraftAccelerationModelTerm
        {
            friend class SphericalHarmonicTerm;

        public:
            //constructor
            GravityTerm(SpacecraftAccelerationModel * acceleration_model_in, body * body_in);

            //methods
            inline body * getBody() const { return this->my_body;  }

            virtual void computeAccelerationTerm() override;
            virtual void computeAccelerationTerm(const bool & generate_derivatives) override;
            virtual void populateInstrumentationFile(std::ofstream & acceleration_model_file) override;
            GravityTerm* clone() const override { return new GravityTerm(*this); }
            virtual void computePointMassGravityTimeDerivatives();
            
        protected:
            void computeScBodyCBtriangle(const bool & generate_derivatives);
            void computePointMassGravityAcceleration(const bool & generate_derivatives);
            void computePointMassGravityDerivatives();
            virtual void computeFrameDragAcceleration();
            
            // perturbing gravitating body
            body * my_body;

            // position vector of the spacecraft relative to the gravitating body
            math::Matrix<doubleType> r_body2sc;
            doubleType r_body2sc_norm;

            // velocity vector of the spacecraft relative to the gravitating body
            math::Matrix<doubleType> v_body2sc;
            doubleType v_body2sc_norm;

            // position vector of the perturbing gravitating body relative to the central body
            math::Matrix<doubleType> r_cb2body;
            doubleType r_cb2body_norm;

            // velocity vector of the perturbing gravitating body relative to the central body
            math::Matrix<doubleType> v_cb2body;
            doubleType v_cb2body_norm;

            // acceleration on the spacecraft due to the gravitating body
            doubleType a_body_on_sc_norm;

            // acceleration on the central body due to the gravitating body
            math::Matrix<doubleType> a_body_on_cb;
            doubleType a_body_on_cb_norm;

            // partial derivative of the acceleration vector with respect to the
            // spacecraft's position vector
            math::Matrix<double> da_cb2scdr_cb2body;

            // explicit partial of the gravitating body position vector relative to the central body
            // w.r.t. the current epoch
            math::Matrix<double> dr_cb2bodydcurrent_epoch;

        };
    }
}

#endif