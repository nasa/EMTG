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

//base class for acceleration model terms

#ifndef SPACECRAFT_ACCELERATION_MODEL_TERM_H
#define SPACECRAFT_ACCELERATION_MODEL_TERM_H

#include <iostream>
#include <string>
#include <vector>

#include "SpacecraftAccelerationModel.h"
#include "doubleType.h"
#include "EMTG_Matrix.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        class SpacecraftAccelerationModel;
        class SpacecraftAccelerationModelTerm
        {
        public:
            // constructor           
            SpacecraftAccelerationModelTerm(SpacecraftAccelerationModel * parent_acceleration_model);

            virtual ~SpacecraftAccelerationModelTerm();

            // methods
            virtual void computeAccelerationTerm() = 0;
            virtual void computeAccelerationTerm(const bool & generate_derivatives) = 0;
            virtual SpacecraftAccelerationModelTerm* clone() const = 0;
            virtual void populateInstrumentationFile(std::ofstream & acceleration_model_file) = 0;

            inline math::Matrix<doubleType> getTermAcceleration() const { return this->term_acceleration; }

            //fields
        protected:

            // acceleration vector contribution of this term
            math::Matrix<doubleType> term_acceleration;

            // pointer to the acceleration model that owns this term
            SpacecraftAccelerationModel * acceleration_model;

        };

        inline SpacecraftAccelerationModelTerm * new_clone(SpacecraftAccelerationModelTerm const & other)
        {
            return other.clone();
        }

    } // end Astrodynamics namespace
} // end EMTG namespace

#endif