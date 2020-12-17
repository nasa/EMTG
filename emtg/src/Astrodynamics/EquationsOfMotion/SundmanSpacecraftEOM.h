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


#ifndef SUNDMAN_SPACECRAFT_EOM_H
#define SUNDMAN_SPACECRAFT_EOM_H

#include "AccelerationModel.h"
#include "EMTG_Matrix.h"
#include "Integrand.h"
#include "SpacecraftAccelerationModel.h"
#include "TimeDomainSpacecraftEOM.h"

namespace EMTG {
    namespace Astrodynamics {

        class SundmanSpacecraftEOM : public TimeDomainSpacecraftEOM
        {

        public:

            // constructor
			SundmanSpacecraftEOM();

            // destructor
            ~SundmanSpacecraftEOM();

            virtual void evaluate(const math::Matrix<doubleType> & state,
                                  math::Matrix<doubleType> & state_dot,
                                  const bool & generate_derivatives);

            virtual void evaluate(const math::Matrix<doubleType> & state,
                                  math::Matrix<doubleType> & state_dot,
                                  const math::Matrix<doubleType> & control,
                                  const bool & generate_derivatives);
            
			inline void setSundmanC(const double & c_in) { this->sundman_constant = c_in; }
			inline void setSundmanGamma(const double & gamma_in) { this->sundman_gamma = gamma_in; }

        private:

			doubleType sundman_constant;
			double sundman_gamma;
            doubleType sundman_transform;

            math::Matrix<doubleType> specific_angular_momentum;
            math::Matrix<doubleType> velocity;
            doubleType mass_flow_rate;
            math::Matrix<double> dsundman_constantdState;

            void configureAccelerationModel(const math::Matrix<doubleType> & state);

            void computeAcceleration(const math::Matrix<doubleType> & state,
                const bool & STM_needed);
            void computeAcceleration(const math::Matrix<doubleType> & state,
                const math::Matrix<doubleType> & control,
                const bool & STM_needed);

            void computeSundmanStatePropagationMatrix();
            void computeSundmanConstantTA(const math::Matrix<doubleType> & state);
            void computeSundmanConstantGA(const math::Matrix<doubleType> & state);
            void computeSundmanConstantScaled(const math::Matrix<doubleType> & state);
        };
    } // end Astrodynamics namespace
} // end EMTG namespace

#endif // SUNDMAN_SPACECRAFT_EOM_H