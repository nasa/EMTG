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


#ifndef TIME_DOMAIN_SPACECRAFT_EOM_H
#define TIME_DOMAIN_SPACECRAFT_EOM_H

#include "AccelerationModel.h"
#include "EMTG_Matrix.h"
#include "Integrand.h"
#include "SpacecraftAccelerationModel.h"

namespace EMTG {
    namespace Astrodynamics {

        class TimeDomainSpacecraftEOM : public Integration::Integrand
        {

        public:

            // constructor
            TimeDomainSpacecraftEOM();

            // destructor
            ~TimeDomainSpacecraftEOM();

            virtual void evaluate(const math::Matrix<doubleType> & state,
                                  math::Matrix<doubleType> & state_dot,
                                  const bool & generate_derivatives);

            virtual void evaluate(const math::Matrix<doubleType> & state,
                                  math::Matrix<doubleType> & state_dot,
                                  const math::Matrix<doubleType> & control,
                                  const bool & generate_derivatives);

            inline void setSpacecraftAccelerationModel(SpacecraftAccelerationModel * spacecraft_acceleration_model_in) 
            { 
                this->spacecraft_acceleration_model = spacecraft_acceleration_model_in; 
                STM.resize(this->spacecraft_acceleration_model->getSTMsize(), this->spacecraft_acceleration_model->getSTMsize(), 0.0);
                STM_dot.resize(this->spacecraft_acceleration_model->getSTMsize(), this->spacecraft_acceleration_model->getSTMsize(), 0.0);
            }

		protected:
			SpacecraftAccelerationModel * spacecraft_acceleration_model;
			math::Matrix<doubleType> acceleration;
			math::Matrix<doubleType> STM;
			math::Matrix<doubleType> STM_dot;
            doubleType current_epoch;

            void configureAccelerationModel(const math::Matrix<doubleType> & state);

			void computeAcceleration(const math::Matrix<doubleType> & state,
				                     const bool & STM_needed);
            void computeAcceleration(const math::Matrix<doubleType> & state,
                                     const math::Matrix<doubleType> & control,
				                     const bool & STM_needed);
			void ballisticEOM(const math::Matrix<doubleType> & state, math::Matrix<doubleType> & state_dot);
			void propulsionEOM(math::Matrix<doubleType> & state_dot, const math::Matrix<doubleType> & control);
			
        };
    } // end Astrodynamics namespace
} // end EMTG namespace

#endif // TIME_DOMAIN_SPACECRAFT_EOM_H