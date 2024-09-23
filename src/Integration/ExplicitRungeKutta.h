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

// Templated version of the Dormand-Prince (DOPRI) 8th order 13 step algorithm
// DOPRI constants are from Numerical Recipies
// Donald Ellison 11/15/2016


#ifndef RUNGEKUTTA8_H
#define RUNGEKUTTA8_H

#include "doubleType.h"
#include "EMTG_enums.h"
#include "FBLT_EOM.h"
#include "Integrand.h"
#include "IntegrationScheme.h"
#include "RungeKuttaTableau.h"

namespace EMTG {
    namespace Integration
    {
        class ExplicitRungeKutta : public IntegrationScheme
        {

        public:

            // constructors
            ExplicitRungeKutta(Integrand * integrand_in,
                               const IntegrationCoefficientsType & RK_type,
                               const size_t & ns_in,
                               const size_t & num_prop_var_deriv_states);

            // destructor
            ~ExplicitRungeKutta();

            virtual void step(const math::Matrix<doubleType> & state_left,
                              const math::Matrix<double> & STM_left,
                              math::Matrix<doubleType> & state_right,
                              math::Matrix<double> & STM_right,
                              const doubleType & step_size,
                              const double& dstep_sizedProp_var,
                              const bool & needSTM);

            virtual void step(const math::Matrix<doubleType> & state_left,
                              const math::Matrix<double> & STM_left,
                              math::Matrix<doubleType> & state_right,
                              math::Matrix<double> & STM_right,
                              const math::Matrix <doubleType> & control,
                              const doubleType & step_size,
                              const double& dstep_sizedProp_var,
                              const bool & needSTM);

            virtual void errorControlledStep(const math::Matrix<doubleType> & state_left,
                                             const math::Matrix<double> & STM_left,
                                             math::Matrix<doubleType> & state_right,
                                             math::Matrix<double> & STM_right,
                                             const math::Matrix <doubleType> & control,
                                             const doubleType & step_size,
                                             const double & dstep_sizedProp_var,
                                             const bool & needSTM,
                                             doubleType & error,
                                             math::Matrix<double> & error_scaling_factors);

        protected:
            // methods

            void evaluateIntegrand(const math::Matrix<doubleType> & state_left,
                                   const size_t & stage_index,
                                   const doubleType & step_size,
                                   const bool & needSTM);

            void evaluateIntegrand(const math::Matrix<doubleType> & state_left,
                                   const math::Matrix <doubleType> & control,
                                   const size_t & stage_index,
                                   const doubleType & step_size,
                                   const bool & needSTM);

            void stateUpdate(const math::Matrix<doubleType> & state_left,
                             const math::Matrix<double> & coefficients,
                             const size_t & stage_index,
                             const doubleType & step_size);

            void stmUpdate(const math::Matrix<double> & STM_left,
                           const math::Matrix<double> & coefficients,
                           const size_t & stage_index,
                           const doubleType & step_size);

            void stageLoop(const math::Matrix<doubleType> & state_left,
                           math::Matrix<doubleType> & state_right,
                           const doubleType & step_size);

            void stageLoop(const math::Matrix<doubleType> & state_left,
                           math::Matrix<doubleType> & state_right,
                           const math::Matrix <doubleType> & control,
                           const doubleType & step_size);

            void stageLoop(const math::Matrix<doubleType> & state_left,
                           math::Matrix<doubleType> & state_right,
                           const math::Matrix<double> & STM_left,
                           math::Matrix<double> & STM_right,
                           const doubleType & step_size);

            void stageLoop(const math::Matrix<doubleType> & state_left,
                           math::Matrix<doubleType> & state_right,
                           const math::Matrix<double> & STM_left,
                           math::Matrix<double> & STM_right,
                           const math::Matrix <doubleType> & control,
                           const doubleType & step_size);

            virtual void computeError(doubleType & error, math::Matrix<double> & error_scaling_factors);

            // fields
        private:
            size_t num_stages;

            EMTG::math::Matrix <double> STM, STM_stage, fx, dstepdState, grad_vec;
            math::Matrix<doubleType> f, y, x_left, x_right;

            RungeKuttaTableau * RK_tableau;

            // commented-out coefficients are equal to zero
            // create local arrays
            math::Matrix<double> A;
            math::Matrix<double> bUpper;
            math::Matrix<double> bLower;
            math::Matrix<double> c;

            math::Matrix<doubleType> gradient_bin;
            std::vector< math::Matrix<double> > STM_bin;

        }; // end rk8 class definition

    } // end integration namespace
} // end EMTG namespace

#endif