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


#ifndef INTEGRAND_H
#define INTEGRAND_H

#include "doubleType.h"
#include "EMTG_Matrix.h"

namespace EMTG {
    namespace Integration {

        class Integrand
        {
        public:
            Integrand();
            virtual ~Integrand();

            virtual void evaluate(const math::Matrix<doubleType> & state,
                                  math::Matrix<doubleType> & state_dot,
                                  const bool & generate_derivatives) = 0;

            virtual void evaluate(const math::Matrix<doubleType> & state,
                                  math::Matrix<doubleType> & state_dot,
                                  const math::Matrix<doubleType> & control,
                                  const bool & generate_derivatives) { throw std::runtime_error("The control overload for Integrand::evaluate has not been implemented!!"); };

            inline void setCurrentIndependentVariable(const doubleType & current_independent_variable_in) 
            { 
                this->current_independent_variable = current_independent_variable_in; 
            }

            inline EMTG::math::Matrix<double> getStatePropMat() const { return this->state_propagation_matrix; }

        protected:
            doubleType current_independent_variable;
            math::Matrix<double> state_propagation_matrix;
        };

    } // end namespace Integration
} // end namespace EMTG

#endif