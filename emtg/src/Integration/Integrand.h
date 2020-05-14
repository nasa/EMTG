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
                                  const math::Matrix<double> & dstate_dProp_vars,
                                  math::Matrix<doubleType> & state_dot,
                                  math::Matrix<double> & dstate_dotdProp_vars,
                                  const bool & generate_derivatives) = 0;

            virtual void evaluate(const math::Matrix<doubleType> & state,
                                  const math::Matrix<double> & dstate_dProp_vars,
                                  math::Matrix<doubleType> & state_dot,
                                  math::Matrix<double> & dstate_dotdProp_vars,
                                  const math::Matrix<doubleType> & control,
                                  const bool & generate_derivatives) = 0;

            inline void setCurrentIndependentVariable(const doubleType & current_independent_variable_in) 
            { 
                this->current_independent_variable = current_independent_variable_in; 
            }
            inline void setdCurrentIndVardPropVar(const double & dcurrent_indvar_dProp_var_in) 
            { 
                this->dcurrent_indvar_dProp_var = dcurrent_indvar_dProp_var_in; 
            }
            inline void setdCurrentIndVardPropVarPrevious(const double & dcurrent_indvar_dProp_var_previous_in)
            {
                this->dcurrent_indvar_dProp_var_previous = dcurrent_indvar_dProp_var_previous_in;
            }

        protected:
            doubleType current_independent_variable;
            doubleType current_epoch;
            
            // Partial of the current independent variable w.r.t. previous propagation variables (flight times)
            double dcurrent_indvar_dProp_var_previous;

            // Partial of the current independent variable w.r.t. the current propagation variable (flight times or total angle)
            double dcurrent_indvar_dProp_var;
        };

    } // end namespace Integration
} // end namespace EMTG

#endif