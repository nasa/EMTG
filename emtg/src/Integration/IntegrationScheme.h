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


#ifndef INTEGRATION_SCHEME_H
#define INTEGRATION_SCHEME_H

#include "doubleType.h"
#include "EMTG_Matrix.h"
#include "Integrand.h"

namespace EMTG {
    namespace Integration {

        class IntegrationScheme
        {
        public:
            IntegrationScheme();
            IntegrationScheme(Integrand * integrand_in);
            virtual ~IntegrationScheme() = default;

            inline void setLeftHandIndependentVariablePtr(doubleType & left_hand_independent_variable_in) 
            { 
                this->left_hand_independent_variable = &left_hand_independent_variable_in; 
            }
            inline void setLeftHandEpochPtr(doubleType & left_hand_epoch_in)
            {
                this->left_hand_epoch = &left_hand_epoch_in;
            }
            inline void setdLeftHandIndVardPropVarPtr(double & dleft_hand_indvar_dProp_var_in) 
            { 
                this->dleft_hand_ind_var_dProp_var = &dleft_hand_indvar_dProp_var_in; 
            }
            inline void setdLeftHandIndVardPropVarPreviousPtr(double & dleft_hand_indvar_dProp_var_previous_in)
            {
                this->dleft_hand_ind_var_dProp_var_previous = &dleft_hand_indvar_dProp_var_previous_in;
            }

            inline void setIntegrand(Integrand * integrand_in) { this->integrand = integrand_in; };

            // set the number of states to integrate over, overriding the constructor parmeter
            inline void setNumStatesToIntegratePtr(size_t & num_states_in) { this->num_states_to_integrate = &num_states_in; };

            virtual void step(const math::Matrix<doubleType> & state_left,
                              const math::Matrix<double> & dstate_leftdProp_vars,
                              math::Matrix<doubleType> & state_right,
                              math::Matrix<double> & dstate_rightdProp_vars,
                              const doubleType & step_size,
                              const double & dstep_sizedProp_var,
                              const bool & needSTM) = 0;

            virtual void step_w_control(const math::Matrix<doubleType> & state_left,
                                        const math::Matrix<double> & dstate_leftdProp_vars,
                                        math::Matrix<doubleType> & state_right,
                                        math::Matrix<double> & dstate_rightdProp_vars,
                                        const math::Matrix <doubleType> & control,
                                        const doubleType & step_size,
                                        const double & dstep_sizedProp_var,
                                        const bool & needSTM) = 0;

            virtual void errorControlledStep(const math::Matrix<doubleType> & state_left,
                                             const math::Matrix<double> & dstate_leftdProp_vars,
                                             math::Matrix<doubleType> & state_right,
                                             math::Matrix<double> & dstate_rightdProp_vars,
                                             const math::Matrix <doubleType> & control,
                                             const doubleType & step_size,
                                             const double & dstep_sizedProp_var,
                                             const bool & needSTM,
                                             doubleType & error,
                                             math::Matrix<double> & error_scaling_factors) {};

            virtual void computeError(doubleType & error, math::Matrix<double> & error_scaling_factors) {};

        protected:
            size_t * num_states_to_integrate;
            Integrand * integrand;
            doubleType * left_hand_independent_variable;
            doubleType * left_hand_epoch;

            // Partial of the current propagation step size w.r.t. the current propagation variable (flight times or total angle)
            double dstep_sizedProp_var;

            // Partial of the step left-hand independent variable w.r.t. previous propagation variables (flight times)
            double * dleft_hand_ind_var_dProp_var_previous;

            // Partial of the step left-hand independent variable w.r.t. the current propagation variable (flight times or total angle)
            double * dleft_hand_ind_var_dProp_var;
        };

    } // end namespace Integration
} // end namespace EMTG

#endif // INTEGRATION_SCHEME_H