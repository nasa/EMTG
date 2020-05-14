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

#ifndef INTEGRATED_PROPAGATOR_H
#define INTEGRATED_PROPAGATOR_H

#include "EMTG_Matrix.h"
#include "Integrand.h"
#include "IntegrationScheme.h"
#include "PropagatorBase.h"
#include "universe.h"

namespace EMTG {
    namespace Astrodynamics {

        class IntegratedPropagator : public PropagatorBase
        {
        public:
            // constructors
            IntegratedPropagator(universe & myUniverse,
                                 const size_t & STM_rows_in,
                                 const size_t & STM_columns_in,
                                 const size_t & STM_start_index_in);

            // methods
            virtual IntegratedPropagator* clone() const = 0;

            // set the integrand
            inline void setIntegrand(Integration::Integrand * integrand_in) { this->integrand = integrand_in; };

            inline void setdStepSizedPropVar(const double & dstep_sizedProp_var_in) { this->dstep_sizedProp_var = dstep_sizedProp_var_in; }
            inline void setdCurrentIndVardPropVar(const double & dcurrent_ind_vardProp_var_in) { this->dcurrent_ind_vardProp_var = dcurrent_ind_vardProp_var_in; }
            inline void setdCurrentIndVardPropVarPrevious(const double & dcurrent_ind_vardProp_var_previous_in)
            {
                this->dcurrent_ind_vardProp_var_previous = dcurrent_ind_vardProp_var_previous_in;
            }

            // set the integration_scheme and link the propagator and integration_scheme current_independent_variable pointers
            inline void setIntegrationScheme(Integration::IntegrationScheme * integration_scheme_in) 
            { 
                this->integration_scheme = integration_scheme_in; 
            };

            // set the actual number of states that you want the integration_scheme to integrate
            // as the integration_scheme's state vector is state+STM
            inline void setNumStatesToIntegrate(const size_t & num_states_in) 
            { 
                this->num_states_to_integrate = num_states_in; 
                this->integration_scheme->setNumStatesToIntegratePtr(this->num_states_to_integrate);
            }

            inline size_t getSTMrowDim() const { return this->num_STM_row; }
            inline size_t getSTMcolDim() const { return this->num_STM_col; }
            inline size_t getSTMstartIndex() const { return this->STM_start_index; }

            virtual void propagate(const doubleType & propagation_span, const bool & needSTM) = 0;
            virtual void propagate(const doubleType & propagation_span, const math::Matrix <doubleType> & control, const bool & needSTM) = 0;

        protected:

            void unpackStates(const math::Matrix<doubleType> & state_left, const bool & STM_needed);
            void packStates(math::Matrix<doubleType> & state_right, 
                            const math::Matrix<double> & dstate_rightdProp_vars,
                            math::Matrix<double> & STM, 
                            const bool & STM_needed);

            // fields

            Integration::Integrand * integrand;
            Integration::IntegrationScheme * integration_scheme;

            // actual number of states that you want the integrtion_scheme to integrate, overriding any of its defaults
            size_t num_states_to_integrate;

            // STM row dimension
            size_t num_STM_row;

            // STM column dimension
            size_t num_STM_col;

            // index in the state vector where STM entries begin
            size_t STM_start_index;

            // dummy control vector for when we call the propagator without control
            math::Matrix <doubleType> zero_control;

            // create augmented state vectors to feed to the integration_scheme
            math::Matrix<doubleType> state_left_augmented;
            math::Matrix<doubleType> state_right_augmented;

            // Partial of the current propagation step size w.r.t. the current propagation variable (flight times or total angle)
            double dstep_sizedProp_var;

        };

        inline IntegratedPropagator * new_clone(IntegratedPropagator const & other)
        {
            return other.clone();
        }

    } // end namespace Astrodynamics
} // end namespace EMTG

#endif