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

#include "IntegratedFixedStepPropagator.h"
#include "missionoptions.h"
#include "universe.h"

namespace EMTG {
    namespace Astrodynamics {

        // constructors

        IntegratedFixedStepPropagator::IntegratedFixedStepPropagator(universe & myUniverse,
                                                                     const size_t & STM_rows_in,
                                                                     const size_t & STM_columns_in,
                                                                     const size_t & STM_start_index_in) :
                                                                     IntegratedPropagator(myUniverse, 
                                                                                          STM_rows_in, 
                                                                                          STM_columns_in, 
                                                                                          STM_start_index_in)

        {

        }

        void IntegratedFixedStepPropagator::propagate_setup(const math::Matrix<doubleType> & state_left,
                                                            math::Matrix<double> & dstate_leftdProp_vars,
                                                            const bool & STM_needed)
        {
            //configure integration scheme pointers
            this->integration_scheme->setLeftHandIndependentVariablePtr(this->current_independent_variable);
            this->integration_scheme->setdLeftHandIndVardPropVarPtr(this->dcurrent_ind_vardProp_var);
            this->integration_scheme->setdLeftHandIndVardPropVarPreviousPtr(this->dcurrent_ind_vardProp_var_previous);

            // Since we are dealing with a fixed time step, the value of the left hand independent variable is not influenced by the current
            // propagation variable (TOF or total angular anomaly), so its partial is zero
            // The partial of the left hand independent variable w.r.t. previous phase propagation variables or the launch epoch is always 1.
            this->dcurrent_ind_vardProp_var_previous = 1.0;
            this->dcurrent_ind_vardProp_var = 0.0;

            // unpack the external state pointer into the integrator's augmented state
            this->unpackStates(state_left, STM_needed);

            dstate_leftdProp_vars(this->index_of_epoch_in_state_vec, 0) = 1.0;
        }

        void IntegratedFixedStepPropagator::propagate_teardown(const math::Matrix<doubleType> & state_left, 
                                                               math::Matrix<doubleType> & state_right,
                                                               const math::Matrix<double> & dstate_rightdProp_vars,
                                                               math::Matrix<double> & STM,
                                                               const doubleType & propagation_span,
                                                               const bool & STM_needed)
        {
            // pack the augmented state back into the external state pointer
            this->packStates(state_right, dstate_rightdProp_vars, STM, STM_needed);

            if ((state_right(7) < state_left(7) && propagation_span > 0.0)
                || (state_right(7) > state_left(7) && propagation_span < 0.0))
            {
                std::cout << "Oh noes!!...mismatch in left epoch vs. right epoch" << std::endl;
            }
        }

        // methods
        void IntegratedFixedStepPropagator::propagate(const doubleType & propagation_span, 
                                                      const bool & STM_needed)
        {
            math::Matrix<doubleType> & state_left = *this->StateLeftPointer;
            math::Matrix<doubleType> & state_right = *this->StateRightPointer;
            math::Matrix<double> & dstate_leftdProp_vars = *this->dStatedIndependentVariablePointer;
            math::Matrix<double> dstate_rightdProp_vars;
            math::Matrix<double> & STM = *this->STMpointer;

            this->propagate_setup(state_left, dstate_leftdProp_vars, STM_needed);

            //TODO: let user set a global integration step size for the mission AND allow the user to override it at the journey level
            doubleType propagation_span_remaining = propagation_span;
            doubleType integration_step_size;
            if (propagation_span_remaining > 0.0)
            {
                while (propagation_span_remaining > 0.0)
                {
                    if (propagation_span_remaining >= this->PropagationStepSize)
                    {
                        integration_step_size = this->PropagationStepSize;
                        this->dstep_sizedProp_var = 0.0;
                    }
                    else
                    {
                        integration_step_size = propagation_span_remaining;
                        this->dstep_sizedProp_var = *this->boundary_target_dstep_sizedProp_var;
                    }

                    this->integration_scheme->step(this->state_left_augmented,
                                                   dstate_leftdProp_vars,
                                                   this->state_right_augmented,
                                                   dstate_rightdProp_vars,
                                                   integration_step_size,
                                                   this->dstep_sizedProp_var,
                                                   STM_needed);

                    this->state_left_augmented = this->state_right_augmented;
                    dstate_leftdProp_vars = dstate_rightdProp_vars;
                    propagation_span_remaining -= integration_step_size;
                    this->current_independent_variable += integration_step_size;

                    try
                    {
                        this->current_epoch = this->state_left_augmented(this->index_of_epoch_in_state_vec);
                    }
                    catch (const std::out_of_range & e)
                    {
                        std::cout << "You did not set the index where epoch is stored in the state vector!! " << e.what() << std::endl;
                    }

                } // end propagation loop
            }//end forward propagation
            else //backward propagation
            {
                while (propagation_span_remaining < 0.0)
                {
                    if (-propagation_span_remaining >= this->PropagationStepSize)
                    {
                        integration_step_size = -this->PropagationStepSize;
                        this->dstep_sizedProp_var = 0.0;
                    }
                    else
                    {
                        integration_step_size = propagation_span_remaining;
                        this->dstep_sizedProp_var = -*this->boundary_target_dstep_sizedProp_var;
                    }

                    this->integration_scheme->step(this->state_left_augmented,
                                                   dstate_leftdProp_vars,
                                                   this->state_right_augmented,
                                                   dstate_rightdProp_vars,
                                                   integration_step_size,
                                                   this->dstep_sizedProp_var,
                                                   STM_needed);

                    this->state_left_augmented = this->state_right_augmented;
                    dstate_leftdProp_vars = dstate_rightdProp_vars;
                    propagation_span_remaining -= integration_step_size;
                    this->current_independent_variable += integration_step_size;

                    try
                    {
                        this->current_epoch = this->state_left_augmented(this->index_of_epoch_in_state_vec);
                    }
                    catch (const std::out_of_range & e)
                    {
                        std::cout << "You did not set the index where epoch is stored in the state vector!! " << e.what() << std::endl;
                    }

                } // end propagation loop
            }//end backward propagation

            this->propagate_teardown(state_left, state_right, dstate_rightdProp_vars, STM, propagation_span, STM_needed);
        }

        void IntegratedFixedStepPropagator::propagate(const doubleType & propagation_span, 
                                                      const math::Matrix <doubleType> & control, 
                                                      const bool & STM_needed)
        {
            math::Matrix<doubleType> & state_left = *this->StateLeftPointer;
            math::Matrix<doubleType> & state_right = *this->StateRightPointer;
            math::Matrix<double> & dstate_leftdProp_vars = *this->dStatedIndependentVariablePointer;
            math::Matrix<double> dstate_rightdProp_vars;
            math::Matrix<double> & STM = *this->STMpointer;
            
            this->propagate_setup(state_left, dstate_leftdProp_vars, STM_needed);

            //TODO: let user set a global integration step size for the mission AND allow the user to override it at the journey level
            doubleType propagation_span_remaining = propagation_span;
            doubleType integration_step_size;
            if (propagation_span_remaining > 0.0)
            {
                while (propagation_span_remaining > 0.0)
                {
                    if (propagation_span_remaining >= this->PropagationStepSize)
                    {
                        integration_step_size = this->PropagationStepSize;
                        this->dstep_sizedProp_var = 0.0;
                    }
                    else
                    {
                        integration_step_size = propagation_span_remaining;
                        this->dstep_sizedProp_var = *this->boundary_target_dstep_sizedProp_var;
                    }

                    this->integration_scheme->step_w_control(this->state_left_augmented,
                                                             dstate_leftdProp_vars,
                                                             this->state_right_augmented,
                                                             dstate_rightdProp_vars,
                                                             control,
                                                             integration_step_size,
                                                             this->dstep_sizedProp_var,
                                                             STM_needed);

                    this->state_left_augmented = this->state_right_augmented;
                    dstate_leftdProp_vars = dstate_rightdProp_vars;
                    propagation_span_remaining -= integration_step_size;
                    this->current_independent_variable += integration_step_size;

                    try
                    {
                        this->current_epoch = this->state_left_augmented(this->index_of_epoch_in_state_vec);
                    }
                    catch (const std::out_of_range & e)
                    {
                        std::cout << "You did not set the index where epoch is stored in the state vector!! " << e.what() << std::endl;
                    }

                } // end propagation loop
            }//end forward propagation
            else //backward propagation
            {
                while (propagation_span_remaining < 0.0)
                {
                    if (-propagation_span_remaining >= this->PropagationStepSize)
                    {
                        integration_step_size = -this->PropagationStepSize;
                        this->dstep_sizedProp_var = 0.0;
                    }
                    else
                    {
                        integration_step_size = propagation_span_remaining;
                        this->dstep_sizedProp_var = -*this->boundary_target_dstep_sizedProp_var;
                    }

                    this->integration_scheme->step_w_control(this->state_left_augmented,
                                                             dstate_leftdProp_vars,
                                                             this->state_right_augmented,
                                                             dstate_rightdProp_vars,
                                                             control,
                                                             integration_step_size,
                                                             this->dstep_sizedProp_var,
                                                             STM_needed);

                    this->state_left_augmented = this->state_right_augmented;
                    dstate_leftdProp_vars = dstate_rightdProp_vars;
                    propagation_span_remaining -= integration_step_size;
                    this->current_independent_variable += integration_step_size;

                    try
                    {
                        this->current_epoch = this->state_left_augmented(this->index_of_epoch_in_state_vec);
                    }
                    catch (const std::out_of_range & e)
                    {
                        std::cout << "You did not set the index where epoch is stored in the state vector!! " << e.what() << std::endl;
                    }

                } // end propagation loop
            }//end backward propagation

            this->propagate_teardown(state_left, state_right, dstate_rightdProp_vars, STM, propagation_span, STM_needed);
        }

    } // end namespace Astrodynamics
} // end namespace EMTG