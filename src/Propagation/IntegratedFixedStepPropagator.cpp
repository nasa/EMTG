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

        IntegratedFixedStepPropagator::IntegratedFixedStepPropagator(const size_t & numStates_in,
                                                                     const size_t & STM_size_in) :
                                                                     IntegratedPropagator(numStates_in, 
                                                                                          STM_size_in)
        {
            this->propagation_history.reserve(1000);
        }

        // methods
        void IntegratedFixedStepPropagator::propagate(const doubleType & propagation_span, 
                                                      const bool & STM_needed)
        {
            math::Matrix<doubleType> & state_left = *this->StateLeftPointer;
            math::Matrix<doubleType> & state_right = *this->StateRightPointer;
            math::Matrix<double> & STM_ptr = *this->STMpointer;

            this->propagatorSetup(state_left, STM_ptr, STM_needed);

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

                    this->integration_scheme->step(this->state_left,
                                                   this->STM_left,
                                                   this->state_right,
                                                   this->STM_right,
                                                   integration_step_size,
                                                   this->dstep_sizedProp_var,
                                                   STM_needed);

                    if (this->store_propagation_history)
                    {
                        this->propagation_history.push_back(this->state_right(this->index_of_epoch_in_state_vec) _GETVALUE - this->state_left(this->index_of_epoch_in_state_vec) _GETVALUE);
                    }

                    this->state_left = this->state_right;
                    this->STM_left = this->STM_right;
                    propagation_span_remaining -= integration_step_size;
                    this->current_independent_variable += integration_step_size;

                    try
                    {
                        this->current_epoch = this->state_left(this->index_of_epoch_in_state_vec);
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

                    this->integration_scheme->step(this->state_left,
                                                   this->STM_left,
                                                   this->state_right,
                                                   this->STM_right,
                                                   integration_step_size,
                                                   this->dstep_sizedProp_var,
                                                   STM_needed);

                    if (this->store_propagation_history)
                    {
                        this->propagation_history.push_back(this->state_right(this->index_of_epoch_in_state_vec) _GETVALUE - this->state_left(this->index_of_epoch_in_state_vec) _GETVALUE);
                    }

                    this->state_left = this->state_right;
                    this->STM_left = this->STM_right;
                    propagation_span_remaining -= integration_step_size;
                    this->current_independent_variable += integration_step_size;

                    try
                    {
                        this->current_epoch = this->state_left(this->index_of_epoch_in_state_vec);
                    }
                    catch (const std::out_of_range & e)
                    {
                        std::cout << "You did not set the index where epoch is stored in the state vector!! " << e.what() << std::endl;
                    }

                } // end propagation loop
            }//end backward propagation

            this->propagatorTeardown(state_left, state_right, STM_ptr, propagation_span);
        }

        void IntegratedFixedStepPropagator::propagate(const doubleType & propagation_span, 
                                                      const math::Matrix <doubleType> & control, 
                                                      const bool & STM_needed)
        {
            math::Matrix<doubleType> & state_left = *this->StateLeftPointer;
            math::Matrix<doubleType> & state_right = *this->StateRightPointer;
            math::Matrix<double> & STM_ptr = *this->STMpointer;
            
            this->propagatorSetup(state_left, STM_left, STM_needed);

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

                    this->integration_scheme->step(this->state_left,
                                                   this->STM_left,
                                                   this->state_right,
                                                   this->STM_right,
                                                   control,
                                                   integration_step_size,
                                                   this->dstep_sizedProp_var,
                                                   STM_needed);

                    if (this->store_propagation_history)
                    {
                        this->propagation_history.push_back(this->state_right(this->index_of_epoch_in_state_vec) _GETVALUE - this->state_left(this->index_of_epoch_in_state_vec) _GETVALUE);
                    }

                    this->state_left = this->state_right;
                    this->STM_left = this->STM_right;
                    propagation_span_remaining -= integration_step_size;
                    this->current_independent_variable += integration_step_size;

                    try
                    {
                        this->current_epoch = this->state_left(this->index_of_epoch_in_state_vec);
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

                    this->integration_scheme->step(this->state_left,
                                                   this->STM_left,
                                                   this->state_right,
                                                   this->STM_right,
                                                   control,
                                                   integration_step_size,
                                                   this->dstep_sizedProp_var,
                                                   STM_needed);

                    if (this->store_propagation_history)
                    {
                        this->propagation_history.push_back(this->state_right(this->index_of_epoch_in_state_vec) _GETVALUE - this->state_left(this->index_of_epoch_in_state_vec) _GETVALUE);
                    }

                    this->state_left = this->state_right;
                    this->STM_left = this->STM_right;
                    propagation_span_remaining -= integration_step_size;
                    this->current_independent_variable += integration_step_size;

                    try
                    {
                        this->current_epoch = this->state_left(this->index_of_epoch_in_state_vec);
                    }
                    catch (const std::out_of_range & e)
                    {
                        std::cout << "You did not set the index where epoch is stored in the state vector!! " << e.what() << std::endl;
                    }

                } // end propagation loop
            }//end backward propagation

            this->propagatorTeardown(state_left, state_right, STM_ptr, propagation_span);
        }

    } // end namespace Astrodynamics
} // end namespace EMTG