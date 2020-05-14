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

#include "IntegratedAdaptiveStepPropagator.h"
#include "missionoptions.h"
#include "universe.h"

namespace EMTG {
    namespace Astrodynamics {

        // constructors

        IntegratedAdaptiveStepPropagator::IntegratedAdaptiveStepPropagator(universe & myUniverse,
            const size_t & STM_rows_in,
            const size_t & STM_columns_in,
            const size_t & STM_start_index_in) :
            IntegratedPropagator(myUniverse,
                STM_rows_in,
                STM_columns_in,
                STM_start_index_in)

        {

        }

        // methods
        void IntegratedAdaptiveStepPropagator::propagate(const doubleType & propagation_span, const bool & STM_needed)
        {
            //TODO: probably won't have to do this...just call an appropriately overloaded IntegrationScheme step method
            this->propagate(propagation_span, this->zero_control, STM_needed);
        }

        void IntegratedAdaptiveStepPropagator::propagate(const doubleType & propagation_span, const math::Matrix <doubleType> & control, const bool & STM_needed)
        {
            math::Matrix<doubleType> & state_left = *this->StateLeftPointer;
            math::Matrix<doubleType> & state_right = *this->StateRightPointer;
            math::Matrix<double> & dstate_leftdProp_vars = *this->dStatedIndependentVariablePointer;
            math::Matrix<double> dstate_rightdProp_vars;
            math::Matrix<double> & STM = *this->STMpointer;

            //configure integration scheme pointers
            this->integration_scheme->setLeftHandIndependentVariablePtr(this->current_epoch);
            this->integration_scheme->setdLeftHandIndVardPropVarPtr(this->dcurrent_ind_vardProp_var);
            this->integration_scheme->setdLeftHandIndVardPropVarPreviousPtr(this->dcurrent_ind_vardProp_var_previous);

            this->dcurrent_ind_vardProp_var_previous = 1.0;
            this->dcurrent_ind_vardProp_var = 0.0;

            // unpack the external state pointer into the integrator's augmented state
            this->unpackStates(state_left, STM_needed);

            this->propagateAdaptiveStep(propagation_span, control, dstate_leftdProp_vars, dstate_rightdProp_vars, STM_needed);

            // If we got here, we made it all of the way through the whole propagation_span
            // pack the augmented state back into the external state pointer
            this->packStates(state_right, dstate_rightdProp_vars, STM, STM_needed);

            // now compute the partials of the state w.r.t. the propagation variables
            if (STM_needed)
            {
                this->computePropVarPartials(propagation_span, control, state_left, state_right, dstate_leftdProp_vars, dstate_rightdProp_vars);

                // the left hand container is also the output container
                dstate_leftdProp_vars = dstate_rightdProp_vars;
            }

        } // end propagate method

        void IntegratedAdaptiveStepPropagator::computePropVarPartials(const doubleType & propagation_span, 
                                                                      const math::Matrix <doubleType> & control, 
                                                                      math::Matrix<doubleType> & state_left,
                                                                      math::Matrix<doubleType> & state_right,
                                                                      math::Matrix<double> & dstate_leftdProp_vars,
                                                                      math::Matrix<double> & dstate_rightdProp_vars)
        {
                math::Matrix<doubleType> states_perturbed_foward;
                math::Matrix<doubleType> states_perturbed_backward;
                //double central_difference_interval = 6.7e-05;
                double central_difference_interval = 10.0;
                const double one_over_two_step = 1.0 / (2.0 * central_difference_interval);

                // store the original value of the independent variable
                doubleType original_independent_variable = this->current_epoch;

                // reset the left hand state container
                // also ensure that STM_needed is set to false so that we don't needlessly compute STM entries
                // during the finite differencing
                this->unpackStates(state_left, false);

                // forward perturb the current independent variable
                this->current_epoch += central_difference_interval;
                this->propagateAdaptiveStep(propagation_span, control, dstate_leftdProp_vars, dstate_rightdProp_vars, false);
                states_perturbed_foward = this->state_right_augmented;

                // backward perturb the current independent variable
                this->unpackStates(state_left, false);
                this->current_epoch = original_independent_variable - central_difference_interval;
                this->propagateAdaptiveStep(propagation_span, control, dstate_leftdProp_vars, dstate_rightdProp_vars, false);
                states_perturbed_backward = this->state_right_augmented;

                for (size_t k = 0; k < this->STM_start_index; ++k)
                {
                    dstate_rightdProp_vars(k, 0) = ((states_perturbed_foward(k) - states_perturbed_backward(k)) * one_over_two_step) _GETVALUE;
                }

                // clean up the perturbation of the current independent variable
                this->current_epoch = original_independent_variable;

                // forward perturb the propagation span
                this->unpackStates(state_left, false);
                this->propagateAdaptiveStep(propagation_span + central_difference_interval, control, dstate_leftdProp_vars, dstate_rightdProp_vars, false);
                states_perturbed_foward = this->state_right_augmented;

                // backward perturb the propagation span
                this->unpackStates(state_left, false);
                this->propagateAdaptiveStep(propagation_span - central_difference_interval, control, dstate_leftdProp_vars, dstate_rightdProp_vars, false);
                states_perturbed_backward = this->state_right_augmented;


                for (size_t k = 0; k < this->STM_start_index; ++k)
                {
                    dstate_rightdProp_vars(k, 1) = ((states_perturbed_foward(k) - states_perturbed_backward(k)) * one_over_two_step) _GETVALUE;

                    // the propagator is only aware of the local propagation span that is passed to it
                    // therefore we must divide by the propagation variable modifier (e.g. 1/N for Sims-Flanagan and FBLT)
                    // difference is (TOF/N + dTOF)    vs.     (TOF + dTOF) / N
                    // if the propagation is backwards, then the sign is flipped
                    dstate_rightdProp_vars(k, 1) *= propagation_span < 0.0 ? -this->boundary_target_dstep_sizedProp_var : this->boundary_target_dstep_sizedProp_var;
                }

        } // end computePropVarPartials method

        void IntegratedAdaptiveStepPropagator::propagateAdaptiveStep(const doubleType & propagation_span, 
                                                                     const math::Matrix <doubleType> & control, 
                                                                     math::Matrix<double> & dstate_leftdProp_vars,
                                                                     math::Matrix<double> & dstate_rightdProp_vars,
                                                                     const bool & STM_needed)
        {


            doubleType accumulatedH = 0.0;
            double daccumulatedHdTOF = 0.0;
            doubleType effectiveH = this->PropagationStepSize;
            doubleType nextStep = effectiveH;

            // TODO: must finite difference for this
            // double deffectiveHdTOF = 0.0;

            doubleType adaptive_step_error = 1.0e-20;
            doubleType worst_error_state = 0.0;

            doubleType t_left_step = this->current_epoch;

            // If the user specifies an integration step size that is too big, then
            // cap it at the propagation length
            if (effectiveH > propagation_span)
            {
                effectiveH = propagation_span;
            }

            bool last_step = false; // at the beginning of a propagation segment we are on the FIRST substep

            // loop until we get all the way through the full propagation_span 
            do
            {
                // take a trial step
                do
                { // cycle until the trial RK step achieves sufficient accuracy

                    effectiveH = nextStep;
                    // Take the trial RK step
                    this->integration_scheme->errorControlledStep(this->state_left_augmented,
                                                                  dstate_leftdProp_vars,
                                                                  this->state_right_augmented,
                                                                  dstate_rightdProp_vars,
                                                                  control,
                                                                  effectiveH,
                                                                  this->dstep_sizedProp_var,
                                                                  STM_needed,
                                                                  adaptive_step_error,
                                                                  this->error_scaling_factors);


                    if (!last_step)
                    {
                        // no error!  give it a real value so we don't divide by zero.
                        if (adaptive_step_error == 0.0)
                        {
                            adaptive_step_error = 1e-15; //Almost zero!
                        }


                        // if we rejected the last sub-step (i.e. the error was too large) shorten the time step
                        if (adaptive_step_error >= this->integrator_tolerance)
                        {
                            //effectiveH = 0.98*effectiveH*pow(this->myOptions->integrator_tolerance / adaptive_step_error, 0.17);
                            nextStep = 0.98 * effectiveH * pow(this->integrator_tolerance / adaptive_step_error, 0.17);
                        }
                        else // make the sub-step a bit longer to save computation time
                        {
                            //effectiveH = 1.01*effectiveH*pow(this->myOptions->integrator_tolerance / adaptive_step_error, 0.18);
                            nextStep = 1.01 * effectiveH * pow(this->integrator_tolerance / adaptive_step_error, 0.18);

                            // if our increased step kicks us too long, make it shorter anyway and just run it
                            if (fabs(propagation_span - accumulatedH) < fabs(nextStep) && !last_step)
                            {
                                nextStep = propagation_span - accumulatedH;
                            }

                        }

                        if (fabs(propagation_span - accumulatedH) < fabs(nextStep) && !last_step)
                        {
                            nextStep = propagation_span - accumulatedH;
                        }

                        //if we make the time step too small, kill the integration - h is too small
                        if (fabs(nextStep) < 1e-13)
                        {
                            throw std::runtime_error("rk7813M: H Got too Small. The integrator has Alexed. Aborting.");
                        }

                    }
                    else if (adaptive_step_error > this->integrator_tolerance && last_step)
                    {
                        // we got here because we thought it was the last substep, and upon calculation it 
                        // was too big and not precise enough so we have to shrink it
                        last_step = false; // not last step after all
                                           //effectiveH = 0.98*effectiveH*pow(this->myOptions->integrator_tolerance / adaptive_step_error, 0.17);
                        nextStep = 0.98 * effectiveH * pow(this->integrator_tolerance / adaptive_step_error, 0.17);
                    }

                } while (adaptive_step_error > this->integrator_tolerance);

                // if we got here, then the trial substep was accurate enough; it becomes the new left
                this->state_left_augmented = this->state_right_augmented;
                dstate_leftdProp_vars = dstate_rightdProp_vars;

                // keep track of our progress through the full RK step
                accumulatedH += effectiveH;

                // move the left hand epoch for the next substep forward to the correct value
                this->current_epoch += effectiveH;

                // if our next step will push us over, reduce it down to be as small as necessary to hit target exactly
                if (fabs(propagation_span - accumulatedH) < fabs(nextStep) && fabs(propagation_span - accumulatedH) > 0 && !last_step)
                {
                    nextStep = propagation_span - accumulatedH;
                    last_step = true; //assume that the next substep will be the last substep now
                }

            } while (fabs(accumulatedH) < fabs(propagation_span));
        } // end propagateAdaptiveStep method


    } // end namespace Astrodynamics
} // end namespace EMTG