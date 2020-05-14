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

#include "Integrand.h"
#include "IntegrationScheme.h"
#include "IntegratedPropagator.h"
#include "PropagatorBase.h"
#include "universe.h"

namespace EMTG {
    namespace Astrodynamics {

        // constructors
        IntegratedPropagator::IntegratedPropagator(Astrodynamics::universe & myUniverse,
                                                   const size_t & STM_rows_in,
                                                   const size_t & STM_columns_in,
                                                   const size_t & STM_start_index_in) :
                                                   PropagatorBase(myUniverse),
                                                   num_STM_row(STM_rows_in),
                                                   num_STM_col(STM_columns_in),
                                                   STM_start_index(STM_start_index_in)
        {
            this->zero_control.resize(4, 1, 0.0);
            this->state_left_augmented.resize(STM_start_index_in + STM_rows_in * STM_columns_in, 1, 0.0);
            this->state_right_augmented.resize(STM_start_index_in + STM_rows_in * STM_columns_in, 1, 0.0);
        }

       void  IntegratedPropagator::unpackStates(const math::Matrix<doubleType> & state_left, const bool & STM_needed)
       {
           // clear the state
           this->state_left_augmented.assign_zeros();

           // clear dStatedIndependentVariable
           this->dStatedIndependentVariablePointer->assign_zeros();

           // insert true state into augmented state
           size_t index = 0;
           for (size_t k = 0; k < this->STM_start_index; ++k)
           {
               // TODO: FIX THIS MONSTRO-CITY: The incoming state vector has the current epoch as its
               // 8th entry. We want to skip that as we don't have current epoch dynamics
               // to integrate
               //if (k == 7)
               //{
               //    ++index;
               //}

               this->state_left_augmented(k) = state_left(index++);
           }

           if (STM_needed)
           {
               // make sure we integrate true states + STM entries
               this->setNumStatesToIntegrate(this->STM_start_index + this->num_STM_row * this->num_STM_col);

               // set STM to identity
               for (size_t k = this->STM_start_index; k < this->state_left_augmented.get_n(); k = k + this->num_STM_row + 1)
               {
                   this->state_left_augmented(k) = 1.0;
               }
           }
           else
           {
               this->setNumStatesToIntegrate(this->STM_start_index);
           }
       }

       void  IntegratedPropagator::packStates(math::Matrix<doubleType> & state_right, 
                                              const math::Matrix<double> & dstate_rightdProp_vars,
                                              math::Matrix<double> & STM, 
                                              const bool & STM_needed)
       {
           size_t index = 0;
           for (size_t k = 0; k < this->STM_start_index; ++k)
           {
               state_right(index++) = this->state_right_augmented(k);
           }

           if (STM_needed)
           {
               size_t STM_entry_index = this->STM_start_index;
               for (size_t i = 0; i < this->num_STM_row; ++i)
               {
                   for (size_t j = 0; j < this->num_STM_col; ++j)
                   {
                       STM(i, j) = (this->state_right_augmented(STM_entry_index++)) _GETVALUE;
                   }
               }

               // insert the current epoch and propagation variable derivatives into the STM
               size_t prop_var_slot = STM.get_m() - 1; // propagation variable slot is always the last STM column
               for (size_t i = 0; i < 10; ++i)
               {
                   // previous flight times (current epoch)
                   STM(i, this->index_of_epoch_in_state_vec) = dstate_rightdProp_vars(i, 0);
                   // current phase propagation variable
                   STM(i, prop_var_slot) = dstate_rightdProp_vars(i, 1);
               }
               STM(prop_var_slot, prop_var_slot) = 1.0;
           }
       }

    } // end namespace Astrodynamics
} // end namespace EMTG