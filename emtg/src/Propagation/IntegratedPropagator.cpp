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
        IntegratedPropagator::IntegratedPropagator(const size_t & numStates_in,
                                                   const size_t & STM_size_in) :
                                                   STM_size(STM_size_in)
        {
            this->numStates = numStates_in;
            this->state_left.resize(this->numStates, 1, 0.0);
            this->state_right.resize(this->numStates, 1, 0.0);
            this->STM_left.resize(STM_size_in, STM_size_in, 0.0);
            this->STM_right.resize(STM_size_in, STM_size_in, 0.0);
        }

       void  IntegratedPropagator::propagatorSetup(const math::Matrix<doubleType> & state_left, 
                                                math::Matrix<double> & STM,           
                                                const bool & STM_needed)
       {

           //configure integration scheme pointers
           this->integration_scheme->setLeftHandIndependentVariablePtr(this->current_independent_variable);

           // clear the state
           this->state_left.assign_zeros();

           size_t index = 0;
           for (size_t k = 0; k < this->numStates; ++k)
           {
               this->state_left(k) = state_left(index++);
           }
           
           this->STM_left = STM;
           if (STM_needed)
           {               
               this->STM_left.construct_identity_matrix();
           }
       }

       void IntegratedPropagator::propagatorTeardown(const math::Matrix<doubleType> & state_left,
                                                     math::Matrix<doubleType> & state_right,
                                                     math::Matrix<double> & STM_ptr,
                                                     const doubleType & propagation_span)
       {
           // pack the augmented state back into the external state pointer
           state_right = this->state_right;

           STM_ptr = this->STM_right;

           if ((state_right(7) < state_left(7) && propagation_span > 0.0)
               || (state_right(7) > state_left(7) && propagation_span < 0.0))
           {
               std::cout << "Oh noes!!...mismatch in left epoch vs. right epoch" << std::endl;
           }
       }

    } // end namespace Astrodynamics
} // end namespace EMTG