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

//base class for propagators

#pragma once

#include "doubleType.h"

#include <string>
#include <vector>

#include <EMTG_Matrix.h>

namespace EMTG
{
    namespace Astrodynamics
    {
        class PropagatorBase
        {
        public:
            //constructors
            PropagatorBase();
			
			virtual ~PropagatorBase() = default;

            //clone
            virtual PropagatorBase* clone() const = 0;

            //methods
            virtual void propagate(const doubleType & PropagationTime, const bool & needSTM) = 0;
            virtual void propagate(const doubleType & propagation_span, const math::Matrix <doubleType> & control, const bool & needSTM) {};

            inline void setStateLeft(math::Matrix <doubleType> & StateLeft) { this->StateLeftPointer = &StateLeft; };
            inline void setStateRight(math::Matrix <doubleType> & StateRight) { this->StateRightPointer = &StateRight; };
            inline void setSTM(math::Matrix<double>& STM) 
            { 
                this->STMpointer = &STM;
                this->prop_var_column_index = STM.get_m() - 1; // propagation variable slot is always the last STM column
            };
            inline void setdStatedIndependentVariable(math::Matrix<double> & dStatedIndependentVariable) { this->dStatedIndependentVariablePointer = &dStatedIndependentVariable; };
            inline void setCurrentEpoch(const doubleType & current_epoch_in) { this->current_epoch = current_epoch_in; };
            inline void setCurrentIndependentVariable(const doubleType & current_indvar_in) { this->current_independent_variable = current_indvar_in; };
            inline void setIndexOfEpochInStateVec(const size_t & index_in) { this->index_of_epoch_in_state_vec = index_in;  };
            inline void setPropVarSTMcolumnIndex(const size_t & index_in) { };
            inline void setPropagationStepSize(const double & PropagationStepSize) { this->PropagationStepSize = PropagationStepSize; };
            inline void setStorePropagationHistory(const bool & switch_in) { this->store_propagation_history = switch_in; };
            inline void clearPropagationHistory() { this->propagation_history.clear(); };
            virtual std::vector<double> getPropagationHistory() const = 0;
            

            //fields
        protected:

            size_t numStates;

            bool store_propagation_history;
            std::vector<double> propagation_history;

            math::Matrix<double>* STMpointer;
            math::Matrix<doubleType>* StateLeftPointer;
            math::Matrix<doubleType>* StateRightPointer;
            math::Matrix<double>* dStatedIndependentVariablePointer;

            doubleType current_epoch;
            doubleType current_independent_variable;
            size_t index_of_epoch_in_state_vec;
            size_t prop_var_column_index;
            
            // Partial of the current independent variable w.r.t. previous propagation variables (flight times)
            double dcurrent_ind_vardProp_var;

            // Partial of the current independent variable w.r.t. the current propagation variable (flight times or anomaly)
            double dcurrent_ind_vardProp_var_previous;

            //propagation step size
            double PropagationStepSize;
        };

        inline PropagatorBase * new_clone(PropagatorBase const & other)
        {
            return other.clone();
        }
    }//end namespace Astrodynamics
}//end namespace EMTG