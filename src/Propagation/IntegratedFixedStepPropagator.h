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

#ifndef INTEGRATED_FIXED_STEP_PROPAGATOR_H
#define INTEGRATED_FIXED_STEP_PROPAGATOR_H

#include "missionoptions.h"
#include "IntegratedPropagator.h"
#include "universe.h"

namespace EMTG {
    namespace Astrodynamics {

        class IntegratedFixedStepPropagator : public IntegratedPropagator
        {
        public:
            // constructors
            IntegratedFixedStepPropagator(const size_t & numStates_in,
                                          const size_t & STM_size_in);

            //clone
            virtual IntegratedFixedStepPropagator* clone() const { return new IntegratedFixedStepPropagator(*this); }

            // methods
            inline void setBoundaryTarget_dStepSizedPropVar(double * boundary_target_dstep_sizedProp_var_in) 
            { 
                this->boundary_target_dstep_sizedProp_var = boundary_target_dstep_sizedProp_var_in;
            }
            virtual void propagate(const doubleType & propagation_span, const bool & needSTM);
            virtual void propagate(const doubleType & propagation_span, const math::Matrix <doubleType> & control, const bool & needSTM);

            // fields
           
        private:
            double* boundary_target_dstep_sizedProp_var;

        };

    } // end namespace Astrodynamics
} // end namespace EMTG

#endif