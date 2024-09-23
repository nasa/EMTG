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

//base class for Kepler-type propagators

#pragma once

#include "PropagatorBase.h"
#include "missionoptions.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        class KeplerPropagator : public PropagatorBase
        {
        public:
            //constructors
            KeplerPropagator() {};
            KeplerPropagator(const size_t& numStates);

            //methods
            virtual void propagate(const doubleType& PropagationTime, const bool& needSTM) = 0;
            inline void setCentralBodyGM(const double& mu) { this->mu = mu; };
            virtual std::vector<double> getPropagationHistory() const = 0;

            void set_dPropagationTime_dIndependentVariable(double * dPropagationTime_dIndependentVariable)
            {
                this->dPropagationTime_dIndependentVariable = dPropagationTime_dIndependentVariable;
            }

        protected:
            //fields
            double mu;
            double* dPropagationTime_dIndependentVariable;
        };
    }//end namespace Astrodynamics
}//end namespace EMTG