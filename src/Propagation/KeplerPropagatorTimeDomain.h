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

//Kepler propagator in the time domain

#pragma once

#include "missionoptions.h"
#include "KeplerPropagator.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        class KeplerPropagatorTimeDomain : public KeplerPropagator
        {
        public:
            // constructors
            KeplerPropagatorTimeDomain() {};
            KeplerPropagatorTimeDomain(const size_t& numStates);

            //clone
            virtual KeplerPropagatorTimeDomain* clone() const { return new KeplerPropagatorTimeDomain(*this); }

            // methods
            void propagate(const doubleType& PropagationTime, const bool& needSTM);

            inline double getF() const { return this->F; }
            inline double getG() const { return this->G; }
            inline double getFt() const { return this->Ft; }
            inline double getGt() const { return this->Gt; }
            inline double getFtt() const { return this->Ftt; }
            inline double getGtt() const { return this->Gtt; }
            virtual std::vector<double> getPropagationHistory() const override { throw std::runtime_error("KeplerPropagatorTimeDomain::getPropagationHistory not implemented"); };

            // fields
            double F, G, Ft, Gt, Ftt, Gtt;
        };
    }//end namespace Astrodynamics
}//end namespace EMTG