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

//propagator factory
//Jacob Englander 1-2-2018

#pragma once


#include "doubleType.h"
#include "EMTG_enums.h"

#include "missionoptions.h"

#include "PropagatorBase.h"

#include "IntegratedPropagator.h"

namespace EMTG
{
    namespace Astrodynamics
    {

        PropagatorBase* CreatePropagator(missionoptions* myOptions,
                                         Astrodynamics::universe* myUniverse,
                                         const size_t& num_states,
                                         math::Matrix <doubleType> & StateLeft,
                                         math::Matrix <doubleType> & StateRight,
                                         math::Matrix <double> & STM,
                                         math::Matrix <double> & dStatedIndependentVariable,
                                         double* dPropagationTime_dIndependentVariable);

        PropagatorBase* CreatePropagator(missionoptions* myOptions,
                                         Astrodynamics::universe* myUniverse,
                                         const size_t & numStates_in,
                                         const size_t & STM_size_in,
                                         math::Matrix <doubleType> & StateLeft,
                                         math::Matrix <doubleType> & StateRight,
                                         math::Matrix <double> & STM,
                                         math::Matrix <double> & dStatedIndependentVariable,
                                         Integration::Integrand * Integrand,
                                         Integration::IntegrationScheme * IntegrationScheme,
                                         double* BoundaryTarget_dStepSizePropVar,
                                         const double PropagationStepSize);

        PropagatorBase* CreateSundmanPropagator(missionoptions* myOptions,
                                                Astrodynamics::universe* myUniverse,
                                                const size_t & numStates,
                                                const size_t & STM_size_in,
                                                math::Matrix <doubleType> & StateLeft,
                                                math::Matrix <doubleType> & StateRight,
                                                math::Matrix <double> & STM,
                                                math::Matrix <double> & dStatedIndependentVariable,
                                                Integration::Integrand * Integrand,
                                                Integration::IntegrationScheme * IntegrationScheme,
                                                double* BoundaryTarget_dStepSizePropVar,
                                                const double PropagationStepSize);

                                         
    }//end namespace HardwareModels
}//end namespace EMTG