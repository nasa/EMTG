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

#include "PropagatorFactory.h"

#include "IntegratedAdaptiveStepPropagator.h"
#include "IntegratedFixedStepPropagator.h"
#include "KeplerPropagatorTimeDomain.h"

#include <exception>

namespace EMTG
{
    namespace Astrodynamics
    {
        //for non-integrators
        PropagatorBase* CreatePropagator(missionoptions* myOptions,
                                         Astrodynamics::universe* myUniverse,
                                         const size_t& num_states,
                                         math::Matrix <doubleType> & StateLeft,
                                         math::Matrix <doubleType> & StateRight,
                                         math::Matrix <double> & STM,
                                         math::Matrix <double> & dStatedIndependentVariable,
                                         double* dPropagationTime_dIndependentVariable)
        {
            KeplerPropagatorTimeDomain* myPropagator = new KeplerPropagatorTimeDomain(num_states);
            myPropagator->setStateLeft(StateLeft);
            myPropagator->setStateRight(StateRight);
            myPropagator->setSTM(STM);
            myPropagator->setdStatedIndependentVariable(dStatedIndependentVariable);
            myPropagator->set_dPropagationTime_dIndependentVariable(dPropagationTime_dIndependentVariable);
            myPropagator->setCentralBodyGM(myUniverse->mu);

            return myPropagator;
        }

        //for integrators
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
                                         const double PropagationStepSize)
        {
            if (myOptions->integratorType == IntegratorType::rk8_fixed)
            {
                IntegratedFixedStepPropagator* myPropagator = new IntegratedFixedStepPropagator(numStates_in, STM_size_in);
                myPropagator->setStateLeft(StateLeft);
                myPropagator->setStateRight(StateRight);
                myPropagator->setSTM(STM);
                myPropagator->setdStatedIndependentVariable(dStatedIndependentVariable);
                myPropagator->setIntegrand(Integrand);
                myPropagator->setIntegrationScheme(IntegrationScheme);
                myPropagator->setBoundaryTarget_dStepSizedPropVar(BoundaryTarget_dStepSizePropVar);
                myPropagator->setPropagationStepSize(PropagationStepSize);


                return myPropagator;
            }
            else if (myOptions->integratorType == IntegratorType::rk7813m_adaptive)
            {
                IntegratedAdaptiveStepPropagator* myPropagator = new IntegratedAdaptiveStepPropagator(numStates_in, STM_size_in);
                myPropagator->setStateLeft(StateLeft);
                myPropagator->setStateRight(StateRight);
                myPropagator->setSTM(STM);
                myPropagator->setdStatedIndependentVariable(dStatedIndependentVariable);
                myPropagator->setIntegrand(Integrand);
                myPropagator->setIntegrationScheme(IntegrationScheme);
                myPropagator->setBoundaryTarget_dStepSizedPropVar(BoundaryTarget_dStepSizePropVar);
                myPropagator->setPropagationStepSize(PropagationStepSize);
                myPropagator->setTolerance(myOptions->integrator_tolerance);

                math::Matrix<double> error_scaling_factors(numStates_in + STM_size_in * STM_size_in, 1, 0.0);

                size_t index = 0;
                // position and velocity error scaling
                for (size_t k = 0; k < 3; ++k)
                {
                    error_scaling_factors(index) = 1.0 / myUniverse->LU;
                    error_scaling_factors(index + 3) = myUniverse->TU / myUniverse->LU;
                    ++index;
                }
                index += 3;

                // mass error scaling
                for (size_t k = 0; k < 3; ++k)
                {
                    error_scaling_factors(index++) = 1.0 / myOptions->maximum_mass;
                }

                // STM position rows error scaling
                for (size_t k = 0; k < 3; ++k)
                {
                    // STM R~ error scaling not required
                    error_scaling_factors(index++) = 1.0;
                    error_scaling_factors(index++) = 1.0;
                    error_scaling_factors(index++) = 1.0;

                    // STM R error scaling
                    error_scaling_factors(index++) = 1.0 / myUniverse->TU;
                    error_scaling_factors(index++) = 1.0 / myUniverse->TU;
                    error_scaling_factors(index++) = 1.0 / myUniverse->TU;

                    // mass and tanks
                    error_scaling_factors(index++) = myOptions->maximum_mass / myUniverse->LU;
                    error_scaling_factors(index++) = myOptions->maximum_mass / myUniverse->LU;
                    error_scaling_factors(index++) = myOptions->maximum_mass / myUniverse->LU;

                    // control
                    error_scaling_factors(index++) = 1.0 / myUniverse->LU;
                    error_scaling_factors(index++) = 1.0 / myUniverse->LU;
                    error_scaling_factors(index++) = 1.0 / myUniverse->LU;
                }

                // STM velocity rows error scaling
                for (size_t k = 0; k < 3; ++k)
                {
                    // STM V~ scaling
                    error_scaling_factors(index++) = myUniverse->TU;
                    error_scaling_factors(index++) = myUniverse->TU;
                    error_scaling_factors(index++) = myUniverse->TU;

                    // STM V error scaling not required
                    error_scaling_factors(index++) = 1.0;
                    error_scaling_factors(index++) = 1.0;
                    error_scaling_factors(index++) = 1.0;
                    
                    // mass and tanks
                    error_scaling_factors(index++) = myOptions->maximum_mass * myUniverse->TU / myUniverse->LU;
                    error_scaling_factors(index++) = myOptions->maximum_mass * myUniverse->TU / myUniverse->LU;
                    error_scaling_factors(index++) = myOptions->maximum_mass * myUniverse->TU / myUniverse->LU;

                    // control
                    error_scaling_factors(index++) = myUniverse->TU / myUniverse->LU;
                    error_scaling_factors(index++) = myUniverse->TU / myUniverse->LU;
                    error_scaling_factors(index++) = myUniverse->TU / myUniverse->LU;
                }

                // mass row scaling
                error_scaling_factors(index++) = myUniverse->LU / myOptions->maximum_mass;
                error_scaling_factors(index++) = myUniverse->LU / myOptions->maximum_mass;
                error_scaling_factors(index++) = myUniverse->LU / myOptions->maximum_mass;
                error_scaling_factors(index++) = myUniverse->LU / (myOptions->maximum_mass * myUniverse->TU);
                error_scaling_factors(index++) = myUniverse->LU / (myOptions->maximum_mass * myUniverse->TU);
                error_scaling_factors(index++) = myUniverse->LU / (myOptions->maximum_mass * myUniverse->TU);
                error_scaling_factors(index++) = 1.0;
                error_scaling_factors(index++) = 1.0;
                error_scaling_factors(index++) = 1.0;
                error_scaling_factors(index++) = 1.0 / myOptions->maximum_mass;
                error_scaling_factors(index++) = 1.0 / myOptions->maximum_mass;
                error_scaling_factors(index++) = 1.0 / myOptions->maximum_mass;

                // control rows do not require scaling as control partials are all zero 
                // (they are decision variables) 

                myPropagator->setErrorScalingFactors(error_scaling_factors);

                return myPropagator;
            }
            else
            {
            throw std::invalid_argument("Integrator type not implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            return NULL;
        }
    }//end namespace HardwareModels
}//end namespace EMTG