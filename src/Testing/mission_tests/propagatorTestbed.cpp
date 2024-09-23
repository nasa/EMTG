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

//event testbed class

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>

#include "KeplerPropagatorTimeDomain.h"
#include "PropagatorFactory.h"

#include "propagatorTestbed.h"

void propagatorTestbed(EMTG::missionoptions& options, 
                       std::vector< EMTG::Astrodynamics::universe > TheUniverse,
                       EMTG::HardwareModels::Spacecraft& mySpacecraft,
                       std::mt19937 RNG,
                       std::uniform_real_distribution<> UniformDouble)
{
    //Kepler propagation
    size_t DerivIndex = 0;
    size_t num_states = 10;
    size_t FirstState = DerivIndex;
    size_t stateIndex_phase_flight_time = 10;

    EMTG::math::Matrix<GSAD::adouble> StateBeforePropagation(num_states, 1, 0.0);
    EMTG::math::Matrix<GSAD::adouble> StateAfterPropagation(num_states, 1, 0.0);
    EMTG::math::Matrix<GSAD::adouble> StateAfterImpulse(num_states, 1, 0.0);
    EMTG::math::Matrix<GSAD::adouble> StateAfterSecondPropagation(num_states, 1, 0.0);
    EMTG::math::Matrix<double> STM(6, EMTG::math::identity);
    EMTG::math::Matrix<double> dStatedIndependentVariable(6, 1, 0.0);
    EMTG::math::Matrix<double> ForwardMTM(11, EMTG::math::identity);
    EMTG::math::Matrix<double> STM2(6, EMTG::math::identity);
    EMTG::math::Matrix<double> dStatedIndependentVariable2(6, 1, 0.0);


    std::vector<double> stateScales({ TheUniverse[0].LU,
        TheUniverse[0].LU,
        TheUniverse[0].LU,
        TheUniverse[0].LU / TheUniverse[0].TU,
        TheUniverse[0].LU / TheUniverse[0].TU,
        TheUniverse[0].LU / TheUniverse[0].TU,
        options.maximum_mass,
        TheUniverse[0].TU,
        options.maximum_mass,
        options.maximum_mass });

    for (size_t stateIndex = 0; stateIndex < num_states; ++stateIndex)
    {
        StateBeforePropagation(stateIndex) = UniformDouble(RNG) * stateScales[stateIndex];
        StateBeforePropagation(stateIndex).setDerivative(DerivIndex++, 1.0);
    }

    GSAD::adouble PhaseFlightTime = UniformDouble(RNG) * 1.0e+7;
    PhaseFlightTime.setDerivative(DerivIndex++, 1.0);

    std::vector<double> dPropagationTimedIndependentVariable({ 0.25, 0.5 });
    std::vector<GSAD::adouble> PropagationStepLength(2);
    PropagationStepLength[0] = PhaseFlightTime * dPropagationTimedIndependentVariable[0];
    PropagationStepLength[1] = PhaseFlightTime * dPropagationTimedIndependentVariable[1];


    EMTG::math::Matrix<GSAD::adouble> ControlVector(3, 1, 0.0);
    for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
    {
        ControlVector(controlIndex) = UniformDouble(RNG);
    }

    GSAD::adouble throttle = ControlVector.norm();
    double dPropagationTime_dIndependentVariable = 1.0;
    EMTG::Astrodynamics::PropagatorBase* myKeplerPropagatorTimeDomain = EMTG::Astrodynamics::CreatePropagator(&options,
        &TheUniverse[0],
        6,
        StateBeforePropagation,
        StateAfterPropagation,
        STM,
        dStatedIndependentVariable,
        &dPropagationTimedIndependentVariable[0]);

    myKeplerPropagatorTimeDomain->propagate(PropagationStepLength[0], true);
    StateAfterPropagation(6) = StateBeforePropagation(6);
    StateAfterPropagation(7) = StateBeforePropagation(7) + PropagationStepLength[0];
    StateAfterPropagation(8) = StateBeforePropagation(8);
    StateAfterPropagation(9) = StateBeforePropagation(9);

    EMTG::math::Matrix<double> AugmentedSTM(11, EMTG::math::identity);

    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 6; ++j)
            AugmentedSTM(i, j) = STM(i, j);
        AugmentedSTM(i, stateIndex_phase_flight_time) = dStatedIndependentVariable(i);
    }
    AugmentedSTM(7, stateIndex_phase_flight_time) = dPropagationTimedIndependentVariable[0];

    //print STM entries
    std::ofstream AugmentedSTMout("tests/KeplerAugmentedSTMout.csv", std::ios::trunc);
    AugmentedSTMout << "i, j, Analytical, Algorithmic, Absolute error, Relative error, an/al, al/an" << std::endl;
    AugmentedSTMout.precision(15);
    for (size_t i = 0; i < 10; ++i)
    {
        for (size_t j = 0; j < 11; ++j)
        {
            double STMalgorithmic = StateAfterPropagation(i).getDerivative(j + FirstState);
            double STManalytical = AugmentedSTM(i, j);
            double abserror = STManalytical - STMalgorithmic;
            double relerror = abserror / STMalgorithmic;
            double an_al = STManalytical / STMalgorithmic;
            double al_an = STMalgorithmic / STManalytical;
            AugmentedSTMout << i << "," << j << "," << STManalytical << "," << STMalgorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
        }
    }
    AugmentedSTMout.close();

    delete myKeplerPropagatorTimeDomain;
}//end main