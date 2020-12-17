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

#include "ForwardBoundedImpulseManeuver.h"
#include "BackwardBoundedImpulseManeuver.h"

#include "SpacecraftAccelerationModel.h"

#include "MTMtestbed.h"

void MTMtestbed(EMTG::missionoptions& options,
    std::vector< EMTG::Astrodynamics::universe > TheUniverse,
    EMTG::HardwareModels::Spacecraft& mySpacecraft,
    std::mt19937 RNG,
    std::uniform_real_distribution<> UniformDouble)
{
    mySpacecraft.setActiveStage(0);

    size_t DerivIndex = 0;
    size_t num_states = 10;
    size_t FirstState = DerivIndex;

    EMTG::math::Matrix<GSAD::adouble> StateLeft(num_states, 1, 0.0);
    EMTG::math::Matrix<GSAD::adouble> StateRight(num_states, 1, 0.0);
    EMTG::math::Matrix<double> ForwardMTM(11, EMTG::math::identity);
    EMTG::math::Matrix<double> BackwardMTM(11, EMTG::math::identity);

    std::vector<double> stateScales({ TheUniverse[0].LU,
        TheUniverse[0].LU,
        TheUniverse[0].LU,
        TheUniverse[0].LU / TheUniverse[0].TU,
        TheUniverse[0].LU / TheUniverse[0].TU,
        TheUniverse[0].LU / TheUniverse[0].TU,
        options.maximum_mass,
        TheUniverse[0].TU,
        options.maximum_mass,
        options.maximum_mass,
        TheUniverse[0].TU});
    TheUniverse[0].set_X_scale_factors(&stateScales);
    std::vector<std::string> Xdescriptions({ "x", "y", "z", "xdot", "ydot", "zdot", "mass", "epoch","hydrazine","xenon","propagation time" });

    size_t numStatesToPropagate = 10;
    EMTG::Astrodynamics::SpacecraftAccelerationModel mySpacecraftAccelerationModel(&options,
        &options.Journeys[0],
        &TheUniverse[0],
        &Xdescriptions,
        &mySpacecraft,
        14); // STM size
    mySpacecraftAccelerationModel.setDutyCycle(options.engine_duty_cycle);

    for (size_t stateIndex = 0; stateIndex < num_states; ++stateIndex)
    {
        if (stateIndex == 7)
            StateLeft(stateIndex) = 53800.0 * 86400.0 + UniformDouble(RNG) * stateScales[stateIndex];
        else
            StateLeft(stateIndex) = UniformDouble(RNG) * stateScales[stateIndex];
        StateLeft(stateIndex).setDerivative(DerivIndex++, stateScales[stateIndex]);
    }

    GSAD::adouble LaunchDate = StateLeft(7);

    GSAD::adouble PhaseFlightTime = UniformDouble(RNG) * 1.0e+7;
    PhaseFlightTime.setDerivative(DerivIndex++, stateScales[10]);

    GSAD::adouble ThrustStepLength = PhaseFlightTime / options.num_timesteps;

    //advance time to half of the thrust step
    StateLeft(7) += ThrustStepLength / 2;

    EMTG::math::Matrix<GSAD::adouble> ControlVector(3, 1, 0.0);
    for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
    {
        ControlVector(controlIndex) = UniformDouble(RNG);
        ControlVector(controlIndex).setDerivative(DerivIndex++, 1.0);
    }

    GSAD::adouble throttle = ControlVector.norm();

    EMTG::Phases::ForwardBoundedImpulseManeuver myForwardBoundedImpulseManeuver(0,//journeyIndex
        0,//phaseIndex,
        0,//stageIndex,
        &TheUniverse[0],
        &mySpacecraft,
        &options);

    EMTG::math::Matrix<GSAD::adouble> DeltaV(3, 1, 0.0);

    GSAD::adouble max_thrust, max_mass_flow_rate, Isp, power, active_power;
    size_t number_of_active_engines;
    int ThrottleLevel;
    std::string ThrottleLevelString("banana");
    double dThrustStepLength_dPropagationVariable = 1.0 / options.num_timesteps;
    double dImpulseEpoch_dPropagationVariable = (1 - 0.5) / options.num_timesteps;

    myForwardBoundedImpulseManeuver.set_ThrustStepLength(ThrustStepLength);
    myForwardBoundedImpulseManeuver.set_dThrustStepLength_dPropagationVariable(dThrustStepLength_dPropagationVariable);
    myForwardBoundedImpulseManeuver.set_dImpulseEpoch_dPropagationVariable_pointer(dImpulseEpoch_dPropagationVariable);
    myForwardBoundedImpulseManeuver.set_LaunchDate(LaunchDate);
    myForwardBoundedImpulseManeuver.set_max_thrust(max_thrust);
    myForwardBoundedImpulseManeuver.set_max_mass_flow_rate(max_mass_flow_rate);
    myForwardBoundedImpulseManeuver.set_Isp(Isp);
    myForwardBoundedImpulseManeuver.set_power(power);
    myForwardBoundedImpulseManeuver.set_active_power(active_power);
    myForwardBoundedImpulseManeuver.set_number_of_active_engines(number_of_active_engines);
    myForwardBoundedImpulseManeuver.set_ThrottleLevel(ThrottleLevel);
    myForwardBoundedImpulseManeuver.set_ThrottleLevelString(ThrottleLevelString);
    myForwardBoundedImpulseManeuver.set_DeltaV(DeltaV);
    myForwardBoundedImpulseManeuver.set_spacecraft_state_minus(StateLeft);
    myForwardBoundedImpulseManeuver.set_spacecraft_state_plus(StateRight);
    myForwardBoundedImpulseManeuver.set_MTM(ForwardMTM);
    myForwardBoundedImpulseManeuver.set_DutyCycle(options.engine_duty_cycle);
    myForwardBoundedImpulseManeuver.set_Control(ControlVector);
    myForwardBoundedImpulseManeuver.set_Throttle(throttle);
    myForwardBoundedImpulseManeuver.set_AccelerationModel(&mySpacecraftAccelerationModel);

    myForwardBoundedImpulseManeuver.process_maneuver(true);


    //print MTM entries
    std::ofstream ForwardMTMout("tests/ForwardBoundedImpulseManeuver.csv", std::ios::trunc);
    ForwardMTMout << "i, j, Analytical, Algorithmic, Absolute error, Relative error, an/al, al/an" << std::endl;
    for (size_t i = 0; i < num_states; ++i)
    {
        for (size_t j = 0; j < num_states + 1; ++j)
        {
            double MTMalgorithmic = StateRight(i).getDerivative(j + FirstState) / stateScales[j];
            double MTManalytical = ForwardMTM(i, j);
            double abserror = MTManalytical - MTMalgorithmic;
            double relerror = abserror / MTMalgorithmic;
            double an_al = MTManalytical / MTMalgorithmic;
            double al_an = MTMalgorithmic / MTManalytical;
            ForwardMTMout << i << "," << j << "," << MTManalytical << "," << MTMalgorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
        }
    }
    ForwardMTMout.close();

    std::ofstream ForwardControlDerivatives("tests/ForwardControlDerivatives.csv", std::ios::trunc);
    ForwardControlDerivatives << "stateIndex, vIndex, Analytical, Algorithmic, Absolute error, Relative error, an/al, al/an" << std::endl;
    EMTG::math::Matrix<double> Forward_dMassAfterManeuver_dThrottleComponents = myForwardBoundedImpulseManeuver.get_dMassAfterManeuver_dThrottleComponents();
    EMTG::math::Matrix<double> Forward_dChemicalFuel_dThrottleComponents = myForwardBoundedImpulseManeuver.get_dChemicalFuel_dThrottleComponents();
    EMTG::math::Matrix<double> Forward_dElectricPropellant_dThrottleComponents = myForwardBoundedImpulseManeuver.get_dElectricPropellant_dThrottleComponents();

    {
        //mass wrt control
        size_t stateIndex = 6;
        for (size_t vIndex = 0; vIndex < 3; ++vIndex)
        {
            size_t GSADindex = num_states + 1 + vIndex;
            size_t Xindex = 1 + vIndex;
            double ControlDerivativealgorithmic = StateRight(stateIndex).getDerivative(GSADindex);
            double ControlDerivativeanalytical = Forward_dMassAfterManeuver_dThrottleComponents(vIndex);
            double abserror = ControlDerivativeanalytical - ControlDerivativealgorithmic;
            double relerror = abserror / ControlDerivativealgorithmic;
            double an_al = ControlDerivativeanalytical / ControlDerivativealgorithmic;
            double al_an = ControlDerivativealgorithmic / ControlDerivativeanalytical;
            ForwardControlDerivatives << stateIndex << "," << vIndex << "," << ControlDerivativeanalytical << "," << ControlDerivativealgorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
        }

        //fuel wrt control
        stateIndex = 8;
        for (size_t vIndex = 0; vIndex < 3; ++vIndex)
        {
            size_t GSADindex = num_states + 1 + vIndex;
            size_t Xindex = 1 + vIndex;
            double ControlDerivativealgorithmic = StateRight(stateIndex).getDerivative(GSADindex);
            double ControlDerivativeanalytical = Forward_dChemicalFuel_dThrottleComponents(vIndex);
            double abserror = ControlDerivativeanalytical - ControlDerivativealgorithmic;
            double relerror = abserror / ControlDerivativealgorithmic;
            double an_al = ControlDerivativeanalytical / ControlDerivativealgorithmic;
            double al_an = ControlDerivativealgorithmic / ControlDerivativeanalytical;
            ForwardControlDerivatives << stateIndex << "," << vIndex << "," << ControlDerivativeanalytical << "," << ControlDerivativealgorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
        }

        //electric propellant wrt control
        stateIndex = 9;
        for (size_t vIndex = 0; vIndex < 3; ++vIndex)
        {
            size_t GSADindex = num_states + 1 + vIndex;
            size_t Xindex = 1 + vIndex;
            double ControlDerivativealgorithmic = StateRight(stateIndex).getDerivative(GSADindex);
            double ControlDerivativeanalytical = Forward_dElectricPropellant_dThrottleComponents(vIndex);
            double abserror = ControlDerivativeanalytical - ControlDerivativealgorithmic;
            double relerror = abserror / ControlDerivativealgorithmic;
            double an_al = ControlDerivativeanalytical / ControlDerivativealgorithmic;
            double al_an = ControlDerivativealgorithmic / ControlDerivativeanalytical;
            ForwardControlDerivatives << stateIndex << "," << vIndex << "," << ControlDerivativeanalytical << "," << ControlDerivativealgorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
        }
    }
    ForwardControlDerivatives.close();

    EMTG::Phases::BackwardBoundedImpulseManeuver myBackwardBoundedImpulseManeuver(0,//journeyIndex
        0,//phaseIndex,
        0,//stageIndex,
        &TheUniverse[0],
        &mySpacecraft,
        &options);

    myBackwardBoundedImpulseManeuver.set_ThrustStepLength(ThrustStepLength);
    myBackwardBoundedImpulseManeuver.set_dThrustStepLength_dPropagationVariable(dThrustStepLength_dPropagationVariable);
    myBackwardBoundedImpulseManeuver.set_dImpulseEpoch_dPropagationVariable_pointer(dImpulseEpoch_dPropagationVariable);
    myBackwardBoundedImpulseManeuver.set_LaunchDate(LaunchDate);
    myBackwardBoundedImpulseManeuver.set_max_thrust(max_thrust);
    myBackwardBoundedImpulseManeuver.set_max_mass_flow_rate(max_mass_flow_rate);
    myBackwardBoundedImpulseManeuver.set_Isp(Isp);
    myBackwardBoundedImpulseManeuver.set_power(power);
    myBackwardBoundedImpulseManeuver.set_active_power(active_power);
    myBackwardBoundedImpulseManeuver.set_number_of_active_engines(number_of_active_engines);
    myBackwardBoundedImpulseManeuver.set_ThrottleLevel(ThrottleLevel);
    myBackwardBoundedImpulseManeuver.set_ThrottleLevelString(ThrottleLevelString);
    myBackwardBoundedImpulseManeuver.set_DeltaV(DeltaV);
    myBackwardBoundedImpulseManeuver.set_spacecraft_state_minus(StateRight);
    myBackwardBoundedImpulseManeuver.set_spacecraft_state_plus(StateLeft);
    myBackwardBoundedImpulseManeuver.set_MTM(BackwardMTM);
    myBackwardBoundedImpulseManeuver.set_DutyCycle(options.engine_duty_cycle);
    myBackwardBoundedImpulseManeuver.set_Control(ControlVector);
    myBackwardBoundedImpulseManeuver.set_Throttle(throttle);
    myBackwardBoundedImpulseManeuver.set_AccelerationModel(&mySpacecraftAccelerationModel);

    myBackwardBoundedImpulseManeuver.process_maneuver(true);


    //print STM entries
    std::ofstream BackwardMTMout("tests/BackwardBoundedImpulseManeuver.csv", std::ios::trunc);
    BackwardMTMout << "i, j, Analytical, Algorithmic, Absolute error, Relative error, an/al, al/an" << std::endl;
    for (size_t i = 0; i < num_states; ++i)
    {
        for (size_t j = 0; j < num_states + 1; ++j)
        {
            double MTMalgorithmic = StateRight(i).getDerivative(j + FirstState) / stateScales[j];
            double MTManalytical = BackwardMTM(i, j);
            double abserror = MTManalytical - MTMalgorithmic;
            double relerror = abserror / MTMalgorithmic;
            double an_al = MTManalytical / MTMalgorithmic;
            double al_an = MTMalgorithmic / MTManalytical;
            BackwardMTMout << i << "," << j << "," << MTManalytical << "," << MTMalgorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
        }
    }
    BackwardMTMout.close();

    std::ofstream BackwardControlDerivatives("tests/BackwardControlDerivatives.csv", std::ios::trunc);
    BackwardControlDerivatives << "stateIndex, vIndex, Analytical, Algorithmic, Absolute error, Relative error, an/al, al/an" << std::endl;
    EMTG::math::Matrix<double> Backward_dMassBeforeManeuver_dThrottleComponents = myBackwardBoundedImpulseManeuver.get_dMassBeforeManeuver_dThrottleComponents();
    EMTG::math::Matrix<double> Backward_dChemicalFuel_dThrottleComponents = myBackwardBoundedImpulseManeuver.get_dChemicalFuel_dThrottleComponents();
    EMTG::math::Matrix<double> Backward_dElectricPropellant_dThrottleComponents = myBackwardBoundedImpulseManeuver.get_dElectricPropellant_dThrottleComponents();

    {
        //mass wrt control
        size_t stateIndex = 6;
        for (size_t vIndex = 0; vIndex < 3; ++vIndex)
        {
            size_t GSADindex = num_states + 1 + vIndex;
            size_t Xindex = 1 + vIndex;
            double ControlDerivativealgorithmic = StateRight(stateIndex).getDerivative(GSADindex);
            double ControlDerivativeanalytical = Backward_dMassBeforeManeuver_dThrottleComponents(vIndex);
            double abserror = ControlDerivativeanalytical - ControlDerivativealgorithmic;
            double relerror = abserror / ControlDerivativealgorithmic;
            double an_al = ControlDerivativeanalytical / ControlDerivativealgorithmic;
            double al_an = ControlDerivativealgorithmic / ControlDerivativeanalytical;
            BackwardControlDerivatives << stateIndex << "," << vIndex << "," << ControlDerivativeanalytical << "," << ControlDerivativealgorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
        }

        //fuel wrt control
        stateIndex = 8;
        for (size_t vIndex = 0; vIndex < 3; ++vIndex)
        {
            size_t GSADindex = num_states + 1 + vIndex;
            size_t Xindex = 1 + vIndex;
            double ControlDerivativealgorithmic = StateRight(stateIndex).getDerivative(GSADindex);
            double ControlDerivativeanalytical = Backward_dChemicalFuel_dThrottleComponents(vIndex);
            double abserror = ControlDerivativeanalytical - ControlDerivativealgorithmic;
            double relerror = abserror / ControlDerivativealgorithmic;
            double an_al = ControlDerivativeanalytical / ControlDerivativealgorithmic;
            double al_an = ControlDerivativealgorithmic / ControlDerivativeanalytical;
            BackwardControlDerivatives << stateIndex << "," << vIndex << "," << ControlDerivativeanalytical << "," << ControlDerivativealgorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
        }

        //electric propellant wrt control
        stateIndex = 9;
        for (size_t vIndex = 0; vIndex < 3; ++vIndex)
        {
            size_t GSADindex = num_states + 1 + vIndex;
            size_t Xindex = 1 + vIndex;
            double ControlDerivativealgorithmic = StateRight(stateIndex).getDerivative(GSADindex);
            double ControlDerivativeanalytical = Backward_dElectricPropellant_dThrottleComponents(vIndex);
            double abserror = ControlDerivativeanalytical - ControlDerivativealgorithmic;
            double relerror = abserror / ControlDerivativealgorithmic;
            double an_al = ControlDerivativeanalytical / ControlDerivativealgorithmic;
            double al_an = ControlDerivativealgorithmic / ControlDerivativeanalytical;
            BackwardControlDerivatives << stateIndex << "," << vIndex << "," << ControlDerivativeanalytical << "," << ControlDerivativealgorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
        }
    }
    BackwardControlDerivatives.close();
}//end main