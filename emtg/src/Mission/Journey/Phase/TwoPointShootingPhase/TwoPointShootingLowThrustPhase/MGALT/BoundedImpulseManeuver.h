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

//maneuver model for MGALT
//Jacob Englander 6/25/2017

#pragma once

#include "doubleType.h"

#include "Spacecraft.h"
#include "EMTG_Matrix.h"
#include "universe.h"
#include "missionoptions.h"
#include "SpacecraftAccelerationModel.h"




namespace EMTG 
{ 
    namespace Phases
    {
        class BoundedImpulseManeuver
        {
        public:
            //constructor
            BoundedImpulseManeuver();
            BoundedImpulseManeuver(const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);


            //set
            inline void set_ThrustStepLength(doubleType& ThrustStepLength) { this->ThrustStepLength_pointer = &ThrustStepLength; }
            inline void set_dThrustStepLength_dPropagationVariable(double& dThrustStepLength_dPropagationVariable) { this->dThrustStepLength_dPropagationVariable_pointer = &dThrustStepLength_dPropagationVariable; }
            inline void set_dImpulseEpoch_dPropagationVariable_pointer(double &dImpulseEpoch_dPropagationVariable_pointer) { this->dImpulseEpoch_dPropagationVariable_pointer = &dImpulseEpoch_dPropagationVariable_pointer; }
            inline void set_DutyCycle(const double& DutyCycle) { this->DutyCycle = DutyCycle; }
            inline void set_Control(math::Matrix<doubleType>& Control) { this->Control_pointer = &Control; }
            inline void set_Throttle(doubleType& Throttle) { this->Throttle_pointer = &Throttle; }
            inline void set_LaunchDate(doubleType& LaunchDate) { this->LaunchDate_pointer = &LaunchDate; }
            inline void set_max_thrust(doubleType& max_thrust) { this->max_thrust_pointer = &max_thrust; }
            inline void set_max_mass_flow_rate(doubleType& max_mass_flow_rate) { this->max_mass_flow_rate_pointer = &max_mass_flow_rate; }
            inline void set_Isp(doubleType& Isp) { this->Isp_pointer = &Isp; }
            inline void set_power(doubleType& power) { this->power_pointer = &power; }
            inline void set_active_power(doubleType& active_power) { this->active_power_pointer = &active_power; }
            inline void set_number_of_active_engines(size_t& number_of_active_engines) { this->number_of_active_engines_pointer = &number_of_active_engines; }
            inline void set_ThrottleLevel(int& ThrottleLevel) { this->ThrottleLevel_pointer = &ThrottleLevel; }
            inline void set_ThrottleLevelString(std::string& ThrottleLevelString) { this->ThrottleLevelString_pointer = &ThrottleLevelString; }
            inline void set_DeltaV(math::Matrix<doubleType>& DeltaV) { this->DeltaV_pointer = &DeltaV; }
            inline void set_spacecraft_state_minus(math::Matrix<doubleType>& spacecraft_state_minus) { this->spacecraft_state_minus_pointer = &spacecraft_state_minus; }
            inline void set_spacecraft_state_plus(math::Matrix<doubleType>& spacecraft_state_plus) { this->spacecraft_state_plus_pointer = &spacecraft_state_plus; }
            inline void set_MTM(math::Matrix<double>& MTM) { this->MTM_pointer = &MTM; }
            inline void set_AccelerationModel(Astrodynamics::SpacecraftAccelerationModel* phaseSpacecraftAccelerationModel) { this->mySpacecraftAccelerationModel = phaseSpacecraftAccelerationModel; }

            inline doubleType get_dvmax() const { return this->dvmax; }
            inline math::Matrix<double> get_dChemicalFuel_dThrottleComponents() const { return this->dChemicalFuel_dThrottleComponents; }
            inline math::Matrix<double> get_dElectricPropellant_dThrottleComponents() const { return this->dElectricPropellant_dThrottleComponents; }

            //methods
            virtual void process_maneuver(const bool& needMTM) = 0;
            void compute_thruster_performance(const bool& needDerivatives);


            //fields

            //pointers to things in the phase
        protected:
            size_t journeyIndex;
            size_t phaseIndex;
            size_t stageIndex;
            double DutyCycle;
            JourneyOptions* myJourneyOptions;
            Astrodynamics::universe* myUniverse;
            HardwareModels::Spacecraft* mySpacecraft;
            missionoptions* myOptions;

            bool hasPerturbations;

            Astrodynamics::SpacecraftAccelerationModel* mySpacecraftAccelerationModel;

            doubleType* LaunchDate_pointer;
            doubleType* max_thrust_pointer;
            doubleType* max_mass_flow_rate_pointer;
            doubleType* Isp_pointer;
            doubleType* power_pointer;
            doubleType* active_power_pointer;
            size_t* number_of_active_engines_pointer;
            int* ThrottleLevel_pointer;
            std::string* ThrottleLevelString_pointer;
            math::Matrix<doubleType>* Control_pointer;
            doubleType* Throttle_pointer;
            math::Matrix<doubleType>* DeltaV_pointer;
            math::Matrix<doubleType>* spacecraft_state_minus_pointer;
            math::Matrix<doubleType>* spacecraft_state_plus_pointer;
            math::Matrix<double>* MTM_pointer;
            doubleType* ThrustStepLength_pointer;
            double* dThrustStepLength_dPropagationVariable_pointer;
            double* dImpulseEpoch_dPropagationVariable_pointer;

            doubleType dvmax;

            doubleType electric_propellant_used;
            doubleType ACS_fuel_used;


            math::Matrix<double> dChemicalFuel_dThrottleComponents;
            math::Matrix<double> dElectricPropellant_dThrottleComponents;

            //internal things
            bool SunIsCentralBody;
            doubleType ManeuverEpoch;
            double dTdP;
            double dmdotdP;
            double dT_du_command;
            double dmdot_du_command;
            double dPdr_sun;
            double dPdt;

            math::Matrix<doubleType> NaturalAcceleration; //3x1
            math::Matrix<doubleType> dNaturalAcceleration_dt;
            math::Matrix<doubleType> dNaturalAcceleration_dState; //3x7
            math::Matrix<doubleType> NaturalDeltaV; //3x1
            math::Matrix<doubleType> dNaturalDeltaV_dPreviousTime;
            math::Matrix<doubleType> dNaturalDeltaV_dPropagationTime;
            math::Matrix<doubleType> dNaturalDeltaV_dState; //3x7

            math::Matrix<double> dR_sc_Sun_dR_CB_sun;
            math::Matrix<double> dR_sc_Sun_dR_sc_CB;
            math::Matrix<double> dP_dr_sc_CB;
            math::Matrix<double> dP_dr_CB_sun;

            math::Matrix<doubleType> R_sc_Sun;
            math::Matrix<doubleType> R_CB_Sun;
            math::Matrix<double> dR_CB_Sun_dt;
            math::Matrix<doubleType> R_sc_CB;
            math::Matrix<doubleType>* state_at_natural_perturbation;
        };
    } //close namespace Astrodynamics
} //close namespace EMTG