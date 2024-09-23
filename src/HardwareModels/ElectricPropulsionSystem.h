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

//electric propulsion system class
//Jacob Englander 11/2/2016

#pragma once

#include "PropulsionSystem.h"
#include "ThrottleTable.h"

#include <iostream>

namespace EMTG
{
    namespace HardwareModels
    {
        class ElectricPropulsionSystem : public PropulsionSystem
        {
        public:
            //constructor
            ElectricPropulsionSystem() {};
            ElectricPropulsionSystem(const PropulsionSystemOptions& propulsionsystemoptions);

            //destructor
            ~ElectricPropulsionSystem() {};

            //get
            inline doubleType getMassFlowRate() const { return this->MassFlowRate; }
            inline double getdTdP() const { return this->dTdP; }
            inline double getdTdu_command() const { return this->dTdu_command; }
            inline double getdMassFlowRatedP() const { return this->dMassFlowRatedP; }
            inline double getdMassFlowRatedu_command() const { return this->dMassFlowRatedu_command; }
            inline double getdIspdP() const { return this->dIspdP; }
            inline double getdIspdu_command() const { return this->dIspdP; }
            inline doubleType getActivePower() const { return this->ActivePower; }
            inline size_t getNumberOfActiveThrusters() const { return (size_t)round(this->smoothed_number_of_active_thrusters _GETVALUE); }
            inline size_t getThrottleLevel() const { return this->ThrottleLevel; }
            inline std::string getThrottleLevelString() const { return this->ThrottleLevelString; }
            

            //evaluate
            void computeSystemMass();
            void computeThrusterPerformance(const doubleType& InputPower, const double& DutyCycle, const ThrottleLogic& myThrottleLogic);
            void computeThrusterPerformance(const doubleType& InputPower, const double& DutyCycle, const ThrottleLogic& myThrottleLogic, const doubleType& u_command);

        protected:
            //methods
            void compute_number_of_active_thrusters(const doubleType& InputPower, const ThrottleLogic& myThrottleLogic);

            //fields
            doubleType MassFlowRate;
            double dTdP;
            double dTdu_command;
            double dMassFlowRatedP;
            double dMassFlowRatedu_command;
            double dIspdP;
            double dIspdu_command;
            doubleType ActivePower;
            doubleType smoothed_number_of_active_thrusters;
            double dNdP;
            doubleType PowerPerThruster;
            double dPowerPerThruster_dP;
            ThrottleTable myThrottleTable;
            int ThrottleLevel;
            std::string ThrottleLevelString;
        };
    }//end namespace HardwareModels
}//end namespace EMTG