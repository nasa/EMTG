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

//throttle setting class for EMTG
//overloaded for use with algorithmic differentiation
//Jacob Englander 5-14-2016

#pragma once

#include "doubleType.h"

#include "EMTG_math.h"
#include <limits>
#include <iostream>
#include <fstream>

namespace EMTG
{
    namespace HardwareModels
    {
        class ThrottleSetting 
        {
        public:
            //default constructor - doesn't do anything, never called
            ThrottleSetting() {};
            //constructor for use with incoming data
            ThrottleSetting(const std::string& ThrottleLevel,
                const double& Mdot,
                const double& BeamCurrent,
                const double& BeamVoltage,
                const double& Thrust,
                const double& Isp,
                const double& efficiency,
                const double& ThrusterPower,
                const double& Sharpness);

            //destructor, doesn't need to do anything
            ~ThrottleSetting() {};

            //get functions
            std::string getThrottleLevelString() const { return this->ThrottleLevelString; }
            double getMdot() const;
            double getThrust() const;
            double getIsp() const;
            double getEfficiency() const;
            double getThrusterPower() const;
            double getBeamCurrent() const;
            double getBeamVoltage() const;
            int getThrottleLevel() const;
            doubleType getH() const;
            doubleType getHeavisideThrust() const;
            doubleType getHeavisideMdot() const;
            doubleType getHeavisideIsp() const;
            doubleType getHeavisideThrottleLevel() const;
            doubleType getHeavisideActivePower() const;
            doubleType getHeavisideBeamCurrent() const;
            doubleType getHeavisideBeamVoltage() const;
            double getHeavisideThrustDerivative() const;
            double getHeavisideMdotDerivative() const;
            double getHeavisideIspDerivative() const;
            double getHeavisideBeamVoltageDerivative() const;
            double getHeavisideThrustCommandedMdotDerivative() const;
            double getHeavisideMdotCommandedMdotDerivative() const;
            double getHeavisideIspCommandedMdotDerivative() const;
            double getHeavisideThrustCommandedVoltageDerivative() const;
            double getHeavisideMdotCommandedVoltageDerivative() const;
            double getHeavisideIspCommandedVoltageDerivative() const;
            double getHeavisideThrustInputPowerDerivative() const;
            double getHeavisideMdotInputPowerDerivative() const;
            double getHeavisideIspInputPowerDerivative() const;

            //set functions
            void setHeavisideThrustHeight(const double& inputHeavisideThrustHeight);
            void setHeavisideMdotHeight(const double& inputHeavisideMdotHeight);
            void setHeavisideIspHeight(const double& inputHeavisideIspHeight);
            void setHeavisideThrottleLevelHeight(const int& inputThrottleLevelHeight);
            void setHeavisideActivePowerHeight(const double& inputHeavisideActivePowerHeight);
            void setHeavisideBeamCurrentHeight(const double& inputHeavisideBeamCurrentHeight);
            void setHeavisideBeamVoltageHeight(const double& inputHeavisideBeamVoltageHeight);
            void setBeamVoltageOn(const double& inputBeamVoltageOn);
            void setBeamVoltageOff(const double& inputBeamVoltageOff);
            void setMdotOn(const double& inputMdotOn);
            void setMdotOff(const double& inputMdotOff);
            void setPowerOff(const double& inputPowerOff);
            void setPowerOn(const double& inputPowerOn);
            void setSharpness(const double& inputSharpness);
            double getSharpness() const;
            //void setHeavisideThrustHeight_mdot_axis(const double& inputHeavisideThrustHeight_mdot_axis);
            //void setHeavisideMdotHeight_mdot_axis(const double& inputHeavisideMdotHeight_mdot_axis);
            //void setHeavisideIspHeight_mdot_axis(const double& inputHeavisideIspHeight_mdot_axis);
            //void setHeavisideActivePowerHeight_mdot_axis(const double& inputHeavisideActivePowerHeight_mdot_axis);
            //void setHeavisideThrottleLevelHeight_mdot_axis(const double& inputHeavisideThrottleLevelHeight_mdot_axis);
            //void setHeavisideThrustHeight_voltage_axis(const double& inputHeavisideThrustHeight_voltage_axis);
            //void setHeavisideMdotHeight_voltage_axis(const double& inputHeavisideMdotHeight_voltage_axis);
            //void setHeavisideIspHeight_voltage_axis(const double& inputHeavisideIspHeight_voltage_axis);
            //void setHeavisideActivePowerHeight_voltage_axis(const double& inputHeavisideActivePowerHeight_voltage_axis);
            //void setHeavisideThrottleLevelHeight_voltage_axis(const double& inputHeavisideThrottleLevelHeight_voltage_axis);

            //function to calculate Heaviside values and derivatives
            void calculateHeavisideValues1D(const doubleType& inputPower);
            void calculateHeavisideValues2D_voltage(const doubleType& commanded_voltage, const doubleType& input_power, const doubleType& max_voltage, const doubleType& max_power);
            void calculateHeavisideValues2D_mdot(const doubleType& commanded_mdot, const doubleType& input_power, const doubleType& max_mdot, const doubleType& max_power);

            //print function
            void printThrottleSetting(std::ofstream& ThrusterOutputFile, const bool& PrintTransitions) const;

            // //comparator
            // bool operator<(ThrottleSetting& RHS)
            // {
            //     double A = this->ThrusterPower;
            //     double B = RHS.getThrusterPower();
            //     return A < B ? true : false;
            //     //return this->ThrusterPower < RHS.getThrusterPower() ? true : false;
            // };

        private:

            //fields for storing input data
            double Thrust, Mdot, Isp, efficiency, ThrusterPower, BeamCurrent, BeamVoltage;
            int ThrottleLevel;
            std::string ThrottleLevelString;

            //fields for Heaviside function
            doubleType Hswitch;
            double Sharpness;
            double HeavisideThrustHeight, HeavisideMdotHeight, HeavisideIspHeight, HeavisideActivePowerHeight, HeavisideBeamCurrentHeight, HeavisideBeamVoltageHeight;
            int HeavisideThrottleLevelHeight;
            doubleType H;
            double dHdP;
            doubleType HThrust, Hmdot, HIsp, HthrottleLevel, Hactivepower, HBeamCurrent, HBeamVoltage;
            double dHthrustdP, dHmdotdP, dHIspdP, dHBeamVoltagedP;

            //2D Heaviside
            doubleType H_commanded_mdot, H_commanded_voltage;
            double Voltage_on, Voltage_off, mdot_on, mdot_off;
            double dH_commanded_mdot_dcommanded_mdot, dH_commanded_voltage_dcommanded_voltage;
            double dH_dcommanded_mdot, dH_dcommanded_voltage;
            double dHthrust_dcommanded_mdot, dHthrust_dcommanded_voltage, dHmdot_dcommanded_mdot, dHmdot_dcommanded_voltage, dHIsp_dcommanded_mdot, dHIsp_dcommanded_voltage;
            
            doubleType mdot_on_switch,mdot_off_switch,voltage_on_switch,voltage_off_switch;
            doubleType expfun_mdot_on,expfun_mdot_off,expfun_voltage_on,expfun_voltage_off;
            doubleType H_commanded_Mdot_on,H_commanded_Mdot_off,H_commanded_Voltage_on,H_commanded_Voltage_off;
            
            doubleType power_on_switch,H_Power_on,H_Power_off,H_power, power_off_switch;
            doubleType expfun_power_on, expfun_power_off;        
            double power_off, power_on;
            double dHthrust_dinput_power,dHIsp_dinput_power,dHmdot_dinput_power;
            double dH_power_on_dinput_power,dH_power_off_dinput_power,dH_power_dinput_power;
            double dH_dinput_power;
                        
        };//end class ThrottleSetting declaration
    }//end namespace HardwareModels
}//end namespace EMTG