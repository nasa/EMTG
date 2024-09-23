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


#include "ThrottleSetting.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace EMTG
{
    namespace HardwareModels
    {
        //constructor for use with incoming data
        ThrottleSetting::ThrottleSetting(const std::string& ThrottleLevel,
            const double& Mdot,
            const double& BeamCurrent,
            const double& BeamVoltage,
            const double& Thrust,
            const double& Isp,
            const double& efficiency,
            const double& ThrusterPower,
            const double& Sharpness)
        {
            this->ThrottleLevelString = ThrottleLevel;
            try
            {
                this->ThrottleLevel = boost::lexical_cast<int>(boost::erase_first_copy(ThrottleLevel, "TL"));
            }
            catch (boost::bad_lexical_cast)
            {
                this->ThrottleLevel = 0; //this is better than a crash and we should not need this field for anything anyway
                //it only gets used for writing ephemeris files and even then we also check thrust values
            }
            this->Thrust = Thrust;
            this->Mdot = Mdot;
            this->BeamCurrent = BeamCurrent;
            this->BeamVoltage = BeamVoltage;
            this->Isp = Isp;
            this->efficiency = efficiency;
            this->ThrusterPower = ThrusterPower;
            this->Sharpness = Sharpness;
        }//end constructor

        //get functions
        double ThrottleSetting::getMdot() const
        {
            return this->Mdot;
        }

        double ThrottleSetting::getThrust() const
        {
            return this->Thrust;
        }

        double ThrottleSetting::getIsp() const
        {
            return this->Isp;
        }

        double ThrottleSetting::getEfficiency() const
        {
            return this->efficiency;
        }

        double ThrottleSetting::getThrusterPower() const
        {
            return this->ThrusterPower;
        }

        double ThrottleSetting::getBeamCurrent() const
        {
            return this->BeamCurrent;
        }

        double ThrottleSetting::getBeamVoltage() const
        {
            return this->BeamVoltage;
        }

        int ThrottleSetting::getThrottleLevel() const
        {
            return this->ThrottleLevel;
        }

        doubleType ThrottleSetting::getH() const
        {
            return this->H;
        }
        
        doubleType ThrottleSetting::getHeavisideThrust() const
        {
            return this->HThrust;
        }

        doubleType ThrottleSetting::getHeavisideMdot() const
        {
            return this->Hmdot;
        }

        doubleType ThrottleSetting::getHeavisideIsp() const
        {
            return this->HIsp;
        }


        doubleType ThrottleSetting::getHeavisideThrottleLevel() const
        {
            return this->HthrottleLevel;
        }

        doubleType ThrottleSetting::getHeavisideActivePower() const
        {
            return this->Hactivepower;
        }

        doubleType ThrottleSetting::getHeavisideBeamCurrent() const
        {
            return this->HBeamCurrent;
        }

        doubleType ThrottleSetting::getHeavisideBeamVoltage() const
        {
            return this->HBeamVoltage;
        }

        double ThrottleSetting::getHeavisideThrustDerivative() const
        {
            return this->dHthrustdP;
        }

        double ThrottleSetting::getHeavisideMdotDerivative() const
        {
            return this->dHmdotdP;
        }

        double ThrottleSetting::getHeavisideIspDerivative() const
        {
            return this->dHIspdP;
        }

        double ThrottleSetting::getHeavisideBeamVoltageDerivative() const
        {
            return this->dHBeamVoltagedP;
        }

        double ThrottleSetting::getHeavisideThrustCommandedMdotDerivative() const
        {
            return this->dHthrust_dcommanded_mdot;
        }

        double ThrottleSetting::getHeavisideMdotCommandedMdotDerivative() const
        {
            return this->dHmdot_dcommanded_mdot;
        }

        double ThrottleSetting::getHeavisideIspCommandedMdotDerivative() const
        {
            return this->dHIsp_dcommanded_mdot;
        }

        double ThrottleSetting::getHeavisideThrustCommandedVoltageDerivative() const
        {
            return this->dHthrust_dcommanded_voltage;
        }

        double ThrottleSetting::getHeavisideMdotCommandedVoltageDerivative() const
        {
            return this->dHmdot_dcommanded_voltage;
        }

        double ThrottleSetting::getHeavisideIspCommandedVoltageDerivative() const
        {
            return this->dHIsp_dcommanded_voltage;
        }
        
        double ThrottleSetting::getHeavisideIspInputPowerDerivative() const
        {
            return this->dHIsp_dinput_power;
        }
        double ThrottleSetting::getHeavisideMdotInputPowerDerivative() const
        {
            return this->dHmdot_dinput_power;
        }
        double ThrottleSetting::getHeavisideThrustInputPowerDerivative() const
        {
            return this->dHthrust_dinput_power;
        }

        //set functions
        void ThrottleSetting::setHeavisideThrustHeight(const double& inputHeavisideThrustHeight)
        {
            this->HeavisideThrustHeight = inputHeavisideThrustHeight;
        }

        void ThrottleSetting::setHeavisideMdotHeight(const double& inputHeavisideMdotHeight)
        {
            this->HeavisideMdotHeight = inputHeavisideMdotHeight;
        }

        void ThrottleSetting::setHeavisideIspHeight(const double& inputHeavisideIspHeight)
        {
            this->HeavisideIspHeight = inputHeavisideIspHeight;
        }

        void ThrottleSetting::setHeavisideThrottleLevelHeight(const int& inputThrottleLevelHeight)
        {
            this->HeavisideThrottleLevelHeight = inputThrottleLevelHeight;
        }

        void ThrottleSetting::setHeavisideActivePowerHeight(const double& inputHeavisideActivePowerHeight)
        {
            this->HeavisideActivePowerHeight = inputHeavisideActivePowerHeight;
        }

        void ThrottleSetting::setHeavisideBeamCurrentHeight(const double& inputHeavisideBeamCurrentHeight)
        {
            this->HeavisideBeamCurrentHeight = inputHeavisideBeamCurrentHeight;
        }

        void ThrottleSetting::setHeavisideBeamVoltageHeight(const double& inputHeavisideBeamVoltageHeight)
        {
            this->HeavisideBeamVoltageHeight = inputHeavisideBeamVoltageHeight;
        }

        void ThrottleSetting::setBeamVoltageOn(const double& inputBeamVoltageOn)
        {
            this->Voltage_on = inputBeamVoltageOn;
        }

        void ThrottleSetting::setBeamVoltageOff(const double& inputBeamVoltageOff)
        {
            this->Voltage_off = inputBeamVoltageOff;
        }

        void ThrottleSetting::setMdotOn(const double& inputMdotOn)
        {
            this->mdot_on = inputMdotOn;
        }

        void ThrottleSetting::setMdotOff(const double& inputMdotOff)
        {
            this->mdot_off = inputMdotOff;
        }
        
        void ThrottleSetting::setPowerOff(const double& inputPowerOff)
        {
            this->power_off = inputPowerOff;
        }
        
        void ThrottleSetting::setPowerOn(const double& inputPowerOn)
        {
            this->power_on = inputPowerOn;
        }
        
        void ThrottleSetting::setSharpness(const double& inputSharpness)
        {
            this->Sharpness = inputSharpness;
        }

        double ThrottleSetting::getSharpness() const
        {
            return this->Sharpness;
        }

        //function to calculate Heaviside values and derivatives
        void ThrottleSetting::calculateHeavisideValues1D(const doubleType& inputPower)
        {
            //compute the Heaviside function and its derivative
            doubleType Hswitch = math::absclip(-this->Sharpness * (inputPower - this->ThrusterPower), 600.0);
            doubleType expfun = exp(Hswitch);
            this->H = 1.0 / (1 + expfun);
            this->dHdP = this->Sharpness * (expfun * H * H)_GETVALUE;

            //Heaviside outputs
            this->HThrust = this->HeavisideThrustHeight * H;
            this->Hmdot = this->HeavisideMdotHeight * H;
            this->HIsp = this->HeavisideIspHeight * H;
            this->HthrottleLevel = this->HeavisideThrottleLevelHeight * H;
            this->Hactivepower = this->HeavisideActivePowerHeight * H;
            this->HBeamCurrent = this->HeavisideBeamCurrentHeight * H;
            this->HBeamVoltage = this->HeavisideBeamVoltageHeight * H;

            //Heaviside output derivatives
            this->dHthrustdP = this->HeavisideThrustHeight * dHdP;
            this->dHmdotdP = this->HeavisideMdotHeight * dHdP;
            this->dHIspdP = this->HeavisideIspHeight * dHdP;
            this->dHBeamVoltagedP = this->HeavisideBeamVoltageHeight * dHdP;
        }

        //function to calculate 2D Heaviside values and derivatives
        void ThrottleSetting::calculateHeavisideValues2D_mdot(const doubleType& commanded_mdot, const doubleType& input_power, const doubleType& max_mdot, const doubleType& max_power)
        {            
            //we assume in this function that we have already computed commanded_mdot and commanded_mdot by evaluating the upper and lower bounds on both mdot and V w.r.t. power
            this->mdot_on_switch        = math::absclip(-this->Sharpness * (commanded_mdot - this->mdot_on)/max_mdot, 600.0);
            this->mdot_off_switch       = math::absclip(-this->Sharpness * (commanded_mdot - this->mdot_off)/max_mdot, 600.0);
            this->power_on_switch       = math::absclip(-this->Sharpness * (input_power - this->power_on)/max_power, 600.0);
            this->power_off_switch      = math::absclip(-this->Sharpness * (input_power - this->power_off)/max_power, 600.0);

            this->expfun_mdot_on        = exp(this->mdot_on_switch);
            this->expfun_mdot_off       = exp(this->mdot_off_switch);
            this->expfun_power_on       = exp(this->power_on_switch);
            this->expfun_power_off      = exp(this->power_off_switch);

            this->H_commanded_Mdot_on  = 1.0 / (1 + this->expfun_mdot_on);
            this->H_commanded_Mdot_off = 1.0 / (1 + this->expfun_mdot_off);
            this->H_Power_on           = 1.0 / (1 + this->expfun_power_on);
            this->H_Power_off          = 1.0 / (1 + this->expfun_power_off);
            
            this->H_commanded_mdot    = this->H_commanded_Mdot_on - this->H_commanded_Mdot_off;
            this->H_power             = this->H_Power_on - this->H_Power_off;
            this->H                   = this->H_commanded_mdot * this->H_power;

            this->dH_commanded_mdot_dcommanded_mdot = (this->Sharpness / max_mdot * max_mdot * 1.1 * (this->expfun_mdot_on  * this->H_commanded_Mdot_on * this->H_commanded_Mdot_on - this->expfun_mdot_off  * this->H_commanded_Mdot_off * this->H_commanded_Mdot_off)) _GETVALUE;
            this->dH_power_dinput_power             = (this->Sharpness / max_power                 * (this->expfun_power_on * this->H_Power_on          * this->H_Power_on          - this->expfun_power_off * this->H_Power_off          * this->H_Power_off)) _GETVALUE;
            
            this->dH_dcommanded_mdot = (this->H_power) _GETVALUE          * this->dH_commanded_mdot_dcommanded_mdot;
            this->dH_dinput_power    = (this->H_commanded_mdot) _GETVALUE * this->dH_power_dinput_power;

            ////Heaviside outputs
            this->HThrust        = this->Thrust * this->H;
            this->Hmdot          = this->Mdot * this->H;
            this->HIsp           = this->Isp * this->H;
            this->HthrottleLevel = this->ThrottleLevel * this->H;
            this->Hactivepower   = this->ThrusterPower * this->H;
            this->HBeamVoltage   = this->BeamVoltage * this->H;

            ////Heaviside output derivatives
            ////these outputs do NOT include the derivatives with respect to control, which are chained in at the ThrottleTable level
            this->dHthrust_dcommanded_mdot = this->Thrust * this->dH_dcommanded_mdot;
            this->dHthrust_dinput_power    = this->Thrust * this->dH_dinput_power;
            this->dHmdot_dcommanded_mdot   = this->Mdot * this->dH_dcommanded_mdot;
            this->dHmdot_dinput_power      = this->Mdot * this->dH_dinput_power;
            this->dHIsp_dcommanded_mdot    = this->Isp * this->dH_dcommanded_mdot;
            this->dHIsp_dinput_power       = this->Isp * this->dH_dinput_power;
        }
        
        //function to calculate 2D Heaviside values and derivatives
        void ThrottleSetting::calculateHeavisideValues2D_voltage(const doubleType& commanded_voltage, const doubleType& input_power, const doubleType& max_voltage, const doubleType& max_power)
        {            
            //we assume in this function that we have already computed commanded_mdot and commanded_voltage by evaluating the upper and lower bounds on both mdot and V w.r.t. power
            this->voltage_on_switch        = math::absclip(-this->Sharpness * (commanded_voltage - this->Voltage_on)/max_voltage, 600.0);
            this->voltage_off_switch       = math::absclip(-this->Sharpness * (commanded_voltage - this->Voltage_off)/max_voltage, 600.0);
            this->power_on_switch          = math::absclip(-this->Sharpness * (input_power - this->power_on)/max_power, 600.0);
            this->power_off_switch         = math::absclip(-this->Sharpness * (input_power - this->power_off)/max_power, 600.0);

            this->expfun_voltage_on        = exp(this->voltage_on_switch);
            this->expfun_voltage_off       = exp(this->voltage_off_switch);
            this->expfun_power_on          = exp(this->power_on_switch);
            this->expfun_power_off         = exp(this->power_off_switch);

            this->H_commanded_Voltage_on  = 1.0 / (1 + this->expfun_voltage_on);
            this->H_commanded_Voltage_off = 1.0 / (1 + this->expfun_voltage_off);
            this->H_Power_on              = 1.0 / (1 + this->expfun_power_on);
            this->H_Power_off             = 1.0 / (1 + this->expfun_power_off);
            
            this->H_commanded_voltage = this->H_commanded_Voltage_on - this->H_commanded_Voltage_off;
            this->H_power             = this->H_Power_on - this->H_Power_off;
            this->H                   = this->H_commanded_voltage * this->H_power;

            this->dH_commanded_voltage_dcommanded_voltage = (this->Sharpness / max_voltage * max_voltage * 1.1 * (this->expfun_voltage_on        * this->H_commanded_Voltage_on * this->H_commanded_Voltage_on - this->expfun_voltage_off * this->H_commanded_Voltage_off * this->H_commanded_Voltage_off)) _GETVALUE;
            this->dH_power_dinput_power                   = (this->Sharpness / max_power                       * (this->expfun_power_on          * this->H_Power_on             * this->H_Power_on             - this->expfun_power_off   * this->H_Power_off             * this->H_Power_off)) _GETVALUE;
            
            this->dH_dcommanded_voltage = (this->H_power) _GETVALUE             * this->dH_commanded_voltage_dcommanded_voltage;
            this->dH_dinput_power       = (this->H_commanded_voltage) _GETVALUE * this->dH_power_dinput_power;

            ////Heaviside outputs
            this->HThrust        = this->Thrust * this->H;
            this->Hmdot          = this->Mdot * this->H;
            this->HIsp           = this->Isp * this->H;
            this->HthrottleLevel = this->ThrottleLevel * this->H;
            this->Hactivepower   = this->ThrusterPower * this->H;
            this->HBeamVoltage   = this->BeamVoltage * this->H;

            ////Heaviside output derivatives
            ////these outputs do NOT include the derivatives with respect to control, which are chained in at the ThrottleTable level
            this->dHthrust_dcommanded_voltage = this->Thrust * this->dH_dcommanded_voltage;
            this->dHthrust_dinput_power       = this->Thrust * this->dH_dinput_power;
            this->dHmdot_dcommanded_voltage   = this->Mdot * this->dH_dcommanded_voltage;
            this->dHmdot_dinput_power         = this->Mdot * this->dH_dinput_power;
            this->dHIsp_dcommanded_voltage    = this->Isp * this->dH_dcommanded_voltage;
            this->dHIsp_dinput_power          = this->Isp * this->dH_dinput_power;
        }

        //print function
        void ThrottleSetting::printThrottleSetting(std::ofstream& ThrusterOutputFile, const bool& PrintTransitions) const
        {
            //#Throttle level, Mass flow rate(mg / s), Beam Current(A), Beam Voltage(V), Thrust(mN), Isp(s), Efficiency, Thruster input power(kW)
            ThrusterOutputFile << this->ThrottleLevelString << ",";
            ThrusterOutputFile << this->Mdot << ",";
            ThrusterOutputFile << this->BeamCurrent << ",";
            ThrusterOutputFile << this->BeamVoltage << ",";
            ThrusterOutputFile << this->Thrust << ",";
            ThrusterOutputFile << this->Isp << ",";
            ThrusterOutputFile << this->efficiency << ",";
            ThrusterOutputFile << this->ThrusterPower;

            if (PrintTransitions)
            {
                ThrusterOutputFile << ",";
                ThrusterOutputFile << this->mdot_on << ",";
                ThrusterOutputFile << this->mdot_off << ",";
                ThrusterOutputFile << this->Voltage_on << ",";
                ThrusterOutputFile << this->Voltage_off;
            }
            
            ThrusterOutputFile << std::endl;
        }
    }//end namespace HardwareModels
}//end namespace EMTG
