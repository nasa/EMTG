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

//throttle table class for EMTG
//overloaded for use with algorithmic differentiation
//Jacob Englander 5-14-2016


#pragma once

#include <string>
#include <vector>

#include "doubleType.h"

#include "ThrottleSetting.h"

#include <boost/algorithm/string.hpp>

#include <fstream>
#include <iostream>

namespace EMTG
{
    namespace HardwareModels
    {
        enum ThrottleSetDef { HighThrust, HighIsp, LowMdot, HighMdot, HighEfficiency, FullSet };

        class ThrottleTable
        {
        public:
            //default constructor - doesn't do anything, never called
            ThrottleTable() : control_type_2d(false) { };
            //constructor for use with incoming data
            ThrottleTable(const std::string& inputfilename, const double& Sharpness);

            //destructor, doesn't need to do anything
            ~ThrottleTable() {};

            //function to parse a throttle table file
            void ParseThrottleTableFile(const std::string& inputfilename);

            //function to create max-thrust and max-Isp sets
            void CreateThrottleSets();

            //function to find nearest neighbor to the left and nearest neighbor below each throttle setting in mdot and voltage, and then set the appropriate mdot and voltage heights
            //void FindBottomLeftNeighbors();

            //function to find the transition voltage and mdot values for each throttle setting
            void FindTransitionValues();
            
            //function to set 2d control variable
            void set2DControlType(const bool & control_type);

            //function to calculate thruster performance at a given power level
            void CalculateThrusterPerformance1D(const doubleType& inputPower, const ThrottleSetDef& PreferredThrottleSet);
            void CalculateThrusterPerformance2D(const doubleType& inputPower, const doubleType& u_command);

            //get functions
            double getPPUminpower() const;
            double getPPUmaxpower() const;
            double getThrottleSetmaxpower(const ThrottleSetDef& PreferredThrottleSet);
            double getThrottleSetminpower(const ThrottleSetDef& PreferredThrottleSet);
            double getPPUefficiency() const;
            doubleType getHeavisideThrust() const;
            doubleType getHeavisideMdot() const;
            doubleType getHeavisideIsp() const;
            doubleType getHeavisideActivePower() const;
            double getHeavisideThrustDerivative() const;
            double getHeavisideMdotDerivative() const;
            double getHeavisideIspDerivative() const;
            int getThrottleLevel() const;
            std::string getThrottleLevelString() const { return this->ThrottleLevelString; }
            double getdHthrust_dcommand() const;
            double getdHmdot_dcommand() const;
            double getdHIsp_dcommand() const;

            //polynomial methods
            void CalculateThrusterPerformance1D_Polynomial(const doubleType& inputPower, const ThrottleSetDef& PreferredThrottleSet);
            void CalculateThrusterPerformance2D_Polynomial(const doubleType& inputPower, const doubleType& u_Mdot);
            doubleType get1D_poly_Thrust() const;
            doubleType get1D_poly_Mdot() const;
            doubleType get2D_poly_Thrust() const;
            doubleType get2D_poly_Mdot() const;
            double get2D_poly_dThrust_dP() const;
            double get2D_poly_dThrust_duMdot() const;
            double get2D_poly_dMdot_dP() const;
            double get2D_poly_dMdot_duMdot() const;
            double get1D_poly_dThrust_dP() const;
            double get1D_poly_dMdot_dP() const;
            double getd2HighThrust_1D_Mdot_dP2() const;
            double getd2HighIsp_1D_Mdot_dP2() const;

            //print function
            void printThrottleSets(const std::string& ThrottleOutputFileName) const;

        private:
            //vectors of throttle settings
            std::vector< ThrottleSetting > ThrottleSettings;
            std::vector< ThrottleSetting > HighThrustSet;
            std::vector< ThrottleSetting > HighIspSet;
            std::vector< ThrottleSetting > LowMdotSet;
            std::vector< ThrottleSetting > HighMdotSet;
            std::vector< ThrottleSetting > HighEfficiencySet;
            std::vector< ThrottleSetting > FullSet;
            
            std::vector<double> Mdot_choices;
            std::vector<double> Voltage_choices;

            //PPU properties
            double PPUefficiency;
            double PPUminpower;
            double PPUmaxpower;
            double Sharpness;
            
            int control_type_2d;

            //thruster performance characteristics
            std::string ThrottleLevelString;
            doubleType ThrusterInputPower;
            doubleType HThrust, Hmdot, HIsp, HThrottleLevel, HactivePower, HBeamCurrent, HBeamVoltage;
            double Mdot_max,Voltage_max,Power_max;
            double dHthrustdP, dHmdotdP, dHIspdP, dHBeamVoltagedP;
            double dHthrust_dcommand, dHmdot_dcommand, dHIsp_dcommand;

            //polynomials
            double high_thrust_Thrust_coeff[5];
            double high_thrust_Mdot_coeff[5];
            double high_Isp_Thrust_coeff[5];
            double high_Isp_Mdot_coeff[5];
            double poly2D_coefficients[5][5];
            doubleType polyThrust_1D;
            doubleType polyMdot_1D;
            double poly_dThrust_dP_1D;
            double poly_dMdotdP_1D;
            doubleType polyThrust_2D;
            doubleType polyMdot_2D;
            double poly_dThrust_dP_2D;
            double poly_dThrust_du_Mdot_2D;
            double poly_dMdot_dP_2D;
            double poly_dMdot_du_Mdot_2D;
            double d2poly_dMdot_dP2_2D;
            double d2HighThrust_1D_Mdot_dP2;
            double d2HighIsp_1D_Mdot_dP2;
            double dHthrust_du_command;
            double dHmdot_du_command;
            double dHIsp_du_command;
        };//end class ThrottleTable
    }//end namespace PropulsionSystem
}//end namespace EMTG