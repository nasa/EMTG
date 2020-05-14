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

//options class for power systems
//Jacob Englander 12-30-2016

#ifndef EMTG_POWERSYSTEMOPTIONS
#define EMTG_POWERSYSTEMOPTIONS

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "EMTG_enums.h"
#include "file_utilities.h"

namespace EMTG
{
    namespace HardwareModels
    {
        class PowerSystemOptions
        {
        public:
            //constructor
            PowerSystemOptions();
            PowerSystemOptions(const std::string& linestring);

            //destructor
            ~PowerSystemOptions() {};

            //methods
            void parse_input_line(const std::string& linestring);
            void write_output_line(std::ofstream& outputfile);

            //get/set
            inline std::string getName() const {return this->name; }
            inline SpacecraftBusPowerType getBusPowerType() const { return this->Spacecraft_Bus_Power_Type; }
            inline SpacecraftPowerSupplyType getPowerSupplyType() const { return this->Spacecraft_Power_Supply_Type; }
            inline SpacecraftPowerSupplyCurveType getPowerSupplyCurveType() const { return this->Spacecraft_Power_Supply_Curve_Type; }
            inline double getDecayRate() const { return this->decay_rate; }
            inline double getPowerSystemDecayRefEpoch() const { return this->power_system_decay_reference_epoch; };
            inline double getP0() const { return this->P0; }
            inline double getPowerMargin() const { return this->PowerMargin; }
            inline std::vector<double> getGammaVector() const { return this->gamma; }
            inline std::vector<double> getBusCoefficientVector() const { return this->bus_coefficients; }
            inline double getmass_per_kw() const { return this->mass_per_kW; }
            inline double getGamma(const size_t& p) const {return this->gamma[p]; }
            inline double getBusCoefficient(const size_t& p) const { return this->bus_coefficients[p]; }

            void setName(const std::string& name) { this->name = name; }
            inline void setP0(const double& newP0) { this->P0 = newP0; }
            void setBusPowerType(const SpacecraftBusPowerType& input);
            void setPowerSupplyType(const SpacecraftPowerSupplyType& input);
            void setPowerSupplyCurveType(const SpacecraftPowerSupplyCurveType& input);
            void setDecayRate(const double& input);
            void setPowerSystemDecayRefEpoch(const double & input);
            void setPowerMargin(const double& input);
            void setGammaVector(const std::vector<double>& input);
            void setBusCoefficientVector(const std::vector<double>& input);
            void setmass_per_kw(const double& input);
            void setGamma(const size_t& p, const double& input);
            void setBusCoefficient(const size_t& p, const double& input);

            //fields
        protected:
            std::string name;
            SpacecraftBusPowerType Spacecraft_Bus_Power_Type;
            SpacecraftPowerSupplyType Spacecraft_Power_Supply_Type;
            SpacecraftPowerSupplyCurveType Spacecraft_Power_Supply_Curve_Type;
            double P0;
            double decay_rate;
            double power_system_decay_reference_epoch;
            double PowerMargin;
            std::vector<double> bus_coefficients;
            std::vector<double> gamma;
            double mass_per_kW;
        };

        class PowerSystemOptionsLibrary
        {
        public:
            //constructor
            PowerSystemOptionsLibrary();
            PowerSystemOptionsLibrary(const std::string& filename);
            PowerSystemOptionsLibrary(std::ifstream& inputfile);
            PowerSystemOptionsLibrary(const std::vector<PowerSystemOptions>& PowerSystemOptionsVector);

            //destructor
            ~PowerSystemOptionsLibrary() {};//destructor does not need to do anything

            //methods
            void parse_input_file(const std::string& filename);
            void parse_input_file(std::ifstream& inputfile);
            void write_output_file(const std::string& filename, const bool& newFile = true);


            //get functions
            PowerSystemOptions getPowerSystem(const std::string& name);
            SpacecraftBusPowerType getBusPowerType(const std::string& name);
            SpacecraftPowerSupplyType getPowerSupplyType(const std::string& name);
            SpacecraftPowerSupplyCurveType getPowerSupplyCurveType(const std::string& name);
            double getDecayRate(const std::string& name);
            double getPowerSystemDecayRefEpoch(const std::string & name);
            double getP0(const std::string& name);
            double getPowerMargin(const std::string& name);
            std::vector<double> getGamma(const std::string& name);
            std::vector<double> getBusCoefficients(const std::string& name);
            double getmass_per_kw(const std::string& name);
			std::vector<std::string> getAllKeys();

        private:
            //private methods
            size_t getIndex(const std::string& name);

            //fields
            std::vector<PowerSystemOptions> PowerSystems;
            std::string currentKey;
            size_t currentIndex;
        };
    }//end namespace HardwareModels
}//end namespace EMTG

#endif //EMTG_POWERSYSTEMOPTIONS