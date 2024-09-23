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
//Jacob Englander 11-16-2016

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "PowerSystemOptions.h"

#include <boost/algorithm/string.hpp>

namespace EMTG
{
    namespace HardwareModels
    {
        PowerSystemOptions::PowerSystemOptions() : 
            name("default_power_system"),
            Spacecraft_Bus_Power_Type(SpacecraftBusPowerType::TypeA_Quadratic),
            Spacecraft_Power_Supply_Curve_Type(SpacecraftPowerSupplyCurveType::Sauer),
            Spacecraft_Power_Supply_Type(SpacecraftPowerSupplyType::Solar),
            P0(10.0),
            mass_per_kW(10.0),
            decay_rate(0.01),
            PowerMargin(0.0),
            bus_coefficients(std::vector<double>(3, 0.0)),
            gamma(std::vector<double>(7, 0.0))
        {
            this->gamma[0] = 1.0;
        }

        PowerSystemOptions::PowerSystemOptions(const std::string& linestring)
        {
            this->parse_input_line(linestring);
        }

        void PowerSystemOptions::parse_input_line(const std::string& linestring)
        {
            std::vector<std::string> linecell;
            boost::split(linecell, linestring, boost::is_any_of(" ,"), boost::token_compress_on);

            size_t cellIndex = 0;

            this->name = linecell[cellIndex++];
            this->Spacecraft_Power_Supply_Type = (SpacecraftPowerSupplyType)std::stoi(linecell[cellIndex++]);
            this->Spacecraft_Power_Supply_Curve_Type = (SpacecraftPowerSupplyCurveType)std::stoi(linecell[cellIndex++]);
            this->Spacecraft_Bus_Power_Type = (SpacecraftBusPowerType)std::stoi(linecell[cellIndex++]);
            this->P0 = std::stod(linecell[cellIndex++]);
            this->mass_per_kW = std::stod(linecell[cellIndex++]);
            this->decay_rate = std::stod(linecell[cellIndex++]);

            this->gamma.clear();
            for (size_t k = 0; k < 7; ++k)
                this->gamma.push_back(std::stod(linecell[cellIndex++]));

            this->bus_coefficients.clear();
            for (size_t k = 0; k < 3; ++k)
                this->bus_coefficients.push_back(std::stod(linecell[cellIndex++]));

        }//parse_input_line

        void PowerSystemOptions::write_output_line(std::ofstream& outputfile)
        {
            outputfile.precision(20);
            outputfile << this->name << " " << this->Spacecraft_Power_Supply_Type << " " << this->Spacecraft_Power_Supply_Curve_Type << " " << this->Spacecraft_Bus_Power_Type;
            outputfile << " " << this->P0 << " " << this->mass_per_kW << " " << this->decay_rate;

            for (size_t k = 0; k < 7; ++k)
                outputfile << " " << this->gamma[k];

            for (size_t k = 0; k < 3; ++k)
                outputfile << " " << this->bus_coefficients[k];

            outputfile << std::endl;

        }//write_output_line
        
        void PowerSystemOptions::setBusPowerType(const SpacecraftBusPowerType& input)
        {
            this->Spacecraft_Bus_Power_Type = input;
        }
        void PowerSystemOptions::setPowerSupplyType(const SpacecraftPowerSupplyType& input)
        {
            this->Spacecraft_Power_Supply_Type = input;
        }
        void PowerSystemOptions::setPowerSupplyCurveType(const SpacecraftPowerSupplyCurveType& input)
        {
            this->Spacecraft_Power_Supply_Curve_Type = input;
        }
        void PowerSystemOptions::setDecayRate(const double& input)
        {
            this->decay_rate = input;
        }
        void PowerSystemOptions::setPowerSystemDecayRefEpoch(const double& input)
        {
            this->power_system_decay_reference_epoch = input;
        }
        void PowerSystemOptions::setPowerMargin(const double& input)
        {
            this->PowerMargin = input;
        }
        void PowerSystemOptions::setGammaVector(const std::vector<double>& input)
        {
            this->gamma = input;
        }
        void PowerSystemOptions::setBusCoefficientVector(const std::vector<double>& input)
        {
            this->bus_coefficients = input;
        }
        void PowerSystemOptions::setmass_per_kw(const double& input)
        {
            this->mass_per_kW = input;
        }
        void PowerSystemOptions::setGamma(const size_t& p, const double& input)
        {
            if (p < this->gamma.size() - 1)
                this->gamma[p] = input;
            else
                throw std::invalid_argument("Array coefficient (gamma) " + std::to_string(p) + " of PowerSystem " + this->name + " does not exist. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }
        void PowerSystemOptions::setBusCoefficient(const size_t& p, const double& input)
        {
            if (p < this->bus_coefficients.size() - 1)
                this->bus_coefficients[p] = input;
            else
                throw std::invalid_argument("Bus power coefficient " + std::to_string(p) + " of PowerSystem " + this->name + " does not exist. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }


        //***************PowerSystemOptionsLibrary

        PowerSystemOptionsLibrary::PowerSystemOptionsLibrary()
        {
            this->PowerSystems.push_back(PowerSystemOptions());
            this->currentKey = this->PowerSystems.front().getName();
            this->currentIndex = 0;
        }

        PowerSystemOptionsLibrary::PowerSystemOptionsLibrary(const std::string& filename)
        {
            this->parse_input_file(filename);
        }

        PowerSystemOptionsLibrary::PowerSystemOptionsLibrary(std::ifstream& inputfile)
        {
            this->parse_input_file(inputfile);
        }

        PowerSystemOptionsLibrary::PowerSystemOptionsLibrary(const std::vector<PowerSystemOptions>& PowerSystemOptionsVector)
        {
            this->PowerSystems = PowerSystemOptionsVector;
            this->currentKey = this->PowerSystems.front().getName();
            this->currentIndex = 0;
        }

        //methods
        void PowerSystemOptionsLibrary::parse_input_file(const std::string& filename)
        {
            std::ifstream inputfile(filename);

            if (!inputfile.is_open())
            {
                throw std::invalid_argument("Cannot find " + filename + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            this->parse_input_file(inputfile);

            inputfile.close();
        }

        void PowerSystemOptionsLibrary::parse_input_file(std::ifstream& inputfile)
        {
            std::string line;
            this->PowerSystems.clear();

            while (EMTG::file_utilities::safeGetline(inputfile, line))
            {
                if (line.size() > 0)
                {
                    if (!(line.front() == *"#"))
                        this->PowerSystems.push_back(line);
                    else if (line == "#EndHardwareBlock")
                        return;
                }
            }
        }
        void PowerSystemOptionsLibrary::write_output_file(const std::string& filename, const bool& newFile)
        {
            std::ofstream outputfile;
            if (newFile)
                outputfile.open(filename, std::ios::trunc);
            else
            {
                outputfile.open(filename, std::ios::app);
            }

            outputfile << "#EMTG power system specification (space or comma delimited)" << std::endl;
            outputfile << "#Spacecraft_Power_Supply_Type choices are (0: Constant, 1: Solar)" << std::endl;
            outputfile << "#Spacecraft_Power_Supply_Curve_Type choices are (0: Sauer, 1: Polynomial)" << std::endl;
            outputfile << "#Spacecraft_Bus_Power_Type choices are (0: TypeA_Quadratic, 1: TypeB_Conditional)" << std::endl;
            outputfile << "#name Spacecraft_Power_Supply_Type Spacecraft_Power_Supply_Curve_Type Spacecraft_Bus_Power_Type P0 mass_per_kW(kg) decay_rate gamma[0:6] BusPower[0:2]" << std::endl;
            outputfile << "#" << std::endl;

            for (size_t Index = 0; Index < this->PowerSystems.size(); ++Index)
                this->PowerSystems[Index].write_output_line(outputfile);


            outputfile << "#EndHardwareBlock" << std::endl;
            outputfile << std::endl;

            outputfile.close();
        }


        //get functions
        PowerSystemOptions PowerSystemOptionsLibrary::getPowerSystem(const std::string& name)
        {
            size_t Index = this->getIndex(name);

            return this->PowerSystems[Index];
        }

        SpacecraftBusPowerType PowerSystemOptionsLibrary::getBusPowerType(const std::string& name)
        {
            return this->PowerSystems[this->getIndex(name)].getBusPowerType();
        }

        SpacecraftPowerSupplyCurveType PowerSystemOptionsLibrary::getPowerSupplyCurveType(const std::string& name)
        {
            return this->PowerSystems[this->getIndex(name)].getPowerSupplyCurveType();
        }

        SpacecraftPowerSupplyType PowerSystemOptionsLibrary::getPowerSupplyType(const std::string& name)
        {
            return this->PowerSystems[this->getIndex(name)].getPowerSupplyType();
        }

        double PowerSystemOptionsLibrary::getDecayRate(const std::string& name)
        {
            return this->PowerSystems[this->getIndex(name)].getDecayRate();
        }
        double PowerSystemOptionsLibrary::getPowerSystemDecayRefEpoch(const std::string& name)
        {
            return this->PowerSystems[this->getIndex(name)].getPowerSystemDecayRefEpoch();
        }
        double PowerSystemOptionsLibrary::getP0(const std::string& name)
        {
            return this->PowerSystems[this->getIndex(name)].getP0();
        }
        double PowerSystemOptionsLibrary::getPowerMargin(const std::string& name)
        {
            return this->PowerSystems[this->getIndex(name)].getPowerMargin();
        }
        std::vector<double> PowerSystemOptionsLibrary::getGamma(const std::string& name)
        {
            return this->PowerSystems[this->getIndex(name)].getGammaVector();
        }
        std::vector<double> PowerSystemOptionsLibrary::getBusCoefficients(const std::string& name)
        {
            return this->PowerSystems[this->getIndex(name)].getBusCoefficientVector();
        }
        double PowerSystemOptionsLibrary::getmass_per_kw(const std::string& name)
        {
            return this->PowerSystems[this->getIndex(name)].getmass_per_kw();
        }

        size_t PowerSystemOptionsLibrary::getIndex(const std::string& name)
        {
            if (name == this->currentKey)
            {
                return this->currentIndex;
            }

            for (size_t Index = 0; Index < this->PowerSystems.size(); ++Index)
                if (this->PowerSystems[Index].getName() == name)
                {
                    this->currentKey = name;
                    this->currentIndex = Index;
                    return this->currentIndex;
                }

            throw std::invalid_argument("PowerSystemOptionsLibrary::Power system '" + name + "' not found. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }
		
		std::vector<std::string> PowerSystemOptionsLibrary::getAllKeys() 
		{
			std::vector<std::string> allKeys;
            for (int idx = 0; idx < this->PowerSystems.size(); idx++)
            {
                allKeys.push_back(this->PowerSystems[idx].getName());
            }
			return allKeys;
		}
    }//end namespace HardwareModels
}//end namespace EMTG