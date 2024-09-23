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

//options class for launch vehicles
//Jacob Englander 12-28-2016

#include <string>
#include <vector>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "LaunchVehicleOptions.h"
#include "EMTG_math.h"

namespace EMTG
{
    namespace HardwareModels
    {
        //LaunchVehicleOptions methods
        LaunchVehicleOptions::LaunchVehicleOptions() :
                                                    name("default"),
                                                    ModelType(0),
                                                    DLA_lowerbound(-math::PI),
                                                    DLA_upperbound(math::PI),
                                                    C3_lowerbound(0.0),
                                                    C3_upperbound(1000.0),
                                                    coefficients(std::vector<double>(1, 1.0))
        {}

        LaunchVehicleOptions::LaunchVehicleOptions(const std::string& linestring)
        {
            this->parse_input_line(linestring);
        }

        void LaunchVehicleOptions::parse_input_line(const std::string& linestring)
        {
            std::vector<std::string> linecell;
            boost::split(linecell, linestring, boost::is_any_of(" ,"), boost::token_compress_on);
            
            size_t cellIndex = 0;

            this->name = linecell[cellIndex++];
            this->ModelType = std::stoi(linecell[cellIndex++]);
            this->DLA_lowerbound = std::stod(linecell[cellIndex++]) * math::deg2rad;
            this->DLA_upperbound = std::stod(linecell[cellIndex++]) * math::deg2rad;
            this->C3_lowerbound = std::stod(linecell[cellIndex++]);
            this->C3_upperbound = std::stod(linecell[cellIndex++]);
            this->AdapterMass = std::stod(linecell[cellIndex++]);

            this->coefficients.clear();
            while (cellIndex < linecell.size())
                this->coefficients.push_back(std::stod(linecell[cellIndex++]));
        }

        void LaunchVehicleOptions::write_output_line(std::ofstream& outputfile)
        {
            outputfile.precision(20);
            outputfile << this->name << " " << this->ModelType << " " << this->DLA_lowerbound * 180.0 / math::PI << " " << this->DLA_upperbound * 180.0 / math::PI;
            outputfile << " " << this->C3_lowerbound << " " << this->C3_upperbound << " " << this->AdapterMass;
            
            for (size_t k = 0; k < this->coefficients.size(); ++k)
                outputfile << " " << this->coefficients[k];

            outputfile << std::endl;
        }

        std::string LaunchVehicleOptions::getName() const
        {
            return this->name;
        }

        size_t LaunchVehicleOptions::getModelType() const
        {
            return this->ModelType;
        }

        double LaunchVehicleOptions::getDLA_lowerbound() const
        {
            return this->DLA_lowerbound;
        }

        double LaunchVehicleOptions::getDLA_upperbound() const
        {
            return this->DLA_upperbound;
        }

        double LaunchVehicleOptions::getC3_lowerbound() const
        {
            return this->C3_lowerbound;
        }

        double LaunchVehicleOptions::getC3_upperbound() const
        {
            return this->C3_upperbound;
        }

        std::vector<double> LaunchVehicleOptions::getCoefficients() const
        {
            return this->coefficients;
        }
        
        double LaunchVehicleOptions::getAdapterMass() const
        {
            return this->AdapterMass;
        }

        void LaunchVehicleOptions::setName(const std::string& value)
        {
            this->name = value;
        }
        void LaunchVehicleOptions::setModelType(const size_t& value)
        {
            this->ModelType = value;
        }
        void LaunchVehicleOptions::setDLA_lowerbound(const double& value)
        {
            this->DLA_lowerbound = value;
        }
        void LaunchVehicleOptions::setDLA_upperbound(const double& value)
        {
            this->DLA_upperbound = value;
        }
        void LaunchVehicleOptions::setC3_lowerbound(const double& value)
        {
            this->C3_lowerbound = value;
        }
        void LaunchVehicleOptions::setC3_upperbound(const double& value)
        {
            this->C3_upperbound = value;
        }
        void LaunchVehicleOptions::setCoefficients(const std::vector<double> Coefficients)
        {
            this->coefficients = Coefficients;
        }
        void LaunchVehicleOptions::setCoefficient(const size_t& p, const double& value)
        {
            if (p < this->coefficients.size() - 1)
                this->coefficients[p] = value;
            else
                throw std::invalid_argument("Coefficient " + std::to_string(p) + " of LaunchVehicle " + this->name + " does not exist. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }
        void LaunchVehicleOptions::setAdapterMass(const double& value)
        {
            this->AdapterMass = value;
        }

        //LaunchVehicleOptionsLibrary methods

        LaunchVehicleOptionsLibrary::LaunchVehicleOptionsLibrary()
        {
            this->LaunchVehicles.push_back(LaunchVehicleOptions());

            this->currentKey = this->LaunchVehicles.front().getName();
            this->currentIndex = 0;
        }

        LaunchVehicleOptionsLibrary::LaunchVehicleOptionsLibrary(const std::string& filename)
        {
            this->parse_input_file(filename);
        }

        LaunchVehicleOptionsLibrary::LaunchVehicleOptionsLibrary(std::ifstream& inputfile)
        {
            this->parse_input_file(inputfile);
        }

        LaunchVehicleOptionsLibrary::LaunchVehicleOptionsLibrary(const std::vector<LaunchVehicleOptions>& LaunchVehicleOptionsVector)
        {
            this->LaunchVehicles = LaunchVehicleOptionsVector;

            this->currentKey = this->LaunchVehicles.front().getName();
            this->currentIndex = 0;
        }

        void LaunchVehicleOptionsLibrary::parse_input_file(const std::string& filename)
        {
            std::ifstream inputfile(filename);

            if (!inputfile.is_open())
            {
                throw std::invalid_argument("Cannot find " + filename + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            this->parse_input_file(inputfile);

            inputfile.close();
        }

        void LaunchVehicleOptionsLibrary::parse_input_file(std::ifstream& inputfile)
        {
            std::string line;
            this->LaunchVehicles.clear();

            while (EMTG::file_utilities::safeGetline(inputfile, line))
            {
                if (line.size() > 0)
                {
                    if (!(line.front() == *"#"))
                        this->LaunchVehicles.push_back(line);
                    else if (line == "#EndHardwareBlock")
                        return;
                }
            }
        }

        void LaunchVehicleOptionsLibrary::write_output_file(const std::string& filename, const bool& newFile)
        {
            std::ofstream outputfile;
            if (newFile)
                 outputfile.open(filename, std::ios::trunc);
            else
            {
                outputfile.open(filename, std::ios::app);
            }

            outputfile << "#EMTG launch vehicle library (space delimited)" << std::endl;
            outputfile << "#tokens are:" << std::endl;
            outputfile << "#name ModelType DLA_lowerbound(degrees) DLA_upperbound(degrees) AdapterMass C3_lowerbound C3_upperbound coefficient[0] coefficient[1] ..." << std::endl;
            outputfile << "#ModelType = 0 is the standard polynomial curve, with the ith coefficient corresponding to the ith power of C3" << std::endl;
            outputfile << "#Additional ModelTypes may be implemented later" << std::endl;
            outputfile << "#" << std::endl;

            for (size_t Index = 0; Index < this->LaunchVehicles.size(); ++Index)
                this->LaunchVehicles[Index].write_output_line(outputfile);

            outputfile << "#EndHardwareBlock" << std::endl;
            outputfile << std::endl;
            outputfile.close();
        }

        LaunchVehicleOptions LaunchVehicleOptionsLibrary::getLaunchVehicle(const std::string& name)
        {
            size_t Index = this->getIndex(name);

            return this->LaunchVehicles[Index];
        }

        size_t LaunchVehicleOptionsLibrary::getIndex(const std::string& name)
        {
            if (name == this->currentKey)
            {
                return this->currentIndex;
            }

            for (size_t Index = 0; Index < this->LaunchVehicles.size(); ++Index)
                if (this->LaunchVehicles[Index].getName() == name)
                {
                    this->currentKey = name;
                    this->currentIndex = Index;
                    return this->currentIndex;
                }

            throw std::invalid_argument("LaunchVehicleOptionsLibrary::Launch vehicle '" + name + "' not found. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }

        size_t LaunchVehicleOptionsLibrary::getModelType(const std::string& name)
        {
            size_t Index = this->getIndex(name);

            return this->LaunchVehicles[Index].getModelType();
        }

        double LaunchVehicleOptionsLibrary::getDLA_lowerbound(const std::string& name)
        {
            size_t Index = this->getIndex(name);

            return this->LaunchVehicles[Index].getDLA_lowerbound();
        }

        double LaunchVehicleOptionsLibrary::getDLA_upperbound(const std::string& name)
        {
            size_t Index = this->getIndex(name);

            return this->LaunchVehicles[Index].getDLA_upperbound();
        }

        double LaunchVehicleOptionsLibrary::getC3_lowerbound(const std::string& name)
        {
            size_t Index = this->getIndex(name);

            return this->LaunchVehicles[Index].getC3_lowerbound();
        }

        double LaunchVehicleOptionsLibrary::getC3_upperbound(const std::string& name)
        {
            size_t Index = this->getIndex(name);

            return this->LaunchVehicles[Index].getC3_upperbound();
        }

        double LaunchVehicleOptionsLibrary::getAdapterMass(const std::string& name)
        {
            size_t Index = this->getIndex(name);

            return this->LaunchVehicles[Index].getAdapterMass();
        }

        std::vector<double> LaunchVehicleOptionsLibrary::getCoefficients(const std::string& name)
        {
            size_t Index = this->getIndex(name);

            return this->LaunchVehicles[Index].getCoefficients();
        }
    }//end namespace HardwareModels
}//end namespace EMTG