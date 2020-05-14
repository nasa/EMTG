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

#pragma once

#include <string>
#include <fstream>
#include <vector>

#include "file_utilities.h"

namespace EMTG
{
    namespace HardwareModels
    {
        class LaunchVehicleOptions
        {
        public:
            //constructor
            LaunchVehicleOptions();
            LaunchVehicleOptions(const std::string& linestring);

            //destructor
            ~LaunchVehicleOptions() {};//destructor does not need to do anything

            //methods
            void parse_input_line(const std::string& linestring);
            void write_output_line(std::ofstream& outputfile);

            //get/set
            std::string getName() const;
            size_t getModelType() const;
            double getDLA_lowerbound() const;
            double getDLA_upperbound() const;
            double getC3_lowerbound() const;
            double getC3_upperbound() const;
            double getAdapterMass() const;
            std::vector<double> getCoefficients() const;

            inline size_t getNumberOfCoefficients() const
            {
                return this->coefficients.size();
            }

            inline double getCoefficient(const size_t& p) const
            {
                return this->coefficients[p];
            }

            void setName(const std::string& value);
            void setModelType(const size_t& value);
            void setDLA_lowerbound(const double& value);
            void setDLA_upperbound(const double& value);
            void setC3_lowerbound(const double& value);
            void setC3_upperbound(const double& value);
            void setCoefficients(const std::vector<double> Coefficients);
            void setCoefficient(const size_t& p, const double& value);
            void setAdapterMass(const double& value);

        private:
            //fields
            std::string name;
            size_t ModelType;
            double DLA_lowerbound;
            double DLA_upperbound;
            double C3_lowerbound;
            double C3_upperbound;
            double AdapterMass;
            std::vector<double> coefficients;
        };

        class LaunchVehicleOptionsLibrary
        {
        public:
            //constructor
            LaunchVehicleOptionsLibrary();
            LaunchVehicleOptionsLibrary(const std::string& filename);
            LaunchVehicleOptionsLibrary(std::ifstream& inputfile);
            LaunchVehicleOptionsLibrary(const std::vector<LaunchVehicleOptions>& LaunchVehicleOptionsVector);

            //destructor
            ~LaunchVehicleOptionsLibrary() {};//destructor does not need to do anything

            //methods
            void parse_input_file(const std::string& filename);
            void parse_input_file(std::ifstream& inputfile);
            void write_output_file(const std::string& filename, const bool& newFile = true);
            
            
            //get functions
            LaunchVehicleOptions getLaunchVehicle(const std::string& name);
            size_t getModelType(const std::string& name);
            double getDLA_lowerbound(const std::string& name);
            double getDLA_upperbound(const std::string& name);
            double getC3_lowerbound(const std::string& name);
            double getC3_upperbound(const std::string& name);
            double getAdapterMass(const std::string& name);
            std::vector<double> getCoefficients(const std::string& name);

            
        private:
            //private methods
            size_t getIndex(const std::string& name);
            
            //fields
            std::vector<LaunchVehicleOptions> LaunchVehicles;
            std::string currentKey;
            size_t currentIndex;
        };
    }//end namespace HardwareModels
}//end namespace EMTG