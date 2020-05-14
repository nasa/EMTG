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

//options class for propulsion system
//Jacob Englander 12-30-2016

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "PropulsionSystemOptions.h"

#include <boost/algorithm/string.hpp>

namespace EMTG
{
    namespace HardwareModels
    {
        //***************PropulsionSystemOptions
        PropulsionSystemOptions::PropulsionSystemOptions() :
            name("default"),
            ThrusterMode(SpacecraftThrusterMode::ConstantThrustIsp),
            Pmin(0.0),
            Pmax(100000.0),
            ThrustCoefficients(std::vector<double>(7, 0.0)),
            MassFlowCoefficients(std::vector<double>(7, 0.0)),
            ThrottleTableFile("none"),
            ConstantThrust(0.10),
            ConstantIsp(3000.0),
            MinimumOrMonopropIsp(220.0),
            FixedEfficiency(0.70),
            MassPerString(30.0),
            NumberOfStrings(1),
            MixtureRatio(1.0),
            g0(9.80665),
            ThrustScaleFactor(1.0)
        {}

        PropulsionSystemOptions::PropulsionSystemOptions(const std::string& linestring) :
            g0(9.80665)
        {
            this->parse_input_line(linestring);
        }

        void PropulsionSystemOptions::parse_input_line(const std::string& linestring)
        {
            std::vector<std::string> linecell;
            boost::split(linecell, linestring, boost::is_any_of(" ,"), boost::token_compress_on);

            size_t cellIndex = 0;

            this->name = linecell[cellIndex++];
            this->ThrusterMode = (SpacecraftThrusterMode)std::stoi(linecell[cellIndex++]);
            this->ThrottleTableFile = linecell[cellIndex++];
            this->MassPerString = std::stod(linecell[cellIndex++]);
            this->NumberOfStrings = std::stoi(linecell[cellIndex++]);
            this->Pmin = std::stod(linecell[cellIndex++]);
            this->Pmax = std::stod(linecell[cellIndex++]);
            this->ConstantThrust = std::stod(linecell[cellIndex++]);
            this->ConstantIsp = std::stod(linecell[cellIndex++]);
            this->MinimumOrMonopropIsp = std::stod(linecell[cellIndex++]);
            this->FixedEfficiency = std::stod(linecell[cellIndex++]);
            this->MixtureRatio = std::stod(linecell[cellIndex++]);
            this->ThrustScaleFactor = std::stod(linecell[cellIndex++]);

            this->ThrustCoefficients.clear();
            for (size_t k = 0; k < 7; ++k)
                this->ThrustCoefficients.push_back(std::stod(linecell[cellIndex++]));

            this->MassFlowCoefficients.clear();
            for (size_t k = 0; k < 7; ++k)
                this->MassFlowCoefficients.push_back(std::stod(linecell[cellIndex++]));

        }//parse_input_line

        void PropulsionSystemOptions::write_output_line(std::ofstream& outputfile)
        {
            outputfile.precision(20);
            outputfile << this->name << " " << this->ThrusterMode << " " << this->ThrottleTableFile << " " << this->MassPerString << " " << this->NumberOfStrings;
            outputfile << " " << this->Pmin << " " << this->Pmax << " " << this->ConstantThrust << " " << this->ConstantIsp << " " << this->MinimumOrMonopropIsp;
            outputfile << " " << this->FixedEfficiency << " " << this->MixtureRatio << " " << this->ThrustScaleFactor;

            for (size_t k = 0; k < 7; ++k)
                outputfile << " " << this->ThrustCoefficients[k];

            for (size_t k = 0; k < 7; ++k)
                outputfile << " " << this->MassFlowCoefficients[k];

            outputfile << std::endl;

        }//write_output_line

        void PropulsionSystemOptions::setName(const std::string& input)
        {
            this->name = input;
        }
        void PropulsionSystemOptions::setNumberOfStrings(const size_t& input)
        {
            this->NumberOfStrings = input;
        }
        void PropulsionSystemOptions::setThrustScaleFactor(const double& input)
        {
            this->ThrustScaleFactor = input;
        }
        void PropulsionSystemOptions::setThrusterMode(const SpacecraftThrusterMode& input)
        {
            this->ThrusterMode = input;
        }
        void PropulsionSystemOptions::setConstantThrust(const double& input)
        {
            this->ConstantThrust = input;
        }
        void PropulsionSystemOptions::setConstantIsp(const double& input)
        {
            this->ConstantIsp = input;
        }
        void PropulsionSystemOptions::setMinimumOrMonopropIsp(const double& input)
        {
            this->MinimumOrMonopropIsp = input;
        }
        void PropulsionSystemOptions::setFixedEfficiency(const double& input)
        {
            this->FixedEfficiency = input;
        }
        void PropulsionSystemOptions::setPmin(const double& input)
        {
            this->Pmin = input;
        }
        void PropulsionSystemOptions::setPmax(const double& input)
        {
            this->Pmax = input;
        }
        void PropulsionSystemOptions::setThrottleTableFile(const std::string& input)
        {
            this->ThrottleTableFile = input;
        }
        void PropulsionSystemOptions::setThrustCoefficients(const std::vector<double>& input)
        {
            this->ThrustCoefficients = input;
        }
        void PropulsionSystemOptions::setMassFlowCoefficients(const std::vector<double>& input)
        {
            this->MassFlowCoefficients = input;
        }
        void PropulsionSystemOptions::setMassPerString(const double& input)
        {
            this->MassPerString = input;
        }
        void PropulsionSystemOptions::setMixtureRatio(const double& input)
        {
            this->MixtureRatio = input;
        }
        void PropulsionSystemOptions::setg0(const double& input)
        {
            this->g0 = input;
        }
		void PropulsionSystemOptions::setSharpness(const double& input)
		{
			this->Sharpness = input;
		}
        void PropulsionSystemOptions::setThrustCoefficient(const size_t& p, const double& input)
        {
            if (p < this->ThrustCoefficients.size() - 1)
                this->ThrustCoefficients[p] = input;
            else
                throw std::invalid_argument("Thrust coefficient " + std::to_string(p) + " of PropulsionSystem " + this->name + " does not exist. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }
        void PropulsionSystemOptions::setMassFlowCoefficient(const size_t& p, const double& input)
        {
            if (p < this->MassFlowCoefficients.size() - 1)
                this->MassFlowCoefficients[p] = input;
            else
                throw std::invalid_argument("Mass flow coefficient " + std::to_string(p) + " of PropulsionSystem " + this->name + " does not exist. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }

        //***************PropulsionSystemOptionsLibrary
        PropulsionSystemOptionsLibrary::PropulsionSystemOptionsLibrary()
        {
            this->PropulsionSystems.push_back(PropulsionSystemOptions());
            this->currentKey = this->PropulsionSystems.front().getName();
            this->currentIndex = 0;
        }

        PropulsionSystemOptionsLibrary::PropulsionSystemOptionsLibrary(const std::string& filename)
        {
            this->parse_input_file(filename);
        }

        PropulsionSystemOptionsLibrary::PropulsionSystemOptionsLibrary(std::ifstream& inputfile)
        {
            this->parse_input_file(inputfile);
        }

        PropulsionSystemOptionsLibrary::PropulsionSystemOptionsLibrary(const std::vector<PropulsionSystemOptions>& PropulsionSystemOptionsVector)
        {
            this->PropulsionSystems = PropulsionSystemOptionsVector;
            this->currentKey = this->PropulsionSystems.front().getName();
            this->currentIndex = 0;
        }

        //methods
        void PropulsionSystemOptionsLibrary::parse_input_file(const std::string& filename)
        {
            std::ifstream inputfile(filename);

            if (!inputfile.is_open())
            {
                throw std::invalid_argument("Cannot find " + filename + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            this->parse_input_file(inputfile);

            inputfile.close();
        }
         
        void PropulsionSystemOptionsLibrary::parse_input_file(std::ifstream& inputfile)
        {
            std::string line;
            this->PropulsionSystems.clear();

            while (EMTG::file_utilities::safeGetline(inputfile, line))
            {
                if (line.size() > 0)
                {
                    if (!(line.front() == *"#"))
                        this->PropulsionSystems.push_back(line);
                    else if (line == "#EndHardwareBlock")
                        return;
                }
            }
        }

        void PropulsionSystemOptionsLibrary::write_output_file(const std::string& filename, const bool& newFile)
        {
            std::ofstream outputfile;
            if (newFile)
                outputfile.open(filename, std::ios::trunc);
            else
            {
                outputfile.open(filename, std::ios::app);
            }

            outputfile << "#EMTG propulsion system specification (space or comma delimited)" << std::endl;
            outputfile << "#SpacecraftThrusterMode choices are (0: ConstantThrustIsp, 1: FixedEfficiencyCSI, 2: FixedEfficiencyVSI, 3: Poly1D, 4: SteppedHThrust1D, 5: SteppedLMdot1D, 6: SteppedHIsp1D, 7: SteppedHefficiency1D, 8: SteppedFullSet1D, 8: Stepped2D, 9: Poly2D)" << std::endl;
            outputfile << "#name ThrusterMode ThrottleTableFile MassPerString(kg) NumberOfStrings Pmin(kW) Pmax(kW) ConstantThrust(N) ConstantIsp(s) MinimumIsp/MonopropIsp(s) FixedEfficiency MixtureRatio ThrustScaleFactor ThrustCoefficients(mN)[0:6] MassFlowCoefficients(mg/s)[0:6]" << std::endl;
            outputfile << "#polynomial coefficients are in order from P^0 to P^5" << std::endl;
            outputfile << "#" << std::endl;

            for (size_t Index = 0; Index < this->PropulsionSystems.size(); ++Index)
                this->PropulsionSystems[Index].write_output_line(outputfile);
            
            
            outputfile << "#EndHardwareBlock" << std::endl;
            outputfile << std::endl;

            outputfile.close();
        }

        PropulsionSystemOptions PropulsionSystemOptionsLibrary::getPropulsionSystem(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)];
        }

        SpacecraftThrusterMode PropulsionSystemOptionsLibrary::getThrusterMode(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getThrusterMode();
        }

        double PropulsionSystemOptionsLibrary::getConstantThrust(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getConstantThrust();
        }

        double PropulsionSystemOptionsLibrary::getConstantIsp(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getConstantIsp();
        }

        double PropulsionSystemOptionsLibrary::getMinimumOrMonopropIsp(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getMinimumOrMonopropIsp();
        }

        double PropulsionSystemOptionsLibrary::getFixedEfficiency(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getFixedEfficiency();
        }

        double PropulsionSystemOptionsLibrary::getPmin(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getPmin();
        }

        double PropulsionSystemOptionsLibrary::getPmax(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getPmax();
        }

        std::string PropulsionSystemOptionsLibrary::getThrottleTableFile(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getThrottleTableFile();
        }

        std::vector<double> PropulsionSystemOptionsLibrary::getThrustCoefficients(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getThrustCoefficients();
        }

        std::vector<double> PropulsionSystemOptionsLibrary::getMassFlowCoefficients(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getMassFlowCoefficients();
        }

        double PropulsionSystemOptionsLibrary::getMassPerString(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getMassPerString();
        }

        size_t PropulsionSystemOptionsLibrary::getNumberOfStrings(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getNumberOfStrings();
        }

        double PropulsionSystemOptionsLibrary::getMixtureRatio(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getMixtureRatio();
        }

        double PropulsionSystemOptionsLibrary::getThrustScaleFactor(const std::string& name)
        {
            return this->PropulsionSystems[this->getIndex(name)].getThrustScaleFactor();
        }

        size_t PropulsionSystemOptionsLibrary::getIndex(const std::string& name)
        {
            if (name == this->currentKey)
            {
                return this->currentIndex;
            }

            for (size_t Index = 0; Index < this->PropulsionSystems.size(); ++Index)
                if (this->PropulsionSystems[Index].getName() == name)
                {
                    this->currentKey = name;
                    this->currentIndex = Index;
                    return this->currentIndex;
                }

            throw std::invalid_argument("PropulsionSystemOptionsLibrary::Propulsion system '" + name + "' not found. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }

    }//end namespace HardwareModels
}//end namespace EMTG