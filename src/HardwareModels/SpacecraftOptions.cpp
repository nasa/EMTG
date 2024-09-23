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

//options class for spacecraft
//Jacob Englander 10-28-2016

#include <string>
#include <vector>

#include "SpacecraftOptions.h"

namespace EMTG
{
    namespace HardwareModels
    {
        //constructors
        SpacecraftOptions::SpacecraftOptions() :
            name("default_spacecraft"),
            number_of_stages(1),
            GlobalElectricPropellantTankCapacity(1000.0),
            GlobalChemicalFuelTankCapacity(1000.0),
            GlobalChemicalOxidizerTankCapacity(1000.0),
            GlobalDryMassBounds(std::vector<double>({ 1000.0, 1500.0 })),
            EnableGlobalElectricPropellantTankConstraint(false),
            EnableGlobalChemicalPropellantTankConstraint(false),
            EnableGlobalDryMassConstraint(false),
            Stages(std::vector<StageOptions>(1))
        {}

        SpacecraftOptions::SpacecraftOptions(const std::string& SpacecraftOptionsFileName) :
            SpacecraftOptions()
        {
            this->number_of_stages = 0;
            this->parse_input_file(SpacecraftOptionsFileName);
        }

        SpacecraftOptions::SpacecraftOptions(const std::vector<StageOptions>& stages) :
            Stages(stages),
            number_of_stages(stages.size()),
            GlobalElectricPropellantTankCapacity(1000.0),
            GlobalChemicalFuelTankCapacity(1000.0),
            GlobalChemicalOxidizerTankCapacity(1000.0),
            GlobalDryMassBounds(std::vector<double>({ 1000.0, 1500.0 })),
            EnableGlobalElectricPropellantTankConstraint(false),
            EnableGlobalChemicalPropellantTankConstraint(false),
            EnableGlobalDryMassConstraint(false)
        {}

        //methods
        void SpacecraftOptions::parse_input_file(const std::string& filename)
        {
            this->Stages.clear();

            std::ifstream inputfile(filename);

            if (!inputfile.is_open())
            {
                throw std::invalid_argument("Cannot find " + filename + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            
            std::string line;


            while (EMTG::file_utilities::safeGetline(inputfile, line))
            {
                if (line.size() > 0)
                {
                    if (line.front() == *"#")
                    {
                        if (line == "#BeginStageBlock")
                        {
                            StageOptions tempStageOptions(inputfile);
                            this->Stages.push_back(tempStageOptions);
                        }
                        else if (line == "#EndSpacecraftBlock")
                        {
                            inputfile.close();

                            this->number_of_stages = this->Stages.size();
                            return;
                        }
                    }
                    else
                    {
                        std::vector<std::string> linecell;
                        boost::split(linecell, line, boost::is_any_of(" "), boost::token_compress_on);
                        if (linecell.front() == "name")
                        {
                            this->setName(linecell[1]);
                        }
                        else if (linecell.front() == "GlobalElectricPropellantTankCapacity")
                        {
                            this->setGlobalElectricPropellantTankCapacity(std::stod(linecell[1]));
                        }
                        else if (linecell.front() == "GlobalFuelTankCapacity")
                        {
                            this->setGlobalFuelTankCapacity(std::stod(linecell[1]));
                        }
                        else if (linecell.front() == "GlobalOxidizerTankCapacity")
                        {
                            this->setGlobalOxidizerTankCapacity(std::stod(linecell[1]));
                        }
                        else if (linecell.front() == "GlobalDryMassBounds")
                        {
                            this->setGlobalDryMassBounds(std::vector<double>(std::stod(linecell[1]), std::stod(linecell[2])));
                        }
                        else if (linecell.front() == "EnableGlobalElectricPropellantTankConstraint")
                        {
                            this->setEnableGlobalElectricPropellantTankConstraint(std::stoi(linecell[1]));
                        }
                        else if (linecell.front() == "EnableGlobalChemicalPropellantTankConstraint")
                        {
                            this->setEnableGlobalChemicalPropellantTankConstraint(std::stoi(linecell[1]));
                        }
						else if (linecell.front() == "EnableGlobalDryMassConstraint")
						{
							this->setEnableGlobalDryMassConstraint(std::stoi(linecell[1]));
						}
                        else
                        {
                            throw std::invalid_argument("SpacecraftOptions::invalid entry in input file. " + linecell.front() + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                        }
                    }
                }
            }
        }

        void SpacecraftOptions::write_output_file(const std::string& filename)
        {
            std::ofstream outputfile;
            outputfile.open(filename, std::ios::trunc);
            outputfile << "#BeginSpacecraftBlock" << std::endl;
            outputfile << "name " << this->name << std::endl;
            outputfile << "EnableGlobalElectricPropellantTankConstraint " << this->getEnableGlobalElectricPropellantTankConstraint() << std::endl;
            outputfile << "EnableGlobalChemicalPropellantTankConstraint " << this->getEnableGlobalChemicalPropellantTankConstraint() << std::endl;
            outputfile << "EnableGlobalDryMassConstraint " << this->getEnableGlobalDryMassConstraint() << std::endl;
            outputfile << "GlobalElectricPropellantTankCapacity " << this->GlobalElectricPropellantTankCapacity << std::endl;
            outputfile << "GlobalFuelTankCapacity " << this->GlobalChemicalFuelTankCapacity << std::endl;
            outputfile << "GlobalOxidizerTankCapacity " << this->GlobalChemicalOxidizerTankCapacity << std::endl;
            outputfile << "GlobalDryMassBounds " << this->getGlobalDryMassBounds()[0] << " " << this->getGlobalDryMassBounds()[1] << std::endl;
            outputfile.close();
            
            for (size_t StageIndex = 0; StageIndex < this->number_of_stages; ++StageIndex)
            {
                this->Stages[StageIndex].write_output_file(filename, false);
            }

            outputfile.open(filename, std::ios::app);
            outputfile << "#EndSpacecraftBlock" << std::endl;
            outputfile.close();
        }


        StageOptions SpacecraftOptions::getStageOptions(const size_t& StageIndex) const
        {
            if (StageIndex > this->number_of_stages)
            {
                throw std::invalid_argument("SpacecraftOptions::StageIndex " + std::to_string(StageIndex)
                    + " exceeds number of stages " + std::to_string(this->number_of_stages) + "! Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            return this->Stages[StageIndex];
        }
    }//end namespace HardwareModels
}//end namespace EMTG