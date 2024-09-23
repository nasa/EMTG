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

//options class for stage
//Jacob Englander 1-4-2017

#include <string>
#include <vector>

#include "StageOptions.h"
#include "PowerSystemOptions.h"
#include "PropulsionSystemOptions.h"

namespace EMTG
{
    namespace HardwareModels
    {
        //constructors
        StageOptions::StageOptions() :
            name("default_stage"),
            BaseDryMass(0.0),
            AdapterMass(0.0),
            EnableChemicalPropellantTankConstraint(false),
            EnableElectricPropellantTankConstraint(false),
            EnableDryMassConstraint(false),
            ElectricPropellantTankCapacity(1000.0),
            ChemicalFuelTankCapacity(1000.0),
            ChemicalOxidizerTankCapacity(1000.0),
            myPowerSystemOptionsName("default_power_system"),
            myElectricPropulsionSystemOptionsName("default_electric_propulsion_system"),
            myChemicalPropulsionSystemOptionsName("default_chemical_propulsion_system"),
            myThrottleLogic(EMTG::ThrottleLogic::MinThrusters),
            ThrottleSharpness(100.0)
        {
            this->myPowerSystemOptions.setName(this->myPowerSystemOptionsName);
            this->myElectricPropulsionSystemOptions.setName(this->myElectricPropulsionSystemOptionsName);
            this->myChemicalPropulsionSystemOptions.setName(this->myChemicalPropulsionSystemOptionsName);
        };

        StageOptions::StageOptions(const std::string& StageOptionsFileName) :
            StageOptions()
        {
            this->parse_input_file(StageOptionsFileName);
        };

        StageOptions::StageOptions(std::ifstream& StageOptionsStream) :
            StageOptions()
        {
            this->parse_input_file(StageOptionsStream);
        }

        StageOptions::StageOptions(const PropulsionSystemOptions& ChemicalPropulsionSystemOptionsInput,
                                   const PropulsionSystemOptions& ElectricPropulsionSystemOptionsInput,
                                   const PowerSystemOptions& PowerSystemOptionsInput) :
            StageOptions()
        {
            this->myPowerSystemOptions = PowerSystemOptionsInput;
            this->myElectricPropulsionSystemOptions = ElectricPropulsionSystemOptionsInput;
            this->myChemicalPropulsionSystemOptions = ChemicalPropulsionSystemOptionsInput;

            this->myPowerSystemOptionsName = this->myPowerSystemOptions.getName();
            this->myElectricPropulsionSystemOptionsName = this->myElectricPropulsionSystemOptions.getName();
            this->myChemicalPropulsionSystemOptionsName = this->myChemicalPropulsionSystemOptions.getName();
        };

        //file i/o
        void StageOptions::parse_input_file(const std::string& filename)
        {
            std::ifstream inputfile(filename);

            if (!inputfile.is_open())
            {
                throw std::invalid_argument("Cannot find " + filename + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            this->parse_input_file(inputfile);

            inputfile.close();
        }

        void StageOptions::parse_input_file(std::ifstream& inputfile)
        {
            std::string line;
            
            while (EMTG::file_utilities::safeGetline(inputfile, line))
            {
                if (line.size() > 0)
                {
                    if (line.front() == *"#")
                    {
                        if (line == "#BeginStagePowerLibraryBlock")
                        {
                            PowerSystemOptionsLibrary tempPowerSystemOptionsLibrary(inputfile);
                            this->myPowerSystemOptions = tempPowerSystemOptionsLibrary.getPowerSystem(this->myPowerSystemOptionsName);
							
							std::vector<std::string> allKeys = tempPowerSystemOptionsLibrary.getAllKeys();
							for (int power_idx = 0; power_idx < allKeys.size(); power_idx++)
							{
								this->allPowerSystemOptions.push_back(tempPowerSystemOptionsLibrary.getPowerSystem(allKeys[power_idx]));
							}
                        }
                        else if (line == "#BeginStagePropulsionLibraryBlock")
                        {
                            PropulsionSystemOptionsLibrary tempPropulsionSystemOptionsLibrary(inputfile);
                            this->myChemicalPropulsionSystemOptions = tempPropulsionSystemOptionsLibrary.getPropulsionSystem(this->myChemicalPropulsionSystemOptionsName);
                            this->myElectricPropulsionSystemOptions = tempPropulsionSystemOptionsLibrary.getPropulsionSystem(this->myElectricPropulsionSystemOptionsName);
							this->myElectricPropulsionSystemOptions.setSharpness(this->getThrottleSharpness());
                        }
                        else if (line == "#EndStageBlock")
                            return;
                    }
                    else
                    {
                        std::vector<std::string> linecell;
                        boost::split(linecell, line, boost::is_any_of(" "), boost::token_compress_on);
                        
                        if (linecell.front() == "name")
                        {
                            this->setName(linecell[1]);
                        }
                        else if (linecell.front() == "BaseDryMass")
                        {
                            this->setBaseDryMass(std::stod(linecell[1]));
                        }
                        else if (linecell.front() == "AdapterMass")
                        {
                            this->setAdapterMass(std::stod(linecell[1]));
                        }
                        else if (linecell.front() == "ElectricPropellantTankCapacity")
                        {
                            this->setElectricPropellantTankCapacity(std::stod(linecell[1]));
                        }
                        else if (linecell.front() == "ChemicalFuelTankCapacity")
                        {
                            this->setChemicalFuelTankCapacity(std::stod(linecell[1]));
                        }
                        else if (linecell.front() == "ChemicalOxidizerTankCapacity")
                        {
                            this->setChemicalOxidizerTankCapacity(std::stod(linecell[1]));
                        }
                        else if (linecell.front() == "EnableElectricPropellantTankConstraint")
                        {
                            this->setEnableElectricPropellantTankConstraint(std::stoi(linecell[1]));
                        }
                        else if (linecell.front() == "EnableChemicalPropellantTankConstraint")
                        {
                            this->setEnableChemicalPropellantTankConstraint(std::stoi(linecell[1]));
                        }
                        else if (linecell.front() == "EnableDryMassConstraint")
                        {
                            this->setEnableDryMassConstraint(std::stoi(linecell[1]));
                        }
                        else if (linecell.front() == "ThrottleLogic")
                        {
                            this->setThrottleLogic((EMTG::ThrottleLogic) std::stoi(linecell[1]));
                        }
                        else if (linecell.front() == "ThrottleSharpness")
                        {
                            this->setThrottleSharpness(std::stod(linecell[1]));
                        }
                        else if (linecell.front() == "PowerSystem")
                        {
                            this->setmyPowerSystemOptionsName(linecell[1]);
                        }
                        else if (linecell.front() == "ElectricPropulsionSystem")
                        {
                            this->setmyElectricPropulsionSystemOptionsName(linecell[1]);
                        }
                        else if (linecell.front() == "ChemicalPropulsionSystem")
                        {
                            this->setmyChemicalPropulsionSystemOptionsName(linecell[1]);
                        }
                        else
                        {
                            throw std::invalid_argument("StageOptions::invalid entry in input file. " + linecell.front() + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                        }
                    }
                }
            }
        }

        void StageOptions::write_output_file(const std::string& filename, const bool& newFile)
        {
            std::ofstream outputfile;
            if (newFile)
                outputfile.open(filename, std::ios::trunc);
            else
            {
                outputfile.open(filename, std::ios::app);
                outputfile << std::endl;
                outputfile << std::endl;
            }

            //do stuff
            outputfile << "#BeginStageBlock" << std::endl;
            outputfile << "name " << this->getName() << std::endl;
            outputfile << "BaseDryMass " << this->getBaseDryMass() << std::endl;
            outputfile << "AdapterMass " << this->getAdapterMass() << std::endl;
            outputfile << "EnableElectricPropellantTankConstraint " << this->getEnableElectricPropellantTankConstraint() << std::endl;
            outputfile << "EnableChemicalPropellantTankConstraint " << this->getEnableChemicalPropellantTankConstraint() << std::endl;
            outputfile << "EnableDryMassConstraint " << this->getEnableDryMassConstraint() << std::endl;
            outputfile << "ElectricPropellantTankCapacity " << this->getElectricPropellantTankCapacity() << std::endl;
            outputfile << "ChemicalFuelTankCapacity " << this->getChemicalFuelTankCapacity() << std::endl;
            outputfile << "ChemicalOxidizerTankCapacity " << this->getChemicalOxidizerTankCapacity() << std::endl;
            outputfile << "ThrottleLogic " << this->getThrottleLogic() <<  std::endl;
            outputfile << "ThrottleSharpness " << this->getThrottleSharpness() << std::endl;
            outputfile << "PowerSystem " << this->getmyPowerSystemOptionsName() << std::endl;
            outputfile << "ElectricPropulsionSystem " << this->getmyElectricPropulsionSystemOptionsName() << std::endl;
            outputfile << "ChemicalPropulsionSystem " << this->getmyChemicalPropulsionSystemOptionsName() << std::endl;
            outputfile << std::endl;

            outputfile << "#BeginStagePowerLibraryBlock" << std::endl;
            outputfile.close();
			if (this->allPowerSystemOptions.size() == 0)
				this->allPowerSystemOptions.push_back(this->myPowerSystemOptions);
            PowerSystemOptionsLibrary tempPowerSystemOptionsLibrary(this->allPowerSystemOptions);
            tempPowerSystemOptionsLibrary.write_output_file(filename, false);

            outputfile.open(filename, std::ios::app);
            outputfile << "#BeginStagePropulsionLibraryBlock" << std::endl;
            outputfile.close();
            std::vector<PropulsionSystemOptions> PropulsionSystemOptionsVector;
            PropulsionSystemOptionsVector.push_back(this->myElectricPropulsionSystemOptions);
            PropulsionSystemOptionsVector.push_back(this->myChemicalPropulsionSystemOptions);
            PropulsionSystemOptionsLibrary tempPropulsionSystemOptionsLibrary(PropulsionSystemOptionsVector);
            tempPropulsionSystemOptionsLibrary.write_output_file(filename, false);

            outputfile.open(filename, std::ios::app);
            outputfile << std::endl;
            outputfile << "#EndStageBlock" << std::endl;
            outputfile.close();
        }
    }//end namespace HardwareModels
}//end namespace EMTG