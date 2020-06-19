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

#ifndef EMTG_STAGEOPTIONS
#define EMTG_STAGEOPTIONS

#include "ElectricPropulsionSystem.h"
#include "ChemicalPropulsionSystem.h"
#include "PowerSystemOptions.h"
#include "HardwareBase.h"

#include <string>
#include <vector>

#include "file_utilities.h"

namespace EMTG
{
    namespace HardwareModels
    {
        class StageOptions
        {
        public:
            //constructor
            StageOptions();//default constructor
            StageOptions(const std::string& StageOptionsFileName);
            StageOptions(std::ifstream& StageOptionsStream);
            StageOptions(const PropulsionSystemOptions& ChemicalPropulsionSystemOptionsInput,
                         const PropulsionSystemOptions& ElectricPropulsionSystemOptionsInput,
                         const PowerSystemOptions& PowerSystemOptionsInput);

            //destructor
            ~StageOptions() {};
                
            //get/set
            void setName(const std::string& name) { this->name = name; }
            void setBaseDryMass(const double& BaseDryMassInput) { this->BaseDryMass = BaseDryMassInput; }
            void setAdapterMass(const double& AdapterMassInput) { this->AdapterMass = AdapterMassInput; }
            void setElectricPropellantTankCapacity(const double& ElectricPropellantTankCapacityInput) { this->ElectricPropellantTankCapacity = ElectricPropellantTankCapacityInput; }
            void setChemicalFuelTankCapacity(const double& ChemicalFuelTankCapacityInput) { this->ChemicalFuelTankCapacity = ChemicalFuelTankCapacityInput; }
            void setChemicalOxidizerTankCapacity(const double& ChemicalOxidizerTankCapacityInput) { this->ChemicalOxidizerTankCapacity = ChemicalOxidizerTankCapacityInput; }
            void setEnableElectricPropellantTankConstraint(const bool& EnableElectricPropellantTankConstraint) { this->EnableElectricPropellantTankConstraint = EnableElectricPropellantTankConstraint; }
            void setEnableChemicalPropellantTankConstraint(const bool& EnableChemicalPropellantTankConstraint) { this->EnableChemicalPropellantTankConstraint = EnableChemicalPropellantTankConstraint; }
            void setEnableDryMassConstraint(const bool& EnableDryMassConstraint) { this->EnableDryMassConstraint = EnableDryMassConstraint; }
            void setmyPowerSystemOptionsName(const std::string& myPowerSystemOptionsName) { this->myPowerSystemOptionsName = myPowerSystemOptionsName; }
            void setmyElectricPropulsionSystemOptionsName(const std::string& myElectricPropulsionSystemOptionsName) { this->myElectricPropulsionSystemOptionsName = myElectricPropulsionSystemOptionsName; }
            void setmyChemicalPropulsionSystemOptionsName(const std::string& myChemicalPropulsionSystemOptionsName) { this->myChemicalPropulsionSystemOptionsName = myChemicalPropulsionSystemOptionsName; }
            void setThrottleLogic(const EMTG::ThrottleLogic& throttlelogic) { this->myThrottleLogic = throttlelogic; }
            void setThrottleSharpness(const double& sharpness) { this->ThrottleSharpness = sharpness; }
			void setElectricPropulsionSystemOptions(const PropulsionSystemOptions& myElectricPropulsionSystemOptions) {this->myElectricPropulsionSystemOptions = myElectricPropulsionSystemOptions; }
            void setChemicalPropulsionSystemOptions(const PropulsionSystemOptions& myChemicalPropulsionSystemOptions) { this->myChemicalPropulsionSystemOptions = myChemicalPropulsionSystemOptions; }
            void setPowerSystemOptions(const PowerSystemOptions& myPowerSystemOptions) { this->myPowerSystemOptions = myPowerSystemOptions; }


            std::string getName() const { return this->name; };
            double getBaseDryMass() const { return this->BaseDryMass; }
            double getAdapterMass() const { return this->AdapterMass; }
            double getElectricPropellantTankCapacity() const { return this->ElectricPropellantTankCapacity; }
            double getChemicalFuelTankCapacity() const { return this->ChemicalFuelTankCapacity; }
            double getChemicalOxidizerTankCapacity() const { return this->ChemicalOxidizerTankCapacity; }
            bool getEnableElectricPropellantTankConstraint() const { return this->EnableElectricPropellantTankConstraint; }
            bool getEnableChemicalPropellantTankConstraint() const { return this->EnableChemicalPropellantTankConstraint; }
            bool getEnableDryMassConstraint() const { return this->EnableDryMassConstraint; }
            std::string getmyChemicalPropulsionSystemOptionsName() const { return this->myChemicalPropulsionSystemOptionsName; }
            std::string getmyElectricPropulsionSystemOptionsName() const { return this->myElectricPropulsionSystemOptionsName; }
            std::string getmyPowerSystemOptionsName() const { return this->myPowerSystemOptionsName; }

            PropulsionSystemOptions getElectricPropulsionSystemOptions() const { return this->myElectricPropulsionSystemOptions; }
            PropulsionSystemOptions getChemicalPropulsionSystemOptions() const { return this->myChemicalPropulsionSystemOptions; }
            PowerSystemOptions getPowerSystemOptions() const { return this->myPowerSystemOptions; }

            EMTG::ThrottleLogic getThrottleLogic() const { return this->myThrottleLogic; }
            double getThrottleSharpness() const { return this->ThrottleSharpness; }

            //methods
            void parse_input_file(const std::string& filename);
            void parse_input_file(std::ifstream& inputfile);
            void write_output_file(const std::string& filename, const bool& newFile = true);

        protected:
            //methods

            //fields
            std::string name;
            double BaseDryMass;//mass before subsystems are added
            double AdapterMass;//jettisoned before the stage does its business
            double ElectricPropellantTankCapacity;
            double ChemicalFuelTankCapacity;
            double ChemicalOxidizerTankCapacity;
            std::string myPowerSystemOptionsName;
            std::string myElectricPropulsionSystemOptionsName;
            std::string myChemicalPropulsionSystemOptionsName;
            PowerSystemOptions myPowerSystemOptions;
            PropulsionSystemOptions myElectricPropulsionSystemOptions;
            PropulsionSystemOptions myChemicalPropulsionSystemOptions;
            bool EnableElectricPropellantTankConstraint;
            bool EnableChemicalPropellantTankConstraint;
            bool EnableDryMassConstraint;
			
			std::vector<PowerSystemOptions> allPowerSystemOptions;

            EMTG::ThrottleLogic myThrottleLogic;
            double ThrottleSharpness;
        };
    }//end namespace HardwareModels
}//end namespace EMTG

#endif //EMTG_STAGEOPTIONS