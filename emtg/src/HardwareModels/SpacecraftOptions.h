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

#pragma once

#include <string>
#include <vector>
#include "StageOptions.h"
#include "file_utilities.h"

namespace EMTG
{
    namespace HardwareModels
    {
        class SpacecraftOptions
        {
        public:
            //constructor
            SpacecraftOptions();
            SpacecraftOptions(const std::string& SpacecraftOptionsFileName);
            SpacecraftOptions(const std::vector<StageOptions>& stages);
            
            //destructor
            ~SpacecraftOptions() {};

            //methods
            void parse_input_file(const std::string& filename);
            void write_output_file(const std::string& filename);
            inline void remove_stage(const size_t &StageIndex)
            {
                this->Stages.erase(this->Stages.begin() + StageIndex);
                this->number_of_stages--;
            }
            
            //get/set
            void setName(const std::string& name) { this->name = name; }
            void setGlobalElectricPropellantTankCapacity(const double& GlobalElectricPropellantTankCapacity) { this->GlobalElectricPropellantTankCapacity = GlobalElectricPropellantTankCapacity; }
            void setGlobalFuelTankCapacity(const double& GlobalFuelTankCapacity) { this->GlobalChemicalFuelTankCapacity = GlobalFuelTankCapacity; }
            void setGlobalOxidizerTankCapacity(const double& GlobalOxidizerTankCapacity) {this->GlobalChemicalOxidizerTankCapacity = GlobalOxidizerTankCapacity; }
            void setGlobalDryMassBounds(const std::vector<double>& GlobalDryMassBounds) { this->GlobalDryMassBounds = GlobalDryMassBounds; }
            void setEnableGlobalElectricPropellantTankConstraint(const bool& EnableGlobalElectricPropellantTankConstraint) { this->EnableGlobalElectricPropellantTankConstraint = EnableGlobalElectricPropellantTankConstraint; }
            void setEnableGlobalChemicalPropellantTankConstraint(const bool& EnableGlobalChemicalPropellantTankConstraint) { this->EnableGlobalChemicalPropellantTankConstraint = EnableGlobalChemicalPropellantTankConstraint; }
            void setEnableGlobalDryMassConstraint(const bool& EnableGlobalDryMassConstraint) { this->EnableGlobalDryMassConstraint = EnableGlobalDryMassConstraint; }
            void setStageOptions(const size_t& StageIndex, const StageOptions& StageOptsIn) { this->Stages[StageIndex] = StageOptsIn;}

            std::string getName() const { return this->name; }
            size_t getNumberOfStages() const { return this->number_of_stages; }
            double getGlobalElectricPropellantTankCapacity() const { return this->GlobalElectricPropellantTankCapacity; }
            double getGlobalChemicalFuelTankCapacity() const { return this->GlobalChemicalFuelTankCapacity; }
            double getGlobalChemicalOxidizerTankCapacity() const { return this->GlobalChemicalOxidizerTankCapacity; }
            std::vector<double> getGlobalDryMassBounds() const { return this->GlobalDryMassBounds; }
            bool getEnableGlobalElectricPropellantTankConstraint() const { return this->EnableGlobalElectricPropellantTankConstraint; }
            bool getEnableGlobalChemicalPropellantTankConstraint() const { return this->EnableGlobalChemicalPropellantTankConstraint; }
            bool getEnableGlobalDryMassConstraint() const { return this->EnableGlobalDryMassConstraint; }
            StageOptions getStageOptions(const size_t& StageIndex) const;

        protected:
            //fields
            std::string name;
            size_t number_of_stages;
            std::vector<StageOptions> Stages;
            double GlobalElectricPropellantTankCapacity;
            double GlobalChemicalFuelTankCapacity;
            double GlobalChemicalOxidizerTankCapacity;
            bool EnableGlobalElectricPropellantTankConstraint;
            bool EnableGlobalChemicalPropellantTankConstraint;
            bool EnableGlobalDryMassConstraint;
            std::vector<double> GlobalDryMassBounds;
        };
    }//end namespace HardwareModels
}//end namespace EMTG