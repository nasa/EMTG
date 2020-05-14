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

//container class for spacecraft
//Jacob Englander 10-28-2016


#include "Spacecraft.h"

namespace EMTG
{
    namespace HardwareModels
    {
        //constructor
        Spacecraft::Spacecraft()
        {
            //create a dummy spacecraft options and intialize from it
            this->initialize(SpacecraftOptions());
        }
        
        Spacecraft::Spacecraft(const SpacecraftOptions& spacecraftoptions)
        {
            this->initialize(spacecraftoptions);
        }


        Spacecraft::Spacecraft(const std::string& spacecraftoptionsfilename)
        {
            SpacecraftOptions tempSpacecraftOptions(spacecraftoptionsfilename);
            this->initialize(tempSpacecraftOptions);
        }

        //initialize method
        void Spacecraft::initialize(const SpacecraftOptions& spacecraftoptions)
        {
            this->mySpacecraftOptions = spacecraftoptions;

            this->number_of_stages = spacecraftoptions.getNumberOfStages();

            for (size_t stageIndex = 0; stageIndex < this->number_of_stages; ++stageIndex)
            {
                this->Stages.push_back(Stage(spacecraftoptions.getStageOptions(stageIndex)));
            }

            this->resetStaging();

            this->name = spacecraftoptions.getName();
        }

        //reset staging
        void Spacecraft::resetStaging()
        {
            this->ActiveStageIndex = 0;
            this->ActiveStage = &(this->Stages[this->ActiveStageIndex]);

            this->computeCurrentDryMass();
        }

        //perform staging event
        void Spacecraft::performStagingEvent()
        {
            //Check to see if there are more stages left. If not, ignore the staging event request.
            if (this->ActiveStageIndex < this->number_of_stages - 1)
            {
                ++this->ActiveStageIndex;
                this->ActiveStage = &this->Stages[this->ActiveStageIndex];

                this->computeCurrentDryMass();
            }
        }

        //set stage
        void Spacecraft::setActiveStage(const size_t& stageIndex)
        {
            this->ActiveStageIndex = stageIndex;

            this->ActiveStage = &this->Stages[this->ActiveStageIndex];

            this->computeCurrentDryMass();
        }

        //compute current dry mass
        void Spacecraft::computeCurrentDryMass()
        {
            //the current dry mass is the sum of all dry masses in the current stage and later, plus the sum of all adapter masses AFTER this stage
            this->CurrentDryMass = 0.0;

            for (size_t StageIndex = this->ActiveStageIndex; StageIndex < this->number_of_stages; ++StageIndex)
            {
                this->Stages[StageIndex].computeDryMass();
                this->CurrentDryMass += this->Stages[StageIndex].getDryMass();

                if (StageIndex > this->ActiveStageIndex)
                    this->CurrentDryMass += this->Stages[StageIndex].getAdapterMass();
            }
        }
        
        //compute power state of the current stage
        void Spacecraft::computePowerState(const doubleType& r_AU, const doubleType& current_epoch)
        {
            this->ActiveStage->computePowerState(r_AU, current_epoch);
        }

        //compute propulsion state of the current stage
        void Spacecraft::computeElectricPropulsionPerformance(const double& DutyCycle)
        {
            this->ActiveStage->computeElectricPropulsionPerformance(DutyCycle);
        }

        void Spacecraft::computeElectricPropulsionPerformance(const double& DutyCycle, const doubleType& u_command)
        {
            this->ActiveStage->computeElectricPropulsionPerformance(DutyCycle, u_command);
        }

        void Spacecraft::computeChemicalPropulsionPerformance(const doubleType& deltav,
            const doubleType& mass_at_maneuver,
            const bool& ForwardManeuverFlag,
            const PropulsionSystemChoice& ThrusterType)
        {
            this->ActiveStage->computeChemicalPropulsionPerformance(deltav, mass_at_maneuver, ForwardManeuverFlag, ThrusterType);
        }

        //EMTG-specific propellant tank things
        void Spacecraft::computePropellantState(const std::vector<doubleType>& X, const double& ElectricPropellantMargin, const double& ChemicalPropellantMargin)
        {
            this->computeChemicalFuelState(X, ChemicalPropellantMargin);
            this->computeChemicalOxidizerState(X, ChemicalPropellantMargin);
            this->computeElectricPropellantState(X, ElectricPropellantMargin);
        }

        void Spacecraft::computeChemicalFuelState(const std::vector<doubleType>& X, const double& PercentMargin)
        {
            this->ChemicalFuelUsed = 0.0;
            this->ChemicalFuelMargin = 0.0;

            for (size_t StageIndex = 0; StageIndex < this->number_of_stages; ++StageIndex)
            {
                this->Stages[StageIndex].computeChemicalFuelState(X, PercentMargin);
                this->ChemicalFuelUsed += this->Stages[StageIndex].getChemicalFuelUsed();
                this->ChemicalFuelMargin += this->Stages[StageIndex].getChemicalFuelMargin();
            }
        }

        void Spacecraft:: computeChemicalOxidizerState(const std::vector<doubleType>& X, const double& PercentMargin)
        {
            this->ChemicalOxidizerUsed = 0.0;
            this->ChemicalOxidizerMargin = 0.0;

            for (size_t StageIndex = 0; StageIndex < this->number_of_stages; ++StageIndex)
            {
                this->Stages[StageIndex].computeChemicalOxidizerState(X, PercentMargin);
                this->ChemicalOxidizerUsed += this->Stages[StageIndex].getChemicalOxidizerUsed();
                this->ChemicalOxidizerMargin += this->Stages[StageIndex].getChemicalOxidizerMargin();
            }
        }

        void Spacecraft:: computeElectricPropellantState(const std::vector<doubleType>& X, const double& PercentMargin)
        {
            this->ElectricPropellantUsed = 0.0;
            this->ElectricPropellantMargin = 0.0;

            for (size_t StageIndex = 0; StageIndex < this->number_of_stages; ++StageIndex)
            {
                this->Stages[StageIndex].computeElectricPropellantState(X, PercentMargin);
                this->ElectricPropellantUsed += this->Stages[StageIndex].getElectricPropellantUsed();
                this->ElectricPropellantMargin += this->Stages[StageIndex].getElectricPropellantMargin();
            }
        }

        void Spacecraft::output_mass_information(std::ofstream& outputfile)
        {
            outputfile << std::endl;

            for (size_t stageIndex = 0; stageIndex < this->number_of_stages; ++stageIndex)
            {
                this->Stages[stageIndex].output_mass_information(outputfile);
            }

            outputfile << std::endl;
            outputfile << "Spacecraft: Electric propellant used (kg): " << this->ElectricPropellantUsed _GETVALUE << std::endl;
            outputfile << "Spacecraft: Chemical fuel used (kg): " << this->ChemicalFuelUsed _GETVALUE << std::endl;
            outputfile << "Spacecraft: Chemical oxidizer used (kg): " << this->ChemicalOxidizerUsed _GETVALUE << std::endl;
            outputfile << "Spacecraft: Electric propellant margin (kg): " << this->ElectricPropellantMargin _GETVALUE << std::endl;
            outputfile << "Spacecraft: Chemical fuel margin (kg): " << this->ChemicalFuelMargin _GETVALUE << std::endl;
            outputfile << "Spacecraft: Chemical oxidizer margin (kg): " << this->ChemicalOxidizerMargin _GETVALUE << std::endl;
            outputfile << "Spacecraft: Total electric propellant (kg): " << (this->ElectricPropellantUsed + this->ElectricPropellantMargin) _GETVALUE << std::endl;
            outputfile << "Spacecraft: Total chemical fuel (kg): " << (this->ChemicalFuelUsed + this->ChemicalFuelMargin) _GETVALUE << std::endl;
            outputfile << "Spacecraft: Total chemical oxidizer (kg): " << (this->ChemicalOxidizerUsed + this->ChemicalOxidizerMargin) _GETVALUE << std::endl;
        }
    }//end namespace HardwareModels
}//end namespace EMTG