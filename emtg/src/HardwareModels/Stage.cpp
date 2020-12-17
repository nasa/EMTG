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

//container class for stage
//Jacob Englander 10-28-2016


#include "Stage.h"

namespace EMTG
{
    namespace HardwareModels
    {
        //constructor
        Stage::Stage(const StageOptions& stageoptions)
        {
            this->initialize(stageoptions);
        }
        
        //initialize method
        void Stage::initialize(const StageOptions& stageoptions)
        {
            this->myStageOptions = stageoptions;
            this->MyPowerSystem = PowerSystem(stageoptions.getPowerSystemOptions());
            this->MyElectricPropulsionSystem = ElectricPropulsionSystem(stageoptions.getElectricPropulsionSystemOptions());
            this->MyChemicalPropulsionSystem = ChemicalPropulsionSystem(stageoptions.getChemicalPropulsionSystemOptions());

            this->name = stageoptions.getName();


            this->ElectricPropellantUsed = 0.0;
            this->ChemicalFuelUsed = 0.0;
            this->ChemicalOxidizerUsed = 0.0;

            this->ElectricPropellantMargin = 0.0;
            this->ChemicalFuelMargin = 0.0;
            this->ChemicalOxidizerMargin = 0.0;

            this->ProducedPower = 0.0;
            this->BusPower = 0.0;
            this->AvailablePower = 0.0;
            this->ActivePower = 0.0;
        }

        //compute dry mass
        void Stage::computeDryMass()
        {
            this->MyChemicalPropulsionSystem.computeSystemMass();
            this->MyElectricPropulsionSystem.computeSystemMass();
            this->MyPowerSystem.computeSystemMass();

            this->DryMass = this->myStageOptions.getBaseDryMass()
                + this->MyChemicalPropulsionSystem.getSystemMass()
                + this->MyElectricPropulsionSystem.getSystemMass()
                + this->MyPowerSystem.getSystemMass();
        }

        void Stage::computePowerState(const doubleType& r_AU, const doubleType& current_epoch)
        {
            this->MyPowerSystem.evaluate_available_power(r_AU, current_epoch);
            this->ProducedPower = this->MyPowerSystem.getProducedPower();
            this->BusPower = this->MyPowerSystem.getBusPower();
            this->AvailablePower = this->MyPowerSystem.getAvailablePower();
        }

        void Stage::computeElectricPropulsionPerformance(const double& DutyCycle)
        {
            this->MyElectricPropulsionSystem.computeThrusterPerformance(this->AvailablePower, DutyCycle, this->myStageOptions.getThrottleLogic());
            this->ActivePower = this->MyElectricPropulsionSystem.getActivePower();
        }

        void Stage::computeElectricPropulsionPerformance(const double& DutyCycle, const doubleType& u_command)
        {
            this->MyElectricPropulsionSystem.computeThrusterPerformance(this->AvailablePower, DutyCycle, this->myStageOptions.getThrottleLogic(), u_command);
            this->ActivePower = this->MyElectricPropulsionSystem.getActivePower();
        }

        void Stage::computeChemicalPropulsionPerformance(const doubleType& deltav,
            const doubleType& mass_at_maneuver,
            const bool& ForwardManeuverFlag,
            const PropulsionSystemChoice& ThrusterType)
        {
            this->MyChemicalPropulsionSystem.computeThrusterPerformance(deltav, mass_at_maneuver, ForwardManeuverFlag, ThrusterType);
        }

        //EMTG-specific propellant tank things
        void Stage::setChemicalFuelTank_Xscaleranges(const std::vector<double>& Xupperbounds, const std::vector<double>& Xlowerbounds)
        {
            this->ChemicalFuelTank_Xscaleranges.clear();
            for (size_t i = 0; i < this->ChemicalFuelTank_Xindices.size(); ++i)
            {
                ChemicalFuelTank_Xscaleranges.push_back(Xupperbounds[this->ChemicalFuelTank_Xindices[i]] - Xlowerbounds[this->ChemicalFuelTank_Xindices[i]]);
            }
        }

        void Stage::setChemicalOxidizerTank_Xscaleranges(const std::vector<double>& Xupperbounds, const std::vector<double>& Xlowerbounds)
        {
            this->ChemicalOxidizerTank_Xscaleranges.clear();
            for (size_t i = 0; i < this->ChemicalOxidizerTank_Xindices.size(); ++i)
            {
                ChemicalOxidizerTank_Xscaleranges.push_back(Xupperbounds[this->ChemicalOxidizerTank_Xindices[i]] - Xlowerbounds[this->ChemicalOxidizerTank_Xindices[i]]);
            }
        }

        void Stage::setElectricPropellantTank_Xscaleranges(const std::vector<double>& Xupperbounds, const std::vector<double>& Xlowerbounds)
        {
            this->ElectricPropellantTank_Xscaleranges.clear();
            for (size_t i = 0; i < this->ElectricPropellantTank_Xindices.size(); ++i)
            {
                ElectricPropellantTank_Xscaleranges.push_back(Xupperbounds[this->ElectricPropellantTank_Xindices[i]] - Xlowerbounds[this->ElectricPropellantTank_Xindices[i]]);
            }
        }


        void Stage::appendChemicalFuelTank_Xscaleranges(const double& Xupperbounds, const double& Xlowerbounds)
        {
            this->ChemicalFuelTank_Xscaleranges.push_back(Xupperbounds - Xlowerbounds);
        }

        void Stage::appendChemicalOxidizerTank_Xscaleranges(const double& Xupperbounds, const double& Xlowerbounds)
        {
            this->ChemicalOxidizerTank_Xscaleranges.push_back(Xupperbounds - Xlowerbounds);
        }

        void Stage::appendElectricPropellantTank_Xscaleranges(const double& Xupperbounds, const double& Xlowerbounds)
        {
            this->ElectricPropellantTank_Xscaleranges.push_back(Xupperbounds - Xlowerbounds);
        }


        void Stage::computeChemicalFuelState(const std::vector<doubleType>& X, const double& PercentMargin)
        {
            this->ChemicalFuelUsed = 0.0;
            for (size_t i = 0; i < this->ChemicalFuelTank_Xindices.size(); ++i)
            {
                this->ChemicalFuelUsed += X[this->ChemicalFuelTank_Xindices[i]];
            }
            this->ChemicalFuelMargin = this->ChemicalFuelUsed * PercentMargin;
        }

        void Stage::computeChemicalOxidizerState(const std::vector<doubleType>& X, const double& PercentMargin)
        {
            this->ChemicalOxidizerUsed = 0.0;
            for (size_t i = 0; i < this->ChemicalOxidizerTank_Xindices.size(); ++i)
            {
                this->ChemicalOxidizerUsed += X[this->ChemicalOxidizerTank_Xindices[i]];
            }
            this->ChemicalOxidizerMargin = this->ChemicalOxidizerUsed * PercentMargin;
        }

        void Stage::computeElectricPropellantState(const std::vector<doubleType>& X, const double& PercentMargin)
        {
            this->ElectricPropellantUsed = 0.0;
            for (size_t i = 0; i < this->ElectricPropellantTank_Xindices.size(); ++i)
            {
                this->ElectricPropellantUsed += X[this->ElectricPropellantTank_Xindices[i]];
            }
            this->ElectricPropellantMargin = this->ElectricPropellantUsed * PercentMargin;
        }

        void Stage::computeStageRequiredFinalMass()
        {
            this->StageRequiredFinalMass = this->DryMass + this->ElectricPropellantMargin + this->ChemicalFuelMargin + this->ChemicalOxidizerMargin;
        }

        void Stage::populateDryMassDerivatives(std::vector<double>& G, const double& ElectricPropellantMargin, const double& ChemicalPropellantMargin)
        {
            size_t Gindex = 0;
            //with respect to virtual electric tanks
            for (size_t i = 0; i < this->ElectricPropellantTank_Xindices.size(); ++i)
            {
                G[this->DryMass_Gindices[Gindex++]] = ElectricPropellantMargin * this->ElectricPropellantTank_Xscaleranges[i] / this->DryMass;
            }

            //with respect to virtual chemical fuel tanks
            for (size_t i = 0; i < this->ChemicalFuelTank_Xindices.size(); ++i)
            {
                G[this->DryMass_Gindices[Gindex++]] = ChemicalPropellantMargin * this->ChemicalFuelTank_Xscaleranges[i] / this->DryMass;
            }

            //with respect to virtual chemical oxidizer tanks
            for (size_t i = 0; i < this->ChemicalOxidizerTank_Xindices.size(); ++i)
            {
                G[this->DryMass_Gindices[Gindex++]] = ChemicalPropellantMargin * this->ChemicalOxidizerTank_Xscaleranges[i] / this->DryMass;
            }

            //with respect to final mass
            G[this->DryMass_Gindices[Gindex++]] = -this->Xscale_StageFinalMass / this->DryMass;
        }

        void Stage::output_mass_information(std::ofstream& outputfile)
        {
            outputfile << std::endl;
            outputfile << this->name << ": Electric propellant used (kg): " << this->ElectricPropellantUsed _GETVALUE << std::endl;
            outputfile << this->name << ": Chemical fuel used (kg): " << this->ChemicalFuelUsed _GETVALUE << std::endl;
            outputfile << this->name << ": Chemical oxidizer used (kg): " << this->ChemicalOxidizerUsed _GETVALUE << std::endl;
            outputfile << this->name << ": Electric propellant margin (kg): " << this->ElectricPropellantMargin _GETVALUE << std::endl;
            outputfile << this->name << ": Chemical fuel margin (kg): " << this->ChemicalFuelMargin _GETVALUE << std::endl;
            outputfile << this->name << ": Chemical oxidizer margin (kg): " << this->ChemicalOxidizerMargin _GETVALUE << std::endl;
            outputfile << this->name << ": Total electric propellant (kg): " << (this->ElectricPropellantUsed + this->ElectricPropellantMargin) _GETVALUE << std::endl;
            outputfile << this->name << ": Total chemical fuel (kg): " << (this->ChemicalFuelUsed + this->ChemicalFuelMargin) _GETVALUE << std::endl;
            outputfile << this->name << ": Total chemical oxidizer (kg): " << (this->ChemicalOxidizerUsed + this->ChemicalOxidizerMargin) _GETVALUE << std::endl;
        }
    }//end namespace HardwareModels
}//end namespace EMTG