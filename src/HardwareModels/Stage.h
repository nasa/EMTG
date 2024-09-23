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


#pragma once

#include <string>
#include <vector>

#include "PowerSystem.h"
#include "ElectricPropulsionSystem.h"
#include "ChemicalPropulsionSystem.h"
#include "HardwareBase.h"
#include "StageOptions.h"

namespace EMTG
{
    namespace HardwareModels
    {
        class Stage : public HardwareBase
        {
        public:
            //constructor
            Stage(const StageOptions& stageoptions);

            //destructor
            ~Stage() {};

            //methods
            void initialize(const StageOptions& stageoptions);
            void computeDryMass();
            void computePowerState(const doubleType& r_AU, const doubleType& current_epoch);
            void computeElectricPropulsionPerformance(const double& DutyCycle);
            void computeElectricPropulsionPerformance(const double& DutyCycle, const doubleType& u_command);
            void computeChemicalPropulsionPerformance(const doubleType& deltav,
                                                      const doubleType& mass_at_maneuver,
                                                      const bool& ForwardManeuverFlag,
                                                      const PropulsionSystemChoice& ThrusterType);

            //output function
            void output_mass_information(std::ofstream& outputfile);

            //get/set
            inline StageOptions& getOptions() { return this->myStageOptions; }

            inline doubleType getEPthrust() const { return this->MyElectricPropulsionSystem.getThrust(); }
            inline doubleType getEPIsp() const { return this->MyElectricPropulsionSystem.getIsp(); }
            inline doubleType getEPMassFlowRate() const { return this->MyElectricPropulsionSystem.getMassFlowRate(); }
            inline double getEPdTdP() const { return this->MyElectricPropulsionSystem.getdTdP(); }
            inline double getEPdTdu_command() const { return this->MyElectricPropulsionSystem.getdTdu_command(); }
            inline double getEPdMassFlowRatedP() const { return this->MyElectricPropulsionSystem.getdMassFlowRatedP(); }
            inline double getEPdMassFlowRatedu_command() const { return this->MyElectricPropulsionSystem.getdMassFlowRatedu_command(); }
            inline double getEPdIspdP() const { return this->MyElectricPropulsionSystem.getdIspdP(); }
            inline double getEPdIspdu_command() const { return this->MyElectricPropulsionSystem.getdIspdu_command(); }
            inline doubleType getEPActivePower() const { return this->MyElectricPropulsionSystem.getActivePower(); }
            inline size_t getEPNumberOfActiveThrusters() const { return this->MyElectricPropulsionSystem.getNumberOfActiveThrusters(); }
            inline size_t getEPThrottleLevel() const { return this->MyElectricPropulsionSystem.getThrottleLevel(); }
            inline std::string getEPThrottleLevelString() const { return this->MyElectricPropulsionSystem.getThrottleLevelString(); }

            inline doubleType getchemthrust() const { return this->MyChemicalPropulsionSystem.getThrust(); }
            inline double getg0() const { return this->MyChemicalPropulsionSystem.getg0(); }

            inline double getBipropIsp() const { return this->MyChemicalPropulsionSystem.getBipropIsp(); }
            inline double getMonopropIsp() const { return this->MyChemicalPropulsionSystem.getMonopropIsp(); }
            inline doubleType getFuelConsumedThisManeuver() const { return this->MyChemicalPropulsionSystem.getFuelConsumedThisManeuver(); }
            inline double getdFuelConsumedThisManeuver_ddeltav() const { return this->MyChemicalPropulsionSystem.getdFuelConsumedThisManeuver_ddeltav(); }
            inline doubleType getOxidizerConsumedThisManeuver() const { return this->MyChemicalPropulsionSystem.getOxidizerConsumedThisManeuver(); }
            inline double getdOxidizerConsumedThisManeuver_ddeltav() const { return this->MyChemicalPropulsionSystem.getdOxidizerConsumedThisManeuver_ddeltav(); }
            inline double getdSpacecraftMassThisManeuver_ddeltav() const { return this->MyChemicalPropulsionSystem.getdSpacecraftMassThisManeuver_ddeltav(); }
            inline double get_dFuelConsumedThisManeuver_dMassAtManeuver() const { return this->MyChemicalPropulsionSystem.get_dFuelConsumedThisManeuver_dMassAtManeuver(); }
            inline double get_dOxidizerConsumedThisManeuver_dMassAtManeuver() const { return this->MyChemicalPropulsionSystem.get_dOxidizerConsumedThisManeuver_dMassAtManeuver(); }
            inline double get_dMassAfterManeuver_dMassAtManeuver() const { return this->MyChemicalPropulsionSystem.get_dMassAfterManeuver_dMassAtManeuver(); }

            inline double getAdapterMass() const { return this->myStageOptions.getAdapterMass(); }
            inline double getDryMass() const { return this->DryMass; }

            inline doubleType getProducedPower() const { return this->ProducedPower; }
            inline doubleType getBusPower() const { return this->BusPower; }
            inline doubleType getAvailablePower() const { return this->AvailablePower; }
            inline double getdPdr() const { return this->MyPowerSystem.getdPdr(); }
            inline double getdPdt() const { return this->MyPowerSystem.getdPdt(); }


            //EMTG-specific propellant tank things
            //****************************************
            void setChemicalFuelTank_Findex(const size_t& ChemicalFuelTank_Findex) { this->ChemicalFuelTank_Findex = ChemicalFuelTank_Findex; }
            void setChemicalOxidizerTank_Findex(const size_t& ChemicalOxidizerTank_Findex) { this->ChemicalOxidizerTank_Findex = ChemicalOxidizerTank_Findex; }
            void setElectricPropellantTank_Findex(const size_t& ElectricPropellantTank_Findex) { this->ElectricPropellantTank_Findex = ElectricPropellantTank_Findex; }
            void setChemicalFuelTank_Xindices(const std::vector<size_t>& ChemicalFuelTank_Xindices) { this->ChemicalFuelTank_Xindices = ChemicalFuelTank_Xindices; }
            void setChemicalOxidizerTank_Xindices(const std::vector<size_t>& ChemicalOxidizerTank_Xindices) { this->ChemicalOxidizerTank_Xindices = ChemicalOxidizerTank_Xindices; }
            void setElectricPropellantTank_Xindices(const std::vector<size_t>& ElectricPropellantTank_Xindices) { this->ElectricPropellantTank_Xindices = ElectricPropellantTank_Xindices; }

            void setChemicalFuelTank_Xscaleranges(const std::vector<double>& Xupperbounds, const std::vector<double>& Xlowerbounds);
            void setChemicalOxidizerTank_Xscaleranges(const std::vector<double>& Xupperbounds, const std::vector<double>& Xlowerbounds);
            void setElectricPropellantTank_Xscaleranges(const std::vector<double>& Xupperbounds, const std::vector<double>& Xlowerbounds);
            inline std::vector<size_t> getChemicalFuelTank_Xindices() const { return this->ChemicalFuelTank_Xindices; }
            inline std::vector<size_t> getChemicalOxidizerTank_Xindices() const { return this->ChemicalOxidizerTank_Xindices; }
            inline std::vector<size_t> getElectricPropellantTank_Xindices() const { return this->ElectricPropellantTank_Xindices; }

            void appendChemicalFuelTank_Xindices(const size_t& ChemicalFuelTank_Xindices) { this->ChemicalFuelTank_Xindices.push_back(ChemicalFuelTank_Xindices); }
            void appendChemicalOxidizerTank_Xindices(const size_t& ChemicalOxidizerTank_Xindices) { this->ChemicalOxidizerTank_Xindices.push_back(ChemicalOxidizerTank_Xindices); }
            void appendElectricPropellantTank_Xindices(const size_t& ElectricPropellantTank_Xindices) { this->ElectricPropellantTank_Xindices.push_back(ElectricPropellantTank_Xindices); }

            void appendChemicalFuelTank_Xscaleranges(const double& Xupperbounds, const double& Xlowerbounds);
            void appendChemicalOxidizerTank_Xscaleranges(const double& Xupperbounds, const double& Xlowerbounds);
            void appendElectricPropellantTank_Xscaleranges(const double& Xupperbounds, const double& Xlowerbounds);

            void computeChemicalFuelState(const std::vector<doubleType>& X, const double& PercentMargin);
            void computeChemicalOxidizerState(const std::vector<doubleType>& X, const double& PercentMargin);
            void computeElectricPropellantState(const std::vector<doubleType>& X, const double& PercentMargin);

            doubleType getChemicalFuelUsed() const { return this->ChemicalFuelUsed; }
            doubleType getChemicalOxidizerUsed() const { return this->ChemicalOxidizerUsed; }
            doubleType getElectricPropellantUsed() const { return this->ElectricPropellantUsed; }
            doubleType getChemicalFuelMargin() const { return this->ChemicalFuelMargin; }
            doubleType getChemicalOxidizerMargin() const { return this->ChemicalOxidizerMargin; }
            doubleType getElectricPropellantMargin() const { return this->ElectricPropellantMargin; }

            //EMTG-specific dry mass things
            void setDryMass_Findex(const size_t& DryMass_Findex) { this->DryMass_Findex = DryMass_Findex; }
            void setXindex_StageFinalMass(const size_t& Xindex_StageFinalMass) { this->Xindex_StageFinalMass = Xindex_StageFinalMass; }
            size_t getXindex_StageFinalMass() const { return this->Xindex_StageFinalMass; }
            void setXscale_StageFinalMass(const double& Xscale_StageFinalMass) { this->Xscale_StageFinalMass = Xscale_StageFinalMass; }
            doubleType getStageRequiredFinalMass() const { return this->StageRequiredFinalMass; }
            void computeStageRequiredFinalMass();
            void populateDryMassDerivatives(std::vector<double>& G, const double& PercentElectricPropellantMargin, const double& PercentChemicalPropellantMargin);
            void appendDryMass_Gindices(const size_t& DryMass_Gindices) { this->DryMass_Gindices.push_back(DryMass_Gindices); }

        protected:     
            //fields
            StageOptions myStageOptions;
            ElectricPropulsionSystem MyElectricPropulsionSystem;
            ChemicalPropulsionSystem MyChemicalPropulsionSystem;
            PowerSystem MyPowerSystem;
            double DryMass;
            doubleType ProducedPower;
            doubleType BusPower;
            doubleType AvailablePower;
            doubleType ActivePower;


            doubleType ChemicalFuelUsed;
            doubleType ChemicalOxidizerUsed;
            doubleType ElectricPropellantUsed;
            doubleType ChemicalFuelMargin;
            doubleType ChemicalOxidizerMargin;
            doubleType ElectricPropellantMargin;

            //EMTG-specific propellant tank things
            bool hasElectricPropellantTankConstraint;
            bool hasChemicalFuelTankConstraint;
            bool hasChemicalOxidizerTankConstraint;
            size_t ChemicalFuelTank_Findex;
            size_t ChemicalOxidizerTank_Findex;
            size_t ElectricPropellantTank_Findex;
            std::vector<size_t> ChemicalFuelTank_Xindices;
            std::vector<size_t> ChemicalOxidizerTank_Xindices;
            std::vector<size_t> ElectricPropellantTank_Xindices;
            std::vector<double> ChemicalFuelTank_Xscaleranges;
            std::vector<double> ChemicalOxidizerTank_Xscaleranges;
            std::vector<double> ElectricPropellantTank_Xscaleranges;

            //EMTG-specific dry mass things
            doubleType StageRequiredFinalMass;
            size_t DryMass_Findex;
            size_t Xindex_StageFinalMass;
            double Xscale_StageFinalMass;
            std::vector<size_t> DryMass_Gindices;
        };
    }//end namespace HardwareModels
}//end namespace EMTG