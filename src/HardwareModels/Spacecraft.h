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


#pragma once

#include <string>
#include <vector>

#include "Stage.h"
#include "SpacecraftOptions.h"
#include "HardwareBase.h"

namespace EMTG
{
    namespace HardwareModels
    {
        class Spacecraft : public HardwareBase
        {
        public:
            //constructor
            Spacecraft();
            Spacecraft(const SpacecraftOptions& spacecraftoptions);
            Spacecraft(const std::string& spacecraftoptionsfilename);

            //methods
            void initialize(const SpacecraftOptions& spacecraftoptions);
            void computeCurrentDryMass();
            void computePowerState(const doubleType& r_AU, const doubleType& current_epoch);
            void computeElectricPropulsionPerformance(const double& DutyCycle);
            void computeElectricPropulsionPerformance(const double& DutyCycle, const doubleType& u_command);
            void computeChemicalPropulsionPerformance(const doubleType& deltav,
                                                      const doubleType& mass_at_maneuver,
                                                      const bool& ForwardManeuverFlag,
                                                      const PropulsionSystemChoice& ThrusterType);
            void resetStaging();
            void performStagingEvent();
            void setActiveStage(const size_t& stageIndex);
            inline size_t getNumberOfStages() const { return this->number_of_stages; }

            //output function
            void output_mass_information(std::ofstream& outputfile);
            
            inline SpacecraftOptions& getSpacecraftOptions() { return this->mySpacecraftOptions; }
            inline StageOptions& getCurrentStageOptions() { return this->ActiveStage->getOptions(); }
            inline StageOptions& getStageOptions(const size_t& stageIndex) { return this->Stages[stageIndex].getOptions(); }

            inline doubleType getEPthrust() const { return this->ActiveStage->getEPthrust(); }
            inline doubleType getEPIsp() const { return this->ActiveStage->getEPIsp(); }
            inline doubleType getEPMassFlowRate() const { return this->ActiveStage->getEPMassFlowRate(); }
            inline double getEPdTdP() const { return this->ActiveStage->getEPdTdP(); }
            inline double getEPdTdu_command() const { return this->ActiveStage->getEPdTdu_command(); }
            inline double getEPdMassFlowRatedP() const { return this->ActiveStage->getEPdMassFlowRatedP(); }
            inline double getEPdMassFlowRatedu_command() const { return this->ActiveStage->getEPdMassFlowRatedu_command(); }
            inline double getEPdIspdP() const { return this->ActiveStage->getEPdIspdP(); }
            inline double getEPdIspdu_command() const { return this->ActiveStage->getEPdIspdu_command(); }
            inline doubleType getEPActivePower() const { return this->ActiveStage->getEPActivePower(); }
            inline size_t getEPNumberOfActiveThrusters() const { return this->ActiveStage->getEPNumberOfActiveThrusters(); }
            inline size_t getEPThrottleLevel() const { return this->ActiveStage->getEPThrottleLevel(); }
            inline std::string getEPThrottleLevelString() const { return this->ActiveStage->getEPThrottleLevelString(); }
            inline doubleType getchemthrust() const { return this->ActiveStage->getchemthrust(); }
            inline double getBipropIsp() const { return this->ActiveStage->getBipropIsp(); }
            inline double getMonopropIsp() const { return this->ActiveStage->getMonopropIsp(); }
            inline double getg0() const { return this->ActiveStage->getg0(); }

            inline doubleType getChemFuelConsumedThisManeuver() const { return this->ActiveStage->getFuelConsumedThisManeuver(); }
            inline double getdFuelConsumedThisManeuver_ddeltav() const { return this->ActiveStage->getdFuelConsumedThisManeuver_ddeltav(); }
            inline doubleType getChemOxidizerConsumedThisManeuver() const { return this->ActiveStage->getOxidizerConsumedThisManeuver(); }
            inline double getdOxidizerConsumedThisManeuver_ddeltav() const { return this->ActiveStage->getdOxidizerConsumedThisManeuver_ddeltav(); }
            inline double getdSpacecraftMassThisManeuver_ddeltav() const { return this->ActiveStage->getdSpacecraftMassThisManeuver_ddeltav(); }
            inline double get_dFuelConsumedThisManeuver_dMassAtManeuver() const { return this->ActiveStage->get_dFuelConsumedThisManeuver_dMassAtManeuver(); }
            inline double get_dOxidizerConsumedThisManeuver_dMassAtManeuver() const { return this->ActiveStage->get_dOxidizerConsumedThisManeuver_dMassAtManeuver(); }
            inline double get_dMassAfterManeuver_dMassAtManeuver() const { return this->ActiveStage->get_dMassAfterManeuver_dMassAtManeuver(); }
            inline doubleType getChemicalFuelMargin() const { return this->ChemicalFuelMargin; }
            inline doubleType getChemicalOxidizerMargin() const { return this->ChemicalOxidizerMargin; }
            inline doubleType getElectricPropellantMargin() const { return this->ElectricPropellantMargin; }

            inline double getCurrentDryMass() const { return this->CurrentDryMass; }


            inline doubleType getProducedPower() const { return this->ActiveStage->getProducedPower(); }
            inline doubleType getBusPower() const { return this->ActiveStage->getBusPower(); }
            inline doubleType getAvailablePower() const { return this->ActiveStage->getAvailablePower(); }
            inline double getdPdr() const { return this->ActiveStage->getdPdr(); }
            inline double getdPdt() const { return this->ActiveStage->getdPdt(); }


            //EMTG-specific propellant tank things
            void setChemicalFuelTank_Findex(const size_t& stageIndex, const size_t& ChemicalFuelTank_Findex) { this->Stages[stageIndex].setChemicalFuelTank_Findex(ChemicalFuelTank_Findex); }
            void setChemicalOxidizerTank_Findex(const size_t& stageIndex, const size_t& ChemicalOxidizerTank_Findex) { this->Stages[stageIndex].setChemicalOxidizerTank_Findex(ChemicalOxidizerTank_Findex); }
            void setElectricPropellantTank_Findex(const size_t& stageIndex, const size_t& ElectricPropellantTank_Findex) { this->Stages[stageIndex].setElectricPropellantTank_Findex(ElectricPropellantTank_Findex); }
            void setChemicalFuelTank_Xindices(const size_t& stageIndex, const std::vector<size_t>& ChemicalFuelTank_Xindices) { this->Stages[stageIndex].setChemicalFuelTank_Xindices(ChemicalFuelTank_Xindices); }
            void setChemicalOxidizerTank_Xindices(const size_t& stageIndex, const std::vector<size_t>& ChemicalOxidizerTank_Xindices) { this->Stages[stageIndex].setChemicalOxidizerTank_Xindices(ChemicalOxidizerTank_Xindices); }
            void setElectricPropellantTank_Xindices(const size_t& stageIndex, const std::vector<size_t>& ElectricPropellantTank_Xindices) { this->Stages[stageIndex].setElectricPropellantTank_Xindices(ElectricPropellantTank_Xindices); }
            void setChemicalFuelTank_Xscaleranges(const size_t& stageIndex, const std::vector<double>& Xupperbounds, const std::vector<double>& Xlowerbounds) { this->Stages[stageIndex].setChemicalFuelTank_Xscaleranges(Xupperbounds, Xlowerbounds); }
            void setChemicalOxidizerTank_Xscaleranges(const size_t& stageIndex, const std::vector<double>& Xupperbounds, const std::vector<double>& Xlowerbounds) { this->Stages[stageIndex].setChemicalOxidizerTank_Xscaleranges(Xupperbounds, Xlowerbounds); }
            void setElectricPropellantTank_Xscaleranges(const size_t& stageIndex, const std::vector<double>& Xupperbounds, const std::vector<double>& Xlowerbounds) { this->Stages[stageIndex].setElectricPropellantTank_Xscaleranges(Xupperbounds, Xlowerbounds); }

            inline std::vector<size_t> getChemicalFuelTank_Xindices(const size_t& stageIndex) const { return this->Stages[stageIndex].getChemicalFuelTank_Xindices(); }
            inline std::vector<size_t> getChemicalOxidizerTank_Xindices(const size_t& stageIndex) const { return this->Stages[stageIndex].getChemicalOxidizerTank_Xindices(); }
            inline std::vector<size_t> getElectricPropellantTank_Xindices(const size_t& stageIndex) const { return this->Stages[stageIndex].getElectricPropellantTank_Xindices(); }

            void setGlobalChemicalFuelTank_Findex(const size_t& ChemicalFuelTank_Findex) { this->GlobalChemicalFuelTank_Findex = ChemicalFuelTank_Findex; }
            void setGlobalChemicalOxidizerTank_Findex(const size_t& ChemicalOxidizerTank_Findex) { this->GlobalChemicalOxidizerTank_Findex = ChemicalOxidizerTank_Findex; }
            void setGlobalElectricPropellantTank_Findex(const size_t& ElectricPropellantTank_Findex) { this->GlobalElectricPropellantTank_Findex = ElectricPropellantTank_Findex; }
            void setGlobalDryMassConstraint_Findex(const size_t& GlobalDryMassConstraint_Findex) { this->GlobalDryMassConstraint_Findex = GlobalDryMassConstraint_Findex; }

            void appendChemicalFuelTank_Xindices(const size_t& stageIndex, const size_t& ChemicalFuelTank_Xindices) { this->Stages[stageIndex].appendChemicalFuelTank_Xindices(ChemicalFuelTank_Xindices); }
            void appendChemicalOxidizerTank_Xindices(const size_t& stageIndex, const size_t& ChemicalOxidizerTank_Xindices) { this->Stages[stageIndex].appendChemicalOxidizerTank_Xindices(ChemicalOxidizerTank_Xindices); }
            void appendElectricPropellantTank_Xindices(const size_t& stageIndex, const size_t& ElectricPropellantTank_Xindices) { this->Stages[stageIndex].appendElectricPropellantTank_Xindices(ElectricPropellantTank_Xindices); }
            void appendChemicalFuelTank_Xscaleranges(const size_t& stageIndex, const double& Xupperbounds, const double& Xlowerbounds) { this->Stages[stageIndex].appendChemicalFuelTank_Xscaleranges(Xupperbounds, Xlowerbounds); }
            void appendChemicalOxidizerTank_Xscaleranges(const size_t& stageIndex, const double& Xupperbounds, const double& Xlowerbounds) { this->Stages[stageIndex].appendChemicalOxidizerTank_Xscaleranges(Xupperbounds, Xlowerbounds); }
            void appendElectricPropellantTank_Xscaleranges(const size_t& stageIndex, const double& Xupperbounds, const double& Xlowerbounds) { this->Stages[stageIndex].appendElectricPropellantTank_Xscaleranges(Xupperbounds, Xlowerbounds); }

            inline std::vector<size_t> getGlobalChemicalFuelTank_Xindices() const { return this->GlobalChemicalFuelTank_Xindices; }
            inline std::vector<size_t> getGlobalChemicalOxidizerTank_Xindices() const { return this->GlobalChemicalOxidizerTank_Xindices; }
            inline std::vector<size_t> getGlobalElectricPropellantTank_Xindices() const { return this->GlobalElectricPropellantTank_Xindices; }
            inline std::vector<size_t> getGlobalDryMassConstraint_Xindices() const { return this->GlobalDryMassConstraint_Xindices; }

            void appendGlobalChemicalFuelTank_Xindices(const size_t& ChemicalFuelTank_Xindices) { this->GlobalChemicalFuelTank_Xindices.push_back(ChemicalFuelTank_Xindices); }
            void appendGlobalChemicalOxidizerTank_Xindices(const size_t& ChemicalOxidizerTank_Xindices) { this->GlobalChemicalOxidizerTank_Xindices.push_back(ChemicalOxidizerTank_Xindices); }
            void appendGlobalElectricPropellantTank_Xindices(const size_t& ElectricPropellantTank_Xindices) { this->GlobalElectricPropellantTank_Xindices.push_back(ElectricPropellantTank_Xindices); }
            void appendGlobalDryMassConstraint_Xindices(const size_t& GlobalDryMassConstraint_Xindices) { this->GlobalDryMassConstraint_Xindices.push_back(GlobalDryMassConstraint_Xindices); }

            void appendGlobalChemicalFuelTank_Xscaleranges(const double& Xupperbounds, const double& Xlowerbounds) { this->GlobalChemicalFuelTank_Xscaleranges.push_back(Xupperbounds - Xlowerbounds); }
            void appendGlobalChemicalOxidizerTank_Xscaleranges(const double& Xupperbounds, const double& Xlowerbounds) { this->GlobalChemicalOxidizerTank_Xscaleranges.push_back(Xupperbounds - Xlowerbounds); }
            void appendGlobalElectricPropellantTank_Xscaleranges(const double& Xupperbounds, const double& Xlowerbounds) { this->GlobalElectricPropellantTank_Xscaleranges.push_back(Xupperbounds - Xlowerbounds); }
            void appendGlobalDryMassConstraint_Xscaleranges(const double& GlobalDryMassConstraint_Xscaleranges) { this->GlobalDryMassConstraint_Xscaleranges.push_back(GlobalDryMassConstraint_Xscaleranges); }

            void computePropellantState(const std::vector<doubleType>& X, const double& ElectricPropellantMargin, const double& ChemicalPropellantMargin);
            void computeChemicalFuelState(const std::vector<doubleType>& X, const double& PercentMargin);
            void computeChemicalOxidizerState(const std::vector<doubleType>& X, const double& PercentMargin);
            void computeElectricPropellantState(const std::vector<doubleType>& X, const double& PercentMargin);

            doubleType getChemicalFuelUsed(const size_t& StageIndex) const { return this->Stages[StageIndex].getChemicalFuelUsed(); }
            doubleType getChemicalOxidizerUsed(const size_t& StageIndex) const { return this->Stages[StageIndex].getChemicalOxidizerUsed(); }
            doubleType getElectricPropellantUsed(const size_t& StageIndex) const { return this->Stages[StageIndex].getElectricPropellantUsed(); }
            
            doubleType getGlobalChemicalFuelUsed() const { return this->ChemicalFuelUsed; }
            doubleType getGlobalChemicalOxidizerUsed() const { return this->ChemicalOxidizerUsed; }
            doubleType getGlobalElectricPropellantUsed() const { return this->ElectricPropellantUsed; }


            //EMTG-specific dry mass things
            void setDryMass_Findex(const size_t& StageIndex, const size_t& DryMass_Findex) { this->Stages[StageIndex].setDryMass_Findex(DryMass_Findex); }
            void setXindex_StageFinalMass(const size_t& StageIndex, const size_t& Xindex_StageFinalMass) { this->Stages[StageIndex].setXindex_StageFinalMass(Xindex_StageFinalMass); }
            size_t getXindex_StageFinalMass(const size_t& StageIndex) const { return this->Stages[StageIndex].getXindex_StageFinalMass(); }
            void setXscale_StageFinalMass(const size_t& StageIndex, const double& Xscale_StageFinalMass) { this->Stages[StageIndex].setXscale_StageFinalMass(Xscale_StageFinalMass); }
            void computeStageRequiredFinalMass(const size_t& StageIndex) { this->Stages[StageIndex].computeStageRequiredFinalMass(); }
            doubleType getStageRequiredFinalMass(const size_t& StageIndex) const { return this->Stages[StageIndex].getStageRequiredFinalMass(); }
            void populateDryMassDerivatives(const size_t& StageIndex, std::vector<double>& G, const double& PercentElectricPropellantMargin, const double& PercentChemicalPropellantMargin) { this->Stages[StageIndex].populateDryMassDerivatives(G, PercentElectricPropellantMargin, PercentChemicalPropellantMargin); }
            void appendDryMass_Gindices(const size_t& StageIndex, const size_t& DryMass_Gindices) { this->Stages[StageIndex].appendDryMass_Gindices(DryMass_Gindices); }



            //fields
        protected:
            SpacecraftOptions mySpacecraftOptions;
            size_t number_of_stages;
            std::vector< Stage > Stages;
            size_t ActiveStageIndex;
            Stage* ActiveStage;

            double CurrentDryMass;
            doubleType ProducedPower;
            doubleType BusPower;
            doubleType AvailablePower;

            doubleType EPthrust;
            doubleType EPIsp;
            doubleType EPMassFlowRate;
            doubleType ChemicalThrust;
            doubleType ChemicalMassFlowRate;

            doubleType ChemicalFuelUsed;
            doubleType ChemicalOxidizerUsed;
            doubleType ElectricPropellantUsed;
            doubleType ChemicalFuelMargin;
            doubleType ChemicalOxidizerMargin;
            doubleType ElectricPropellantMargin;

            //EMTG-specific propellant tank things
            size_t GlobalChemicalFuelTank_Findex;
            size_t GlobalChemicalOxidizerTank_Findex;
            size_t GlobalElectricPropellantTank_Findex;
            std::vector<size_t> GlobalChemicalFuelTank_Xindices;
            std::vector<size_t> GlobalChemicalOxidizerTank_Xindices;
            std::vector<size_t> GlobalElectricPropellantTank_Xindices;
            std::vector<double> GlobalChemicalFuelTank_Xscaleranges;
            std::vector<double> GlobalChemicalOxidizerTank_Xscaleranges;
            std::vector<double> GlobalElectricPropellantTank_Xscaleranges;

            //EMTG-specific global dry mass constraint
            size_t GlobalDryMassConstraint_Findex;
            std::vector<size_t> GlobalDryMassConstraint_Xindices;
            std::vector<double> GlobalDryMassConstraint_Xscaleranges;
        };
    }//end namespace HardwareModels
}//end namespace EMTG