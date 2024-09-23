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

//chemical propulsion system class
//Jacob Englander 11/2/2016

#pragma once

#include "PropulsionSystem.h"

namespace EMTG
{
    namespace HardwareModels
    {
        class ChemicalPropulsionSystem : public PropulsionSystem
        {
        public:
            //constructor
            ChemicalPropulsionSystem() : PropulsionSystem() {};
            ChemicalPropulsionSystem(const PropulsionSystemOptions& propulsionsystemoptions) : PropulsionSystem(propulsionsystemoptions), 
                MonopropIsp(propulsionsystemoptions.getMinimumOrMonopropIsp()),
                BipropIsp(propulsionsystemoptions.getConstantIsp()),
                MixtureRatio(propulsionsystemoptions.getMixtureRatio()),
                Thrust(propulsionsystemoptions.getConstantThrust())
            {};

            //get
            inline doubleType getThrust() const { return this->Thrust; }
            inline double getBipropIsp() const { return this->BipropIsp; }
            inline double getMonopropIsp() const { return this->MonopropIsp; }

            inline doubleType getFuelConsumedThisManeuver() const { return this->FuelConsumedThisManeuver; }
            inline double getdFuelConsumedThisManeuver_ddeltav() const { return this->dFuelConsumedThisManeuver_ddeltav; }
            inline doubleType getOxidizerConsumedThisManeuver() const { return this->OxidizerConsumedThisManeuver; }
            inline double getdOxidizerConsumedThisManeuver_ddeltav() const { return this->dOxidizerConsumedThisManeuver_ddeltav; }
            inline double getdSpacecraftMassThisManeuver_ddeltav() const { return this->dSpacecraftMassThisManeuver_ddeltav; }

            inline double get_dFuelConsumedThisManeuver_dMassAtManeuver() const { return this->dFuelConsumedThisManeuver_dMassAtManeuver; }
            inline double get_dOxidizerConsumedThisManeuver_dMassAtManeuver() const { return this->dOxidizerConsumedThisManeuver_dMassAtManeuver; }
            inline double get_dMassAfterManeuver_dMassAtManeuver() const { return this->dMassAfterManeuver_dMassAtManeuver; } //note "after" means "in the direction of propagation"

            //evaluate
            void computeThrusterPerformance(const doubleType& deltav, 
                                            const doubleType& mass_at_maneuver,
                                            const bool& ForwardManeuverFlag,
                                            const PropulsionSystemChoice& ThrusterType);
            void computeSystemMass();

        protected:

            //fields
            double Thrust;
            double BipropIsp;
            double MonopropIsp;
            double MixtureRatio;
            doubleType FuelConsumedThisManeuver;
            double dFuelConsumedThisManeuver_ddeltav;
            doubleType OxidizerConsumedThisManeuver;
            double dOxidizerConsumedThisManeuver_ddeltav;
            double dSpacecraftMassThisManeuver_ddeltav;
            double dFuelConsumedThisManeuver_dMassAtManeuver;
            double dOxidizerConsumedThisManeuver_dMassAtManeuver;
            double dMassAfterManeuver_dMassAtManeuver; //note "after" means "in the direction of propagation"
        };
    }//end namespace HardwareModels
}//end namespace EMTG