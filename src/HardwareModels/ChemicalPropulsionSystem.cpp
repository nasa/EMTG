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

#include "ChemicalPropulsionSystem.h"

#include <cmath>

namespace EMTG
{
    namespace HardwareModels
    {
        
        //method to compute system mass
        void ChemicalPropulsionSystem::computeSystemMass()
        {
            this->SystemMass = this->myPropulsionSystemOptions.getMassPerString() * this->myPropulsionSystemOptions.getNumberOfStrings();
        }

        void ChemicalPropulsionSystem::computeThrusterPerformance(const doubleType& deltav,
                                                                  const doubleType& mass_at_maneuver,
                                                                  const bool& ForwardManeuverFlag,
                                                                  const PropulsionSystemChoice& ThrusterType)
        {
            //for now we're assuming constant Isp, although it might some day be a function of tank fill?
            //if the biprop flag is set, we use the biprop Isp, otherwise we use the monoprop Isp
            if (ThrusterType == PropulsionSystemChoice::Biprop)
            {
                if (ForwardManeuverFlag)
                {
                    double scale = -1000.0 / (this->g0 * this->BipropIsp);
                    doubleType expfun = exp(deltav * scale);
                    doubleType mass_after_maneuver = mass_at_maneuver * expfun;
                    doubleType biprop_consumed = mass_at_maneuver - mass_after_maneuver;

                    double FuelScale = 1.0 / (1 + this->MixtureRatio);
                    double OxidizerScale = 1.0 / (1 + 1 / this->MixtureRatio);
                    this->FuelConsumedThisManeuver = biprop_consumed * FuelScale;
                    this->OxidizerConsumedThisManeuver = biprop_consumed * OxidizerScale;

                    this->dFuelConsumedThisManeuver_ddeltav = (mass_at_maneuver *  scale * expfun * FuelScale) _GETVALUE;
                    this->dOxidizerConsumedThisManeuver_ddeltav = (mass_at_maneuver * scale * expfun * OxidizerScale)_GETVALUE;
                    this->dSpacecraftMassThisManeuver_ddeltav = (mass_at_maneuver * scale * expfun) _GETVALUE;
                    this->dFuelConsumedThisManeuver_dMassAtManeuver = FuelScale * (1.0 - expfun) _GETVALUE;
                    this->dOxidizerConsumedThisManeuver_dMassAtManeuver = OxidizerScale * (1.0 - expfun) _GETVALUE;
                    this->dMassAfterManeuver_dMassAtManeuver = expfun _GETVALUE;
                }
                else //maneuver is backward in time
                {
                    double scale = 1000.0 / (this->g0 * this->BipropIsp);
                    doubleType expfun = exp(deltav * scale);
                    doubleType mass_before_maneuver = mass_at_maneuver * expfun;
                    doubleType biprop_consumed = mass_before_maneuver - mass_at_maneuver;

                    double FuelScale = 1.0 / (1 + this->MixtureRatio);
                    double OxidizerScale = this->MixtureRatio / (1.0 + this->MixtureRatio);
                    this->FuelConsumedThisManeuver = biprop_consumed * FuelScale;
                    this->OxidizerConsumedThisManeuver = biprop_consumed * OxidizerScale;

                    this->dFuelConsumedThisManeuver_ddeltav = (mass_at_maneuver *  scale * expfun * FuelScale )_GETVALUE;
                    this->dOxidizerConsumedThisManeuver_ddeltav = (mass_at_maneuver * scale * expfun * OxidizerScale )_GETVALUE;
                    this->dSpacecraftMassThisManeuver_ddeltav = (mass_at_maneuver * scale * expfun) _GETVALUE;
                    this->dFuelConsumedThisManeuver_dMassAtManeuver = FuelScale * (expfun - 1.0) _GETVALUE;
                    this->dOxidizerConsumedThisManeuver_dMassAtManeuver = OxidizerScale * (expfun - 1.0) _GETVALUE;
                    this->dMassAfterManeuver_dMassAtManeuver = expfun _GETVALUE;
                }
            }
            else if (ThrusterType == PropulsionSystemChoice::Monoprop)
            {
                if (ForwardManeuverFlag)
                {
                    double scale = -1000.0 / (this->g0 * this->MonopropIsp);
                    doubleType expfun = exp(deltav * scale);
                    doubleType mass_after_maneuver = mass_at_maneuver * expfun;

                    this->FuelConsumedThisManeuver = (mass_at_maneuver - mass_after_maneuver);
                    this->OxidizerConsumedThisManeuver = 0.0;

                    this->dFuelConsumedThisManeuver_ddeltav = (mass_at_maneuver * scale * expfun )_GETVALUE;
                    this->dOxidizerConsumedThisManeuver_ddeltav = 0.0;
                    this->dSpacecraftMassThisManeuver_ddeltav = (mass_at_maneuver * scale * expfun) _GETVALUE;
                    this->dFuelConsumedThisManeuver_dMassAtManeuver = (1.0 - expfun) _GETVALUE;
                    this->dOxidizerConsumedThisManeuver_dMassAtManeuver = 0.0;
                    this->dMassAfterManeuver_dMassAtManeuver = expfun _GETVALUE;
                }
                else //maneuver is backward in time
                {
                    double scale = 1000.0 / (this->g0 * this->MonopropIsp);
                    doubleType expfun = exp(deltav * scale);
                    doubleType mass_before_maneuver = mass_at_maneuver * expfun;

                    this->FuelConsumedThisManeuver = (mass_before_maneuver - mass_at_maneuver);
                    this->OxidizerConsumedThisManeuver = 0.0;

                    this->dFuelConsumedThisManeuver_ddeltav = (mass_at_maneuver * scale * expfun )_GETVALUE;
                    this->dOxidizerConsumedThisManeuver_ddeltav = 0.0;
                    this->dSpacecraftMassThisManeuver_ddeltav = (mass_at_maneuver * scale * expfun) _GETVALUE;
                    this->dFuelConsumedThisManeuver_dMassAtManeuver = (expfun - 1.0) _GETVALUE;
                    this->dOxidizerConsumedThisManeuver_dMassAtManeuver = 0.0;
                    this->dMassAfterManeuver_dMassAtManeuver = expfun _GETVALUE;
                }
            }
            else
            {
                throw std::invalid_argument("Invalid thruster type " + PropulsionSystemChoiceStrings[ThrusterType] + " for chemical maneuver. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__)); 
            }
        }
    }//end namespace HardwareModels
}//end namespace EMTG