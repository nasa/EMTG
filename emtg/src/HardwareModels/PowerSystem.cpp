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

//power system class
//Jacob Englander 12/30/2016

#include <cmath>

#include "PowerSystem.h"

namespace EMTG
{
    namespace HardwareModels
    {
        
        //constructor
        PowerSystem::PowerSystem(const PowerSystemOptions& powersystemoptions)
        {
            this->initialize(powersystemoptions);
        }

        //initialize method
        void PowerSystem::initialize(const PowerSystemOptions& powersystemoptions)
        {
            this->MyPowerSystemOptions = powersystemoptions;
            this->ProducedPower = 0.0;
            this->BusPower = 0.0;
            this->AvailablePower = 0.0;
        }

        //evaluate
        void PowerSystem::evaluate_available_power(const doubleType& r,//in AU
                                                   const doubleType& current_epoch) //in MJD seconds
        {
            doubleType r2 = r * r;

            //radial component
            if (this->MyPowerSystemOptions.getPowerSupplyType() == SpacecraftPowerSupplyType::ConstantPower)
            {
                this->ProducedPower = this->MyPowerSystemOptions.getP0();
                this->dPdr = 0.0;
            }
            else if (this->MyPowerSystemOptions.getPowerSupplyType() == SpacecraftPowerSupplyType::Solar)
            {
                double g0 = this->MyPowerSystemOptions.getGamma(0);
                double g1 = this->MyPowerSystemOptions.getGamma(1);
                double g2 = this->MyPowerSystemOptions.getGamma(2);
                double g3 = this->MyPowerSystemOptions.getGamma(3);
                double g4 = this->MyPowerSystemOptions.getGamma(4);
                double g5 = this->MyPowerSystemOptions.getGamma(5);
                double g6 = this->MyPowerSystemOptions.getGamma(6);

                if (this->MyPowerSystemOptions.getPowerSupplyCurveType() == SpacecraftPowerSupplyCurveType::Sauer)
                {

                    this->ProducedPower = this->MyPowerSystemOptions.getP0() / r2 * ((g0 + g1 / r + g2 / r2) / (1.0 + g3 * r + g4 * r2));

                    doubleType r3 = r2*r;
                    doubleType r4 = r3*r;
                    doubleType r5 = r4*r;

                    this->dPdr = -(this->MyPowerSystemOptions.getP0() * (4 * g2 + 3 * g1*r + g3*(3 * g0*r3 + 4 * g1*r2 + 5 * g2*r) + 2 * g0*r2 + g4*(4 * g0*r4 + 5 * g1*r3 + 6 * g2*r2)) / (r5*(g4*r2 + g3*r + 1)*(g4*r2 + g3*r + 1))) _GETVALUE;
                }
                else if (this->MyPowerSystemOptions.getPowerSupplyCurveType() == SpacecraftPowerSupplyCurveType::Polynomial)
                {
                    doubleType r3 = r2*r;
                    doubleType r4 = r3*r;
                    doubleType r5 = r4*r;
                    doubleType r6 = r5*r;

                    this->ProducedPower = this->MyPowerSystemOptions.getP0() * (g0 / r2 + g1 / r + g2 + g3*r + g4*r2 + g5*r3 + g6*r4);

                    this->dPdr = this->MyPowerSystemOptions.getP0() * (-2.0*g0 / r3 - g1 / r2 + g3 + 2.0*g4*r + 3.0*g5*r2 + 4.0*g6*r3) _GETVALUE;
                }
                
            }
            else
            {
                throw std::invalid_argument("PowerSystem::invalid Spacecraft_Power_Supply_Type " + SpacecraftPowerSupplyCurveTypeStrings[this->MyPowerSystemOptions.getPowerSupplyType()] + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            //time component
            if (this->MyPowerSystemOptions.getDecayRate() > 1.0e-5)
            {
                doubleType decay_coeff = exp(-this->MyPowerSystemOptions.getDecayRate() * (current_epoch - this->MyPowerSystemOptions.getPowerSystemDecayRefEpoch()) / (365.25 * 86400.0));
                this->ProducedPower *= decay_coeff;

                dPdr *= decay_coeff _GETVALUE;
                dPdt = -this->MyPowerSystemOptions.getDecayRate() / (365.25 * 86400.0) * this->ProducedPower _GETVALUE;
                
            }
            else
                dPdt = 0.0;

            //bus power component
            if (this->MyPowerSystemOptions.getBusPowerType() == SpacecraftBusPowerType::TypeA_Quadratic)
            {
                this->BusPower = this->MyPowerSystemOptions.getBusCoefficient(0) + this->MyPowerSystemOptions.getBusCoefficient(1) / r + this->MyPowerSystemOptions.getBusCoefficient(2) / r2;
                this->AvailablePower = this->ProducedPower - this->BusPower;

                this->dPdr -= (-this->MyPowerSystemOptions.getBusCoefficient(1) / r2 - 2.0 * this->MyPowerSystemOptions.getBusCoefficient(2) / (r2*r)) _GETVALUE;
            }
            else if (this->MyPowerSystemOptions.getBusPowerType() == SpacecraftBusPowerType::TypeB_Conditional)
            {
                if (this->ProducedPower > this->MyPowerSystemOptions.getBusCoefficient(0))
                {
                    this->BusPower = this->MyPowerSystemOptions.getBusCoefficient(0);
                    this->AvailablePower = this->ProducedPower - this->BusPower;
                }
                else
                {
                    this->BusPower = this->MyPowerSystemOptions.getBusCoefficient(0) + this->MyPowerSystemOptions.getBusCoefficient(1) * (this->MyPowerSystemOptions.getBusCoefficient(2) - this->ProducedPower);
                    this->AvailablePower = this->ProducedPower - this->BusPower;

                    dPdr -= -this->MyPowerSystemOptions.getBusCoefficient(1) * this->dPdr;
                }
            }
            else
            {
                throw std::invalid_argument("PowerSystem::invalid Spacecraft_Bus_Power_Type " + SpacecraftBusPowerTypeStrings[this->MyPowerSystemOptions.getBusPowerType()] + ".Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            //apply power margin
            this->AvailablePower *= (1.0 - this->MyPowerSystemOptions.getPowerMargin());
            this->dPdr *= (1.0 - this->MyPowerSystemOptions.getPowerMargin());
            this->dPdt *= (1.0 - this->MyPowerSystemOptions.getPowerMargin());
        }//evaluate_available_power

        //compute system mass
        void PowerSystem::computeSystemMass()
        {
            this->SystemMass = this->MyPowerSystemOptions.getP0() * this->MyPowerSystemOptions.getmass_per_kw();
        }

        //get
        doubleType PowerSystem::getProducedPower() const
        {
            return this->ProducedPower;
        }

        doubleType PowerSystem::getBusPower() const
        {
            return this->BusPower;
        }

        doubleType PowerSystem::getAvailablePower() const
        {
            return this->AvailablePower;
        }

        double PowerSystem::getdPdr() const
        {
            return this->dPdr;
        }

        double PowerSystem::getdPdt() const
        {
            return this->dPdt;
        }

        double PowerSystem::getSystemMass() const
        {
            return this->SystemMass;
        }

    }//end namespace HardwareModels
}//end namespace EMTG