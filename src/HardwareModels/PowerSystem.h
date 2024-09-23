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

#pragma once

#include "HardwareBase.h"
#include "PowerSystemOptions.h"

namespace EMTG
{
    namespace HardwareModels
    {
        class PowerSystem : public HardwareBase
        {
        public:
            //constructor
            PowerSystem() : MyPowerSystemOptions(PowerSystemOptions()) {};
            PowerSystem(const PowerSystemOptions& powersystemoptions);

            //destructor

            //methods
            void initialize(const PowerSystemOptions& powersystemoptions);

            //method to evaluate available power, including derivatives
            void evaluate_available_power(const doubleType& r, //in AU
                const doubleType& current_epoch);//in MJD seconds

            //method to compute system mass
            void computeSystemMass();            

            //get
            doubleType getProducedPower() const;
            doubleType getBusPower() const;
            doubleType getAvailablePower() const;
            double getdPdr() const;
            double getdPdt() const;
            double getSystemMass() const;

        private:

            //fields
            PowerSystemOptions MyPowerSystemOptions;
            doubleType ProducedPower;
            doubleType BusPower;
            doubleType AvailablePower;
            double dPdr;
            double dPdt;
            double SystemMass;
        };
    }//end namespace HardwareModels
}//end namespace EMTG