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

//base propulsion system class
//Jacob Englander 11/2/2016

#pragma once

#include "HardwareBase.h"
#include "PropulsionSystemOptions.h"

namespace EMTG
{
    namespace HardwareModels
    {
        class PropulsionSystem : public HardwareBase
        {
        public:
            //constructor
            PropulsionSystem() { this->myPropulsionSystemOptions = PropulsionSystemOptions(); };
            PropulsionSystem(const PropulsionSystemOptions& propulsionsystemoptions) : myPropulsionSystemOptions(propulsionsystemoptions) {};

            //destructor

            //methods
            void initialize(const PropulsionSystemOptions& propulsionsystemoptions) { this->myPropulsionSystemOptions = propulsionsystemoptions; };
            virtual void computeSystemMass() = 0;

            //get
            inline doubleType getIsp() const { return this->Isp; }
            inline doubleType getThrust() const { return this->Thrust; }
            inline double getSystemMass() const { return this->SystemMass; }
            inline double getg0() const { return this->g0; }

        protected:
            //fields common to all propulsion system
            doubleType Isp;
            doubleType Thrust;
            double SystemMass;
            PropulsionSystemOptions myPropulsionSystemOptions;
            double g0 = 9.806649999999999423;
        };
    }//end namespace HardwareModels
}//end namespace EMTG