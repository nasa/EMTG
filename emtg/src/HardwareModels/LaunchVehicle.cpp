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

#include "LaunchVehicle.h"

#include <iostream>

namespace EMTG
{
    namespace HardwareModels
    {
        //constructor
        LaunchVehicle::LaunchVehicle()
        {
            LaunchVehicleOptions dummyLaunchVehicleOptions;
            this->initialize(dummyLaunchVehicleOptions);
        }

        LaunchVehicle::LaunchVehicle(const LaunchVehicleOptions& launchvehicleoptions)
        {
            this->initialize(launchvehicleoptions);
        }

        //initialize method
        void LaunchVehicle::initialize(const LaunchVehicleOptions& launchvehicleoptions)
        {
            this->name = launchvehicleoptions.getName();
            this->myOptions = launchvehicleoptions;
        }

        //method to compute performance and performance derivative
        void LaunchVehicle::computePerformance(const doubleType& C3, const double& LV_margin)
        {
            switch (this->myOptions.getModelType())
            {
                case 0: //polynomial
                {
                    this->DeliveredMass = 0.0;
                    this->dmdC3 = 0.0;
                    doubleType C3power = 1.0;
                    double C3power_derivative = 0.0;

                    for (size_t p = 0; p < this->myOptions.getNumberOfCoefficients(); ++p)
                    {
                        this->DeliveredMass += this->myOptions.getCoefficient(p) * C3power;
                        this->dmdC3 += p * this->myOptions.getCoefficient(p) * C3power_derivative;
                        C3power_derivative = C3power _GETVALUE;
                        C3power *= C3;
                    }
                    break;
                }
                default:
                    throw std::invalid_argument("LaunchVehicle::ModelType " + std::to_string(this->myOptions.getModelType()) + " is not currently implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            this->DeliveredMass *= (1.0 - LV_margin);
            this->DeliveredMass -= this->myOptions.getAdapterMass();
            this->dmdC3 *= (1.0 - LV_margin);
        }
    }//end namespace HardwareModels
}//end namespace EMTG