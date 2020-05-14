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

//Jacob Englander 6/22/2018
//maneuver spec item
//this class is for one item that goes into a maneuver spec line
//
//contains the following information
//<FRAME>,<EPOCH(ET)>,<THRX>,<THRY>,<THRZ>,<THRMAG[N]>,<SMASS[kg]>,<MDOT[kg/s]>,<DUTY>,<FMASS[kg]>,<DV[km/s]>,<DUR[s]>
//where <THRX>, <THRY>, <THRZ> are unit vector components
//


#pragma once

#include "doubleType.h"
#include "EMTG_Matrix.h"

#include <string>
#include <fstream>

namespace EMTG
{
    class maneuver_spec_item
    {
    public:
        //constructors
        maneuver_spec_item();

        maneuver_spec_item(const std::string& frame,
            const doubleType& epoch,
            const math::Matrix<doubleType>& ControlVector,
            const doubleType& StartMass,
            const doubleType& FinalMass,
            const doubleType& ThrustMagnitude,
            const doubleType& MassFlowRate,
            const doubleType& ManeuverDuration,
            const doubleType& EnforcedDutyCycle); //this is the duty cycle enforced by the problem, not the one chosen by the optimizer

        //methods
        void initialize(const std::string& frame,
            const doubleType& epoch,
            const math::Matrix<doubleType>& ControlVector,
            const doubleType& StartMass,
            const doubleType& FinalMass,
            const doubleType& ThrustMagnitude,
            const doubleType& MassFlowRate,
            const doubleType& ManeuverDuration,
            const doubleType& EnforcedDutyCycle);

        void write(std::ofstream& outputfile);

    protected:
        std::string frame;
        doubleType ManeuverStartEpochETseconds;
        doubleType ManeuverStartEpochETJD;
        std::string ManeuverStartEpochETGregorian;
        math::Matrix<doubleType> ControlVector;
        math::Matrix<doubleType> ControlUnitVector;
        doubleType StartMass;
        doubleType FinalMass;
        doubleType ThrustMagnitude;
        doubleType MassFlowRate;
        doubleType ManeuverDuration;
        doubleType EnforcedDutyCycle;
        doubleType TrueDutyCycle; 
        doubleType Deltav;
    };//end class maneuver_spec_item
}//end namespace EMTG
