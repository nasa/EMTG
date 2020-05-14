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
//target spec line
//this class is for one EMTG-MIRAGE target spec line
//
//contains the following information
//<SAMPLE>,<CBODY>,<FRAME>,<EPOCH(ET)>,<X[km]>,<Y[km]>,<Z[km]>,<VX[km/s]>,<VY[km/s]>,<VZ[km/s]>,<B.R[km]>,<B.T [km]>
//


#pragma once

#include "doubleType.h"
#include "EMTG_Matrix.h"

#include <string>
#include <fstream>

namespace EMTG
{
    class target_spec_line
    {
    public:
        //constructors
        target_spec_line() : BdotR(0.0), BdotT(0.0) {};

        target_spec_line(const std::string& eventID,
            const std::string& frame,
            const std::string& CentralBody,
            const doubleType& epoch,
            const math::Matrix<doubleType>& State);

        target_spec_line(const std::string& eventID,
            const std::string& frame,
            const std::string& CentralBody,
            const doubleType& epoch,
            const math::Matrix<doubleType>& State,
            const doubleType& BdotR,
            const doubleType& BdotT);

        //methods
        void initialize(const std::string& eventID,
            const std::string& frame,
            const std::string& CentralBody,
            const doubleType& epoch,
            const math::Matrix<doubleType>& State);

        void initialize(const std::string& eventID,
            const std::string& frame,
            const std::string& CentralBody,
            const doubleType& epoch,
            const math::Matrix<doubleType>& State,
            const doubleType& BdotR,
            const doubleType& BdotT);

        void set_eventID(const std::string& eventID) { this->eventID = eventID; }

        void write(std::ofstream& outputfile);

    protected:
        std::string eventID;
        std::string frame;
        std::string CentralBody;

        doubleType TargetEpochETseconds;
        doubleType TargetEpochETJD;
        std::string TargetEpochETGregorian;
        math::Matrix<doubleType> State;
        doubleType BdotR;
        doubleType BdotT;
    };//end class maneuver_spec_item
}//end namespace EMTG
