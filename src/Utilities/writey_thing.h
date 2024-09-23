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

#pragma once

#include "doubleType.h"

#include <string>
#include <vector>

#include <EMTG_Matrix.h>
#include "missionoptions.h"
#include "universe.h"

#include "maneuver_spec_line.h"
#include "target_spec_line.h"

namespace EMTG
{
    class writey_thing
    {
    public:
        //constructor
        writey_thing() {};

        writey_thing(missionoptions* myOptions, Astrodynamics::universe* myUniverse);

    protected:
        virtual void initialize(missionoptions* myOptions, Astrodynamics::universe* myUniverse);
    
        //write method
        virtual void write_output_line(std::ofstream& outputfile,
            size_t& eventcount,
            const std::string& event_type,
            const std::string& boundary_name,
            const doubleType& timestep_size,
            const doubleType& flyby_altitude,
            const doubleType& BdotR,
            const doubleType& BdotT,
            const doubleType& angle1,
            const doubleType& angle2,
            const doubleType& C3,
            const math::Matrix<doubleType>& state,
            const math::Matrix<doubleType>& dV,
            const math::Matrix<doubleType>& ThrustVector,
            const doubleType& dVmag,
            const doubleType& Thrust,
            const doubleType& Isp,
            const doubleType& AvailPower,
            const doubleType& mdot,
            const int& number_of_active_engines,
            const doubleType& active_power,
            const std::string& ThrottleLevel);

        virtual void write_ephemeris_line(std::ofstream& outputfile,
            const math::Matrix<doubleType>& state);

        virtual void write_ephemeris_line(std::ofstream& outputfile,
            const math::Matrix<doubleType>& state,
            const math::Matrix<doubleType>& ControlVector,
            const doubleType& ThrustMagnitude,
            const doubleType& MassFlowRate,
            const doubleType& Isp,
            const int& NumberOfActiveThrusters,
            const doubleType& ActivePower,
            const std::string& ThrottleLevel);

        virtual void write_ephemeris_state(std::ofstream& outputfile,
            const math::Matrix<doubleType>& state);

        //fields
        missionoptions* myOptions;
        Astrodynamics::universe* myUniverse;
    };//end clss writey_thing
}