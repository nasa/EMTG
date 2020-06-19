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

#include "writey_thing.h"


#include "mjd_to_mdyhms.h"
#include "SpiceUsr.h"

#include <iostream>
#include <fstream>
#include <sstream>

namespace EMTG
{
    //constructor

    writey_thing::writey_thing(missionoptions* myOptions, Astrodynamics::universe* myUniverse)
    {
        this->initialize(myOptions, myUniverse);
    }

    void writey_thing::initialize(missionoptions* myOptions, Astrodynamics::universe* myUniverse)
    {
        this->myOptions = myOptions;
        this->myUniverse = myUniverse;
    }

    //******************************************output methods
    void writey_thing::write_output_line(std::ofstream& outputfile,
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
            const std::string& ThrottleLevel)
    {
        //create a 3-element storage vector that will be used every time something needs to be rotated to the local frame
        math::Matrix<doubleType> display_vector(3, 1, 0.0);
        math::Matrix<doubleType> rot_in_vec(3, 1, 0.0);

        outputfile.width(5); outputfile << eventcount++;
        outputfile.width(3); outputfile << " | ";

        //output the event epoch in both MJD and MM/DD/YYYY
        double current_epoch_MJD = state(7) _GETVALUE / 86400.0;
        outputfile.width(16); outputfile.setf(std::ios::fixed, std::ios::floatfield); outputfile.precision(8); outputfile << current_epoch_MJD + 2400000.5;
        outputfile.width(3); outputfile << " | ";
        int month, day, year, hrs, mins;
        double secs;
        mjd_to_mdyhms(current_epoch_MJD, &month, &day, &year, &hrs, &mins, &secs);
        std::stringstream datestream;
        datestream << month << "/" << day << "/" << year;
        outputfile.width(11); outputfile << datestream.str();
        outputfile.width(3); outputfile << " | ";

        outputfile.width(12); outputfile << event_type;
        outputfile.width(3); outputfile << " | ";

        outputfile.width(25); outputfile << boundary_name;
        outputfile.width(3); outputfile << " | ";

        outputfile.width(15); outputfile << timestep_size _GETVALUE;
        outputfile.width(3); outputfile << " | ";

        //rp, BdotR, and BdotT
        if (event_type == "upwr_flyby"
            || event_type == "pwr_flyby"
            || event_type == "periapse"
            || event_type == "zeroflyby"
            || (event_type == "intercept" && flyby_altitude > 0.0))
        {
            outputfile.width(19); outputfile << flyby_altitude _GETVALUE;
            outputfile.width(3); outputfile << " | ";
            if (BdotR == BdotR)
            {
                outputfile.width(19); outputfile << BdotR _GETVALUE;
            }
            else
            {
                outputfile.width(19); outputfile << "-";
            }
            outputfile.width(3); outputfile << " | ";
            if (BdotT == BdotT)
            {
                outputfile.width(19); outputfile << BdotT _GETVALUE;
            }
            else
            {
                outputfile.width(19); outputfile << "-";
            }
            outputfile.width(3); outputfile << " | ";
        }
        else if (boundary_name == "Hyp-arrival")
        {
            outputfile.width(19); outputfile << state.getSubMatrix1D(0, 2).norm()  _GETVALUE;
            outputfile.width(3); outputfile << " | ";
            if (BdotR == BdotR)
            {
                outputfile.width(19); outputfile << BdotR _GETVALUE;
            }
            else
            {
                outputfile.width(19); outputfile << "-";
            }
            outputfile.width(3); outputfile << " | ";
            if (BdotT == BdotT)
            {
                outputfile.width(19); outputfile << BdotT _GETVALUE;
            }
            else
            {
                outputfile.width(19); outputfile << "-";
            }
            outputfile.width(3); outputfile << " | ";
        }
        else if (event_type.find("spiral") < 1024
            || (event_type == "launch" && boundary_name == "periapse"))
        {
            outputfile.width(19); outputfile << flyby_altitude _GETVALUE;
            outputfile.width(3); outputfile << " | ";
            outputfile.width(19); outputfile << "-";
            outputfile.width(3); outputfile << " | ";
            outputfile.width(19); outputfile << "-";
            outputfile.width(3); outputfile << " | ";
        }
        else
        {
            outputfile.width(19); outputfile << "-";
            outputfile.width(3); outputfile << " | ";
            outputfile.width(19); outputfile << "-";
            outputfile.width(3); outputfile << " | ";
            outputfile.width(19); outputfile << "-";
            outputfile.width(3); outputfile << " | ";
        }

        //control angles
        if (!(event_type == "coast" || event_type == "Scoast" || event_type == "TCM" || event_type == "nav-coast"))
        {
            double print_angle;
            if (angle1 != angle1)
                print_angle = 0.0;
            else
                print_angle = angle1 _GETVALUE;

            outputfile.precision(3); outputfile.width(8); outputfile << fmod(print_angle, 2.0*math::PI) * 180.0 / math::PI;
            outputfile.width(3); outputfile << " | ";

            if (angle2 != angle2)
                print_angle = 0.0;
            else
                print_angle = angle2 _GETVALUE;

            outputfile.precision(3); outputfile.width(8); outputfile << fmod(print_angle, 2.0*math::PI) * 180.0 / math::PI;
            outputfile.width(3); outputfile << " | ";
        }
        else
        {
            outputfile.precision(3); outputfile.width(8); outputfile << 0.0;
            outputfile.width(3); outputfile << " | ";
            outputfile.precision(3); outputfile.width(8); outputfile << 0.0;
            outputfile.width(3); outputfile << " | ";
        }

        //C3
        if (event_type == "upwr_flyby"
            || event_type == "pwr_flyby"
            || event_type == "launch"
            || event_type == "intercept"
            || event_type == "insertion"
            || event_type == "departure"
            || event_type == "interface"
            || event_type == "zeroflyby"
            || event_type == "rendezvous"
            || event_type == "momtransfer")
        {
            outputfile.precision(5);
            outputfile.width(14); outputfile << C3 _GETVALUE;
        }
        else
        {
            outputfile.width(14); outputfile << "-";
        }
        outputfile.width(3); outputfile << " | ";


        //state at event +
        outputfile.precision(8);
        for (size_t k = 0; k < 3; ++k)
            rot_in_vec(k) = state(k);
        this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, rot_in_vec, this->myOptions->output_file_frame, display_vector, state(7));
        for (size_t k = 0; k < 3; ++k)
        {
            outputfile.width(19); outputfile << display_vector(k) _GETVALUE;
            outputfile.width(3); outputfile << " | ";
        }
        outputfile.precision(8);
        for (size_t k = 0; k < 3; ++k)
            rot_in_vec(k) = state(k + 3);
        this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, rot_in_vec, this->myOptions->output_file_frame, display_vector, state(7));
        for (size_t k = 0; k < 3; ++k)
        {
            outputfile.width(19); outputfile << display_vector(k) _GETVALUE;
            outputfile.width(3); outputfile << " | ";
        }


        //deltaV - we can't track deltaV vectors for flybys or TCMs
        if (event_type == "upwr_flyby"
            || event_type == "pwr_flyby"
            || event_type == "LT_rndzvs"
            || event_type == "zeroflyby"
            || event_type == "LT_spiral"
            || event_type == "begin_spiral"
            || event_type == "end_spiral"
            || event_type == "TCM")
        {
            for (int k = 0; k < 3; ++k)
            {
                outputfile.width(19); outputfile << "-";
                outputfile.width(3); outputfile << " | ";
            }
        }
        else
        {
            outputfile.precision(12);
            for (size_t k = 0; k < 3; ++k)
                rot_in_vec(k) = dV(k);
            this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, rot_in_vec, this->myOptions->output_file_frame, display_vector, state(7));
            for (int k = 0; k < 3; ++k)
            {
                outputfile.width(19); outputfile << display_vector(k) _GETVALUE;
                outputfile.width(3); outputfile << " | ";
            }
        }

        //Thrust
        outputfile.precision(12);
        if (event_type == "SFthrust"
            || event_type == "FBLTthrust"
            || event_type == "SSFthrust"
            || event_type == "PSBIthrust"
            || event_type == "PSFBthrust")
        {
            for (size_t k = 0; k < 3; ++k)
                rot_in_vec(k) = ThrustVector(k) _GETVALUE;
            this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, rot_in_vec, this->myOptions->output_file_frame, display_vector, state(7));

            for (int k = 0; k < 3; ++k)
            {
                outputfile.width(19); outputfile << display_vector(k) _GETVALUE;
                outputfile.width(3); outputfile << " | ";
            }
        }
        else
        {
            for (int k = 0; k < 3; ++k)
            {
                outputfile.width(19); outputfile << "-";
                outputfile.width(3); outputfile << " | ";
            }
        }

        //dV magnitude
        outputfile.precision(5);
        outputfile.width(17); outputfile << dVmag _GETVALUE;
        outputfile.width(3); outputfile << " | ";

        //thrust, Isp, power
        outputfile.precision(5);
        if (event_type == "coast"
            || event_type == "Scoast"
            || event_type == "force-coast"
            || event_type == "nav-coast"
            || event_type == "upwr_flyby"
            || event_type == "intercept"
            || event_type == "interface"
            || event_type == "LT_rndzvs"
            || event_type == "match_point"
            || event_type == "match-vinf"
            || event_type == "zeroflyby"
            || event_type == "waiting"
            || event_type == "mission_end")
        {
            outputfile.width(14); outputfile << "-";
            outputfile.width(3); outputfile << " | ";
        }
        else if (Thrust > 1.0e-6)
        {
            outputfile.width(14); outputfile << Thrust _GETVALUE;
            outputfile.width(3); outputfile << " | ";
        }
        else
        {
            outputfile.width(14); outputfile << "impulse";
            outputfile.width(3); outputfile << " | ";
        }

        outputfile.precision(0);
        if (event_type == "coast"
            || event_type == "Scoast"
            || event_type == "force-coast"
            || event_type == "nav-coast"
            || event_type == "upwr_flyby"
            || event_type == "periapse"
            || event_type == "intercept"
            || event_type == "interface"
            || event_type == "LT_rndzvs"
            || event_type == "match_point"
            || event_type == "match-vinf"
            || event_type == "zeroflyby"
            || event_type == "waiting"
            || event_type == "mission_end")
        {
            outputfile.width(14); outputfile << "-";
            outputfile.width(3); outputfile << " | ";
        }
        else if (Isp > 0 && event_type != "launch")
        {
            outputfile.width(14); outputfile << Isp _GETVALUE;
            outputfile.width(3); outputfile << " | ";
        }
        else if (event_type == "launch")
        {
            outputfile.width(14); outputfile << "LV-supplied";
            outputfile.width(3); outputfile << " | ";
        }
        else if (event_type == "departure" || event_type == "separation")
        {
            outputfile.width(14); outputfile << "-";
            outputfile.width(3); outputfile << " | ";
        }
        else
        {
            outputfile.width(14); outputfile << "UNHANDLED EVENT TYPE";
            outputfile.width(3); outputfile << " | ";
        }


        outputfile.precision(5);
        if (AvailPower > 0)
        {
            outputfile.width(14); outputfile << AvailPower _GETVALUE;
            outputfile.width(3); outputfile << " | ";
        }
        else
        {
            outputfile.width(14); outputfile << "-";
            outputfile.width(3); outputfile << " | ";
        }

        if (event_type == "SFthrust"
            || event_type == "FBLTthrust"
            || event_type == "SSFthrust"
            || event_type == "PSBIthrust"
            || event_type == "PSFBthrust"
            || event_type == "begin_spiral"
            || event_type == "end_spiral"
            || event_type == "LT_spiral")
        {
            outputfile.precision(8);
            outputfile.width(19); outputfile << std::scientific << mdot _GETVALUE << std::fixed;

        }
        else
        {
            outputfile.width(19); outputfile << "-";
        }
        outputfile.width(3); outputfile << " | ";

        //mass
        outputfile.width(14); outputfile.precision(4); outputfile << state(6) _GETVALUE;

        outputfile.width(3); outputfile << " | ";

        //number of active engines
        outputfile.width(14);
        if (event_type == "SFthrust"
            || event_type == "FBLTthrust"
            || event_type == "SSFthrust"
            || event_type == "PSBIthrust"
            || event_type == "PSFBthrust"
            || event_type == "begin_spiral"
            || event_type == "end_spiral"
            || event_type == "LT_spiral")
            outputfile << number_of_active_engines;
        else
            outputfile << "-";

        outputfile.width(3); outputfile << " | ";
        outputfile.precision(5);
        if (event_type == "SFthrust"
            || event_type == "FBLTthrust"
            || event_type == "SSFthrust"
            || event_type == "PSBIthrust"
            || event_type == "PSFBthrust"
            || event_type == "begin_spiral"
            || event_type == "end_spiral"
            || event_type == "LT_spiral")
        {
            outputfile.width(14); outputfile << active_power _GETVALUE;
            outputfile.width(3); outputfile << " | ";
        }
        else
        {
            outputfile.width(14); outputfile << 0.0;
            outputfile.width(3); outputfile << " | ";
        }

        if (ThrottleLevel != "none" &&
                (event_type == "SFthrust"
                || event_type == "FBLTthrust"
                || event_type == "SSFthrust"
                || event_type == "PSBIthrust"
                || event_type == "PSFBthrust"
                || event_type == "begin_spiral"
                || event_type == "end_spiral"
                || event_type == "LT_spiral"))
        {
            outputfile.width(14); outputfile << ThrottleLevel;
            outputfile.width(3); outputfile << " | ";
        }
        else
        {
            outputfile.width(14); outputfile << "-";
            outputfile.width(3); outputfile << " | ";
        }

        outputfile << std::endl;
    }//end write_output_line

    void writey_thing::write_ephemeris_line(std::ofstream& outputfile,
        const math::Matrix<doubleType>& state)
    {
        //Step 1: state
        this->write_ephemeris_state(outputfile, state);

        //Step 2: newline
        outputfile << std::endl;
    }


    void writey_thing::write_ephemeris_line(std::ofstream& outputfile,
        const math::Matrix<doubleType>& state,
        const math::Matrix<doubleType>& ControlVector,
        const doubleType& ThrustMagnitude,
        const doubleType& MassFlowRate,
        const doubleType& Isp,
        const int& NumberOfActiveThrusters,
        const doubleType& ActivePower,
        const std::string& ThrottleLevel)
    {
        //Step 1: state
        this->write_ephemeris_state(outputfile, state);

        //Step 2: control vector
        if (this->myOptions->append_control_to_ephemeris_output)
        {
            for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
                outputfile << ", " << ControlVector(controlIndex)_GETVALUE;
        }

        //Step 3: thrust magnitude
        if (this->myOptions->append_thrust_to_ephemeris_output)
            outputfile << ", " << ThrustMagnitude _GETVALUE;

        //Step 4: mass flow rate
        if (this->myOptions->append_mdot_to_ephemeris_output)
            outputfile << ", " << MassFlowRate _GETVALUE;

        //Step 5: Isp
        if (this->myOptions->append_Isp_to_ephemeris_output)
            outputfile << ", " << Isp _GETVALUE;

        //Step 6: number of active thrusters
        if (this->myOptions->append_number_of_active_engines_to_ephemeris_output)
            outputfile << ", " << NumberOfActiveThrusters;

        //Step 7: active power
        if (this->myOptions->append_active_power_to_ephemeris_output)
            outputfile << ", " << ActivePower _GETVALUE;

        //Step 8: throttle level
        if (this->myOptions->append_throttle_level_to_ephemeris_output)
            outputfile << ", " << ThrottleLevel;

        //Step 9: newline
        outputfile << std::endl;
    }

    void writey_thing::write_ephemeris_state(std::ofstream& outputfile,
        const math::Matrix<doubleType>& state)
    {
        SpiceChar epochstring[32];
        timout_c(state(7) _GETVALUE - (51544.5 * 86400.0), "YYYY Mon DD ::TDB HR:MN:SC.######", 32, epochstring);
        outputfile << epochstring;

        //here we assume that the state vector is in sun-centered, Earth equatorial at J2000 coordinates
        outputfile.precision(14);
        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
            outputfile << ", " << state(stateIndex)_GETVALUE;

        if (this->myOptions->append_mass_to_ephemeris_output)
            outputfile << ", " << state(6);
    }
}//end namespace EMTG