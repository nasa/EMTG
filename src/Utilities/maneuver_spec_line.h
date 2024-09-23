//Jacob Englander 6/22/2018
//a box to put maneuver spec items in
#pragma once

#include "maneuver_spec_item.h"

#include <vector>
#include <string>
#include <fstream>

namespace EMTG
{
    class maneuver_spec_line
    {
    public:
        //constructors
        maneuver_spec_line() {};
        maneuver_spec_line(const std::string& eventID) : eventID(eventID) {};

        void append_maneuver_spec_item(const std::string& frame,
            const doubleType& epoch,
            const math::Matrix<doubleType>& ControlVector,
            const doubleType& StartMass,
            const doubleType& FinalMass,
            const doubleType& ThrustMagnitude,
            const doubleType& MassFlowRate,
            const doubleType& ManeuverDuration,
            const doubleType& EnforcedDutyCycle);

        void set_eventID(const std::string& eventID) { this->eventID = eventID; }

        void write(std::ofstream& outputfile);

    protected:
        std::string eventID;
        std::vector<maneuver_spec_item> ManeuverItems;
    };//end class maneuver_spec_line
}//close namespace EMTG