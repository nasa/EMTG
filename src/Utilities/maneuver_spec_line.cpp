//Jacob Englander 6/22/2018
//a box to put maneuver spec items in

#include "maneuver_spec_line.h"

namespace EMTG
{
    void maneuver_spec_line::append_maneuver_spec_item(const std::string& frame,
        const doubleType& epoch,
        const math::Matrix<doubleType>& ControlVector,
        const doubleType& StartMass,
        const doubleType& FinalMass,
        const doubleType& ThrustMagnitude,
        const doubleType& MassFlowRate,
        const doubleType& ManeuverDuration,
        const doubleType& EnforcedDutyCycle)
    {
        this->ManeuverItems.push_back(maneuver_spec_item(frame,
            epoch,
            ControlVector,
            StartMass,
            FinalMass,
            ThrustMagnitude,
            MassFlowRate,
            ManeuverDuration,
            EnforcedDutyCycle));
    }//end append_maneuver_spec_item()
    
    void maneuver_spec_line::write(std::ofstream& outputfile)
    {
        outputfile << this->eventID;
        outputfile << ", " << this->ManeuverItems.size();

        for (maneuver_spec_item thisItem : this->ManeuverItems)
            thisItem.write(outputfile);

        outputfile << std::endl;
    }//end write()
}//close namespace EMTG