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