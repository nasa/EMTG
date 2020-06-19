//Jacob Englander 6/22/2018
//maneuver spec item
//this class is for one item that goes into an EMTG-MIRAGE maneuver spec line
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
