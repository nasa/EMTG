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
