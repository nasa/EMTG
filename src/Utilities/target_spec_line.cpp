//Jacob Englander 6/22/2018
//target spec line
//this class is for one EMTG-MIRAGE target spec line
//
//contains the following information
//<SAMPLE>,<CBODY>,<FRAME>,<EPOCH(ET)>,<X[km]>,<Y[km]>,<Z[km]>,<VX[km/s]>,<VY[km/s]>,<VZ[km/s]>,<B.R[km]>,<B.T [km]>,<EPOCH(ET seconds)>
//

#include "target_spec_line.h"

#include "SpiceUsr.h"

namespace EMTG
{
    target_spec_line::target_spec_line(const std::string& eventID,
        const std::string& frame,
        const std::string& CentralBody,
        const doubleType& epoch,
        const math::Matrix<doubleType>& State)
    {
        this->initialize(eventID,
            frame,
            CentralBody,
            epoch,
            State);
    }

    target_spec_line::target_spec_line(const std::string& eventID,
        const std::string& frame,
        const std::string& CentralBody,
        const doubleType& epoch,
        const math::Matrix<doubleType>& State,
        const doubleType& BdotR,
        const doubleType& BdotT)
    {
        this->initialize(eventID,
            frame,
            CentralBody,
            epoch,
            State,
            BdotR,
            BdotT);
    }

    //methods
    void target_spec_line::initialize(const std::string& eventID,
        const std::string& frame,
        const std::string& CentralBody,
        const doubleType& epoch,
        const math::Matrix<doubleType>& State)
    {
        this->eventID = eventID;
        this->frame = frame;
        this->CentralBody = CentralBody;
        this->TargetEpochETseconds = epoch;
        this->State = State;
        this->BdotR = 0.0;
        this->BdotT = 0.0;

        //process time
        this->TargetEpochETJD = this->TargetEpochETseconds / 86400.0 + 2400000.5;
        char GregorianDate_temp[68];
        timout_c(unitim_c(this->TargetEpochETJD _GETVALUE, "JDTDB", "TDB"), "YYYY Mon DD ::TDB HR:MN:SC.######", 68, GregorianDate_temp);
        this->TargetEpochETGregorian = std::string(GregorianDate_temp);
    }

    void target_spec_line::initialize(const std::string& eventID,
        const std::string& frame,
        const std::string& CentralBody,
        const doubleType& epoch,
        const math::Matrix<doubleType>& State,
        const doubleType& BdotR,
        const doubleType& BdotT)
    {
        this->initialize(eventID,
            frame,
            CentralBody,
            epoch,
            State);

        this->BdotR = BdotR;
        this->BdotT = BdotT;
    }

    void target_spec_line::write(std::ofstream& outputfile)
    {
        outputfile.precision(14);

        outputfile << this->eventID;
        outputfile << ", 1";
        outputfile << ", " << this->CentralBody;
        outputfile << ", " << this->frame;
        outputfile << ", " << this->TargetEpochETGregorian;
        outputfile << ", " << this->State(0);
        outputfile << ", " << this->State(1);
        outputfile << ", " << this->State(2);
        outputfile << ", " << this->State(3);
        outputfile << ", " << this->State(4);
        outputfile << ", " << this->State(5);
        outputfile << ", " << this->State(6);
        outputfile << ", " << this->BdotR;
        outputfile << ", " << this->BdotT;
        outputfile.precision(20);
        outputfile << ", " << this->TargetEpochETseconds;
        outputfile << std::endl;
    }
}//end namespace EMTG
