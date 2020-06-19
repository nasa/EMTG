//Jacob Englander 6/22/2018
//maneuver spec item
//this class is for one item that goes into an EMTG-MIRAGE maneuver spec line
//
//contains the following information
//<FRAME>,<EPOCH(ET)>,<THRX>,<THRY>,<THRZ>,<THRMAG[N]>,<SMASS[kg]>,<MDOT[kg/s]>,<DUTY>,<FMASS[kg]>,<DV[km/s]>,<DUR[s]>,<EPOCH(ET seconds)>
//where <THRX>, <THRY>, <THRZ> are unit vector components
//

#include "maneuver_spec_item.h"
#include "EMTG_math.h"
#include "SpiceUsr.h"

namespace EMTG
{
    maneuver_spec_item::maneuver_spec_item()
    {
        this->initialize("EME2000",
            51544.5,
            math::Matrix<doubleType>(3, 1, 0.0),
            1000.0,
            990.0,
            0.1,
            1.0e-6,
            3600.0,
            1.0);
    }//end default constructor


    maneuver_spec_item::maneuver_spec_item(const std::string& frame,
        const doubleType& epoch,
        const math::Matrix<doubleType>& ControlVector,
        const doubleType& StartMass,
        const doubleType& FinalMass,
        const doubleType& ThrustMagnitude,
        const doubleType& MassFlowRate,
        const doubleType& ManeuverDuration,
        const doubleType& EnforcedDutyCycle) //this is the duty cycle enforced by the problem, not the one chosen by the optimizer
    {
        this->initialize(frame,
            epoch,
            ControlVector,
            StartMass,
            FinalMass,
            ThrustMagnitude,
            MassFlowRate,
            ManeuverDuration,
            EnforcedDutyCycle);
    }//end constructor with inputs

    void maneuver_spec_item::initialize(const std::string& frame,
        const doubleType& epoch,
        const math::Matrix<doubleType>& ControlVector,
        const doubleType& StartMass,
        const doubleType& FinalMass,
        const doubleType& ThrustMagnitude,
        const doubleType& MassFlowRate,
        const doubleType& ManeuverDuration,
        const doubleType& EnforcedDutyCycle)
    {
        this->frame = frame;
        this->ManeuverStartEpochETseconds = epoch;
        this->ControlVector = ControlVector;
        this->StartMass = StartMass;
        this->FinalMass = FinalMass;
        this->ThrustMagnitude = ThrustMagnitude;
        this->ManeuverDuration = ManeuverDuration;
        this->MassFlowRate = MassFlowRate;
        this->EnforcedDutyCycle = EnforcedDutyCycle;

        //process time
        this->ManeuverStartEpochETJD = this->ManeuverStartEpochETseconds / 86400.0 + 2400000.5;
        char GregorianDate_temp[68];
        timout_c(unitim_c(this->ManeuverStartEpochETJD _GETVALUE, "JDTDB", "TDB"), "YYYY Mon DD ::TDB HR:MN:SC.######", 68, GregorianDate_temp);
        this->ManeuverStartEpochETGregorian = std::string(GregorianDate_temp);

        //process duty cycle
        doubleType ControlVectorMagnitude = this->ControlVector.norm();
        this->ControlUnitVector = this->ControlVector / (ControlVectorMagnitude + 1.0e-20);
        this->TrueDutyCycle = this->EnforcedDutyCycle * ControlVectorMagnitude;

        //process delta-v
        const double g0 = 9.80665;
        doubleType Isp = this->ThrustMagnitude / g0 / this->MassFlowRate;
        this->Deltav = -log(this->FinalMass / this->StartMass) * this->ThrustMagnitude / this->MassFlowRate / 1000.0;
    }//end initialize()


    void maneuver_spec_item::write(std::ofstream& outputfile)
    {
        outputfile.precision(14);

        outputfile << ", " << this->frame;
        outputfile << ", " << this->ManeuverStartEpochETGregorian;
        outputfile << ", " << this->ControlUnitVector(0);
        outputfile << ", " << this->ControlUnitVector(1);
        outputfile << ", " << this->ControlUnitVector(2);
        outputfile << ", " << this->ThrustMagnitude;
        outputfile << ", " << this->StartMass;
        outputfile << ", " << this->MassFlowRate;
        outputfile << ", " << this->TrueDutyCycle;
        outputfile << ", " << this->FinalMass;
        outputfile << ", " << this->Deltav;
        outputfile << ", " << this->ManeuverDuration;
        outputfile.precision(20);
        outputfile << ", " << this->ManeuverStartEpochETseconds;
    }//end write()
}//end namespace EMTG