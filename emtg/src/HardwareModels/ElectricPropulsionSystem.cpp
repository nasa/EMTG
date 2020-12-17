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

//electric propulsion system class
//Jacob Englander 11/2/2016

#include "ElectricPropulsionSystem.h"

#include <iostream>

namespace EMTG
{
    namespace HardwareModels
    {        
        //constructor
        ElectricPropulsionSystem::ElectricPropulsionSystem(const PropulsionSystemOptions& propulsionsystemoptions) :
            PropulsionSystem(propulsionsystemoptions),
            ThrottleLevel(0),
            ThrottleLevelString("none")
        {
            if (!(this->myPropulsionSystemOptions.getThrottleTableFile() == "none") && !(this->myPropulsionSystemOptions.getThrottleTableFile() == "MyThrottleTable.throttletable"))
                this->myThrottleTable = ThrottleTable(this->myPropulsionSystemOptions.getThrottleTableFile(), this->myPropulsionSystemOptions.getSharpness());

            this->ActivePower = 0.0;
        }

        //method to compute system mass
        void ElectricPropulsionSystem::computeSystemMass()
        {
            this->SystemMass = this->myPropulsionSystemOptions.getMassPerString() * this->myPropulsionSystemOptions.getNumberOfStrings();
        }

        //evaluate methods
        void ElectricPropulsionSystem::compute_number_of_active_thrusters(const doubleType& InputPower, const ThrottleLogic& myThrottleLogic)
        {
            switch (myThrottleLogic)
            {
                case ThrottleLogic::MaxThrusters:
                {
                    //the idea here it to stack several heaviside step functions, each one with a transition power at n*this->myPropulsionSystemOptions.getPmin()
                    this->smoothed_number_of_active_thrusters = 0.0;
                    this->dNdP = 0.0;
                    for (size_t n = 1; n <= this->myPropulsionSystemOptions.getNumberOfStrings(); ++n)
                    {
                        double Pcrit = n*this->myPropulsionSystemOptions.getPmin();
                        doubleType Pswitch = math::absclip(-this->myPropulsionSystemOptions.getSharpness() * (InputPower - Pcrit), 600.0);
                        doubleType expfun = exp(Pswitch);
                        doubleType ThisThruster = 1.0 / (1.0 + expfun);
                        this->smoothed_number_of_active_thrusters += ThisThruster;

                        double DEN = (expfun * (1.0 / expfun + 1.0) * (1.0 / expfun + 1.0))_GETVALUE;
                        if (DEN > 0.0)
                            this->dNdP += this->myPropulsionSystemOptions.getSharpness() / DEN;

                    }
                    if (this->smoothed_number_of_active_thrusters > 1 || InputPower > this->myPropulsionSystemOptions.getPmax())
                    {
                        if (InputPower / this->smoothed_number_of_active_thrusters > this->myPropulsionSystemOptions.getPmax())
                        {
                            if (this->myPropulsionSystemOptions.getThrusterMode() > 3)
                            {
                                PowerPerThruster = InputPower / this->smoothed_number_of_active_thrusters;
                                dPowerPerThruster_dP = ((this->smoothed_number_of_active_thrusters - InputPower * this->dNdP) / (this->smoothed_number_of_active_thrusters*this->smoothed_number_of_active_thrusters)) _GETVALUE;
                            }
                            else
                            {
                                PowerPerThruster = this->myPropulsionSystemOptions.getPmax();
                                dPowerPerThruster_dP = 0.0;
                            }
                        }
                        else
                        {
                            PowerPerThruster = InputPower / this->smoothed_number_of_active_thrusters;
                            dPowerPerThruster_dP = ((this->smoothed_number_of_active_thrusters - InputPower * this->dNdP) / (this->smoothed_number_of_active_thrusters*this->smoothed_number_of_active_thrusters))_GETVALUE;
                        }
                    }
                    else if (this->smoothed_number_of_active_thrusters < 1 && InputPower < this->myPropulsionSystemOptions.getPmin())
                    {
                        if (this->myPropulsionSystemOptions.getThrusterMode() > 3)
                        {
                            PowerPerThruster = InputPower;
                            dPowerPerThruster_dP = 1.0;
                        }
                        else
                        {
                            PowerPerThruster = this->myPropulsionSystemOptions.getPmin();
                            dPowerPerThruster_dP = 0.0;
                        }
                    }
                    else
                    {
                        PowerPerThruster = InputPower;
                        dPowerPerThruster_dP = 1.0;
                    }
                    break;
                }//end maximum number of thrusters
                case ThrottleLogic::MinThrusters:
                {
                    //the idea here it to stack several heaviside step functions, each one with a transition power at (n-1)*this->myPropulsionSystemOptions.getPmax() except the bottom one whose transition InputPower is 
                    this->smoothed_number_of_active_thrusters = 0.0;
                    this->dNdP = 0.0;
                    for (size_t n = 1; n <= this->myPropulsionSystemOptions.getNumberOfStrings(); ++n)
                    {
                        double Pcrit = n == 1 ? this->myPropulsionSystemOptions.getPmin() : (n - 1)*this->myPropulsionSystemOptions.getPmax();
                        doubleType Pswitch = math::absclip(-this->myPropulsionSystemOptions.getSharpness() * (InputPower - Pcrit), 600.0);
                        doubleType expfun = exp(Pswitch);
                        doubleType ThisThruster = 1.0 / (1.0 + expfun);
                        this->smoothed_number_of_active_thrusters += ThisThruster;

                        double DEN = (expfun * (1.0 / expfun + 1.0) * (1.0 / expfun + 1.0))_GETVALUE;
                        if (DEN > 0.0)
                            this->dNdP += this->myPropulsionSystemOptions.getSharpness() / DEN;
                    }
                    if (this->smoothed_number_of_active_thrusters > 1 || InputPower > this->myPropulsionSystemOptions.getPmax())
                    {
                        if (InputPower / this->smoothed_number_of_active_thrusters > this->myPropulsionSystemOptions.getPmax())
                        {
                            if (this->myPropulsionSystemOptions.getThrusterMode() > 3)
                            {
                                PowerPerThruster = InputPower / this->smoothed_number_of_active_thrusters;
                                dPowerPerThruster_dP = ((this->smoothed_number_of_active_thrusters - InputPower * this->dNdP) / (this->smoothed_number_of_active_thrusters*this->smoothed_number_of_active_thrusters)) _GETVALUE;
                            }
                            else
                            {
                                PowerPerThruster = this->myPropulsionSystemOptions.getPmax();
                                dPowerPerThruster_dP = 0.0;
                            }
                        }
                        else
                        {
                            PowerPerThruster = InputPower / this->smoothed_number_of_active_thrusters;
                            dPowerPerThruster_dP = ((this->smoothed_number_of_active_thrusters - InputPower * this->dNdP) / (this->smoothed_number_of_active_thrusters*this->smoothed_number_of_active_thrusters))_GETVALUE;
                        }
                    }
                    else if (this->smoothed_number_of_active_thrusters < 1 && InputPower < this->myPropulsionSystemOptions.getPmin())
                    {
                        if (this->myPropulsionSystemOptions.getThrusterMode() > 3)
                        {
                            PowerPerThruster = InputPower;
                            dPowerPerThruster_dP = 1.0;
                        }
                        else
                        {
                            PowerPerThruster = this->myPropulsionSystemOptions.getPmin();
                            dPowerPerThruster_dP = 0.0;
                        }
                    }
                    else
                    {
                        PowerPerThruster = InputPower;
                        dPowerPerThruster_dP = 1.0;
                    }

                }//end minimum number of thrusters
            }
        }

        void ElectricPropulsionSystem::computeThrusterPerformance(const doubleType& InputPower, const double& DutyCycle, const ThrottleLogic& myThrottleLogic)
        {
            this->compute_number_of_active_thrusters(InputPower, myThrottleLogic);
            double g0 = this->myPropulsionSystemOptions.getg0();

            switch (this->myPropulsionSystemOptions.getThrusterMode())
            {
                case SpacecraftThrusterMode::ConstantThrustIsp:
                {
                    this->Thrust = this->myPropulsionSystemOptions.getConstantThrust();
                    this->Isp = this->myPropulsionSystemOptions.getConstantIsp();
                    this->MassFlowRate = this->Thrust / this->Isp / g0;
                    this->dTdP = 0.0;
                    this->dMassFlowRatedP = 0.0;
                    this->dIspdP = 0.0;

                    break;
                }
                case SpacecraftThrusterMode::FixedEfficiencyCSI:
                {
                    this->Isp = this->myPropulsionSystemOptions.getConstantIsp();

                    //compute thrust
                    doubleType temp_dTdP = (2000.0 * (this->myPropulsionSystemOptions.getFixedEfficiency() / (this->Isp * g0)));

                    this->dTdP = temp_dTdP _GETVALUE;
                    this->Thrust = dTdP * this->PowerPerThruster;
                    this->ActivePower = this->PowerPerThruster;

                    //compute mass flow rate
                    doubleType temp_dmdotdP = temp_dTdP / this->Isp / g0;
                    this->dMassFlowRatedP = temp_dmdotdP _GETVALUE;
                    this->MassFlowRate = temp_dmdotdP * this->PowerPerThruster;
                    this->dIspdP = 0.0;
                    this->dIspdu_command = 0.0;

                    break;
                }
                case SpacecraftThrusterMode::Poly1D:
                {
                    //now either power the engines, or, if insufficient power, switch off the engines
                    doubleType T, F, P = PowerPerThruster;
                    doubleType Poffset = PowerPerThruster * 1.0e-10;
                    //next compute Thrust in mN
                    doubleType P2 = P * P;
                    doubleType P3 = P2 * P;
                    doubleType P4 = P3 * P;
                    doubleType P5 = P4 * P;
                    doubleType P6 = P5 * P;

                    double at = this->myPropulsionSystemOptions.getThrustCoefficient(0);
                    double bt = this->myPropulsionSystemOptions.getThrustCoefficient(1);
                    double ct = this->myPropulsionSystemOptions.getThrustCoefficient(2);
                    double dt = this->myPropulsionSystemOptions.getThrustCoefficient(3);
                    double et = this->myPropulsionSystemOptions.getThrustCoefficient(4);
                    double ft = this->myPropulsionSystemOptions.getThrustCoefficient(5);
                    double gt = this->myPropulsionSystemOptions.getThrustCoefficient(6);

                    double af = this->myPropulsionSystemOptions.getMassFlowCoefficient(0);
                    double bf = this->myPropulsionSystemOptions.getMassFlowCoefficient(1);
                    double cf = this->myPropulsionSystemOptions.getMassFlowCoefficient(2);
                    double df = this->myPropulsionSystemOptions.getMassFlowCoefficient(3);
                    double ef = this->myPropulsionSystemOptions.getMassFlowCoefficient(4);
                    double ff = this->myPropulsionSystemOptions.getMassFlowCoefficient(5);
                    double gf = this->myPropulsionSystemOptions.getMassFlowCoefficient(6);
                    
                    T = gt*P6 + ft*P5 + et*P4 + dt*P3 + ct*P2 + bt*P + at;

                    //next compute mass flow rate in mg/s
                    F = gf*P6 + ff*P5 + ef*P4 + df*P3 + cf*P2 + bf*P + af;

                    //return Thrust in N (convert from mN) and mass flow rate in kg/s (convert from mg/s)
                    this->Thrust = 1.0e-3 * T * this->smoothed_number_of_active_thrusters;
                    this->MassFlowRate = 1.0e-6 * F * this->smoothed_number_of_active_thrusters;
                    this->Isp = this->Thrust / this->MassFlowRate / this->myPropulsionSystemOptions.getg0();
                    this->ActivePower = this->smoothed_number_of_active_thrusters * P;

                    if (PowerPerThruster <= this->myPropulsionSystemOptions.getPmax())
                    {
                        this->dTdP = ((6 * gt*P5 + 5 * ft*P4 + 4 * et*P3 + 3 * dt*P2 + 2 * ct*P + bt) * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                            + T * dNdP) _GETVALUE * 1.0e-3;
                        this->dMassFlowRatedP = ((6 * gf*P5 + 5 * ff*P4 + 4 * ef*P3 + 3 * df*P2 + 2 * cf*P + bf) * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                            + F * dNdP) _GETVALUE * 1.0e-6;

                        this->dIspdP = 1.0 / this->myPropulsionSystemOptions.getg0() * ((this->dTdP * this->MassFlowRate - this->Thrust * this->dMassFlowRatedP) / (this->MassFlowRate * this->MassFlowRate)) _GETVALUE;
                    }
                    else
                    {
                        this->dTdP = 0.0;
                        this->dMassFlowRatedP = 0.0;
                        this->dIspdP = 0.0;
                    }

                    break;
                }
                case SpacecraftThrusterMode::SteppedHThrust1D:
                {
                    this->myThrottleTable.CalculateThrusterPerformance1D(this->PowerPerThruster, HardwareModels::HighThrust);
                    
                    this->Thrust = 1e-3 * this->myThrottleTable.getHeavisideThrust() * this->smoothed_number_of_active_thrusters;
                    this->MassFlowRate = 1e-6 * this->myThrottleTable.getHeavisideMdot() * this->smoothed_number_of_active_thrusters;
                    this->ThrottleLevel = this->myThrottleTable.getThrottleLevel();
                    this->ThrottleLevelString = this->myThrottleTable.getThrottleLevelString();
                    this->ActivePower = this->myThrottleTable.getHeavisideActivePower() * this->smoothed_number_of_active_thrusters;
                    this->Isp = this->myThrottleTable.getHeavisideIsp();


                    this->dTdP = 1e-3 * (this->myThrottleTable.getHeavisideThrustDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->myThrottleTable.getHeavisideThrust() * this->dNdP) _GETVALUE;
                    this->dMassFlowRatedP = 1e-6 * (this->myThrottleTable.getHeavisideMdotDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->myThrottleTable.getHeavisideMdot() * this->dNdP) _GETVALUE;

                    this->dIspdP = 1.0 / this->myPropulsionSystemOptions.getg0() * ((this->dTdP * this->MassFlowRate - this->Thrust * this->dMassFlowRatedP) / (this->MassFlowRate * this->MassFlowRate)) _GETVALUE;

                    break;
                }
                case SpacecraftThrusterMode::SteppedLMdot1D:
                {
                    this->myThrottleTable.CalculateThrusterPerformance1D(this->PowerPerThruster, HardwareModels::LowMdot);

                    this->Thrust = 1e-3 * this->myThrottleTable.getHeavisideThrust() * this->smoothed_number_of_active_thrusters;
                    this->MassFlowRate = 1e-6 * this->myThrottleTable.getHeavisideMdot() * this->smoothed_number_of_active_thrusters;
                    this->ThrottleLevel = this->myThrottleTable.getThrottleLevel();
                    this->ThrottleLevelString = this->myThrottleTable.getThrottleLevelString();
                    this->ActivePower = this->myThrottleTable.getHeavisideActivePower() * this->smoothed_number_of_active_thrusters;
                    this->Isp = this->myThrottleTable.getHeavisideIsp();


                    this->dTdP = 1e-3 * (this->myThrottleTable.getHeavisideThrustDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->myThrottleTable.getHeavisideThrust() * this->dNdP) _GETVALUE;
                    this->dMassFlowRatedP = 1e-6 * (this->myThrottleTable.getHeavisideMdotDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->myThrottleTable.getHeavisideMdot() * this->dNdP) _GETVALUE;
                    
                    this->dIspdP = 1e-6 * (this->myThrottleTable.getHeavisideIspDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->Isp * this->dNdP) _GETVALUE;

                    break;
                }
                case SpacecraftThrusterMode::SteppedHEfficiency1D:
                {
                    this->myThrottleTable.CalculateThrusterPerformance1D(this->PowerPerThruster, HardwareModels::HighEfficiency);

                    this->Thrust = 1e-3 * this->myThrottleTable.getHeavisideThrust() * this->smoothed_number_of_active_thrusters;
                    this->MassFlowRate = 1e-6 * this->myThrottleTable.getHeavisideMdot() * this->smoothed_number_of_active_thrusters;
                    this->ThrottleLevel = this->myThrottleTable.getThrottleLevel();
                    this->ThrottleLevelString = this->myThrottleTable.getThrottleLevelString();
                    this->ActivePower = this->myThrottleTable.getHeavisideActivePower() * this->smoothed_number_of_active_thrusters;
                    this->Isp = this->myThrottleTable.getHeavisideIsp();


                    this->dTdP = 1e-3 * (this->myThrottleTable.getHeavisideThrustDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->myThrottleTable.getHeavisideThrust() * this->dNdP) _GETVALUE;
                    this->dMassFlowRatedP = 1e-6 * (this->myThrottleTable.getHeavisideMdotDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->myThrottleTable.getHeavisideMdot() * this->dNdP) _GETVALUE;

                    this->dIspdP = 1e-6 * (this->myThrottleTable.getHeavisideIspDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->Isp * this->dNdP) _GETVALUE;
                    break;
                }
                case SpacecraftThrusterMode::SteppedHisp1D:
                {
                    this->myThrottleTable.CalculateThrusterPerformance1D(this->PowerPerThruster, HardwareModels::HighIsp);

                    this->Thrust = 1e-3 * this->myThrottleTable.getHeavisideThrust() * this->smoothed_number_of_active_thrusters;
                    this->MassFlowRate = 1e-6 * this->myThrottleTable.getHeavisideMdot() * this->smoothed_number_of_active_thrusters;
                    this->ThrottleLevel = this->myThrottleTable.getThrottleLevel();
                    this->ThrottleLevelString = this->myThrottleTable.getThrottleLevelString();
                    this->ActivePower = this->myThrottleTable.getHeavisideActivePower() * this->smoothed_number_of_active_thrusters;
                    this->Isp = this->myThrottleTable.getHeavisideIsp();


                    this->dTdP = 1e-3 * (this->myThrottleTable.getHeavisideThrustDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->myThrottleTable.getHeavisideThrust() * this->dNdP) _GETVALUE;
                    this->dMassFlowRatedP = 1e-6 * (this->myThrottleTable.getHeavisideMdotDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->myThrottleTable.getHeavisideMdot() * this->dNdP) _GETVALUE;

                    this->dIspdP = 1e-6 * (this->myThrottleTable.getHeavisideIspDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->Isp * this->dNdP) _GETVALUE;
                    break;
                }
                case SpacecraftThrusterMode::SteppedFullSet1D:
                {
                    this->myThrottleTable.CalculateThrusterPerformance1D(this->PowerPerThruster, HardwareModels::FullSet);

                    this->Thrust = 1e-3 * this->myThrottleTable.getHeavisideThrust() * this->smoothed_number_of_active_thrusters;
                    this->MassFlowRate = 1e-6 * this->myThrottleTable.getHeavisideMdot() * this->smoothed_number_of_active_thrusters;
                    this->ThrottleLevel = this->myThrottleTable.getThrottleLevel();
                    this->ThrottleLevelString = this->myThrottleTable.getThrottleLevelString();
                    this->ActivePower = this->myThrottleTable.getHeavisideActivePower() * this->smoothed_number_of_active_thrusters;
                    this->Isp = this->myThrottleTable.getHeavisideIsp();


                    this->dTdP = 1e-3 * (this->myThrottleTable.getHeavisideThrustDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->myThrottleTable.getHeavisideThrust() * this->dNdP) _GETVALUE;
                    this->dMassFlowRatedP = 1e-6 * (this->myThrottleTable.getHeavisideMdotDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->myThrottleTable.getHeavisideMdot() * this->dNdP) _GETVALUE;

                    this->dIspdP = 1e-6 * (this->myThrottleTable.getHeavisideIspDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP
                        + this->Isp * this->dNdP) _GETVALUE;
                    break;
                }
                default:
                {
                    throw std::invalid_argument("ElectricPropulsionSystem::computeThrusterPerformance() Invalid SpacecraftThrusterMode " + SpacecraftThrusterModeStrings[this->myPropulsionSystemOptions.getThrusterMode()] + " for 1D thruster model. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
            }//end switch

            //this is the 1D version, so u_command has no derivatives
            this->dMassFlowRatedu_command = 0.0;
            this->dTdu_command = 0.0;
            this->dIspdu_command = 0.0;

            //scale by ThrustScaleFactor and duty cycle
            this->Thrust *= this->myPropulsionSystemOptions.getThrustScaleFactor() * DutyCycle;
            this->dTdP *= this->myPropulsionSystemOptions.getThrustScaleFactor() * DutyCycle;
            this->MassFlowRate *= DutyCycle;
            this->dMassFlowRatedP *= DutyCycle;
            //dIspdP does not need to be scaled by DutyCycle

        }//end 1D computeThrusterPerformance

        void ElectricPropulsionSystem::computeThrusterPerformance(const doubleType& InputPower, const double& DutyCycle, const ThrottleLogic& myThrottleLogic, const doubleType& u_command)
        {
            this->compute_number_of_active_thrusters(InputPower, myThrottleLogic);

            switch (this->myPropulsionSystemOptions.getThrusterMode())
            {
                case SpacecraftThrusterMode::FixedEfficiencyVSI:
                {
                    this->Isp = u_command * (this->myPropulsionSystemOptions.getConstantIsp() - this->myPropulsionSystemOptions.getMinimumOrMonopropIsp()) + this->myPropulsionSystemOptions.getMinimumOrMonopropIsp();

                    //compute thrust
                    doubleType temp_dTdP = (2000.0 * (this->myPropulsionSystemOptions.getFixedEfficiency() / (this->Isp * this->g0)));

                    this->dTdP = temp_dTdP _GETVALUE;
                    this->Thrust = dTdP * this->PowerPerThruster;
                    this->ActivePower = this->PowerPerThruster;

                    //compute mass flow rate
                    doubleType temp_dmdotdP = temp_dTdP / this->Isp / this->g0;
                    this->dMassFlowRatedP = temp_dmdotdP _GETVALUE;
                    this->MassFlowRate = temp_dmdotdP * this->PowerPerThruster;

                    //derivative with respect to u_command
                    this->dTdu_command = (-1.0 / (this->Isp) * this->Thrust * (this->myPropulsionSystemOptions.getConstantIsp() - this->myPropulsionSystemOptions.getMinimumOrMonopropIsp()))_GETVALUE;
                    this->dMassFlowRatedu_command = (-2.0 / (this->Isp) * this->MassFlowRate * (this->myPropulsionSystemOptions.getConstantIsp() - this->myPropulsionSystemOptions.getMinimumOrMonopropIsp()))_GETVALUE;
                    
                    this->Thrust *= this->myPropulsionSystemOptions.getThrustScaleFactor() * DutyCycle;
                    this->dTdP *= this->myPropulsionSystemOptions.getThrustScaleFactor() * DutyCycle;
                    this->dTdu_command *= this->myPropulsionSystemOptions.getThrustScaleFactor() * DutyCycle;
                    this->MassFlowRate *= DutyCycle;
                    this->dMassFlowRatedP *= DutyCycle;
                    this->dMassFlowRatedu_command *= DutyCycle;
                    this->dIspdP = 0.0;
                    this->dIspdu_command = (this->myPropulsionSystemOptions.getConstantIsp() - this->myPropulsionSystemOptions.getMinimumOrMonopropIsp()) + this->myPropulsionSystemOptions.getMinimumOrMonopropIsp();

                    break;
                }
                case SpacecraftThrusterMode::Stepped2D:
                {
                    myThrottleTable.CalculateThrusterPerformance2D(this->PowerPerThruster, u_command);

                    this->Thrust = 1e-3 * myThrottleTable.getHeavisideThrust() * this->smoothed_number_of_active_thrusters;
                    this->MassFlowRate = 1e-6 * myThrottleTable.getHeavisideMdot() * this->smoothed_number_of_active_thrusters;
                    this->ThrottleLevel = myThrottleTable.getThrottleLevel();
                    this->ThrottleLevelString = this->myThrottleTable.getThrottleLevelString();
                    this->ActivePower = myThrottleTable.getHeavisideActivePower() * this->smoothed_number_of_active_thrusters;
                    this->Isp = myThrottleTable.getHeavisideIsp();
                    this->dTdP = 1e-3 * (myThrottleTable.getHeavisideThrustDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP + myThrottleTable.getHeavisideThrust() * this->dNdP) _GETVALUE;
                    this->dMassFlowRatedP = 1e-6 * (myThrottleTable.getHeavisideMdotDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP + myThrottleTable.getHeavisideMdot() * this->dNdP) _GETVALUE;
                    this->dTdu_command = 1e-6 * (myThrottleTable.getdHthrust_dcommand() * this->smoothed_number_of_active_thrusters) _GETVALUE;
                    this->dMassFlowRatedu_command = 1e-6 * (myThrottleTable.getdHmdot_dcommand() * this->smoothed_number_of_active_thrusters) _GETVALUE;

                    this->Thrust *= this->myPropulsionSystemOptions.getThrustScaleFactor() * DutyCycle;
                    this->dTdP *= this->myPropulsionSystemOptions.getThrustScaleFactor() * DutyCycle;
                    this->dTdu_command *= this->myPropulsionSystemOptions.getThrustScaleFactor() * DutyCycle;
                    this->MassFlowRate *= DutyCycle;
                    this->dMassFlowRatedP *= DutyCycle;
                    this->dMassFlowRatedu_command *= DutyCycle;
                    this->dIspdP = (myThrottleTable.getHeavisideIspDerivative() * this->smoothed_number_of_active_thrusters * this->dPowerPerThruster_dP + myThrottleTable.getHeavisideIsp() * this->dNdP) _GETVALUE;
                    this->dIspdu_command = (myThrottleTable.getdHIsp_dcommand() * this->smoothed_number_of_active_thrusters) _GETVALUE;

                    break;
                }
                default:
                {
                    //drop to 1D
                    this->computeThrusterPerformance(InputPower, DutyCycle, myThrottleLogic);
                }
            }//end switch

            //scale by ThrustScaleFactor and duty cycle

        }//end 2D computeThrusterPerformance
    }//end namespace HardwareModels
}//end namespace EMTG