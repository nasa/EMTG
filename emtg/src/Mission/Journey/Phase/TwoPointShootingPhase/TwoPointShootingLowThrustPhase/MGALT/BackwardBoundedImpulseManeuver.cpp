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

//force model for EMTGv9
//Jacob Englander 6/25/2017

#include "BackwardBoundedImpulseManeuver.h"

namespace EMTG
{
    namespace Phases
    {
        BackwardBoundedImpulseManeuver::BackwardBoundedImpulseManeuver(const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stageIndex,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            BoundedImpulseManeuver(journeyIndex,
                phaseIndex,
                stageIndex,
                myUniverse,
                mySpacecraft,
                myOptions)
        {
            this->dMassBeforeManeuver_dThrottleComponents.resize(3, 1, 0.0);
        }

        void BackwardBoundedImpulseManeuver::process_maneuver(const bool& needMTM)
        {
            //Step 0: dereference our pointers to local references
            doubleType& max_thrust = *this->max_thrust_pointer;
            doubleType& max_mass_flow_rate = *this->max_mass_flow_rate_pointer;
            doubleType& power = *this->power_pointer;
            doubleType& active_power = *this->active_power_pointer;
            math::Matrix<doubleType>& DeltaV = *this->DeltaV_pointer;
            math::Matrix<doubleType>& spacecraft_state_minus = *this->spacecraft_state_minus_pointer;
            math::Matrix<doubleType>& spacecraft_state_plus = *this->spacecraft_state_plus_pointer;
            math::Matrix<double>& MTM = *this->MTM_pointer;
            doubleType& ThrustStepLength = *this->ThrustStepLength_pointer;
            double& dThrustStepLength_dPropagationVariable = *this->dThrustStepLength_dPropagationVariable_pointer;
            double& dImpulseEpoch_dPropagationVariable = *this->dImpulseEpoch_dPropagationVariable_pointer;
            math::Matrix<doubleType>& Control = *this->Control_pointer;
            doubleType& Throttle = *this->Throttle_pointer;


            //Step 1: extract spacecraft position and epoch
            for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
            {
                this->R_sc_CB(stateIndex) = spacecraft_state_plus(stateIndex);
            }
            this->state_at_natural_perturbation = spacecraft_state_plus_pointer;
            this->ManeuverEpoch = spacecraft_state_plus(7);

            //Step 2: compute thruster performance
            this->compute_thruster_performance(needMTM);

            //Step 3: compute the propellant used
            this->electric_propellant_used = Throttle * max_mass_flow_rate * ThrustStepLength;
            this->ACS_fuel_used = (this->myOptions->trackACS ? this->myOptions->ACS_kg_per_day * ThrustStepLength / 86400.0 : 0.0);

            //Step 4: compute the effective mass, which is the mass BEFORE the maneuver
            spacecraft_state_minus(6) = spacecraft_state_plus(6) + this->electric_propellant_used + this->ACS_fuel_used;
            spacecraft_state_minus(6) = spacecraft_state_minus(6) > 1.0e-3 ? spacecraft_state_minus(6) : 1.0e-3;

            //Step 5: compute the maximum delta-v for this maneuver
            this->dvmax = max_thrust / spacecraft_state_minus(6) * ThrustStepLength;

            //Step 6: compute the actual delta-v
            for (size_t dVindex = 0; dVindex < 3; ++dVindex)
                DeltaV(dVindex) = Control(dVindex) * this->dvmax;

            //Step 7: modify the state vector - do not let the mass drop below 1 gram
            for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
            {
                spacecraft_state_minus(stateIndex) = spacecraft_state_plus(stateIndex);
                spacecraft_state_minus(stateIndex + 3) = spacecraft_state_plus(stateIndex + 3) - DeltaV(stateIndex) - this->NaturalDeltaV(stateIndex);
            }
            spacecraft_state_minus(7) = spacecraft_state_plus(7);

            //Step 8: tank states
            spacecraft_state_minus(8) = spacecraft_state_plus(8) - this->ACS_fuel_used;
            spacecraft_state_minus(9) = spacecraft_state_plus(9) - this->electric_propellant_used;

            //Step 9: build MTM if appropriate
            if (needMTM)
            {
                double P0 = this->myOptions->power_at_1_AU;
                double r = this->R_sc_Sun.norm() _GETVALUE; //this value has been calculated by BoundedImpulseManeuver::compute_thruster_performance()
                double m = spacecraft_state_minus(6) _GETVALUE;
                double dTdt_previous = this->dTdP * this->dPdt;
                double dmdotdt_previous = this->dmdotdP * this->dPdt;
                double dmdt_previous = (Throttle * ThrustStepLength * dmdotdt_previous)_GETVALUE;
                double dT_dTOF = 0;// dTdP * (dPdt * dImpulseEpoch_dPropagationVariable);
                double dmdot_dTOF = 0;// dmdotdP * (dPdt * dImpulseEpoch_dPropagationVariable);
                double dm_dTOF = (Throttle * (max_mass_flow_rate * dThrustStepLength_dPropagationVariable + ThrustStepLength * dmdot_dTOF)) _GETVALUE;
                double dPdP0 = power _GETVALUE / P0;
                double dTdP0 = dTdP * dPdP0;
                double dmdP0 = (Throttle * ThrustStepLength * dmdotdP * dPdP0)_GETVALUE;
                double ddeltamACS_dTOF = (this->myOptions->trackACS ? this->myOptions->ACS_kg_per_day / 86400.0 : 0.0) * dThrustStepLength_dPropagationVariable;

                //Step 9.1 derivatives of delta-V
                for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
                {
                    double c = Control(controlIndex) _GETVALUE;
                    //loop over control components, producing block "A" in the matrix
                    for (size_t Xidx = 0; Xidx < 3; ++Xidx)
                    {
                        MTM(3 + controlIndex, Xidx) = (-c / m * ThrustStepLength * this->dPdr_sun / this->myOptions->AU * this->R_sc_Sun(Xidx) / r
                            * (dTdP - Throttle * max_thrust * ThrustStepLength / m * dmdotdP)) _GETVALUE
                            - this->dNaturalDeltaV_dState(controlIndex, Xidx)_GETVALUE;
                    }

                    //now get derivative with respect to mass
                    MTM(3 + controlIndex, 6) = c * (this->dvmax / m) _GETVALUE
                        - this->dNaturalDeltaV_dState(controlIndex, 6)_GETVALUE;

                    //derivative with respect to previous times
                    MTM(3 + controlIndex, 7) = (-c / (m * m) * ThrustStepLength * (m * dTdt_previous - max_thrust * dmdt_previous) )_GETVALUE
                        - this->dNaturalDeltaV_dPreviousTime(controlIndex)_GETVALUE;

                    //derivative with respect to phase flight time
                    MTM(3 + controlIndex, 10) = (-c * max_thrust / m * dThrustStepLength_dPropagationVariable
                        - c * ThrustStepLength * (1 / (m*m) * (m * dT_dTOF - max_thrust * (dm_dTOF + ddeltamACS_dTOF))))_GETVALUE
                        - this->dNaturalDeltaV_dPropagationTime(controlIndex)_GETVALUE;

                    //with respect to chosen P0
                    //MTM(3 + controlIndex, 9) = (-c * ThrustStepLength / (m * m) * (m * dTdP * dPdP0 - max_thrust * dmdP0))_GETVALUE;

                    //with respect to chosen u_command
                    //MTM(3 + controlIndex, 10) = (-c * ThrustStepLength * (this->dT_du_command / m - Throttle * ThrustStepLength * max_thrust * this->dmdot_du_command / (m * m)))_GETVALUE;
                }

                //Step 9.2 derivatives of mass
                //with respect to position
                for (int Xidx = 0; Xidx < 3; ++Xidx)
                    MTM(6, Xidx) = (Throttle * ThrustStepLength * dmdotdP * dPdr_sun / this->myOptions->AU
                    * this->R_sc_Sun(Xidx) / r) _GETVALUE;

                //with respect to previous time
                MTM(6, 7) = dmdt_previous;

                //with respect to phase flight time
                MTM(6, 10) = dm_dTOF
                                + ddeltamACS_dTOF;

                //with respect to chosen P0
                //MTM(6, 9) = dmdP0;

                //with respect to chosen u_command
                //MTM(6, 10) = (Throttle * ThrustStepLength * this->dmdot_du_command)_GETVALUE;

                //Step 9.3 derivatives of chemical (ACS) fuel tank
                MTM(8, 10) = -ddeltamACS_dTOF;

                //Step 9.4 derivatives of electric tank
                //with respect to position
                for (int Xidx = 0; Xidx < 3; ++Xidx)
                    MTM(9, Xidx) = -(Throttle * ThrustStepLength * dmdotdP * dPdr_sun / this->myOptions->AU
                        * this->R_sc_Sun(Xidx) / r) _GETVALUE;

                //with respect to boundary time
                MTM(9, 7) = -dmdt_previous;

                //with respect to phase flight time
                MTM(9, 10) = -dm_dTOF;

                //with respect to chosen P0
                //MTM(9, 9) = dmdP0;

                //with respect to chosen u_command
                //MTM(9, 10) = (Throttle * ThrustStepLength * dmdot_du_command)_GETVALUE;


                //Step 10: derivatives with respect to control
                for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
                {
                    this->dMassBeforeManeuver_dThrottleComponents(controlIndex) = ((-max_mass_flow_rate * ThrustStepLength)
                        * (Control(controlIndex) / Throttle))_GETVALUE;

                    this->dChemicalFuel_dThrottleComponents(controlIndex) = 0.0;

                    this->dElectricPropellant_dThrottleComponents(controlIndex) = -this->dMassBeforeManeuver_dThrottleComponents(controlIndex);
                }
            }//end derivatives
                
        }//end process_maneuver()
    }//close namespace Astrodynamics
}//close namespace EMTG