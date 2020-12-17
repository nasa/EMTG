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

//MGALT maneuver base class for EMTGv9
//Jacob Englander 6/25/2017

#include "BoundedImpulseManeuver.h"

namespace EMTG
{
    namespace Phases
    {
        BoundedImpulseManeuver::BoundedImpulseManeuver() :
            dR_sc_Sun_dR_CB_sun(3, math::identity),
            dR_sc_Sun_dR_sc_CB(3, math::identity),
            R_sc_Sun(3, 1, 0.0),
            R_CB_Sun(3, 1, 0.0),
            R_sc_CB(3, 1, 0.0),
            dR_CB_Sun_dt(3, 1, 0.0),
            SunIsCentralBody(true),
            electric_propellant_used(0.0),
            ACS_fuel_used(0.0),
            dChemicalFuel_dThrottleComponents(math::Matrix<double>(3, 1, 0.0)),
            dElectricPropellant_dThrottleComponents(math::Matrix<double>(3, 1, 0.0)),
            NaturalAcceleration(3, 1, 0.0),
            dNaturalAcceleration_dt(3, 1, 0.0),
            dNaturalAcceleration_dState(3, 7, 0.0),
            NaturalDeltaV(3, 1, 0.0),
            dNaturalDeltaV_dPreviousTime(3, 1, 0.0),
            dNaturalDeltaV_dPropagationTime(3, 1, 0.0),
            dNaturalDeltaV_dState(3, 7, 0.0)
        {}

        BoundedImpulseManeuver::BoundedImpulseManeuver(const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stageIndex,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            BoundedImpulseManeuver()
        {
            this->journeyIndex = journeyIndex;
            this->phaseIndex = phaseIndex;
            this->stageIndex = stageIndex;
            this->myUniverse = myUniverse;
            this->mySpacecraft = mySpacecraft;
            this->myOptions = myOptions;
            this->myJourneyOptions = &myOptions->Journeys[this->journeyIndex];

            if (this->myUniverse->central_body_SPICE_ID == 10)
            {
                this->SunIsCentralBody = true;
            }
            else
            {
                this->SunIsCentralBody = false;
            }

            if (this->myOptions->perturb_SRP || this->myOptions->perturb_thirdbody || this->myOptions->perturb_J2 ||
				this->myJourneyOptions->perturb_drag)
            {
                this->hasPerturbations = true;
            }
            else
            {
                this->hasPerturbations = false;
            }
        }//end BoundedImpulseManeuver()

        void BoundedImpulseManeuver::compute_thruster_performance(const bool& needDerivatives)
        {
            //Step 0: dereference our pointers to local references
            doubleType& max_thrust = *this->max_thrust_pointer;
            doubleType& max_mass_flow_rate = *this->max_mass_flow_rate_pointer;
            doubleType& Isp = *this->Isp_pointer;
            doubleType& power = *this->power_pointer;
            doubleType& active_power = *this->active_power_pointer;
            size_t& number_of_active_engines = *this->number_of_active_engines_pointer;
            int& ThrottleLevel = *this->ThrottleLevel_pointer;
            std::string& ThrottleLevelString = *this->ThrottleLevelString_pointer;
            math::Matrix<doubleType>& Control = *this->Control_pointer;
            doubleType& LaunchDate = *this->LaunchDate_pointer;
            doubleType& ThrustStepLength = *this->ThrustStepLength_pointer;
            double& dThrustStepLength_dPropagationVariable = *this->dThrustStepLength_dPropagationVariable_pointer;
            double& dImpulseEpoch_dPropagationVariable = *this->dImpulseEpoch_dPropagationVariable_pointer;

            //Step 1: call the acceleration model to get the natural accelerations
            if (this->hasPerturbations)
            {
                this->mySpacecraftAccelerationModel->setEpoch(this->ManeuverEpoch);
                this->mySpacecraftAccelerationModel->computeAcceleration(*this->state_at_natural_perturbation, needDerivatives);

                this->NaturalAcceleration = this->mySpacecraftAccelerationModel->getAccelerationVec();
                math::Matrix<double> fx = this->mySpacecraftAccelerationModel->getfx();
                for (size_t i : { 0, 1, 2 })
                {
                    for (size_t j : {0, 1, 2, 3, 4, 5, 6})
                    {
                        this->dNaturalAcceleration_dState(i, j) = fx(i + 3, j);
                    }
                    this->dNaturalAcceleration_dt(i) = fx(i + 3, 7);
                }

                //Step 2: compute the natural delta-v vector and its partial derivatives with respect to state and time
                this->NaturalDeltaV = this->NaturalAcceleration * ThrustStepLength;

                this->dNaturalDeltaV_dState = this->dNaturalAcceleration_dState * ThrustStepLength _GETVALUE;

                this->dNaturalDeltaV_dPreviousTime = this->dNaturalAcceleration_dt * ThrustStepLength _GETVALUE;

                this->dNaturalDeltaV_dPropagationTime = this->NaturalAcceleration * dThrustStepLength_dPropagationVariable; //no need for the epoch derivative because it gets taken care of by dNaturalDeltaV_dPreviousTime and STM chains
            }

            //Step 3: where am I relative to the sun? - I need this for derivatives
            if (this->SunIsCentralBody)
            {
                this->R_sc_Sun = this->R_sc_CB;
            }
            else
            {
                //where is the central body relative to the sun?
                doubleType central_body_state_and_derivatives[12];
                this->myUniverse->locate_central_body(this->ManeuverEpoch, 
                    central_body_state_and_derivatives,
                    *this->myOptions, 
                    needDerivatives);

                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                }

                if (needDerivatives)
                {
                    for (size_t stateIndex : {0, 1, 2})
                        this->dR_CB_Sun_dt(stateIndex) = central_body_state_and_derivatives[6 + stateIndex]_GETVALUE;
                }

                this->R_sc_Sun = this->R_sc_CB + R_CB_Sun;
            }
            
            //Step 4: call the power model
            doubleType r_sc_sun_AU = this->R_sc_Sun.norm() / this->myOptions->AU;
            this->mySpacecraft->computePowerState(r_sc_sun_AU, this->ManeuverEpoch);

            //Step 5: call the thruster model
            if (Control.get_n() == 4)
            {
                this->mySpacecraft->computeElectricPropulsionPerformance(this->DutyCycle, Control(3));
            }
            else
            {
                this->mySpacecraft->computeElectricPropulsionPerformance(this->DutyCycle);
            }

            //Step 6: store the thruster model outputs
            max_thrust = this->mySpacecraft->getEPthrust() * 1.0e-3;
            max_mass_flow_rate = this->mySpacecraft->getEPMassFlowRate();
            Isp = this->mySpacecraft->getEPIsp();
            power = this->mySpacecraft->getAvailablePower();
            active_power = this->mySpacecraft->getEPActivePower();
            number_of_active_engines = this->mySpacecraft->getEPNumberOfActiveThrusters();
            ThrottleLevel = this->mySpacecraft->getEPThrottleLevel();
            ThrottleLevelString = this->mySpacecraft->getEPThrottleLevelString();

            //Step 7: store the derivatives of the thruster model
            if (needDerivatives)
            {
                double& LU = this->myUniverse->LU;
                double& TU = this->myUniverse->TU;

                this->dTdP = this->mySpacecraft->getEPdTdP() * 1.0e-3;
                this->dmdotdP = this->mySpacecraft->getEPdMassFlowRatedP();
                this->dT_du_command = this->mySpacecraft->getEPdTdu_command() * 1.0e-3;
                this->dmdot_du_command = this->mySpacecraft->getEPdMassFlowRatedu_command();
                this->dPdr_sun = this->mySpacecraft->getdPdr();
                double drdt = ((  this->dR_CB_Sun_dt(0) * this->R_sc_Sun(0)
                                + this->dR_CB_Sun_dt(1) * this->R_sc_Sun(1)
                                + this->dR_CB_Sun_dt(2) * this->R_sc_Sun(2)) / this->R_sc_Sun.norm()) _GETVALUE;
                this->dPdt = this->mySpacecraft->getdPdt()
                    + this->dPdr_sun / this->myOptions->AU * drdt;
            }
        }//end compute_thruster_performance()
    }//close namespace Astrodynamics
}//close namespace EMTG