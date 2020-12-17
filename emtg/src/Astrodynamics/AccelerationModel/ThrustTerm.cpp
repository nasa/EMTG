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

//central force term

#include "SpacecraftAccelerationModel.h"
#include "doubleType.h"
#include "ThrustTerm.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        //constructors
        ThrustTerm::ThrustTerm(SpacecraftAccelerationModel * acceleration_model_in) : 
            SpacecraftAccelerationModelTerm::SpacecraftAccelerationModelTerm(acceleration_model_in)
        {
            this->r_cb2sc.resize(3, 1, 0.0);
            this->v_cb2sc.resize(3, 1, 0.0);

            // by default, the thrust term will assume it is following an input Cartesian unit vector
            this->control_law = ThrustControlLaw::CartesianUnit;
        }

        ThrustTerm::~ThrustTerm() {}

        void ThrustTerm::processControlLogic()
        {
            if (this->control_law == ThrustControlLaw::CartesianUnit)
            {
                // control vector contains {ux, uy, uz, u_command}  
                this->u_x = this->acceleration_model->control(0);
                this->u_y = this->acceleration_model->control(1);
                this->u_z = this->acceleration_model->control(2);
                this->u_command = this->acceleration_model->control.get_n() > 3 ? this->acceleration_model->control(3) : 0.0; //TODO this could be faster
            }
            else
            {
                throw std::invalid_argument("The selected ThrustTerm control law is not implemented.");
            }
        }

        //methods
        void ThrustTerm::computeAccelerationTerm()
        {
            // populate some local spacecraft state vector information
            this->r_cb2sc = this->acceleration_model->r_cb2sc;
            this->r_cb2sc_norm = this->acceleration_model->r_cb2sc_norm;
            this->v_cb2sc = this->acceleration_model->v_cb2sc;
            this->v_cb2sc_norm = this->acceleration_model->v_cb2sc_norm;

            this->processControlLogic();

            HardwareModels::Spacecraft * my_spacecraft = this->acceleration_model->my_spacecraft;

            //***********************************************************
            //
            // Poll the power model
            //
            //***********************************************************

            // compute the spacecraft power state
            my_spacecraft->computePowerState(this->acceleration_model->r_sun2sc_norm / this->acceleration_model->my_options->AU, this->acceleration_model->current_epoch);

            // compute the spacecraft propulsion state
            my_spacecraft->computeElectricPropulsionPerformance(this->acceleration_model->duty_cycle, this->u_command);

            // extract various data
            this->max_thrust = my_spacecraft->getEPthrust() * 1.0e-3; // N to kN conversion
            this->max_mass_flow_rate = my_spacecraft->getEPMassFlowRate(); // includes duty cycle
            this->power = my_spacecraft->getProducedPower();
            this->number_of_active_engines = my_spacecraft->getEPNumberOfActiveThrusters();

            //****************************************
            //
            // Compute thruster acceleration term
            //
            //****************************************

            // modify the thrust by the duty cycle of the engine - this is already done inside the spacecraft model

            // compute the thrust...control vector consists of throttle decision parameters
            doubleType thrust_acceleration_vector[3];
            thrust_acceleration_vector[0] = this->u_x * this->max_thrust / this->acceleration_model->spacecraft_mass;
            thrust_acceleration_vector[1] = this->u_y * this->max_thrust / this->acceleration_model->spacecraft_mass;
            thrust_acceleration_vector[2] = this->u_z * this->max_thrust / this->acceleration_model->spacecraft_mass;
            for (size_t k = 0; k < 3; ++k) 
            {
                this->term_acceleration(k) = thrust_acceleration_vector[k];
                this->acceleration_model->acceleration(k) += this->term_acceleration(k);
            }

        }// end computeAccelerationTerm()

        void ThrustTerm::computeAccelerationTerm(const bool & generate_derivatives)
        {
            this->computeAccelerationTerm();

            HardwareModels::Spacecraft * my_spacecraft = this->acceleration_model->my_spacecraft;

            // extract partial derivative data from the spacecraft / hardware models
            this->dTdP = my_spacecraft->getEPdTdP() * 1.0e-3; // N to kN conversion
            this->dmdotdP = my_spacecraft->getEPdMassFlowRatedP();
            this->dTdu_command = my_spacecraft->getEPdTdu_command() * 1.0e-3; // N to kN conversion
            this->dmdotdu_command = my_spacecraft->getEPdMassFlowRatedu_command();
            this->dPdr_sun = my_spacecraft->getdPdr();
            this->dPdt = my_spacecraft->getdPdt();

            // We pick up this extra AU scale factor in the derivative
            // because the power model takes distance from the Sun in AUs as an input
            this->dPdr_sun = (this->dPdr_sun) / this->acceleration_model->my_options->AU;

            const double epsilon = 1.0e-25; // leak parameter to avoid divide by zero if thruster is off
            this->acceleration_model->control_norm = sqrt(this->u_x * this->u_x + this->u_y * this->u_y + this->u_z * this->u_z + epsilon);
            doubleType one_over_u_norm_plus_epsilon = 1.0 / this->acceleration_model->control_norm;
            double du_normdux = (u_x * one_over_u_norm_plus_epsilon) _GETVALUE;
            double du_normduy = (u_y * one_over_u_norm_plus_epsilon) _GETVALUE;
            double du_normduz = (u_z * one_over_u_norm_plus_epsilon) _GETVALUE;


            // NOTE: dTdP includes the multiplication by the duty cycle!
            doubleType D_dTdP_dPdr_over_msc = (this->dTdP) * (this->dPdr_sun) / this->acceleration_model->spacecraft_mass; 

            ////////////////////////////////////////////////////////////
            // state propagation matrix contributions from the thruster
            ////////////////////////////////////////////////////////////

            // A21 dadr (everything except for the third body terms, they are handled in GravityTerm)
            this->acceleration_model->fx(3, 0) += (u_x * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sc(0)) _GETVALUE;
            this->acceleration_model->fx(3, 1) += (u_x * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sc(1)) _GETVALUE;
            this->acceleration_model->fx(3, 2) += (u_x * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sc(2)) _GETVALUE;
            this->acceleration_model->fx(4, 0) += (u_y * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sc(0)) _GETVALUE;
            this->acceleration_model->fx(4, 1) += (u_y * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sc(1)) _GETVALUE;
            this->acceleration_model->fx(4, 2) += (u_y * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sc(2)) _GETVALUE;
            this->acceleration_model->fx(5, 0) += (u_z * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sc(0)) _GETVALUE;
            this->acceleration_model->fx(5, 1) += (u_z * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sc(1)) _GETVALUE;
            this->acceleration_model->fx(5, 2) += (u_z * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sc(2)) _GETVALUE;

            // A31 dmdotdr
            doubleType u_norm_D_dmdot_dP = this->acceleration_model->control_norm * this->dmdotdP;
            this->acceleration_model->fx(6, 0) += (-u_norm_D_dmdot_dP * this->dPdr_sun * this->acceleration_model->dr_sun2sc_normdr_cb2sc(0)) _GETVALUE;
            this->acceleration_model->fx(6, 1) += (-u_norm_D_dmdot_dP * this->dPdr_sun * this->acceleration_model->dr_sun2sc_normdr_cb2sc(1)) _GETVALUE;
            this->acceleration_model->fx(6, 2) += (-u_norm_D_dmdot_dP * this->dPdr_sun * this->acceleration_model->dr_sun2sc_normdr_cb2sc(2)) _GETVALUE;

            // dVirtualElectricPropellantdr
            this->acceleration_model->fx(9, 0) += (u_norm_D_dmdot_dP * this->dPdr_sun * this->acceleration_model->dr_sun2sc_normdr_cb2sc(0)) _GETVALUE;
            this->acceleration_model->fx(9, 1) += (u_norm_D_dmdot_dP * this->dPdr_sun * this->acceleration_model->dr_sun2sc_normdr_cb2sc(1)) _GETVALUE;
            this->acceleration_model->fx(9, 2) += (u_norm_D_dmdot_dP * this->dPdr_sun * this->acceleration_model->dr_sun2sc_normdr_cb2sc(2)) _GETVALUE;
            
            // A23 dadm
            doubleType D_thrust_over_msc2 = max_thrust / (this->acceleration_model->spacecraft_mass  * this->acceleration_model->spacecraft_mass);
            this->acceleration_model->fx(3, 6) += (-u_x * D_thrust_over_msc2) _GETVALUE;
            this->acceleration_model->fx(4, 6) += (-u_y * D_thrust_over_msc2) _GETVALUE;
            this->acceleration_model->fx(5, 6) += (-u_z * D_thrust_over_msc2) _GETVALUE;


            // explicit time partials
            // we need to account for our motion relative to the sun dr_cb2sc_dt is handled via A21 above
            // dr_cb2sun_dt must be handled here
            this->acceleration_model->fx(3, 7) += (u_x * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sun(0) * this->acceleration_model->dr_cb2sundcurrent_epoch(0)
                                                 + u_x * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sun(1) * this->acceleration_model->dr_cb2sundcurrent_epoch(1)
                                                 + u_x * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sun(2) * this->acceleration_model->dr_cb2sundcurrent_epoch(2)) _GETVALUE;
            this->acceleration_model->fx(4, 7) += (u_y * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sun(0) * this->acceleration_model->dr_cb2sundcurrent_epoch(0)
                                                 + u_y * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sun(1) * this->acceleration_model->dr_cb2sundcurrent_epoch(1)
                                                 + u_y * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sun(2) * this->acceleration_model->dr_cb2sundcurrent_epoch(2)) _GETVALUE;
            this->acceleration_model->fx(5, 7) += (u_z * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sun(0) * this->acceleration_model->dr_cb2sundcurrent_epoch(0)
                                                 + u_z * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sun(1) * this->acceleration_model->dr_cb2sundcurrent_epoch(1)
                                                 + u_z * D_dTdP_dPdr_over_msc * this->acceleration_model->dr_sun2sc_normdr_cb2sun(2) * this->acceleration_model->dr_cb2sundcurrent_epoch(2)) _GETVALUE;

            // for mass flow, again dr_cb2sc_dt is handled in A31
            // dr_cb2sun_dt must be handled here
            this->acceleration_model->fx(6, 7) += -(u_norm_D_dmdot_dP * this->dPdr_sun * ( this->acceleration_model->dr_sun2sc_normdr_cb2sun(0) * this->acceleration_model->dr_cb2sundcurrent_epoch(0)
                                                                                        + this->acceleration_model->dr_sun2sc_normdr_cb2sun(1) * this->acceleration_model->dr_cb2sundcurrent_epoch(1)
                                                                                        + this->acceleration_model->dr_sun2sc_normdr_cb2sun(2) * this->acceleration_model->dr_cb2sundcurrent_epoch(2))) _GETVALUE;
            this->acceleration_model->fx(6, 7) += -(u_norm_D_dmdot_dP * this->dPdt) _GETVALUE;

            // similarly for the tank partials
            this->acceleration_model->fx(9, 7) += (u_norm_D_dmdot_dP * this->dPdr_sun * ( this->acceleration_model->dr_sun2sc_normdr_cb2sun(0) * this->acceleration_model->dr_cb2sundcurrent_epoch(0)
                                                                                       + this->acceleration_model->dr_sun2sc_normdr_cb2sun(1) * this->acceleration_model->dr_cb2sundcurrent_epoch(1)
                                                                                       + this->acceleration_model->dr_sun2sc_normdr_cb2sun(2) * this->acceleration_model->dr_cb2sundcurrent_epoch(2))) _GETVALUE;
            this->acceleration_model->fx(9, 7) += (u_norm_D_dmdot_dP * this->dPdt) _GETVALUE;


            if (this->control_law == ThrustControlLaw::CartesianUnit)
            {
                // A24 dadu
                doubleType D_thrust_over_msc = max_thrust / this->acceleration_model->spacecraft_mass;
                this->acceleration_model->fx(3, 10) += (D_thrust_over_msc)_GETVALUE;
                this->acceleration_model->fx(3, 11) += 0.0;
                this->acceleration_model->fx(3, 12) += 0.0;
                this->acceleration_model->fx(4, 10) += 0.0;
                this->acceleration_model->fx(4, 11) += (D_thrust_over_msc)_GETVALUE;
                this->acceleration_model->fx(4, 12) += 0.0;
                this->acceleration_model->fx(5, 10) += 0.0;
                this->acceleration_model->fx(5, 11) += 0.0;
                this->acceleration_model->fx(5, 12) += (D_thrust_over_msc)_GETVALUE;

                // A34 dmdotdu
                doubleType D_mdot = this->max_mass_flow_rate;
                this->acceleration_model->fx(6, 10) += (-du_normdux * D_mdot) _GETVALUE;
                this->acceleration_model->fx(6, 11) += (-du_normduy * D_mdot) _GETVALUE;
                this->acceleration_model->fx(6, 12) += (-du_normduz * D_mdot) _GETVALUE;

                // dVirtualElectricPropellantdu
                this->acceleration_model->fx(9, 10) += (du_normdux * D_mdot) _GETVALUE;
                this->acceleration_model->fx(9, 11) += (du_normduy * D_mdot) _GETVALUE;
                this->acceleration_model->fx(9, 12) += (du_normduz * D_mdot) _GETVALUE;
            }

        }// end computeAccelerationTerm(bool)

        void ThrustTerm::populateInstrumentationFile(std::ofstream & acceleration_model_file)
        {
            acceleration_model_file << "," << sqrt(this->term_acceleration(0)*this->term_acceleration(0)
                                                 + this->term_acceleration(1)*this->term_acceleration(1)
                                                 + this->term_acceleration(2)*this->term_acceleration(2));
            for (size_t k = 0; k < this->term_acceleration.get_n(); ++k)
            {
                acceleration_model_file << "," << this->term_acceleration(k);
            }
        }
    }//close namespace Astrodynamics
}//close namespace EMTG