// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2014 - 2017 United States Government as represented by the
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
//Donald Ellison August 1st 2015

#include <cmath>
#include <iostream>
#include <fstream>

#include "FBLT_acceleration_model.h"
#include "doubleType.h"
#include "universe.h"
#include "missionoptions.h"
#include "Spacecraft.h"
#include "EMTG_time_utilities.h"
#include "EMTG_math.h"

namespace EMTG {
    namespace Astrodynamics {

        int FBLT_acceleration_model(const EMTG::missionoptions& options,
            Astrodynamics::universe & Universe,
            HardwareModels::Spacecraft * mySpacecraft,
            std::vector <doubleType> & spacecraft_state_relative_to_central_body,
            EMTG::math::Matrix <double> & dspacecraft_state_relative_to_central_bodydTOF,
            const doubleType & epoch_step_left,
            std::vector <double> & depoch_left_segmentdTOF,
            const double & c,
            const doubleType & h,
            const double & dhdTOF,
            const doubleType & launch_epoch,
            const std::vector <doubleType> & control,
            const doubleType & u_command,
            const double& DutyCycle,
            EMTG::math::Matrix <double> & dfdTOF,
            doubleType & max_thrust,
            doubleType & max_mass_flow_rate,
            doubleType & Isp,
            doubleType & power,
            doubleType & active_power,
            int & number_of_active_engines,
            EMTG::math::Matrix<doubleType> & fx,
            std::vector <doubleType> & acceleration_vector,
            const bool & generate_derivatives,
            const int& j,
            const int& p,
            const int& step)
        {

            doubleType epoch = epoch_step_left + c * h;
            doubleType epoch_JD = (epoch / 86400 + 2400000.5) _GETVALUE;
            //clear the acceleration vector otherwise it ramps up absurdly!
            for (size_t k = 0; k < 3; ++k)
                acceleration_vector[k] = 0.0;

            //If we are orbiting the Sun, the vector and derivative math is simplified
            bool isHeliocentric = Universe.central_body_SPICE_ID == 10;

            //derivative of current epoch w.r.t. TOF
            //Not in paper, this is an implementation issue that lets us know
            //EXACTLY what the epoch is from the integrator's point of reference
            static std::vector<double> depochdTOF(2);
            depochdTOF[0] = depoch_left_segmentdTOF[0] + c * dhdTOF;
            depochdTOF[1] = depoch_left_segmentdTOF[1];

            //we need to know how far we are from the Sun
            //for the purposes of the solar electric power model
            static std::vector<doubleType> spacecraft_position_relative_to_sun(3);

            static doubleType spacecraft_distance_from_sun;
            static doubleType spacecraft_distance_from_sun_in_AU;

            static std::vector<doubleType> central_body_state(12, 0.0);
            static EMTG::math::Matrix <double> dcentral_body_state_dTOF(6, 2, 0.0);

            static std::vector<doubleType> sun_state_relative_to_central_body(12, 0.0);
            static EMTG::math::Matrix <doubleType> dsun_state_relative_to_central_body_dTOF(6, 2, 0.0);

            doubleType drdx;
            doubleType drdy;
            doubleType drdz;

            //dadr3B -- sensitivity of s/c acceleration to third body position
            static std::vector <doubleType> dadr3B(9, 0.0);
            static std::vector <doubleType> nbodyTOFterms_x(2, 0.0);
            static std::vector <doubleType> nbodyTOFterms_y(2, 0.0);
            static std::vector <doubleType> nbodyTOFterms_z(2, 0.0);


            //we must re-initialise these to be zero since we declared them as static
            for (size_t i = 0; i < 2; ++i) {
                nbodyTOFterms_x[i] = 0.0;
                nbodyTOFterms_y[i] = 0.0;
                nbodyTOFterms_z[i] = 0.0;
            }


            doubleType spacecraft_distance_from_central_body;
            spacecraft_distance_from_central_body = sqrt(spacecraft_state_relative_to_central_body[0] * spacecraft_state_relative_to_central_body[0] +
                spacecraft_state_relative_to_central_body[1] * spacecraft_state_relative_to_central_body[1] +
                spacecraft_state_relative_to_central_body[2] * spacecraft_state_relative_to_central_body[2]);

            doubleType spacecraft_mass = spacecraft_state_relative_to_central_body[6];

            //If the Sun is NOT the central body
            if (!isHeliocentric) {
                // locate the central body's position w.r.t. the Sun
                // SPICE is a choke point for the AD routine, we must stop AD at this point and pick it back up again on the other side of the spkez call

#ifdef AD_INSTRUMENTATION
                Universe.locate_central_body(epoch, central_body_state.data(), options, true);
                std::vector<size_t> timevars = epoch.getDerivativeIndicies();
                for (size_t i = 0; i < 6; ++i)
                {
                    central_body_state[i].setDerivative(timevars[0], central_body_state[i].getDerivative(timevars[0]) * depochdTOF[1]);
                    central_body_state[i].setDerivative(timevars[1], central_body_state[i].getDerivative(timevars[1]) * depochdTOF[1]);
                    central_body_state[i].setDerivative(timevars[2], central_body_state[i].getDerivative(timevars[2]) * depochdTOF[0]);

                    // this is a hack for now, see the note below in the equivalent third body code
                    if (timevars.back() == 10)
                    {
                        central_body_state[i].setDerivative(timevars.back(), central_body_state[i].getDerivative(timevars.back()) * depochdTOF[0]);
                    }
                }

                // load the central body's position into the adouble vector and create the normalized version while we're at it
                for (size_t i = 0; i < 12; ++i)
                {
                    //    central_body_state[i].setValue(central_body_state_temp[i]);
                    sun_state_relative_to_central_body[i] = -central_body_state[i];
                }


#else
                Universe.locate_central_body(epoch, central_body_state.data(), options, generate_derivatives);
                for (size_t i = 0; i < 12; ++i) {
                    sun_state_relative_to_central_body[i] = -central_body_state[i];
                }
#endif

                // Form the TOF derivatives for the central body
                // These are needed in Eq. 32  and 38 (second derivative of second term for both)
                /*for (size_t i = 0; i < 3; ++i) {
                dcentral_body_state_dTOF(i, 0) = central_body_state[i + 3] * depochdTOF[0];
                dcentral_body_state_dTOF(i, 1) = central_body_state[i + 3] * depochdTOF[1];

                dcentral_body_state_dTOF(i + 3, 0) = central_body_state[i + 9] * depochdTOF[0];
                dcentral_body_state_dTOF(i + 3, 1) = central_body_state[i + 9] * depochdTOF[1];
                }*/


                for (size_t i = 0; i < 3; ++i)
                {
                    dcentral_body_state_dTOF(i, 0) = central_body_state[i + 6]_GETVALUE * depochdTOF[0];
                    dcentral_body_state_dTOF(i, 1) = central_body_state[i + 6]_GETVALUE * depochdTOF[1];

                    dcentral_body_state_dTOF(i + 3, 0) = central_body_state[i + 9]_GETVALUE * depochdTOF[0];
                    dcentral_body_state_dTOF(i + 3, 1) = central_body_state[i + 9]_GETVALUE * depochdTOF[1];
                }


                //The TOF derivatives of the Sun's state are just the negation of the central body TOF derivatives
                //The paper's frame of reference is fixed to the central body, SPICE's is fixed to the Sun
                //To stay consistent with the paper we'll create this extra variable
                for (size_t i = 0; i < 6; ++i) {
                    dsun_state_relative_to_central_body_dTOF(i, 0) = -dcentral_body_state_dTOF(i, 0);
                    dsun_state_relative_to_central_body_dTOF(i, 1) = -dcentral_body_state_dTOF(i, 1);
                }

                // Since we are not orbiting the Sun, our position with respect to it and the central body are distinct
                // This is Equation 5 in the paper although SPICE uses the Sun as its absolute reference frame so we add the two vectors
                spacecraft_position_relative_to_sun[0] = (spacecraft_state_relative_to_central_body[0] - sun_state_relative_to_central_body[0]);
                spacecraft_position_relative_to_sun[1] = (spacecraft_state_relative_to_central_body[1] - sun_state_relative_to_central_body[1]);
                spacecraft_position_relative_to_sun[2] = (spacecraft_state_relative_to_central_body[2] - sun_state_relative_to_central_body[2]);

                // Equation 6
                spacecraft_distance_from_sun = sqrt(spacecraft_position_relative_to_sun[0] * spacecraft_position_relative_to_sun[0] +
                    spacecraft_position_relative_to_sun[1] * spacecraft_position_relative_to_sun[1] +
                    spacecraft_position_relative_to_sun[2] * spacecraft_position_relative_to_sun[2]);

            }
            else {

                // if we are orbiting the Sun, it's quite silly to look up the position of the Sun relative to itself, so don't bother 
                spacecraft_position_relative_to_sun[0] = spacecraft_state_relative_to_central_body[0];
                spacecraft_position_relative_to_sun[1] = spacecraft_state_relative_to_central_body[1];
                spacecraft_position_relative_to_sun[2] = spacecraft_state_relative_to_central_body[2];

                //Equation 6
                spacecraft_distance_from_sun = spacecraft_distance_from_central_body;
            }

#ifdef REPORT_FORCE_MODEL
            std::ofstream outputfile("force_model.txt", ios::trunc);
            outputfile.precision(15);
            outputfile << scientific;
            outputfile << "EMTG model debug file" << endl;
            outputfile << endl;
            outputfile << "Epoch (seconds since J2000): " << epoch << endl;
            outputfile << "Epoch (MJD): " << epoch / 86400.0 << endl;
            outputfile << "Epoch (JD): " << epoch / 86400.0 + 2400000.5 << endl;
            outputfile << "Epoch (UTC): " << time_utilities::convert_ET_to_UTC_string(epoch - (51544.5 * 86400.0)) << endl;
            outputfile << "Central body: " << Universe.central_body_name << " (" << Universe.central_body_SPICE_ID << ")" << endl;
            outputfile << "Central body mu (km^2/s^3): " << Universe.mu << endl;
            outputfile << "Universe length unit (km): " << Universe.LU << endl;
            outputfile << endl;
            outputfile << "-------------Spacecraft State--------------------------------------" << endl;
            outputfile << "Spacecraft position relative to central body (km):" << endl;
            for (size_t k = 0; k < 3; ++k)
                outputfile << " " << spacecraft_state_relative_to_central_body[k];
            outputfile << endl;
            outputfile << "Spacecraft distance from central body (km): " << spacecraft_distance_from_central_body << endl;
            outputfile << "Spacecraft position relative to Sun (km):" << endl;
            for (size_t k = 0; k < 3; ++k)
                outputfile << " " << spacecraft_position_relative_to_sun[k];
            outputfile << endl;
            outputfile << "Spacecraft distance from sun (km): " << spacecraft_distance_from_sun << endl;
            outputfile << "Spacecraft mass (kg): " << spacecraft_mass << endl;
            outputfile << endl;
#endif
            //***********************************************************
            //
            //Call to the power model
            //
            //***********************************************************

            spacecraft_distance_from_sun_in_AU = spacecraft_distance_from_sun / options.AU;

            //compute the spacecraft power state
            mySpacecraft->computePowerState(spacecraft_distance_from_sun_in_AU, epoch_step_left - launch_epoch, options.power_margin);

            //compute the spacecraft propulsion state
            mySpacecraft->computeElectricPropulsionPerformance(DutyCycle, u_command);

            //extract various data
            max_thrust = mySpacecraft->getEPthrust() * 1.0e-3; //N to kN conversion
            max_mass_flow_rate = mySpacecraft->getEPMassFlowRate();
            power = mySpacecraft->getProducedPower();
            number_of_active_engines = mySpacecraft->getEPNumberOfActiveThrusters();
            double dTdP = mySpacecraft->getEPdTdP() * 1.0e-3; //N to kN conversion
            double dmdotdP = mySpacecraft->getEPdMassFlowRatedP();
            double dT_du_command = mySpacecraft->getEPdTdu_command() * 1.0e-3; //N to kN conversion
            double dmdot_du_command = mySpacecraft->getEPdMassFlowRatedu_command();
            double dPdr_sun = mySpacecraft->getdPdr();
            double dPdt = mySpacecraft->getdPdt();

#ifdef REPORT_FORCE_MODEL
            outputfile << "-------------Power and Propulsion----------------------------------" << endl;
            outputfile << "max_thrust (kN): " << max_thrust << endl;
            outputfile << "max_mass_flow_rate (kg/s): " << max_mass_flow_rate << endl;
            outputfile << "Isp (s): " << Isp << endl;
            outputfile << "power (kW): " << power << endl;
            outputfile << "active_power (kW): " << active_power << endl;
            outputfile << "number_of_active_engines: " << number_of_active_engines << endl;
            outputfile << endl;
#endif

            //The power model returns everything in MKS, so we have to convert to normalized units
            //max_thrust = (max_thrust) / MU / LU * TU * TU;
            //max_mass_flow_rate = (max_mass_flow_rate) / MU * TU;
            //we don't need to scale Isp, because none of the EOM's directly depend on it
            //we don't want to scale Isp, because it will mess up the output to file
            //*Isp = (*Isp) / TU;
            //dTdP = (dTdP) / MU / LU * TU * TU;
            //dmdotdP = (dmdotdP) / MU * TU;
            //dTdIsp = (dTdIsp) / MU / LU * TU * TU * TU;
            //dmdotdIsp = (dmdotdIsp) / MU * TU;
            //dPdr_sun = (dPdr_sun) / options.AU * LU;
            dPdr_sun = (dPdr_sun) / options.AU;
            //dPdt = (dPdt) * TU;



            //****************************************
            //
            //Force vector modification: thruster
            //
            //****************************************

            //modify the thrust by the duty cycle of the engine - this is already done inside the spacecraft model


            //compute the thrust...control vector consists of throttle decision parameters for this FBLT segment
            doubleType ThrustVector[3];
            for (size_t k = 0; k < 3; ++k) {
                ThrustVector[k] = control[k] * max_thrust / spacecraft_mass;
                acceleration_vector[k] += ThrustVector[k];
            }

#ifdef REPORT_FORCE_MODEL
            outputfile << "-------------Applied Thrust----------------------------------------" << endl;
            outputfile << "duty cycle: " << mySpacecraft->getDutyCycle() << endl;
            outputfile << "allowable thrust (kN): " << max_thrust << endl;
            outputfile << "control vector:" << endl;
            for (size_t k = 0; k < 3; ++k)
                outputfile << " " << control[k];
            outputfile << endl;
            outputfile << "applied thrust vector (kN):" << endl;
            for (size_t k = 0; k < 3; ++k)
                outputfile << " " << ThrustVector[k];
            outputfile << endl;
            outputfile << "thrust acceleration vector (km/s^2):" << endl;
            for (size_t k = 0; k < 3; ++k)
                outputfile << " " << ThrustVector[k] / spacecraft_mass;
            outputfile << endl;
            outputfile << endl;
#endif

            const double epsilon = 1.0e-25; //avoid divide by zero if thruster is off
            doubleType control_norm = sqrt(control[0] * control[0] + control[1] * control[1] + control[2] * control[2] + epsilon);

            if (generate_derivatives) {

                //gradient of the magnitude of the spacecraft position vector w.r.t. the Sun
                //Equation 48 of conference paper
                doubleType one_over_spacecraft_distance_from_sun = 1.0 / spacecraft_distance_from_sun;
                doubleType one_over_spacecraft_distance_from_central_body = 1.0 / spacecraft_distance_from_central_body;

                //if we are not orbiting the Sun, then this derivative is calculated slightly differently
                if (!isHeliocentric) {
                    drdx = spacecraft_position_relative_to_sun[0] * one_over_spacecraft_distance_from_sun;
                    drdy = spacecraft_position_relative_to_sun[1] * one_over_spacecraft_distance_from_sun;
                    drdz = spacecraft_position_relative_to_sun[2] * one_over_spacecraft_distance_from_sun;
                }
                else {
                    drdx = spacecraft_state_relative_to_central_body[0] * one_over_spacecraft_distance_from_central_body;
                    drdy = spacecraft_state_relative_to_central_body[1] * one_over_spacecraft_distance_from_central_body;
                    drdz = spacecraft_state_relative_to_central_body[2] * one_over_spacecraft_distance_from_central_body;
                }

                //****************************************
                //
                //State Propagation (fx) matrix calculation
                //
                //****************************************

                //gradient of the magnitude of the control vector
                //Equation 54 in paper
                static doubleType one_over_control_norm_plus_epsilon;
                one_over_control_norm_plus_epsilon = 1.0 / control_norm;
                doubleType dcontrol_normdux = control[0] * one_over_control_norm_plus_epsilon;
                doubleType dcontrol_normduy = control[1] * one_over_control_norm_plus_epsilon;
                doubleType dcontrol_normduz = control[2] * one_over_control_norm_plus_epsilon;

                //Top Row Identity
                //See Equation 45
                fx(0, 3) = 1.0;
                fx(1, 4) = 1.0;
                fx(2, 5) = 1.0;

                //A21 dadr (everything except for the third body terms, they are handled below)
                //Equation 46 (terms 1,2 and 5)
                double mu_CB = Universe.mu;
                doubleType spacecraft_distance_from_CB3 = spacecraft_distance_from_central_body * spacecraft_distance_from_central_body * spacecraft_distance_from_central_body;
                doubleType spacecraft_distance_from_CB5 = spacecraft_distance_from_CB3 * spacecraft_distance_from_central_body * spacecraft_distance_from_central_body;
                doubleType three_muCB_over_spacecraft_distance_from_CB5 = 3.0 * mu_CB / spacecraft_distance_from_CB5;
                doubleType muCB_over_spacecraft_distance_from_CB3 = 1.0 * mu_CB / spacecraft_distance_from_CB3;
                doubleType D_dTdP_dPdr_over_msc = (dTdP) * (dPdr_sun) / spacecraft_mass; // NOTE: dTdP includes the multiplication by the duty cycle
                std::vector<doubleType> State = spacecraft_state_relative_to_central_body;

                fx(3, 0) = three_muCB_over_spacecraft_distance_from_CB5 * State[0] * State[0] - muCB_over_spacecraft_distance_from_CB3 + control[0] * D_dTdP_dPdr_over_msc * drdx;
                fx(3, 1) = three_muCB_over_spacecraft_distance_from_CB5 * State[0] * State[1] + control[0] * D_dTdP_dPdr_over_msc * drdy;
                fx(3, 2) = three_muCB_over_spacecraft_distance_from_CB5 * State[0] * State[2] + control[0] * D_dTdP_dPdr_over_msc * drdz;
                fx(4, 0) = three_muCB_over_spacecraft_distance_from_CB5 * State[1] * State[0] + control[1] * D_dTdP_dPdr_over_msc * drdx;
                fx(4, 1) = three_muCB_over_spacecraft_distance_from_CB5 * State[1] * State[1] - muCB_over_spacecraft_distance_from_CB3 + control[1] * D_dTdP_dPdr_over_msc * drdy;
                fx(4, 2) = three_muCB_over_spacecraft_distance_from_CB5 * State[1] * State[2] + control[1] * D_dTdP_dPdr_over_msc * drdz;
                fx(5, 0) = three_muCB_over_spacecraft_distance_from_CB5 * State[2] * State[0] + control[2] * D_dTdP_dPdr_over_msc * drdx;
                fx(5, 1) = three_muCB_over_spacecraft_distance_from_CB5 * State[2] * State[1] + control[2] * D_dTdP_dPdr_over_msc * drdy;
                fx(5, 2) = three_muCB_over_spacecraft_distance_from_CB5 * State[2] * State[2] - muCB_over_spacecraft_distance_from_CB3 + control[2] * D_dTdP_dPdr_over_msc * drdz;


                //A23 dadm
                //Equation 49/50
                doubleType D_thrust_over_msc2 = max_thrust / (spacecraft_mass * spacecraft_mass);

                fx(3, 6) = -control[0] * D_thrust_over_msc2;
                fx(4, 6) = -control[1] * D_thrust_over_msc2;
                fx(5, 6) = -control[2] * D_thrust_over_msc2;


                //A24 dadu
                //Equation 51
                doubleType D_thrust_over_msc = max_thrust / spacecraft_mass;

                fx(3, 7) = D_thrust_over_msc;
                fx(3, 8) = 0.0;
                fx(3, 9) = 0.0;
                fx(4, 7) = 0.0;
                fx(4, 8) = D_thrust_over_msc;
                fx(4, 9) = 0.0;
                fx(5, 7) = 0.0;
                fx(5, 8) = 0.0;
                fx(5, 9) = D_thrust_over_msc;


                //A31 dmdotdr
                //Equation 52
                doubleType control_norm_D_dmdotdP_dPdr = control_norm * (dmdotdP) * (dPdr_sun);

                fx(6, 0) = -control_norm_D_dmdotdP_dPdr * drdx;
                fx(6, 1) = -control_norm_D_dmdotdP_dPdr * drdy;
                fx(6, 2) = -control_norm_D_dmdotdP_dPdr * drdz;


                //A34 dmdotdu
                //Equation 53
                doubleType D_mdot = (max_mass_flow_rate);

                fx(6, 7) = -dcontrol_normdux * D_mdot;
                fx(6, 8) = -dcontrol_normduy * D_mdot;
                fx(6, 9) = -dcontrol_normduz * D_mdot;
            }

            // compute J2 oblateness if desired
            if (options.perturb_J2) {

                //****************************************
                //
                //Force vector modification: central body J2
                //
                //****************************************

                // central body J2 perturbations
                // THE GRAVITY FIELD OF THE SATURNIAN SYSTEM FROM SATELLITE OBSERVATIONS AND SPACECRAFT TRACKING DATA
                // R. A. Jacobson et al. The Astronomical Journal 132:2520-2526, 2006 December.
                // double J2_saturn = 16290.71e-6.0; 

                // Saturn's Gravitational Field, Internal Rotation, and Interior Structure
                // John D. Anderson and Gerald Schubert, Science 07 Sep. 2007
                // Vol. 317, Issue 5843, pp. 1384 - 1387
                // DOI: 10.1126 / science.1144835
                double mu_CB = Universe.mu;
                double J2_CB = Universe.central_body_J2; // (unitless)
                double R_CB = Universe.central_body_radius; // saturn's radius in km
                double R_CB2 = R_CB * R_CB;

                // J2 calculations must be performed in the central-body fixed equatorial frame
                static std::vector <doubleType> spacecraft_position_relative_to_central_body_BCF(3, 0.0);

                static EMTG::math::Matrix <doubleType> temp_matrix_in(3, 1, 0.0);
                static EMTG::math::Matrix <doubleType> temp_matrix_out(3, 1, 0.0);

                // rotate position ICRF->BCF
                for (size_t k = 0; k < 3; ++k) {
                    temp_matrix_in(k, 0) = spacecraft_state_relative_to_central_body[k];
                }
                Universe.LocalFrame.rotate_frame_to_frame(EMTG::ReferenceFrame::ICRF,
                    temp_matrix_in,
                    EMTG::ReferenceFrame::TrueOfDate_BCF,
                    temp_matrix_out,
                    epoch_JD);
                for (size_t k = 0; k < 3; ++k) {
                    spacecraft_position_relative_to_central_body_BCF[k] = temp_matrix_out(k, 0);
                }

                doubleType spacecraft_distance_from_central_body_BCF = sqrt(spacecraft_position_relative_to_central_body_BCF[0] * spacecraft_position_relative_to_central_body_BCF[0] +
                    spacecraft_position_relative_to_central_body_BCF[1] * spacecraft_position_relative_to_central_body_BCF[1] +
                    spacecraft_position_relative_to_central_body_BCF[2] * spacecraft_position_relative_to_central_body_BCF[2]);

                doubleType spacecraft_distance_from_central_body_BCF2 = spacecraft_distance_from_central_body_BCF * spacecraft_distance_from_central_body_BCF;
                doubleType spacecraft_distance_from_central_body_BCF5 = spacecraft_distance_from_central_body_BCF * spacecraft_distance_from_central_body_BCF2 * spacecraft_distance_from_central_body_BCF2;

                // HAVE NOT YET HANDLED ROTATION OF STATE FROM EARTH EQUATORIAL TO CENTRAL-BODY-FIXED EQUATORIAL
                // We need to perform an in-place rotation of the J2000 Earth equatorial state vector to 
                // central-body fixed in order to perform the J2 force calculation
                // There will also be an explicit TOF dependence for the rotation matrix entries that must be accounted for
                //std::vector<doubleType> spacecraft_state_CB_equatorial_fixed = Universe.LocalFrame.
                doubleType z2 = spacecraft_position_relative_to_central_body_BCF[2] * spacecraft_position_relative_to_central_body_BCF[2];
                doubleType zcoeff = 5.0 * z2 / (spacecraft_distance_from_central_body_BCF2);
                doubleType J2coeff = -3.0 * J2_CB * mu_CB * R_CB2 / (2.0 * spacecraft_distance_from_central_body_BCF5);
                doubleType J2_acceleration_BCF[3];
                doubleType J2_acceleration_ICRF[3];

                // compute the J2 acceleration in BCF coordinates
                J2_acceleration_BCF[0] = J2coeff * (1.0 - zcoeff) * spacecraft_position_relative_to_central_body_BCF[0];
                J2_acceleration_BCF[1] = J2coeff * (1.0 - zcoeff) * spacecraft_position_relative_to_central_body_BCF[1];
                J2_acceleration_BCF[2] = J2coeff * (3.0 - zcoeff) * spacecraft_position_relative_to_central_body_BCF[2];

                // rotate the J2 acceleration vector BCF->ICRF
                for (size_t k = 0; k < 3; ++k) {
                    temp_matrix_in(k, 0) = J2_acceleration_BCF[k];
                }
                Universe.LocalFrame.rotate_frame_to_frame(EMTG::ReferenceFrame::TrueOfDate_BCF,
                    temp_matrix_in,
                    EMTG::ReferenceFrame::ICRF,
                    temp_matrix_out,
                    epoch_JD);
                for (size_t k = 0; k < 3; ++k) {
                    J2_acceleration_ICRF[k] = temp_matrix_out(k, 0);
                }
                for (size_t k = 0; k < 3; ++k) {
                    acceleration_vector[k] += J2_acceleration_ICRF[k];
                }


#ifdef REPORT_FORCE_MODEL
                outputfile << "-------------Central Body J2------------------------------" << endl;
                outputfile << "J2 force magnitude: " << EMTG::math::norm(J2_force, 3) << endl;
                outputfile << "J2 force vector (kN):" << endl;
                for (size_t k = 0; k < 3; ++k)
                    outputfile << " " << J2_acceleration_ICRF[k];
                outputfile << endl;
                outputfile << "central body J2 acceleration vector (km/s^2):" << endl;
                for (size_t k = 0; k < 3; ++k)
                    outputfile << " " << J2_force[k] / spacecraft_mass;
                outputfile << endl;
                outputfile << endl;
#endif

                if (generate_derivatives) {
                    std::vector<doubleType> State = spacecraft_position_relative_to_central_body_BCF;
                    doubleType xy = State[0] * State[1];
                    doubleType xz = State[0] * State[2];
                    doubleType yz = State[1] * State[2];
                    doubleType x2 = State[0] * State[0];
                    doubleType y2 = State[1] * State[1];
                    doubleType x4 = x2*x2;
                    doubleType y4 = y2*y2;
                    doubleType z4 = z2*z2;
                    doubleType x2y2 = x2*y2;
                    doubleType x2z2 = x2*z2;
                    doubleType y2z2 = y2*z2;
                    doubleType spacecraft_distance_from_central_body_BCF9 = spacecraft_distance_from_central_body_BCF5 * spacecraft_distance_from_central_body_BCF2 * spacecraft_distance_from_central_body_BCF2;
                    doubleType J2_deriv_coeff = J2_CB * mu_CB * R_CB2 / (2.0 * spacecraft_distance_from_central_body_BCF9);
                    doubleType J2_deriv_coeff3 = J2_deriv_coeff * J2_deriv_coeff * J2_deriv_coeff;
                    // grab the ICRF->BCF and BCF->ICRF frame rotation matrices
                    static EMTG::math::Matrix <doubleType> R_from_ICRF_to_BCF;
                    static EMTG::math::Matrix <doubleType> R_from_BCF_to_ICRF;
                    R_from_ICRF_to_BCF = Universe.LocalFrame.get_R_from_ICRF_to_BCF();
                    R_from_BCF_to_ICRF = Universe.LocalFrame.get_R_from_BCF_to_ICRF();


                    doubleType dax_BCFdx_BCF = 3.0  * J2_deriv_coeff * (4.0 * x4 + 3.0 * x2y2 - 27.0* x2z2 - y4 + 3.0 * y2z2 + 4.0 * z4);
                    doubleType dax_BCFdy_BCF = 15.0 * J2_deriv_coeff * xy * (x2 + y2 - 6.0 * z2);
                    doubleType dax_BCFdz_BCF = 15.0 * J2_deriv_coeff * xz * (3.0 * x2 + 3.0 * y2 - 4.0 * z2);
                    doubleType day_BCFdx_BCF = 15.0 * J2_deriv_coeff * xy * (x2 + y2 - 6.0 * z2);
                    doubleType day_BCFdy_BCF = 3.0  * J2_deriv_coeff * (-1.0 * x4 + 3.0 * x2y2 + 3.0 * x2z2 + 4.0 * y4 - 27.0 * y2z2 + 4.0 * z4);
                    doubleType day_BCFdz_BCF = 15.0 * J2_deriv_coeff * yz * (3.0 * x2 + 3.0 * y2 - 4.0 * z2);
                    doubleType daz_BCFdx_BCF = 15.0 * J2_deriv_coeff * xz * (3.0 * x2 + 3.0 * y2 - 4.0 * z2);
                    doubleType daz_BCFdy_BCF = 15.0 * J2_deriv_coeff * yz * (3.0 * x2 + 3.0 * y2 - 4.0 * z2);
                    doubleType daz_BCFdz_BCF = -3.0  * J2_deriv_coeff * (3.0 * x4 + 6.0 * x2y2 - 24.0 * x2z2 + 3.0 * y4 - 24.0 * y2z2 + 8.0 * z4);

                    doubleType dax_BCFdx_ICRF = dax_BCFdx_BCF * R_from_ICRF_to_BCF(0, 0) + dax_BCFdy_BCF * R_from_ICRF_to_BCF(1, 0) + dax_BCFdz_BCF * R_from_ICRF_to_BCF(2, 0);
                    doubleType dax_BCFdy_ICRF = dax_BCFdx_BCF * R_from_ICRF_to_BCF(0, 1) + dax_BCFdy_BCF * R_from_ICRF_to_BCF(1, 1) + dax_BCFdz_BCF * R_from_ICRF_to_BCF(2, 1);
                    doubleType dax_BCFdz_ICRF = dax_BCFdx_BCF * R_from_ICRF_to_BCF(0, 2) + dax_BCFdy_BCF * R_from_ICRF_to_BCF(1, 2) + dax_BCFdz_BCF * R_from_ICRF_to_BCF(2, 2);
                    doubleType day_BCFdx_ICRF = day_BCFdx_BCF * R_from_ICRF_to_BCF(0, 0) + day_BCFdy_BCF * R_from_ICRF_to_BCF(1, 0) + day_BCFdz_BCF * R_from_ICRF_to_BCF(2, 0);
                    doubleType day_BCFdy_ICRF = day_BCFdx_BCF * R_from_ICRF_to_BCF(0, 1) + day_BCFdy_BCF * R_from_ICRF_to_BCF(1, 1) + day_BCFdz_BCF * R_from_ICRF_to_BCF(2, 1);
                    doubleType day_BCFdz_ICRF = day_BCFdx_BCF * R_from_ICRF_to_BCF(0, 2) + day_BCFdy_BCF * R_from_ICRF_to_BCF(1, 2) + day_BCFdz_BCF * R_from_ICRF_to_BCF(2, 2);
                    doubleType daz_BCFdx_ICRF = daz_BCFdx_BCF * R_from_ICRF_to_BCF(0, 0) + daz_BCFdy_BCF * R_from_ICRF_to_BCF(1, 0) + daz_BCFdz_BCF * R_from_ICRF_to_BCF(2, 0);
                    doubleType daz_BCFdy_ICRF = daz_BCFdx_BCF * R_from_ICRF_to_BCF(0, 1) + daz_BCFdy_BCF * R_from_ICRF_to_BCF(1, 1) + daz_BCFdz_BCF * R_from_ICRF_to_BCF(2, 1);
                    doubleType daz_BCFdz_ICRF = daz_BCFdx_BCF * R_from_ICRF_to_BCF(0, 2) + daz_BCFdy_BCF * R_from_ICRF_to_BCF(1, 2) + daz_BCFdz_BCF * R_from_ICRF_to_BCF(2, 2);

                    doubleType dax_ICRFdx_ICRF = R_from_BCF_to_ICRF(0, 0) * dax_BCFdx_ICRF + R_from_BCF_to_ICRF(0, 1) * day_BCFdx_ICRF + R_from_BCF_to_ICRF(0, 2) * daz_BCFdx_ICRF;
                    doubleType dax_ICRFdy_ICRF = R_from_BCF_to_ICRF(0, 0) * dax_BCFdy_ICRF + R_from_BCF_to_ICRF(0, 1) * day_BCFdy_ICRF + R_from_BCF_to_ICRF(0, 2) * daz_BCFdy_ICRF;
                    doubleType dax_ICRFdz_ICRF = R_from_BCF_to_ICRF(0, 0) * dax_BCFdz_ICRF + R_from_BCF_to_ICRF(0, 1) * day_BCFdz_ICRF + R_from_BCF_to_ICRF(0, 2) * daz_BCFdz_ICRF;
                    doubleType day_ICRFdx_ICRF = R_from_BCF_to_ICRF(1, 0) * dax_BCFdx_ICRF + R_from_BCF_to_ICRF(1, 1) * day_BCFdx_ICRF + R_from_BCF_to_ICRF(1, 2) * daz_BCFdx_ICRF;
                    doubleType day_ICRFdy_ICRF = R_from_BCF_to_ICRF(1, 0) * dax_BCFdy_ICRF + R_from_BCF_to_ICRF(1, 1) * day_BCFdy_ICRF + R_from_BCF_to_ICRF(1, 2) * daz_BCFdy_ICRF;
                    doubleType day_ICRFdz_ICRF = R_from_BCF_to_ICRF(1, 0) * dax_BCFdz_ICRF + R_from_BCF_to_ICRF(1, 1) * day_BCFdz_ICRF + R_from_BCF_to_ICRF(1, 2) * daz_BCFdz_ICRF;
                    doubleType daz_ICRFdx_ICRF = R_from_BCF_to_ICRF(2, 0) * dax_BCFdx_ICRF + R_from_BCF_to_ICRF(2, 1) * day_BCFdx_ICRF + R_from_BCF_to_ICRF(2, 2) * daz_BCFdx_ICRF;
                    doubleType daz_ICRFdy_ICRF = R_from_BCF_to_ICRF(2, 0) * dax_BCFdy_ICRF + R_from_BCF_to_ICRF(2, 1) * day_BCFdy_ICRF + R_from_BCF_to_ICRF(2, 2) * daz_BCFdy_ICRF;
                    doubleType daz_ICRFdz_ICRF = R_from_BCF_to_ICRF(2, 0) * dax_BCFdz_ICRF + R_from_BCF_to_ICRF(2, 1) * day_BCFdz_ICRF + R_from_BCF_to_ICRF(2, 2) * daz_BCFdz_ICRF;

                    // A21 dadr
                    fx(3, 0) += dax_ICRFdx_ICRF;
                    fx(3, 1) += dax_ICRFdy_ICRF;
                    fx(3, 2) += dax_ICRFdz_ICRF;
                    fx(4, 0) += day_ICRFdx_ICRF;
                    fx(4, 1) += day_ICRFdy_ICRF;
                    fx(4, 2) += day_ICRFdz_ICRF;
                    fx(5, 0) += daz_ICRFdx_ICRF;
                    fx(5, 1) += daz_ICRFdy_ICRF;
                    fx(5, 2) += daz_ICRFdz_ICRF;
                }
            }

            // compute SRP perturbations if desired
            if (options.perturb_SRP) {

                //****************************************
                //
                //Force vector modification: SRP
                //
                //****************************************

                static double Cr = options.coefficient_of_reflectivity; //coefficient of reflectivity [0,2]
                static double scArea = options.spacecraft_area / (1000.0 * 1000.0); //spacecraft area exposed to sunlight km^2
                static double K = 1.0; //percentage of Sun seen by spacecraft [0,1]
                static double solarFlux = 1360.0; // W/m^2 = kg/s^3 --> equal to 3.846e26 Watts / (4 * pi * r_km^2)
                const static double c = 299792458.0 / 1000.0; //speed of light in a vacuum km/s

                static double SRPcoeff = Cr * scArea * K * solarFlux / c * (options.AU * options.AU);
                static doubleType FSRP;

                FSRP = SRPcoeff / (spacecraft_distance_from_sun * spacecraft_distance_from_sun); // normalized force

                doubleType SRP_Vector[3];
                for (size_t k = 0; k < 3; ++k) {
                    SRP_Vector[k] = FSRP * spacecraft_position_relative_to_sun[k] / spacecraft_distance_from_sun;
                    acceleration_vector[k] += SRP_Vector[k] / spacecraft_mass;
                }

#ifdef REPORT_FORCE_MODEL
                outputfile << "-------------Solar Radiation Pressure------------------------------" << endl;
                outputfile << "SRP force magnitude: " << FSRP << endl;
                outputfile << "SRP force vector (kN):" << endl;
                for (size_t k = 0; k < 3; ++k)
                    outputfile << " " << SRP_Vector[k];
                outputfile << endl;
                outputfile << "SRP acceleration vector (km/s^2):" << endl;
                for (size_t k = 0; k < 3; ++k)
                    outputfile << " " << SRP_Vector[k] / spacecraft_mass;
                outputfile << endl;
                outputfile << endl;
#endif

                doubleType one_over_r_sun_sc3 = 1.0 / (spacecraft_distance_from_sun * spacecraft_distance_from_sun * spacecraft_distance_from_sun);
                doubleType one_over_r_sun_sc5 = one_over_r_sun_sc3 * 1.0 / (spacecraft_distance_from_sun * spacecraft_distance_from_sun);
                doubleType SRPcoeff_over_msc = SRPcoeff / spacecraft_mass;
                std::vector<doubleType> State = spacecraft_position_relative_to_sun;

                if (generate_derivatives) {
                    // A21 dadr
                    fx(3, 0) += SRPcoeff_over_msc * (-3.0 * State[0] * State[0] * one_over_r_sun_sc5 + one_over_r_sun_sc3);
                    fx(3, 1) += SRPcoeff_over_msc * (-3.0 * State[0] * State[1] * one_over_r_sun_sc5);
                    fx(3, 2) += SRPcoeff_over_msc * (-3.0 * State[0] * State[2] * one_over_r_sun_sc5);
                    fx(4, 0) += SRPcoeff_over_msc * (-3.0 * State[1] * State[0] * one_over_r_sun_sc5);
                    fx(4, 1) += SRPcoeff_over_msc * (-3.0 * State[1] * State[1] * one_over_r_sun_sc5 + one_over_r_sun_sc3);
                    fx(4, 2) += SRPcoeff_over_msc * (-3.0 * State[1] * State[2] * one_over_r_sun_sc5);
                    fx(5, 0) += SRPcoeff_over_msc * (-3.0 * State[2] * State[0] * one_over_r_sun_sc5);
                    fx(5, 1) += SRPcoeff_over_msc * (-3.0 * State[2] * State[1] * one_over_r_sun_sc5);
                    fx(5, 2) += SRPcoeff_over_msc * (-3.0 * State[2] * State[2] * one_over_r_sun_sc5 + one_over_r_sun_sc3);

                    // A23 dadm
                    fx(3, 6) += -FSRP / (spacecraft_mass * spacecraft_mass) * spacecraft_position_relative_to_sun[0] / spacecraft_distance_from_sun;
                    fx(4, 6) += -FSRP / (spacecraft_mass * spacecraft_mass) * spacecraft_position_relative_to_sun[1] / spacecraft_distance_from_sun;
                    fx(5, 6) += -FSRP / (spacecraft_mass * spacecraft_mass) * spacecraft_position_relative_to_sun[2] / spacecraft_distance_from_sun;
                }

            }

            //compute third body perturbations if desired
            if (options.perturb_thirdbody) {
                doubleType Agravity = 0.0;
#ifdef REPORT_FORCE_MODEL
                outputfile << "-------------Third-Body Perturbations------------------------------" << endl;
                outputfile << endl;
#endif

                //perturbations due to the Sun
                //don't bother calculating these if we are orbiting the sun; that would be silly
                //if we are not orbiting the Sun, then we ALWAYS calculate Sun perturbations

                if (!isHeliocentric) {

                    // Force of gravity due to the Sun in normalized units
                    double mu_sun = 1.32712440018e11;
                    Agravity = -(mu_sun) / (spacecraft_distance_from_sun * spacecraft_distance_from_sun);
                    doubleType sun_distance_from_central_body = math::norm(sun_state_relative_to_central_body.data(), 3);
                    doubleType a_sun_on_central_body = -(mu_sun) / (sun_distance_from_central_body*sun_distance_from_central_body);

                    //****************************************
                    //
                    //Force vector modification: Sun's gravity
                    //
                    //****************************************
                    doubleType a3rdbody[3];
                    for (size_t k = 0; k < 3; ++k) {
                        a3rdbody[k] = Agravity * spacecraft_position_relative_to_sun[k] / spacecraft_distance_from_sun + a_sun_on_central_body * sun_state_relative_to_central_body[k] / sun_distance_from_central_body;
                        acceleration_vector[k] += a3rdbody[k];
                    }

#ifdef REPORT_FORCE_MODEL
                    outputfile << "-------------Sun------" << endl;
                    outputfile << "Sun force magnitude (kN): " << fabs(Fgravity) << endl;
                    outputfile << "Sun force vector (kN):" << endl;
                    for (size_t k = 0; k < 3; ++k)
                        outputfile << " " << F3rdbody[k];
                    outputfile << endl;
                    outputfile << "Sun acceleration magnitude (km/s^2): " << fabs(Fgravity / spacecraft_mass) << endl;
                    outputfile << "Sun acceleration vector (km/s^2):" << endl;
                    for (size_t k = 0; k < 3; ++k)
                        outputfile << " " << F3rdbody[k] / spacecraft_mass;
                    outputfile << endl;
                    outputfile << endl;
#endif
                    if (generate_derivatives) {

                        doubleType spacecraft_distance_from_sun3 = spacecraft_distance_from_sun * spacecraft_distance_from_sun * spacecraft_distance_from_sun;
                        doubleType spacecraft_distance_from_sun5 = spacecraft_distance_from_sun3 * spacecraft_distance_from_sun * spacecraft_distance_from_sun;
                        doubleType three_muSun_over_spacecraft_distance_from_sun5 = 3.0 * mu_sun / spacecraft_distance_from_sun5;
                        doubleType muSun_over_spacecraft_distance_from_sun3 = mu_sun / spacecraft_distance_from_sun3;
                        std::vector<doubleType> State = spacecraft_position_relative_to_sun;
                        std::vector<doubleType> StateSun = sun_state_relative_to_central_body;
                        doubleType sun_distance_from_central_body3 = sun_distance_from_central_body * sun_distance_from_central_body * sun_distance_from_central_body;
                        doubleType sun_distance_from_central_body5 = sun_distance_from_central_body3 * sun_distance_from_central_body * sun_distance_from_central_body;
                        doubleType three_muSun_over_sun_distance_from_central_body5 = 3.0 * mu_sun / sun_distance_from_central_body5;
                        doubleType muSun_over_sun_distance_from_central_body3 = mu_sun / sun_distance_from_central_body3;

                        //****************************************
                        //
                        //State Propagation (fx) matrix calculation
                        //
                        //****************************************
                        // Equation 41 (terms 3 and 4)
                        fx(3, 0) += three_muSun_over_spacecraft_distance_from_sun5 * State[0] * State[0] - muSun_over_spacecraft_distance_from_sun3;
                        fx(3, 1) += three_muSun_over_spacecraft_distance_from_sun5 * State[0] * State[1];
                        fx(3, 2) += three_muSun_over_spacecraft_distance_from_sun5 * State[0] * State[2];
                        fx(4, 0) += three_muSun_over_spacecraft_distance_from_sun5 * State[1] * State[0];
                        fx(4, 1) += three_muSun_over_spacecraft_distance_from_sun5 * State[1] * State[1] - muSun_over_spacecraft_distance_from_sun3;
                        fx(4, 2) += three_muSun_over_spacecraft_distance_from_sun5 * State[1] * State[2];
                        fx(5, 0) += three_muSun_over_spacecraft_distance_from_sun5 * State[2] * State[0];
                        fx(5, 1) += three_muSun_over_spacecraft_distance_from_sun5 * State[2] * State[1];
                        fx(5, 2) += three_muSun_over_spacecraft_distance_from_sun5 * State[2] * State[2] - muSun_over_spacecraft_distance_from_sun3;

                        // TODO: add in frame drag derivatives
                        // dadr3B -- sensitivity of s/c acceleration to Sun's position (gravity component)
                        // Equation 62 (first term)
                        // daxdr3B
                        dadr3B[0] = -three_muSun_over_spacecraft_distance_from_sun5 * State[0] * State[0] + muSun_over_spacecraft_distance_from_sun3 + three_muSun_over_sun_distance_from_central_body5 * StateSun[0] * StateSun[0] - muSun_over_sun_distance_from_central_body3;
                        dadr3B[1] = -three_muSun_over_spacecraft_distance_from_sun5 * State[0] * State[1] + three_muSun_over_sun_distance_from_central_body5 * StateSun[0] * StateSun[1];
                        dadr3B[2] = -three_muSun_over_spacecraft_distance_from_sun5 * State[0] * State[2] + three_muSun_over_sun_distance_from_central_body5 * StateSun[0] * StateSun[2];

                        // daydr3B
                        dadr3B[3] = -three_muSun_over_spacecraft_distance_from_sun5 * State[1] * State[0] + three_muSun_over_sun_distance_from_central_body5 * StateSun[1] * StateSun[0];
                        dadr3B[4] = -three_muSun_over_spacecraft_distance_from_sun5 * State[1] * State[1] + muSun_over_spacecraft_distance_from_sun3 + three_muSun_over_sun_distance_from_central_body5 * StateSun[1] * StateSun[1] - muSun_over_sun_distance_from_central_body3;
                        dadr3B[5] = -three_muSun_over_spacecraft_distance_from_sun5 * State[1] * State[2] + three_muSun_over_sun_distance_from_central_body5 * StateSun[1] * StateSun[2];

                        // dazdr3B
                        dadr3B[6] = -three_muSun_over_spacecraft_distance_from_sun5 * State[2] * State[0] + three_muSun_over_sun_distance_from_central_body5 * StateSun[2] * StateSun[0];
                        dadr3B[7] = -three_muSun_over_spacecraft_distance_from_sun5 * State[2] * State[1] + three_muSun_over_sun_distance_from_central_body5 * StateSun[2] * StateSun[1];
                        dadr3B[8] = -three_muSun_over_spacecraft_distance_from_sun5 * State[2] * State[2] + muSun_over_spacecraft_distance_from_sun3 + three_muSun_over_sun_distance_from_central_body5 * StateSun[2] * StateSun[2] - muSun_over_sun_distance_from_central_body3;


                        // Equation 32 (second term for i = Sun)
                        for (size_t p = 0; p < 2; ++p) {
                            nbodyTOFterms_x[p] += dadr3B[0] * dsun_state_relative_to_central_body_dTOF(0, p) + dadr3B[1] * dsun_state_relative_to_central_body_dTOF(1, p) + dadr3B[2] * dsun_state_relative_to_central_body_dTOF(2, p);
                            nbodyTOFterms_y[p] += dadr3B[3] * dsun_state_relative_to_central_body_dTOF(0, p) + dadr3B[4] * dsun_state_relative_to_central_body_dTOF(1, p) + dadr3B[5] * dsun_state_relative_to_central_body_dTOF(2, p);
                            nbodyTOFterms_z[p] += dadr3B[6] * dsun_state_relative_to_central_body_dTOF(0, p) + dadr3B[7] * dsun_state_relative_to_central_body_dTOF(1, p) + dadr3B[8] * dsun_state_relative_to_central_body_dTOF(2, p);
                        }

                    }
                } // end Sun perturbations

                  //perturbations due to bodies in the perturbation list
                for (size_t b = 0; b < Universe.perturbation_menu.size(); ++b) {
                    static std::vector<doubleType> third_body_state_relative_to_central_body(12);
                    static std::vector<doubleType> spacecraft_position_relative_to_third_body(3);
                    doubleType spacecraft_distance_from_third_body;
                    doubleType third_body_distance_from_central_body;

                    static EMTG::math::Matrix <double> dthird_body_state_relative_to_central_bodydTOF(6, 2, 0.0);

                    //extract the third body state vector from the ephemeris
#ifdef AD_INSTRUMENTATION

                    static std::vector<double> third_body_state_relative_to_central_body_temp(12, 0.0);
                    static std::vector<double> dthird_body_state_relative_to_central_body_in_kmdTOF(12, 0.0);

                    Universe.bodies[Universe.perturbation_menu[b]].locate_body(epoch, third_body_state_relative_to_central_body.data(), true, options);

                    std::vector<size_t> timevars = epoch.getDerivativeIndicies();
                    for (size_t i = 0; i < 6; ++i)
                    {
                        // OK we have GOT to find a better way to do this
                        // Depending on whether this is a full or partial RK step, epoch may or may not have a derivative w.r.t. eta
                        third_body_state_relative_to_central_body[i].setDerivative(timevars[0], third_body_state_relative_to_central_body[i].getDerivative(timevars[0]) * depochdTOF[1]);
                        third_body_state_relative_to_central_body[i].setDerivative(timevars[1], third_body_state_relative_to_central_body[i].getDerivative(timevars[1]) * depochdTOF[1]);
                        //third_body_state_relative_to_central_body[i].setDerivative(timevars[2], third_body_state_relative_to_central_body[i].getDerivative(timevars[2]) * depochdTOF[0]);

                        //if (timevars.back() == 10)
                        if (timevars.back() == 2)
                        {
                            third_body_state_relative_to_central_body[i].setDerivative(timevars.back(), third_body_state_relative_to_central_body[i].getDerivative(timevars.back()) * depochdTOF[0]);
                        }
                    }
#else
                    Universe.bodies[Universe.perturbation_menu[b]].locate_body(epoch,
                        third_body_state_relative_to_central_body.data(),
                        generate_derivatives && options.derivative_type > 2,
                        options);
#endif

                    //compute the TOF derivatives
                    /*for (size_t i = 0; i < 3; ++i) {
                    dthird_body_state_relative_to_central_bodydTOF(i, 0) = third_body_state_relative_to_central_body[i + 3] * depochdTOF[0];
                    dthird_body_state_relative_to_central_bodydTOF(i, 1) = third_body_state_relative_to_central_body[i + 3] * depochdTOF[1];

                    dthird_body_state_relative_to_central_bodydTOF(i + 3, 0) = third_body_state_relative_to_central_body[i + 9] * depochdTOF[0];
                    dthird_body_state_relative_to_central_bodydTOF(i + 3, 1) = third_body_state_relative_to_central_body[i + 9] * depochdTOF[1];
                    }*/


                    for (size_t i = 0; i < 3; ++i)
                    {
                        dthird_body_state_relative_to_central_bodydTOF(i, 0) = (third_body_state_relative_to_central_body[i + 6] * depochdTOF[0])_GETVALUE;
                        dthird_body_state_relative_to_central_bodydTOF(i, 1) = (third_body_state_relative_to_central_body[i + 6] * depochdTOF[1])_GETVALUE;

                        dthird_body_state_relative_to_central_bodydTOF(i + 3, 0) = (third_body_state_relative_to_central_body[i + 9] * depochdTOF[0])_GETVALUE;
                        dthird_body_state_relative_to_central_bodydTOF(i + 3, 1) = (third_body_state_relative_to_central_body[i + 9] * depochdTOF[1])_GETVALUE;
                    }

                    // compute the distance from the central body to the perturbing body
                    third_body_distance_from_central_body = sqrt(third_body_state_relative_to_central_body[0] * third_body_state_relative_to_central_body[0] +
                        third_body_state_relative_to_central_body[1] * third_body_state_relative_to_central_body[1] +
                        third_body_state_relative_to_central_body[2] * third_body_state_relative_to_central_body[2]);

                    //form the vector from the third body to the spacecraft and its TOF derivatives in the central body reference frame
                    //See Figure 1 in paper
                    spacecraft_position_relative_to_third_body[0] = spacecraft_state_relative_to_central_body[0] - third_body_state_relative_to_central_body[0];
                    spacecraft_position_relative_to_third_body[1] = spacecraft_state_relative_to_central_body[1] - third_body_state_relative_to_central_body[1];
                    spacecraft_position_relative_to_third_body[2] = spacecraft_state_relative_to_central_body[2] - third_body_state_relative_to_central_body[2];

                    spacecraft_distance_from_third_body = sqrt(spacecraft_position_relative_to_third_body[0] * spacecraft_position_relative_to_third_body[0] +
                        spacecraft_position_relative_to_third_body[1] * spacecraft_position_relative_to_third_body[1] +
                        spacecraft_position_relative_to_third_body[2] * spacecraft_position_relative_to_third_body[2]);

                    //to avoid singularities and to avoid screwing up the flyby model, we will only include a third body perturbations when we are not too close to the body
                    bool include_body = true;
                    bool departing_first_body = (step < (int)options.num_timesteps / 2) &&
                        (Universe.bodies[Universe.perturbation_menu[b]].body_code == options.Journeys[j].sequence[p]
                            || Universe.bodies[Universe.perturbation_menu[b]].spice_ID * 100 + 99 == Universe.bodies[options.Journeys[j].sequence[p] - 1].spice_ID);
                    bool arriving_last_body = (step >= (int)options.num_timesteps / 2 || h < 0.0) &&
                        (Universe.bodies[Universe.perturbation_menu[b]].body_code == options.Journeys[j].sequence[p + 1]
                            || Universe.bodies[Universe.perturbation_menu[b]].spice_ID * 100 + 99 == Universe.bodies[options.Journeys[j].sequence[p + 1] - 1].spice_ID);

                    if (spacecraft_distance_from_third_body < Universe.bodies[Universe.perturbation_menu[b]].r_SOI
                        && (departing_first_body || arriving_last_body))
                        include_body = false;

                    //if (spacecraft_distance_from_third_body < Universe.bodies[Universe.perturbation_menu[b]].r_SOI)
                    //|| spacecraft_distance_from_third_body < Universe.bodies[Universe.perturbation_menu[b]].exclusion_radius)
                    //include_body = false;

                    if (include_body) {
                        Agravity = -(Universe.bodies[Universe.perturbation_menu[b]].mu) / (spacecraft_distance_from_third_body*spacecraft_distance_from_third_body);
                        doubleType r_body_to_central_body = math::norm(third_body_state_relative_to_central_body.data(), 3);
                        doubleType a_body_on_central_body = -(Universe.bodies[Universe.perturbation_menu[b]].mu) / (r_body_to_central_body*r_body_to_central_body);

                        //********************************************************
                        //
                        //Force vector modification: Perturbing gravitational body
                        //
                        //********************************************************
                        doubleType a3rdBody[3];
                        for (size_t k = 0; k < 3; ++k) {
                            a3rdBody[k] = Agravity * spacecraft_position_relative_to_third_body[k] / spacecraft_distance_from_third_body + a_body_on_central_body * third_body_state_relative_to_central_body[k] / r_body_to_central_body;
                            acceleration_vector[k] += a3rdBody[k];
                        }

#ifdef REPORT_FORCE_MODEL
                        outputfile << "-------------" << Universe.bodies[Universe.perturbation_menu[b]].name << " (" << Universe.bodies[Universe.perturbation_menu[b]].spice_ID << ")------" << endl;
                        outputfile << Universe.bodies[Universe.perturbation_menu[b]].name << " mu (km^2/s^3): " << Universe.bodies[Universe.perturbation_menu[b]].mu << endl;
                        outputfile << Universe.bodies[Universe.perturbation_menu[b]].name << " distance from central body (km): " << sqrt(third_body_state_relative_to_central_body[0] * third_body_state_relative_to_central_body[0]
                            + third_body_state_relative_to_central_body[1] * third_body_state_relative_to_central_body[1]
                            + third_body_state_relative_to_central_body[2] * third_body_state_relative_to_central_body[2]) << endl;
                        outputfile << Universe.bodies[Universe.perturbation_menu[b]].name << " distance from central body (LU): " << sqrt(third_body_state_relative_to_central_body[0] * third_body_state_relative_to_central_body[0]
                            + third_body_state_relative_to_central_body[1] * third_body_state_relative_to_central_body[1]
                            + third_body_state_relative_to_central_body[2] * third_body_state_relative_to_central_body[2]) / LU << endl;
                        outputfile << Universe.bodies[Universe.perturbation_menu[b]].name << " state vector relative to central body (km): " << endl;
                        for (size_t k = 0; k < 3; ++k)
                            outputfile << " " << third_body_state_relative_to_central_body[k]_GETVALUE;
                        outputfile << endl;
                        outputfile << "Spacecraft distance from " << Universe.bodies[Universe.perturbation_menu[b]].name << " (km): " << spacecraft_distance_from_third_body << endl;
                        outputfile << "Spacecraft distance from " << Universe.bodies[Universe.perturbation_menu[b]].name << " (LU): " << spacecraft_distance_from_third_body / Universe.LU << endl;
                        outputfile << Universe.bodies[Universe.perturbation_menu[b]].name << " force magnitude (kN): " << fabs(Fgravity) << endl;
                        outputfile << Universe.bodies[Universe.perturbation_menu[b]].name << " force vector (kN):" << endl;
                        for (size_t k = 0; k < 3; ++k)
                            outputfile << " " << F3rdBody[k];
                        outputfile << endl;
                        outputfile << Universe.bodies[Universe.perturbation_menu[b]].name << " acceleration magnitude (km/s^2): " << fabs(Fgravity / spacecraft_mass) << endl;
                        outputfile << Universe.bodies[Universe.perturbation_menu[b]].name << " acceleration vector (km/s^2):" << endl;
                        for (size_t k = 0; k < 3; ++k)
                            outputfile << " " << F3rdBody[k] / spacecraft_mass;
                        outputfile << endl;
                        outputfile << endl;
#endif

                        if (generate_derivatives) {

                            doubleType mu_3B = Universe.bodies[Universe.perturbation_menu[b]].mu;
                            doubleType spacecraft_distance_from_third_body3 = spacecraft_distance_from_third_body * spacecraft_distance_from_third_body * spacecraft_distance_from_third_body;
                            doubleType spacecraft_distance_from_third_body5 = spacecraft_distance_from_third_body3 * spacecraft_distance_from_third_body * spacecraft_distance_from_third_body;
                            doubleType three_mu3B_over_spacecraft_distance_from_third_body5 = 3.0 * mu_3B / spacecraft_distance_from_third_body5;
                            doubleType mu3B_over_spacecraft_distance_from_third_body3 = mu_3B / spacecraft_distance_from_third_body3;
                            std::vector<doubleType> State = spacecraft_position_relative_to_third_body;
                            std::vector<doubleType> State3B = third_body_state_relative_to_central_body;
                            doubleType third_body_distance_from_central_body3 = third_body_distance_from_central_body * third_body_distance_from_central_body * third_body_distance_from_central_body;
                            doubleType third_body_distance_from_central_body5 = third_body_distance_from_central_body3 * third_body_distance_from_central_body * third_body_distance_from_central_body;
                            doubleType three_mu3B_over_third_body_distance_from_central_body5 = 3.0 * mu_3B / third_body_distance_from_central_body5;
                            doubleType mu3B_over_third_body_distance_from_central_body3 = mu_3B / third_body_distance_from_central_body3;


                            //****************************************
                            //
                            //State Propagation (fx) matrix calculation
                            //
                            //****************************************
                            //Equation 41 (terms 3 and 4)
                            fx(3, 0) += three_mu3B_over_spacecraft_distance_from_third_body5 * State[0] * State[0] - mu3B_over_spacecraft_distance_from_third_body3;
                            fx(3, 1) += three_mu3B_over_spacecraft_distance_from_third_body5 * State[0] * State[1];
                            fx(3, 2) += three_mu3B_over_spacecraft_distance_from_third_body5 * State[0] * State[2];
                            fx(4, 0) += three_mu3B_over_spacecraft_distance_from_third_body5 * State[1] * State[0];
                            fx(4, 1) += three_mu3B_over_spacecraft_distance_from_third_body5 * State[1] * State[1] - mu3B_over_spacecraft_distance_from_third_body3;
                            fx(4, 2) += three_mu3B_over_spacecraft_distance_from_third_body5 * State[1] * State[2];
                            fx(5, 0) += three_mu3B_over_spacecraft_distance_from_third_body5 * State[2] * State[0];
                            fx(5, 1) += three_mu3B_over_spacecraft_distance_from_third_body5 * State[2] * State[1];
                            fx(5, 2) += three_mu3B_over_spacecraft_distance_from_third_body5 * State[2] * State[2] - mu3B_over_spacecraft_distance_from_third_body3;


                            // TODO: add frame drag derivatives here
                            //dadr3B -- sensitivity of s/c acceleration to third body position
                            //Equation 62
                            //daxdr3B
                            dadr3B[0] = -three_mu3B_over_spacecraft_distance_from_third_body5 * State[0] * State[0] + mu3B_over_spacecraft_distance_from_third_body3 + three_mu3B_over_third_body_distance_from_central_body5 * State3B[0] * State3B[0] - mu3B_over_third_body_distance_from_central_body3;
                            dadr3B[1] = -three_mu3B_over_spacecraft_distance_from_third_body5 * State[0] * State[1]                                                  + three_mu3B_over_third_body_distance_from_central_body5 * State3B[0] * State3B[1];
                            dadr3B[2] = -three_mu3B_over_spacecraft_distance_from_third_body5 * State[0] * State[2]                                                  + three_mu3B_over_third_body_distance_from_central_body5 * State3B[0] * State3B[2];

                            //daydr3B
                            dadr3B[3] = -three_mu3B_over_spacecraft_distance_from_third_body5 * State[1] * State[0]                                                  + three_mu3B_over_third_body_distance_from_central_body5 * State3B[1] * State3B[0];
                            dadr3B[4] = -three_mu3B_over_spacecraft_distance_from_third_body5 * State[1] * State[1] + mu3B_over_spacecraft_distance_from_third_body3 + three_mu3B_over_third_body_distance_from_central_body5 * State3B[1] * State3B[1] - mu3B_over_third_body_distance_from_central_body3;
                            dadr3B[5] = -three_mu3B_over_spacecraft_distance_from_third_body5 * State[1] * State[2]                                                  + three_mu3B_over_third_body_distance_from_central_body5 * State3B[1] * State3B[2];

                            //dazdr3B
                            dadr3B[6] = -three_mu3B_over_spacecraft_distance_from_third_body5 * State[2] * State[0]                                                  + three_mu3B_over_third_body_distance_from_central_body5 * State3B[2] * State3B[0];
                            dadr3B[7] = -three_mu3B_over_spacecraft_distance_from_third_body5 * State[2] * State[1]                                                  + three_mu3B_over_third_body_distance_from_central_body5 * State3B[2] * State3B[1];
                            dadr3B[8] = -three_mu3B_over_spacecraft_distance_from_third_body5 * State[2] * State[2] + mu3B_over_spacecraft_distance_from_third_body3 + three_mu3B_over_third_body_distance_from_central_body5 * State3B[2] * State3B[2] - mu3B_over_third_body_distance_from_central_body3;


                            //Equation 32 (second term - summation over i bodies, not including Sun)
                            for (size_t p = 0; p < 2; ++p)
                            {
                                nbodyTOFterms_x[p] += dadr3B[0] * dthird_body_state_relative_to_central_bodydTOF(0, p) +
                                    dadr3B[1] * dthird_body_state_relative_to_central_bodydTOF(1, p) +
                                    dadr3B[2] * dthird_body_state_relative_to_central_bodydTOF(2, p);
                                nbodyTOFterms_y[p] += dadr3B[3] * dthird_body_state_relative_to_central_bodydTOF(0, p) +
                                    dadr3B[4] * dthird_body_state_relative_to_central_bodydTOF(1, p) +
                                    dadr3B[5] * dthird_body_state_relative_to_central_bodydTOF(2, p);
                                nbodyTOFterms_z[p] += dadr3B[6] * dthird_body_state_relative_to_central_bodydTOF(0, p) +
                                    dadr3B[7] * dthird_body_state_relative_to_central_bodydTOF(1, p) +
                                    dadr3B[8] * dthird_body_state_relative_to_central_bodydTOF(2, p);
                            }

                        }
                    }

                } //end perturbation body loop
            }//end third body perturbation code

             //***********************************************************
             //
             //Current and previous phase TOF derivative calculation for this DOPRI 7(8) stage
             //
             //***********************************************************
            if (generate_derivatives) {
                static std::vector<doubleType> dTdTOF(2, 0.0);
                static std::vector<doubleType> dmdotdTOF(2, 0.0);

                //Equation 37 (second term)
                double dPdt_launch = 0.0; // only non-zero if we are modelling power model decay
                double dt_launchdTOF = 0.0; // only non-zero if we are modelling power model decay


                static std::vector <doubleType> dr_sundTOF(2, 0.0);
                static std::vector <doubleType> dPdTOF(2, 0.0);


                //we need to compute phase TOF derivatives for each previous phase and the current one
                if (!isHeliocentric) {
                    //Equation 65 in the paper
                    //If we are NOT orbiting the Sun, then the second term in this equation must be calculated
                    for (size_t p = 0; p < 2; ++p) {
                        dr_sundTOF[p] = drdx * dspacecraft_state_relative_to_central_bodydTOF(0, p) +
                            drdy * dspacecraft_state_relative_to_central_bodydTOF(1, p) +
                            drdz * dspacecraft_state_relative_to_central_bodydTOF(2, p) +
                            (drdx)* dsun_state_relative_to_central_body_dTOF(0, p) +
                            (drdy)* dsun_state_relative_to_central_body_dTOF(1, p) +
                            (drdz)* dsun_state_relative_to_central_body_dTOF(2, p);

                        dPdTOF[p] = ((dPdr_sun)* dr_sundTOF[p] + dPdt_launch * dt_launchdTOF);
                        dTdTOF[p] = (dTdP)* dPdTOF[p];
                        dmdotdTOF[p] = (dmdotdP)* dPdTOF[p];
                    }


                    //Sensitivity of s/c acceleration to changes in Sun's position (power component only - gravity taken care of already in perturbation loop)
                    //Equation 62, third term 
                    static std::vector <doubleType> dadsr_sun(9, 0.0);
                    doubleType D_dTdP_dPdr_over_msc = (dTdP) * (dPdr_sun) / spacecraft_mass;
                    doubleType ux_D_dTdP_dPdr_over_msc = control[0] * D_dTdP_dPdr_over_msc;
                    doubleType uy_D_dTdP_dPdr_over_msc = control[1] * D_dTdP_dPdr_over_msc;
                    doubleType uz_D_dTdP_dPdr_over_msc = control[2] * D_dTdP_dPdr_over_msc;

                    dadsr_sun[0] = -ux_D_dTdP_dPdr_over_msc * (drdx);
                    dadsr_sun[1] = -ux_D_dTdP_dPdr_over_msc * (drdy);
                    dadsr_sun[2] = -ux_D_dTdP_dPdr_over_msc * (drdz);
                    dadsr_sun[3] = -uy_D_dTdP_dPdr_over_msc * (drdx);
                    dadsr_sun[4] = -uy_D_dTdP_dPdr_over_msc * (drdy);
                    dadsr_sun[5] = -uy_D_dTdP_dPdr_over_msc * (drdz);
                    dadsr_sun[6] = -uz_D_dTdP_dPdr_over_msc * (drdx);
                    dadsr_sun[7] = -uz_D_dTdP_dPdr_over_msc * (drdy);
                    dadsr_sun[8] = -uz_D_dTdP_dPdr_over_msc * (drdz);

                    //Construct the partial derivatives w.r.t. TOF for the state gradient
                    for (size_t p = 0; p < 2; ++p) {
                        //dvdTOF
                        //Equation 31 (needed for Eq. 29)
                        dfdTOF(0, p) = dspacecraft_state_relative_to_central_bodydTOF(3, p);
                        dfdTOF(1, p) = dspacecraft_state_relative_to_central_bodydTOF(4, p);
                        dfdTOF(2, p) = dspacecraft_state_relative_to_central_bodydTOF(5, p);

                        //dadTOF
                        //Equation 32 (whole thing)
                        dfdTOF(3, p) = (fx(3, 0)*dspacecraft_state_relative_to_central_bodydTOF(0, p) +
                            fx(3, 1)*dspacecraft_state_relative_to_central_bodydTOF(1, p) +
                            fx(3, 2)*dspacecraft_state_relative_to_central_bodydTOF(2, p) +
                            nbodyTOFterms_x[p] +
                            fx(3, 6)*dspacecraft_state_relative_to_central_bodydTOF(6, p) +
                            dadsr_sun[0] * dsun_state_relative_to_central_body_dTOF(0, p) +
                            dadsr_sun[1] * dsun_state_relative_to_central_body_dTOF(1, p) +
                            dadsr_sun[2] * dsun_state_relative_to_central_body_dTOF(2, p))_GETVALUE;

                        dfdTOF(4, p) = (fx(4, 0)*dspacecraft_state_relative_to_central_bodydTOF(0, p) +
                            fx(4, 1)*dspacecraft_state_relative_to_central_bodydTOF(1, p) +
                            fx(4, 2)*dspacecraft_state_relative_to_central_bodydTOF(2, p) +
                            nbodyTOFterms_y[p] +
                            fx(4, 6)*dspacecraft_state_relative_to_central_bodydTOF(6, p) +
                            dadsr_sun[3] * dsun_state_relative_to_central_body_dTOF(0, p) +
                            dadsr_sun[4] * dsun_state_relative_to_central_body_dTOF(1, p) +
                            dadsr_sun[5] * dsun_state_relative_to_central_body_dTOF(2, p))_GETVALUE;

                        dfdTOF(5, p) = (fx(5, 0)*dspacecraft_state_relative_to_central_bodydTOF(0, p) +
                            fx(5, 1)*dspacecraft_state_relative_to_central_bodydTOF(1, p) +
                            fx(5, 2)*dspacecraft_state_relative_to_central_bodydTOF(2, p) +
                            nbodyTOFterms_z[p] +
                            fx(5, 6)*dspacecraft_state_relative_to_central_bodydTOF(6, p) +
                            dadsr_sun[6] * dsun_state_relative_to_central_body_dTOF(0, p) +
                            dadsr_sun[7] * dsun_state_relative_to_central_body_dTOF(1, p) +
                            dadsr_sun[8] * dsun_state_relative_to_central_body_dTOF(2, p))_GETVALUE;

                        dfdTOF(6, p) = (-control_norm * dmdotdTOF[p])_GETVALUE;
                    }
                }
                else {
                    //Equation 65, if we are orbiting the Sun, the second term is zero
                    for (size_t p = 0; p < 2; ++p) {
                        dr_sundTOF[p] = drdx * dspacecraft_state_relative_to_central_bodydTOF(0, p) +
                            drdy * dspacecraft_state_relative_to_central_bodydTOF(1, p) +
                            drdz * dspacecraft_state_relative_to_central_bodydTOF(2, p);

                        dPdTOF[p] = ((dPdr_sun)* dr_sundTOF[p] + dPdt_launch * dt_launchdTOF);
                        dTdTOF[p] = (dTdP)* dPdTOF[p];
                        dmdotdTOF[p] = (dmdotdP)* dPdTOF[p];
                    }

                    for (size_t p = 0; p < 2; ++p) {

                        //dvdTOF
                        //Equation 59
                        dfdTOF(0, p) = dspacecraft_state_relative_to_central_bodydTOF(3, p);
                        dfdTOF(1, p) = dspacecraft_state_relative_to_central_bodydTOF(4, p);
                        dfdTOF(2, p) = dspacecraft_state_relative_to_central_bodydTOF(5, p);

                        //dadTOF
                        //Equation 60 (whole thing)
                        dfdTOF(3, p) = (fx(3, 0)*dspacecraft_state_relative_to_central_bodydTOF(0, p) +
                            fx(3, 1)*dspacecraft_state_relative_to_central_bodydTOF(1, p) +
                            fx(3, 2)*dspacecraft_state_relative_to_central_bodydTOF(2, p) +
                            nbodyTOFterms_x[p] +
                            fx(3, 6)*dspacecraft_state_relative_to_central_bodydTOF(6, p))_GETVALUE;

                        dfdTOF(4, p) = (fx(4, 0)*dspacecraft_state_relative_to_central_bodydTOF(0, p) +
                            fx(4, 1)*dspacecraft_state_relative_to_central_bodydTOF(1, p) +
                            fx(4, 2)*dspacecraft_state_relative_to_central_bodydTOF(2, p) +
                            nbodyTOFterms_y[p] +
                            fx(4, 6)*dspacecraft_state_relative_to_central_bodydTOF(6, p))_GETVALUE;

                        dfdTOF(5, p) = (fx(5, 0)*dspacecraft_state_relative_to_central_bodydTOF(0, p) +
                            fx(5, 1)*dspacecraft_state_relative_to_central_bodydTOF(1, p) +
                            fx(5, 2)*dspacecraft_state_relative_to_central_bodydTOF(2, p) +
                            nbodyTOFterms_z[p] +
                            fx(5, 6)*dspacecraft_state_relative_to_central_bodydTOF(6, p))_GETVALUE;

                        dfdTOF(6, p) = (-control_norm * dmdotdTOF[p])_GETVALUE;

                    }
                }
            }// end TOF derivative code

#ifdef REPORT_FORCE_MODEL
            outputfile.close();
#endif
            // double acceleration_magnitude = (math::norm(acceleration_vector.data(), 3)) _GETVALUE;
            return 0;
        }// end FBLT_acceleration_model
    }
} //close namespace