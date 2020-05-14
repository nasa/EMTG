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

//differential equations of motion for EMTGv9
//Donald Ellison August 1st 2015

#include <cmath>
#include <vector>

#include "missionoptions.h"
#include "EMTG_Matrix.h"
#include "EMTG_math.h"
#include "doubleType.h"
#include "FBLT_acceleration_model.h"
#include "Spacecraft.h"
#include "universe.h"
#include "FBLT_EOM.h"

namespace EMTG {
    namespace Astrodynamics {
        namespace EOM
        {
            //cartesian coordinate equations of motion for a spacecraft with a thrust term

            void EOM_inertial_continuous_thrust(std::vector <doubleType> & spacecraft_state,
                EMTG::math::Matrix <double> & dxdTOF,
                const doubleType & epoch_step_left,
                std::vector <double> & depoch_step_leftdTOF,
                const double & c,
                doubleType & h,
                const double & dhdTOF,
                const doubleType & launch_epoch,
                const std::vector <doubleType> & u,
                const doubleType& u_command,
                const double& DutyCycle,
                std::vector <doubleType> & f, // EOM gradient vector
                EMTG::math::Matrix <double> & dfdTOF,
                doubleType & thrust,
                doubleType & mdot,
                doubleType & Isp,
                doubleType & power,
                doubleType & active_power,
                int & number_of_active_engines,
                int & STMrows,
                int & STMcolumns,
                const missionoptions& options,
                EMTG::Astrodynamics::universe& Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                const bool& GenerateDerivatives,
                const int& j,
                const int& p,
                const int& step,
                const double& mu,
                const double& LU,
                const double& TU,
                const double& MU)
            {
                //assume that all inputs are normalized!

                static std::vector <doubleType> AccelerationVector(3, 0.0);
                static std::vector <doubleType> AccelerationVector_mks(3, 0.0);
                static std::vector<doubleType> spacecraft_state_mks(7, 0.0);
                static math::Matrix<double> dxdTOF_mks(7, 2, 0.0);
                static math::Matrix<double> dfdTOF_mks(7, 2, 0.0);

                //determine the forces acting on the spacecraft 
                //calculate the state propagation matrix fx
                //calculate the TOF derivatives of the EOM
                static EMTG::math::Matrix<doubleType> fx(STMrows, STMcolumns, 0.0);
                static EMTG::math::Matrix<doubleType> fx_mks(STMrows, STMcolumns, 0.0);

                //Force model is specified in MKS
                for (size_t i = 0; i < 3; ++i)
                {
                    spacecraft_state_mks[i] = spacecraft_state[i] * LU;
                    spacecraft_state_mks[i + 3] = spacecraft_state[i + 3] * LU / TU;

                    dxdTOF_mks(i, 0) = dxdTOF(i, 0) * LU;
                    dxdTOF_mks(i, 1) = dxdTOF(i, 1) * LU;

                    dxdTOF_mks(i + 3, 0) = dxdTOF(i + 3, 0) * LU / (TU);
                    dxdTOF_mks(i + 3, 1) = dxdTOF(i + 3, 1) * LU / (TU);
                }
                spacecraft_state_mks[6] = spacecraft_state[6] * MU;

                dxdTOF_mks(6, 0) = dxdTOF(6, 0) * MU;
                dxdTOF_mks(6, 1) = dxdTOF(6, 1) * MU;

                doubleType epoch_step_left_mks = epoch_step_left * TU;
                static std::vector<double> depoch_step_leftdTOF_mks(2, 0.0);
                for (size_t k = 0; k < 2; ++k)
                {
                    depoch_step_leftdTOF_mks[k] = depoch_step_leftdTOF[k] * TU;
                }
                doubleType h_mks = h * TU;
                double dhdTOF_mks = dhdTOF * TU;
                doubleType launch_epoch_mks = launch_epoch * TU;

                // fx and dfdTOF do NOT need to be scaled as they are OUTPUTS of FBLT_acceleration_model only

                EMTG::Astrodynamics::FBLT_acceleration_model(options,
                    Universe,
                    mySpacecraft,
                    spacecraft_state_mks,
                    dxdTOF_mks,
                    epoch_step_left_mks,
                    depoch_step_leftdTOF_mks,
                    c,
                    h_mks,
                    dhdTOF_mks,
                    launch_epoch_mks,
                    u,
                    u_command,
                    DutyCycle,
                    dfdTOF_mks,
                    thrust,
                    mdot,
                    Isp,
                    power,
                    active_power,
                    number_of_active_engines,
                    fx_mks,
                    AccelerationVector_mks,
                    GenerateDerivatives,
                    j,
                    p,
                    step);


                //now scale AccelerationVector, fx, and dfdTOF back to canonical units

                //AccelerationVector
                for (size_t k = 0; k < 3; ++k)
                {
                    AccelerationVector[k] = AccelerationVector_mks[k] * TU * TU / LU;
                }


                //fx
                fx = fx_mks;

                //(scale entries that are not unity or zero)
                for (size_t i = 3; i < 6; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                        fx(i, j) *= TU * TU;

                    for (size_t j = 7; j < 10; ++j)
                        fx(i, j) *= TU * TU / LU;
                    fx(i, 6) *= MU * TU * TU / LU;
                }

                for (size_t j = 0; j < 6; ++j)
                    fx(6, j) *= LU / MU;
                for (size_t j = 0; j < 3; ++j)
                    fx(6, j) *= TU;
                for (size_t j = 7; j < 10; ++j)
                    fx(6, j) *= TU / MU;


                //dfdTOF
                for (size_t i = 0; i < 3; ++i)
                {
                    dfdTOF(i, 0) = dfdTOF_mks(i, 0) * TU / LU;
                    dfdTOF(i, 1) = dfdTOF_mks(i, 1) * TU / LU;
                    dfdTOF(i + 3, 0) = dfdTOF_mks(i + 3, 0) * TU * TU / LU;
                    dfdTOF(i + 3, 1) = dfdTOF_mks(i + 3, 1) * TU * TU / LU;
                }
                dfdTOF(6, 0) = dfdTOF_mks(6, 0) * TU / MU;
                dfdTOF(6, 1) = dfdTOF_mks(6, 1) * TU / MU;



                //some checks to avoid singularities in the EOM
                doubleType r = sqrt(spacecraft_state[0] * spacecraft_state[0] + spacecraft_state[1] * spacecraft_state[1] + spacecraft_state[2] * spacecraft_state[2]); //magnitude of position vector
                r = fabs(r) < EMTG::math::SMALL ? EMTG::math::sgn(r) * EMTG::math::SMALL : r;
                spacecraft_state[6] = spacecraft_state[6] < EMTG::math::SMALL ? EMTG::math::SMALL : spacecraft_state[6];

                //calculate the EOM
                //n-body gravitational forces and thruster effects are 
                //all included in "AccelerationVector"

                //xdot
                f[0] = spacecraft_state[3];

                //ydot
                f[1] = spacecraft_state[4];

                //zdot
                f[2] = spacecraft_state[5];

                //xddot
                f[3] = -Universe.mu*spacecraft_state[0] / (r*r*r) + AccelerationVector[0];

                //yddot
                f[4] = -Universe.mu *spacecraft_state[1] / (r*r*r) + AccelerationVector[1];

                //zddot
                f[5] = -Universe.mu*spacecraft_state[2] / (r*r*r) + AccelerationVector[2];

                //mdot
                //Equation 3 from paper
                f[6] = -sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) * mdot * (TU / MU);

#ifdef REPORT_EOM
                std::vector<std::string> StateDescriptions;
                StateDescriptions.push_back("x");
                StateDescriptions.push_back("y");
                StateDescriptions.push_back("z");
                StateDescriptions.push_back("xdot");
                StateDescriptions.push_back("ydot");
                StateDescriptions.push_back("zdot");
                StateDescriptions.push_back("m");
                StateDescriptions.push_back("u_x");
                StateDescriptions.push_back("u_y");
                StateDescriptions.push_back("u_z");
                //check fx
                math::Matrix<double> Fx_alg(10, 10);
                math::Matrix<double> Fx_analytical(10, 10);
                math::Matrix<double> Fxcheck(10, 10);
                for (size_t i = 0; i < 10; ++i)
                {
                    //state
                    for (size_t j = 0; j < 7; ++j)
                    {
                        Fx_alg(i, j) = f[i].getDerivative(j);
                        Fx_analytical(i, j) = fx(i, j) _GETVALUE;
                    }

                    //control
                    for (size_t j = 7; j < 10; ++j)
                    {
                        Fx_alg(i, j) = f[i].getDerivative(j + 4);
                        Fx_analytical(i, j) = fx(i, j) _GETVALUE;
                    }
                }
                Fx_analytical.print_to_file("fx_analytical.txt", StateDescriptions, StateDescriptions);
                Fx_alg.print_to_file("Fx_alg.txt", StateDescriptions, StateDescriptions);
                Fx_analytical.element_divide(Fx_alg).print_to_file("Fx_analytical_over_algorithmic.txt", StateDescriptions, StateDescriptions);


                //check dxdTOF
                math::Matrix<double> dxdTOF_algorithmic(7, 2);
                for (size_t i = 0; i < 7; ++i)
                {
                    dxdTOF_algorithmic(i, 0) = spacecraft_state[i].getDerivative(9);
                    dxdTOF_algorithmic(i, 1) = spacecraft_state[i].getDerivative(8);
                }
                dxdTOF.print_to_file("dxdTOF_analytical.txt");
                dxdTOF_algorithmic.print_to_file("dxdTOF_state_alg.txt");
                dxdTOF.element_divide(dxdTOF_algorithmic).print_to_file("dxdTOF_analytical_over_algorithmic.txt");

                //check dfdTOF
                math::Matrix<double> dfdTOF_algorithmic(7, 2);
                for (size_t i = 0; i < 7; ++i)
                {
                    dfdTOF_algorithmic(i, 0) = f[i].getDerivative(9);
                    dfdTOF_algorithmic(i, 1) = f[i].getDerivative(8);
                }
                dfdTOF.print_to_file("dfdTOF_analytical.txt");
                dfdTOF_algorithmic.print_to_file("dfdTOF_algorithmic.txt");
                dfdTOF.element_divide(dfdTOF_algorithmic).print_to_file("dfdTOF_analytical_over_algorithmic.txt");

#endif
                size_t n_states_and_controls = 7;
                /*
                if (options.Journeys[j].phase_type_input[p] == EMTG::PhaseType::FBLT || options.mission_type == EMTG::PhaseType::FBLT)
                {
                n_states_and_controls = 10;
                //control vector time rates of change (these are decision variables and are not impacted by dynamics and are constant across an FBLT time step)
                //these are just a place holders in the STM and will never have a diffeqs associated with them
                //Equation 18 from paper
                f[7] = 0.0000000000000000;
                f[8] = 0.0000000000000000;
                f[9] = 0.0000000000000000;

                //TOF time rates of change (this is also a decision variable; it is not impacted by dynamics and is constant across an FBLT time step/phase)
                //this is just a place holder in the STM and will never have a diffeq associated with it
                //f[10] = 0.0000000000000000;
                }
                */

                //*******************
                //
                //Phi dot calculation
                //
                //*******************

                if (GenerateDerivatives)
                {
                    //Form the STM
                    //The STM is comprised of state variables
                    static EMTG::math::Matrix<doubleType> STM(STMrows, STMcolumns, 0.0);
                    static EMTG::math::Matrix<doubleType> STMdot(STMrows, STMcolumns, 0.0);

                    size_t xcount = n_states_and_controls;
                    for (size_t i = 0; i < STMrows; ++i)
                    {
                        for (size_t j = 0; j < STMcolumns; ++j)
                        {
                            STM(i, j) = spacecraft_state[xcount];
                            ++xcount;
                        }
                    }

#ifdef REPORT_EOM
                    //check pre-calc STM
                    math::Matrix<double> STM_precalc_analytical(7, 7);
                    math::Matrix<double> STM_precalc_algorithmic(7, 7);
                    for (size_t i = 0; i < 7; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                        {
                            STM_precalc_algorithmic(i, j) = spacecraft_state[i].getDerivative(j);
                            STM_precalc_analytical(i, j) = STM(i, j) _GETVALUE;
                        }
                    }
                    STM_precalc_analytical.print_to_file("STM_precalc_analytical.txt");
                    STM_precalc_algorithmic.print_to_file("STM_precalc_algorithmic.txt");
                    STM_precalc_analytical.element_divide(STM_precalc_algorithmic).print_to_file("STM_precalc_analytical_over_algorithmic.txt");
#endif

                    //differential equation for STM creation
                    //Equation 14 from paper
                    STMdot = fx * STM;

                    //Package the Phidot entries into the gradient vector f behind the spacecraft's state variables as well as the control gradients which are always zero
                    //(thrust is constant over an FBLT segment)
                    //Equation 15 from paper

                    xcount = n_states_and_controls; //STM entry diffeq's are placed behind the augmented state vector
                    for (size_t i = 0; i < STMrows; ++i)
                    {
                        for (size_t j = 0; j < STMcolumns; ++j)
                        {
                            f[xcount++] = STMdot(i, j);
                        }
                    }
                }//end derivatives

                else
                {

#ifndef AD_INSTRUMENTATION
                    size_t xcount = n_states_and_controls; //STM entry diffeq's are placed behind the augmented state vector
                    for (size_t i = 0; i < STMrows; ++i)
                    {
                        for (size_t j = 0; j < STMcolumns; ++j)
                        {
                            f[xcount++] = 0.0;
                        }
                    }
#endif

                }

            }//end EOM function


        } //end EOM namespace
    } //end Astrodynamics namespace
} //end EMTG namespace