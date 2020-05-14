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

// Templated version of the Dormand-Prince (DOPRI) 8th order 13 step algorithm
// DOPRI constants are from Numerical Recipies
// Donald Ellison 11/15/2016


#ifndef RUNGEKUTTA8_H
#define RUNGEKUTTA8_H


#include "FBLT_EOM.h"
#include "doubleType.h"
#include "Integrand.h"
#include "IntegrationScheme.h"

namespace EMTG {
    namespace Integration
    {
        class RungeKutta8 : public IntegrationScheme
        {

        public:

            // constructors
            RungeKutta8();
            RungeKutta8(Integrand * integrand_in, const size_t & ns_in, const size_t & num_prop_var_deriv_states);

            // destructor
            ~RungeKutta8();

            // methods
            virtual void step(const math::Matrix<doubleType> & state_left,
                              const math::Matrix<double> & dstate_leftdProp_vars,
                              math::Matrix<doubleType> & state_right,
                              math::Matrix<double> & dstate_rightdProp_vars,
                              const doubleType & step_size,
                              const double& dstep_sizedProp_var,
                              const bool & needSTM);

            virtual void step_w_control(const math::Matrix<doubleType> & state_left,
                                        const math::Matrix<double> & dstate_leftdProp_vars,
                                        math::Matrix<doubleType> & state_right,
                                        math::Matrix<double> & dstate_rightdProp_vars,
                                        const math::Matrix <doubleType> & control,
                                        const doubleType & step_size,
                                        const double& dstep_sizedProp_var,
                                        const bool & needSTM);

            virtual void errorControlledStep(const math::Matrix<doubleType> & state_left,
                                             const math::Matrix<double> & dstate_leftdProp_vars,
                                             math::Matrix<doubleType> & state_right,
                                             math::Matrix<double> & dstate_rightdProp_vars,
                                             const math::Matrix <doubleType> & control,
                                             const doubleType & step_size,
                                             const double & dstep_sizedProp_var,
                                             const bool & needSTM,
                                             doubleType & error,
                                             math::Matrix<double> & error_scaling_factors);

            virtual void computeError(doubleType & error, math::Matrix<double> & error_scaling_factors);

            // fields
        private:

            //////////////
            // stage 1
            //////////////
            void stage1(const math::Matrix<doubleType> & state_left,
                        const math::Matrix<double> & dstate_leftdProp_vars,
                        const doubleType & step_size,
                        const bool & needSTM);

            void stage1wSTM(const math::Matrix<doubleType> & state_left,
                            const math::Matrix<double> & dstate_leftdProp_vars,
                            const doubleType & step_size,
                            const bool & needSTM);

            void stage1_w_control(const math::Matrix<doubleType> & state_left,
                                  const math::Matrix<double> & dstate_leftdProp_vars,
                                  const math::Matrix<doubleType> & control,
                                  const doubleType & step_size,
                                  const bool & needSTM);

            void stage1wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                      const math::Matrix<double> & dstate_leftdProp_vars,
                                      const math::Matrix<doubleType> & control,
                                      const doubleType & step_size,
                                      const bool & needSTM);

            void stage1_state_update(const math::Matrix<doubleType> & state_left,
                                     const doubleType & step_size,
                                     const size_t & num_states_to_integrate);

            void stage1_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                                       const doubleType & step_size);

            //////////////
            // stage 2
            //////////////
            void stage2(const math::Matrix<doubleType> & state_left,
                        const math::Matrix<double> & dstate_leftdProp_vars,
                        const doubleType & step_size,
                        const bool & needSTM);

            void stage2wSTM(const math::Matrix<doubleType> & state_left,
                            const math::Matrix<double> & dstate_leftdProp_vars,
                            const doubleType & step_size,
                            const bool & needSTM);

            void stage2_w_control(const math::Matrix<doubleType> & state_left,
                                  const math::Matrix<double> & dstate_leftdProp_vars,
                                  const math::Matrix<doubleType> & control,
                                  const doubleType & step_size,
                                  const bool & needSTM);

            void stage2wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                      const math::Matrix<double> & dstate_leftdProp_vars,
                                      const math::Matrix<doubleType> & control,
                                      const doubleType & step_size,
                                      const bool & needSTM);

            void stage2_state_update(const math::Matrix<doubleType> & state_left,
                                     const doubleType & step_size,
                                     const size_t & num_states_to_integrate);

            void stage2_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                                       const doubleType & step_size);

            //////////////
            // stage 3
            //////////////
            void stage3(const math::Matrix<doubleType> & state_left,
                        const math::Matrix<double> & dstate_leftdProp_vars,
                        const doubleType & step_size,
                        const bool & needSTM);

            void stage3wSTM(const math::Matrix<doubleType> & state_left,
                            const math::Matrix<double> & dstate_leftdProp_vars,
                            const doubleType & step_size,
                            const bool & needSTM);

            void stage3_w_control(const math::Matrix<doubleType> & state_left,
                                  const math::Matrix<double> & dstate_leftdProp_vars,
                                  const math::Matrix<doubleType> & control,
                                  const doubleType & step_size,
                                  const bool & needSTM);

            void stage3wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                      const math::Matrix<double> & dstate_leftdProp_vars,
                                      const math::Matrix<doubleType> & control,
                                      const doubleType & step_size,
                                      const bool & needSTM);

            void stage3_state_update(const math::Matrix<doubleType> & state_left,
                                     const doubleType & step_size,
                                     const size_t & num_states_to_integrate);

            void stage3_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                                       const doubleType & step_size);

            //////////////
            // stage 4
            //////////////
            void stage4(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage4wSTM(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage4_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage4wSTM_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage4_state_update(const math::Matrix<doubleType> & state_left,
                const doubleType & step_size,
                const size_t & num_states_to_integrate);

            void stage4_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size);

            //////////////
            // stage 5
            //////////////
            void stage5(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage5wSTM(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage5_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage5wSTM_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage5_state_update(const math::Matrix<doubleType> & state_left,
                const doubleType & step_size,
                const size_t & num_states_to_integrate);

            void stage5_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size);

            //////////////
            // stage 6
            //////////////
            void stage6(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage6wSTM(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage6_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage6wSTM_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage6_state_update(const math::Matrix<doubleType> & state_left,
                const doubleType & step_size,
                const size_t & num_states_to_integrate);

            void stage6_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size);

            //////////////
            // stage 7
            //////////////
            void stage7(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage7wSTM(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage7_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage7wSTM_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage7_state_update(const math::Matrix<doubleType> & state_left,
                const doubleType & step_size,
                const size_t & num_states_to_integrate);

            void stage7_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size);
            
            //////////////
            // stage 8
            //////////////
            void stage8(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage8wSTM(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage8_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage8wSTM_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage8_state_update(const math::Matrix<doubleType> & state_left,
                const doubleType & step_size,
                const size_t & num_states_to_integrate);

            void stage8_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size);

            //////////////
            // stage 9
            //////////////
            void stage9(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage9wSTM(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage9_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage9wSTM_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage9_state_update(const math::Matrix<doubleType> & state_left,
                const doubleType & step_size,
                const size_t & num_states_to_integrate);

            void stage9_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size);

            //////////////
            // stage 10
            //////////////
            void stage10(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage10wSTM(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage10_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage10wSTM_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage10_state_update(const math::Matrix<doubleType> & state_left,
                const doubleType & step_size,
                const size_t & num_states_to_integrate);

            void stage10_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size);
            
            //////////////
            // stage 11
            //////////////
            void stage11(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage11wSTM(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage11_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage11wSTM_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage11_state_update(const math::Matrix<doubleType> & state_left,
                const doubleType & step_size,
                const size_t & num_states_to_integrate);

            void stage11_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size);

            //////////////
            // stage 12
            //////////////
            void stage12(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage12wSTM(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage12_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage12wSTM_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage12_state_update(const math::Matrix<doubleType> & state_left,
                const doubleType & step_size,
                const size_t & num_states_to_integrate);

            void stage12_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size);

            void stage13eighthOrder(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            //////////////
            // stage 13
            //////////////
            void stage13eighthOrderwSTM(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size,
                const bool & needSTM);

            void stage13eighthOrder_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage13eighthOrderwSTM_w_control(const math::Matrix<doubleType> & state_left,
                const math::Matrix<double> & dstate_leftdProp_vars,
                const math::Matrix<doubleType> & control,
                const doubleType & step_size,
                const bool & needSTM);

            void stage13eighthOrder_state_update(const math::Matrix<doubleType> & state_left,
                const doubleType & step_size,
                const size_t & num_states_to_integrate);

            void stage13eighthOrder_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                const doubleType & step_size);

            void stage13seventhOrder(const math::Matrix<doubleType> & state_left,
                                     const math::Matrix<double> & dstate_leftdProp_vars,
                                     const math::Matrix<doubleType> & control,
                                     const doubleType & step_size,
                                     const bool & needSTM);
            /*
            void stage13seventhOrderwSTM(const math::Matrix<doubleType> & state_left,
                                         const math::Matrix<double> & dstate_leftdProp_vars,
                                         const math::Matrix<doubleType> & control,
                                         const doubleType & step_size,
                                         const bool & needSTM);
            */

            size_t ns;

            size_t num_prop_var_deriv_states;

            math::Matrix<doubleType> f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13,
                                     y, x_seventh, x_eighth, x_left, x_right;

            EMTG::math::Matrix <double> df1dProp_vars, df2dProp_vars, df3dProp_vars,
                                        df4dProp_vars, df5dProp_vars, df6dProp_vars,
                                        df7dProp_vars, df8dProp_vars, df9dProp_vars,
                                        df10dProp_vars, df11dProp_vars, df12dProp_vars, 
                                        df13dProp_vars, dx_eighthdProp_vars,
                                        dydProp_vars;
            
            // commented-out coefficients are equal to zero
            double a21;
            double a31;
            double a32;
            double a41;
            // double a42; 
            double a43;
            double a51;
            // double a52; 
            double a53;
            double a54;
            double a61;
            // double a62; 
            // double a63; 
            double a64;
            double a65;
            double a71;
            // double a72; 
            // double a73; 
            double a74;
            double a75;
            double a76;
            double a81;
            // double a82; 
            // double a83; 
            double a84;
            double a85;
            double a86;
            double a87;
            double a91;
            // double a92; 
            // double a93; 
            double a94;
            double a95;
            double a96;
            double a97;
            double a98;
            double a10_1;
            // double a10_2; 
            // double a10_3; 
            double a10_4;
            double a10_5;
            double a10_6;
            double a10_7;
            double a10_8;
            double a10_9;
            double a11_1;
            // double a11_2; 
            // double a11_3; 
            double a11_4;
            double a11_5;
            double a11_6;
            double a11_7;
            double a11_8;
            double a11_9;
            double a11_10;
            double a12_1;
            // double a12_2; 
            // double a12_3; 
            double a12_4;
            double a12_5;
            double a12_6;
            double a12_7;
            double a12_8;
            double a12_9;
            double a12_10;
            double a12_11;
            double a13_1;
            // double a13_2;
            // double a13_3;
            double a13_4;
            double a13_5;
            double a13_6;
            double a13_7;
            double a13_8;
            double a13_9;
            double a13_10;
            double a13_11;

            double b1lower;
            // double b2lower;
            // double b3lower;
            // double b4lower;
            // double b5lower;
            double b6lower;
            double b7lower;
            double b8lower;
            double b9lower;
            double b10lower; 
            double b11lower; 
            double b12lower; 
            // double b13lower;

            double b1upper;
            // double b2upper;
            // double b3upper;
            // double b4upper;
            // double b5upper;
            double b6upper;
            double b7upper;
            double b8upper;
            double b9upper;
            double b10upper;
            double b11upper;
            double b12upper;
            double b13upper;

            double c2;
            double c3;
            double c4;
            double c5;
            double c6;
            double c7;
            double c8;
            double c9;
            double c10;
            double c11;
            double c12;
            double c13;

        }; // end rk8 class definition

    } // end integration namespace
} // end EMTG namespace

#endif