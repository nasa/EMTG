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

#include "Integrand.h"
#include "RungeKutta8.h"

namespace EMTG {
    namespace Integration
    {

        //standard constructor
        RungeKutta8::RungeKutta8(Integrand * integrand_in, const size_t & ns_in, const size_t & num_prop_var_deriv_states) : IntegrationScheme(integrand_in),
            f1(ns_in, 1, 0.0), f2(ns_in, 1, 0.0), f3(ns_in, 1, 0.0), 
            f4(ns_in, 1, 0.0), f5(ns_in, 1, 0.0), f6(ns_in, 1, 0.0), 
            f7(ns_in, 1, 0.0), f8(ns_in, 1, 0.0), f9(ns_in, 1, 0.0), 
            f10(ns_in, 1, 0.0), f11(ns_in, 1, 0.0), f12(ns_in, 1, 0.0), 
            f13(ns_in, 1, 0.0), y(ns_in, 1, 0.0), x_seventh(ns_in, 1, 0.0), x_eighth(ns_in, 1, 0.0),
            df1dProp_vars(num_prop_var_deriv_states, 2, 0.0), df2dProp_vars(num_prop_var_deriv_states, 2, 0.0), df3dProp_vars(num_prop_var_deriv_states, 2, 0.0),
            df4dProp_vars(num_prop_var_deriv_states, 2, 0.0), df5dProp_vars(num_prop_var_deriv_states, 2, 0.0), df6dProp_vars(num_prop_var_deriv_states, 2, 0.0),
            df7dProp_vars(num_prop_var_deriv_states, 2, 0.0), df8dProp_vars(num_prop_var_deriv_states, 2, 0.0), df9dProp_vars(num_prop_var_deriv_states, 2, 0.0),
            df10dProp_vars(num_prop_var_deriv_states, 2, 0.0), df11dProp_vars(num_prop_var_deriv_states, 2, 0.0), df12dProp_vars(num_prop_var_deriv_states, 2, 0.0),
            df13dProp_vars(num_prop_var_deriv_states, 2, 0.0), dx_eighthdProp_vars(num_prop_var_deriv_states, 2, 0.0),
            dydProp_vars(num_prop_var_deriv_states, 2, 0.0)
        {
            this->ns = ns_in;

            this->a21 = 1.0 / 18.0;
            this->a31 = 1.0 / 48.0;
            this->a32 = 1.0 / 16.0;
            this->a41 = 1.0 / 32.0;
            // this-> a42 = 0.0
            this->a43 = 3.0 / 32.0;
            this->a51 = 5.0 / 16.0;
            // this->a52 = 0.0;
            this->a53 = -75.0 / 64.0;
            this->a54 = 75.0 / 64.0;
            this->a61 = 3.0 / 80.0;
            // this->a62 = 0.0;
            // this->a63 = 0.0;
            this->a64 = 3.0 / 16.0;
            this->a65 = 3.0 / 20.0;
            this->a71 = 29443841.0 / 614563906.0;
            // this->a72 = 0.0;
            // this->a73 = 0.0;
            this->a74 = 77736538.0 / 692538347.0;
            this->a75 = -28693883.0 / 1125000000.0;
            this->a76 = 23124283.0 / 1800000000.0;
            this->a81 = 16016141.0 / 946692911.0;
            // this->a82 = 0.0;
            // this->a83 = 0.0;
            this->a84 = 61564180.0 / 158732637.0;
            this->a85 = 22789713.0 / 633445777.0;
            this->a86 = 545815736.0 / 2771057229.0;
            this->a87 = -180193667.0 / 1043307555.0;
            this->a91 = 39632708.0 / 573591083.0;
            // this->a92 = 0.0;
            // this->a93 = 0.0;
            this->a94 = -433636366.0 / 683701615.0;
            this->a95 = -421739975.0 / 2616292301.0;
            this->a96 = 100302831.0 / 723423059.0;
            this->a97 = 790204164.0 / 839813087.0;
            this->a98 = 800635310.0 / 3783071287.0;
            this->a10_1 = 246121993.0 / 1340847787.0;
            // this->a10_2 = 0.0;
            // this->a10_3 = 0.0;
            this->a10_4 = -37695042795.0 / 15268766246.0;
            this->a10_5 = -309121744.0 / 1061227803.0;
            this->a10_6 = -12992083.0 / 490766935.0;
            this->a10_7 = 6005943493.0 / 2108947869.0;
            this->a10_8 = 393006217.0 / 1396673457.0;
            this->a10_9 = 123872331.0 / 1001029789.0;
            this->a11_1 = -1028468189.0 / 846180014.0;
            // this->a11_2 = 0.0;
            // this->a11_3 = 0.0;
            this->a11_4 = 8478235783.0 / 508512852.0;
            this->a11_5 = 1311729495.0 / 1432422823.0;
            this->a11_6 = -10304129995.0 / 1701304382.0;
            this->a11_7 = -48777925059.0 / 3047939560.0;
            this->a11_8 = 15336726248.0 / 1032824649.0;
            this->a11_9 = -45442868181.0 / 3398467696.0;
            this->a11_10 = 3065993473.0 / 597172653.0;
            this->a12_1 = 185892177.0 / 718116043.0;
            // this->a12_2 = 0.0;
            // this->a12_3 = 0.0;
            this->a12_4 = -3185094517.0 / 667107341.0;
            this->a12_5 = -477755414.0 / 1098053517.0;
            this->a12_6 = -703635378.0 / 230739211.0;
            this->a12_7 = 5731566787.0 / 1027545527.0;
            this->a12_8 = 5232866602.0 / 850066563.0;
            this->a12_9 = -4093664535.0 / 808688257.0;
            this->a12_10 = 3962137247.0 / 1805957418.0;
            this->a12_11 = 65686358.0 / 487910083.0;
            this->a13_1 = 403863854.0 / 491063109.0;
            // this->a13_2 = 0.0;
            // this->a13_3 = 0.0;
            this->a13_4 = -5068492393.0 / 434740067.0;
            this->a13_5 = -411421997.0 / 543043805.0;
            this->a13_6 = 652783627.0 / 914296604.0;
            this->a13_7 = 11173962825.0 / 925320556.0;
            this->a13_8 = -13158990841.0 / 6184727034.0;
            this->a13_9 = 3936647629.0 / 1978049680.0;
            this->a13_10 = -160528059.0 / 685178525.0;
            this->a13_11 = 248638103.0 / 1413531060.0;

            this->b1lower = 13451932.0 / 455176623.0;
            // this->b2lower = 0.0; 
            // this->b3lower = 0.0; 
            // this->b4lower = 0.0; 
            // this->b5lower = 0.0; 
            this->b6lower = -808719846.0 / 976000145.0;
            this->b7lower = 1757004468.0 / 5645159321.0;
            this->b8lower = 656045339.0 / 265891186.0;
            this->b9lower = -3867574721.0 / 1518517206.0;
            this->b10lower = 465885868.0 / 322736535.0;
            this->b11lower = 53011238.0 / 667516719.0;
            this->b12lower = 2.0 / 45.0;
            // this->b13lower = 0;

            this->b1upper = 14005451.0 / 335480064.0;
            // this->b2upper = 0.0;
            // this->b3upper = 0.0;
            // this->b4upper = 0.0;
            // this->b5upper = 0.0;
            this->b6upper = -59238493.0 / 1068277825.0;
            this->b7upper = 181606767.0 / 758867731.0;
            this->b8upper = 561292985.0 / 797845732.0;
            this->b9upper = -1041891430.0 / 1371343529.0;
            this->b10upper = 760417239.0 / 1151165299.0;
            this->b11upper = 118820643.0 / 751138087.0;
            this->b12upper = -528747749.0 / 2220607170.0;
            this->b13upper = 1.0 / 4.0;

            //RK node constants
            //these encode the positions within the RK step of each stage
            this->c2 = 1.0 / 18.0;
            this->c3 = 1.0 / 12.0;
            this->c4 = 1.0 / 8.0;
            this->c5 = 5.0 / 16.0;
            this->c6 = 3.0 / 8.0;
            this->c7 = 59.0 / 400.0;
            this->c8 = 93.0 / 200.0;
            this->c9 = 5490023248.0 / 9719169821.0;
            this->c10 = 13.0 / 20.0;
            this->c11 = 1201146811.0 / 1299019798.0;
            this->c12 = 1.0;
            this->c13 = 1.0;

        }

        //destructor
        RungeKutta8 ::~RungeKutta8() {}

        //////////////
        // stage 1
        //////////////
        void RungeKutta8::stage1(const math::Matrix<doubleType> & state_left, 
                                 const math::Matrix<double> & dstate_leftdProp_vars,
                                 const doubleType & step_size,
                                 const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var);
            this->integrand->evaluate(state_left, dstate_leftdProp_vars, this->f1, this->df1dProp_vars, needSTM);

            this->stage1_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage1wSTM(const math::Matrix<doubleType> & state_left,
                                     const math::Matrix<double> & dstate_leftdProp_vars,
                                     const doubleType & step_size,
                                     const bool & needSTM)
        {
            this->stage1(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage1_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage1_w_control(const math::Matrix<doubleType> & state_left,
                                           const math::Matrix<double> & dstate_leftdProp_vars,
                                           const math::Matrix<doubleType> & control,
                                           const doubleType & step_size,
                                           const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var);
            this->integrand->evaluate(state_left, dstate_leftdProp_vars, this->f1, this->df1dProp_vars, control, needSTM);

            this->stage1_state_update(state_left, step_size, num_states_to_integrate);
            
        }

        void RungeKutta8::stage1wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                               const math::Matrix<double> & dstate_leftdProp_vars,
                                               const math::Matrix<doubleType> & control,
                                               const doubleType & step_size,
                                               const bool & needSTM)
        {
            this->stage1_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage1_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage1_state_update(const math::Matrix<doubleType> & state_left,
                                              const doubleType & step_size,
                                              const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                             this->a21*step_size*this->f1(k);
            }
        }

        void RungeKutta8::stage1_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                                                const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->a21*step_size*this->df1dProp_vars(k, 0))_GETVALUE;
                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->a21*(this->dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))) _GETVALUE;
            }
        }


        //////////////
        // stage 2
        //////////////
        void RungeKutta8::stage2(const math::Matrix<doubleType> & state_left,
                                 const math::Matrix<double> & dstate_leftdProp_vars,
                                 const doubleType & step_size,
                                 const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c2 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c2 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f2, this->df2dProp_vars, needSTM);

            this->stage2_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage2wSTM(const math::Matrix<doubleType> & state_left,
                                     const math::Matrix<double> & dstate_leftdProp_vars,
                                     const doubleType & step_size,
                                     const bool & needSTM)
        {
            this->stage2(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage2_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage2_w_control(const math::Matrix<doubleType> & state_left,
                                           const math::Matrix<double> & dstate_leftdProp_vars,
                                           const math::Matrix<doubleType> & control,
                                           const doubleType & step_size,
                                           const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c2 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c2 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f2, this->df2dProp_vars, control, needSTM);

            this->stage2_state_update(state_left, step_size, num_states_to_integrate);

        }

        void RungeKutta8::stage2wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                               const math::Matrix<double> & dstate_leftdProp_vars,
                                               const math::Matrix<doubleType> & control,
                                               const doubleType & step_size,
                                               const bool & needSTM)
        {
            this->stage2_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage2_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage2_state_update(const math::Matrix<doubleType> & state_left,
                                              const doubleType & step_size,
                                              const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                             this->a31*step_size*this->f1(k) +
                             this->a32*step_size*this->f2(k);
            }
        }

        void RungeKutta8::stage2_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                                                const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->a31*step_size*this->df1dProp_vars(k, 0)
                                         + this->a32*step_size*this->df2dProp_vars(k, 0)) _GETVALUE;
                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->a31*(this->dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))
                                         + this->a32*(this->dstep_sizedProp_var*this->f2(k) + step_size * this->df2dProp_vars(k, 1))) _GETVALUE;
            }
        }


        //////////////
        // stage 3
        //////////////
        void RungeKutta8::stage3(const math::Matrix<doubleType> & state_left,
                                 const math::Matrix<double> & dstate_leftdProp_vars,
                                 const doubleType & step_size,
                                 const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c3 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c3 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f3, this->df3dProp_vars, needSTM);

            this->stage3_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage3wSTM(const math::Matrix<doubleType> & state_left,
                                     const math::Matrix<double> & dstate_leftdProp_vars,
                                     const doubleType & step_size,
                                     const bool & needSTM)
        {
            this->stage3(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage3_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage3_w_control(const math::Matrix<doubleType> & state_left,
                                           const math::Matrix<double> & dstate_leftdProp_vars,
                                           const math::Matrix<doubleType> & control,
                                           const doubleType & step_size,
                                           const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c3 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c3 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f3, this->df3dProp_vars, control, needSTM);

            this->stage3_state_update(state_left, step_size, num_states_to_integrate);

        }

        void RungeKutta8::stage3wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                               const math::Matrix<double> & dstate_leftdProp_vars,
                                               const math::Matrix<doubleType> & control,
                                               const doubleType & step_size,
                                               const bool & needSTM)
        {
            this->stage3_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage3_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage3_state_update(const math::Matrix<doubleType> & state_left,
                                              const doubleType & step_size,
                                              const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                             this->a41*step_size*this->f1(k) +
                             this->a43*step_size*this->f3(k);
            }
        }

        void RungeKutta8::stage3_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                                                const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->a41*step_size*this->df1dProp_vars(k, 0)
                                         + this->a43*step_size*this->df3dProp_vars(k, 0)) _GETVALUE;
                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->a41*(dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))
                                         + this->a43*(dstep_sizedProp_var*this->f3(k) + step_size * this->df3dProp_vars(k, 1))) _GETVALUE;
            }
        }

        //////////////
        // stage 4
        //////////////
        void RungeKutta8::stage4(const math::Matrix<doubleType> & state_left,
                                 const math::Matrix<double> & dstate_leftdProp_vars,
                                 const doubleType & step_size,
                                 const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c4 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c4 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f4, this->df4dProp_vars, needSTM);

            this->stage4_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage4wSTM(const math::Matrix<doubleType> & state_left,
                                     const math::Matrix<double> & dstate_leftdProp_vars,
                                     const doubleType & step_size,
                                     const bool & needSTM)
        {
            this->stage4(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage4_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage4_w_control(const math::Matrix<doubleType> & state_left,
                                           const math::Matrix<double> & dstate_leftdProp_vars,
                                           const math::Matrix<doubleType> & control,
                                           const doubleType & step_size,
                                           const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c4 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c4 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f4, this->df4dProp_vars, control, needSTM);

            this->stage4_state_update(state_left, step_size, num_states_to_integrate);

        }

        void RungeKutta8::stage4wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                               const math::Matrix<double> & dstate_leftdProp_vars,
                                               const math::Matrix<doubleType> & control,
                                               const doubleType & step_size,
                                               const bool & needSTM)
        {
            this->stage4_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage4_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage4_state_update(const math::Matrix<doubleType> & state_left,
                                              const doubleType & step_size,
                                              const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                             this->a51*step_size*this->f1(k) +
                             this->a53*step_size*this->f3(k) +
                             this->a54*step_size*this->f4(k);
            }
        }

        void RungeKutta8::stage4_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                                                const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->a51*(step_size*this->df1dProp_vars(k, 0))
                                         + this->a53*(step_size*this->df3dProp_vars(k, 0))
                                         + this->a54*(step_size*this->df4dProp_vars(k, 0))) _GETVALUE;

                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->a51*(dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))
                                         + this->a53*(dstep_sizedProp_var*this->f3(k) + step_size * this->df3dProp_vars(k, 1))
                                         + this->a54*(dstep_sizedProp_var*this->f4(k) + step_size * this->df4dProp_vars(k, 1))) _GETVALUE;
            }
        }


        //////////////
        // stage 5
        //////////////
        void RungeKutta8::stage5(const math::Matrix<doubleType> & state_left,
                                 const math::Matrix<double> & dstate_leftdProp_vars,
                                 const doubleType & step_size,
                                 const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c5 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c5 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f5, this->df5dProp_vars, needSTM);

            this->stage5_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage5wSTM(const math::Matrix<doubleType> & state_left,
                                     const math::Matrix<double> & dstate_leftdProp_vars,
                                     const doubleType & step_size,
                                     const bool & needSTM)
        {
            this->stage5(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage5_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage5_w_control(const math::Matrix<doubleType> & state_left,
                                           const math::Matrix<double> & dstate_leftdProp_vars,
                                           const math::Matrix<doubleType> & control,
                                           const doubleType & step_size,
                                           const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c5 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c5 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f5, this->df5dProp_vars, control, needSTM);

            this->stage5_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage5wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                               const math::Matrix<double> & dstate_leftdProp_vars,
                                               const math::Matrix<doubleType> & control,
                                               const doubleType & step_size,
                                               const bool & needSTM)
        {
            this->stage5_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage5_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage5_state_update(const math::Matrix<doubleType> & state_left,
                                              const doubleType & step_size,
                                              const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                             this->a61*step_size*this->f1(k) +
                             this->a64*step_size*this->f4(k) +
                             this->a65*step_size*this->f5(k);
            }
        }

        void RungeKutta8::stage5_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
            const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0) 
                                         + this->a61*(step_size*this->df1dProp_vars(k, 0))
                                         + this->a64*(step_size*this->df4dProp_vars(k, 0))
                                         + this->a65*(step_size*this->df5dProp_vars(k, 0))) _GETVALUE;

                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->a61*(dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))
                                         + this->a64*(dstep_sizedProp_var*this->f4(k) + step_size * this->df4dProp_vars(k, 1))
                                         + this->a65*(dstep_sizedProp_var*this->f5(k) + step_size * this->df5dProp_vars(k, 1))) _GETVALUE;
            }
        }

        //////////////
        // stage 6
        //////////////
        void RungeKutta8::stage6(const math::Matrix<doubleType> & state_left,
                                 const math::Matrix<double> & dstate_leftdProp_vars,
                                 const doubleType & step_size,
                                 const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c6 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c6 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f6, this->df6dProp_vars, needSTM);

            this->stage6_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage6wSTM(const math::Matrix<doubleType> & state_left,
                                     const math::Matrix<double> & dstate_leftdProp_vars,
                                     const doubleType & step_size,
                                     const bool & needSTM)
        {
            this->stage6(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage6_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage6_w_control(const math::Matrix<doubleType> & state_left,
                                           const math::Matrix<double> & dstate_leftdProp_vars,
                                           const math::Matrix<doubleType> & control,
                                           const doubleType & step_size,
                                           const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c6 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c6 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f6, this->df6dProp_vars, control, needSTM);

            this->stage6_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage6wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                               const math::Matrix<double> & dstate_leftdProp_vars,
                                               const math::Matrix<doubleType> & control,
                                               const doubleType & step_size,
                                               const bool & needSTM)
        {
            this->stage6_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage6_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage6_state_update(const math::Matrix<doubleType> & state_left,
                                              const doubleType & step_size,
                                              const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                             this->a71*step_size*this->f1(k) +
                             this->a74*step_size*this->f4(k) +
                             this->a75*step_size*this->f5(k) +
                             this->a76*step_size*this->f6(k);
            }
        }

        void RungeKutta8::stage6_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                                                const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->a71*(step_size*this->df1dProp_vars(k, 0))
                                         + this->a74*(step_size*this->df4dProp_vars(k, 0))
                                         + this->a75*(step_size*this->df5dProp_vars(k, 0))
                                         + this->a76*(step_size*this->df6dProp_vars(k, 0))) _GETVALUE;

                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->a71*(dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))
                                         + this->a74*(dstep_sizedProp_var*this->f4(k) + step_size * this->df4dProp_vars(k, 1))
                                         + this->a75*(dstep_sizedProp_var*this->f5(k) + step_size * this->df5dProp_vars(k, 1))
                                         + this->a76*(dstep_sizedProp_var*this->f6(k) + step_size * this->df6dProp_vars(k, 1))) _GETVALUE;
            }
        }


        //////////////
        // stage 7
        //////////////
        void RungeKutta8::stage7(const math::Matrix<doubleType> & state_left,
                                 const math::Matrix<double> & dstate_leftdProp_vars,
                                 const doubleType & step_size,
                                 const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c7 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c7 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f7, this->df7dProp_vars, needSTM);

            this->stage7_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage7wSTM(const math::Matrix<doubleType> & state_left,
                                     const math::Matrix<double> & dstate_leftdProp_vars,
                                     const doubleType & step_size,
                                     const bool & needSTM)
        {
            this->stage7(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage7_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage7_w_control(const math::Matrix<doubleType> & state_left,
                                           const math::Matrix<double> & dstate_leftdProp_vars,
                                           const math::Matrix<doubleType> & control,
                                           const doubleType & step_size,
                                           const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c7 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c7 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f7, this->df7dProp_vars, control, needSTM);

            this->stage7_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage7wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                               const math::Matrix<double> & dstate_leftdProp_vars,
                                               const math::Matrix<doubleType> & control,
                                               const doubleType & step_size,
                                               const bool & needSTM)
        {
            this->stage7_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage7_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage7_state_update(const math::Matrix<doubleType> & state_left,
                                              const doubleType & step_size,
                                              const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                             this->a81*step_size*this->f1(k) +
                             this->a84*step_size*this->f4(k) +
                             this->a85*step_size*this->f5(k) +
                             this->a86*step_size*this->f6(k) +
                             this->a87*step_size*this->f7(k);
            }
        }

        void RungeKutta8::stage7_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
            const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->a81*(step_size*this->df1dProp_vars(k, 0))
                                         + this->a84*(step_size*this->df4dProp_vars(k, 0))
                                         + this->a85*(step_size*this->df5dProp_vars(k, 0))
                                         + this->a86*(step_size*this->df6dProp_vars(k, 0))
                                         + this->a87*(step_size*this->df7dProp_vars(k, 0))) _GETVALUE;

                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->a81*(dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))
                                         + this->a84*(dstep_sizedProp_var*this->f4(k) + step_size * this->df4dProp_vars(k, 1))
                                         + this->a85*(dstep_sizedProp_var*this->f5(k) + step_size * this->df5dProp_vars(k, 1))
                                         + this->a86*(dstep_sizedProp_var*this->f6(k) + step_size * this->df6dProp_vars(k, 1))
                                         + this->a87*(dstep_sizedProp_var*this->f7(k) + step_size * this->df7dProp_vars(k, 1))) _GETVALUE;
            }
        }


        //////////////
        // stage 8
        //////////////
        void RungeKutta8::stage8(const math::Matrix<doubleType> & state_left,
                                 const math::Matrix<double> & dstate_leftdProp_vars,
                                 const doubleType & step_size,
                                 const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c8 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c8 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f8, this->df8dProp_vars, needSTM);

            this->stage8_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage8wSTM(const math::Matrix<doubleType> & state_left,
                                     const math::Matrix<double> & dstate_leftdProp_vars,
                                     const doubleType & step_size,
                                     const bool & needSTM)
        {
            this->stage8(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage8_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage8_w_control(const math::Matrix<doubleType> & state_left,
                                           const math::Matrix<double> & dstate_leftdProp_vars,
                                           const math::Matrix<doubleType> & control,
                                           const doubleType & step_size,
                                           const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c8 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c8 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f8, this->df8dProp_vars, control, needSTM);

            this->stage8_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage8wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                               const math::Matrix<double> & dstate_leftdProp_vars,
                                               const math::Matrix<doubleType> & control,
                                               const doubleType & step_size,
                                               const bool & needSTM)
        {
            this->stage8_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage8_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage8_state_update(const math::Matrix<doubleType> & state_left,
                                              const doubleType & step_size,
                                              const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                             this->a91*step_size*this->f1(k) +
                             this->a94*step_size*this->f4(k) +
                             this->a95*step_size*this->f5(k) +
                             this->a96*step_size*this->f6(k) +
                             this->a97*step_size*this->f7(k) +
                             this->a98*step_size*this->f8(k);
            }
        }

        void RungeKutta8::stage8_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                                                const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->a91*(step_size*this->df1dProp_vars(k, 0))
                                         + this->a94*(step_size*this->df4dProp_vars(k, 0))
                                         + this->a95*(step_size*this->df5dProp_vars(k, 0))
                                         + this->a96*(step_size*this->df6dProp_vars(k, 0))
                                         + this->a97*(step_size*this->df7dProp_vars(k, 0))
                                         + this->a98*(step_size*this->df8dProp_vars(k, 0))) _GETVALUE;

                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->a91*(dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))
                                         + this->a94*(dstep_sizedProp_var*this->f4(k) + step_size * this->df4dProp_vars(k, 1))
                                         + this->a95*(dstep_sizedProp_var*this->f5(k) + step_size * this->df5dProp_vars(k, 1))
                                         + this->a96*(dstep_sizedProp_var*this->f6(k) + step_size * this->df6dProp_vars(k, 1))
                                         + this->a97*(dstep_sizedProp_var*this->f7(k) + step_size * this->df7dProp_vars(k, 1))
                                         + this->a98*(dstep_sizedProp_var*this->f8(k) + step_size * this->df8dProp_vars(k, 1))) _GETVALUE;
            }
        }


        //////////////
        // stage 9
        //////////////
        void RungeKutta8::stage9(const math::Matrix<doubleType> & state_left,
                                 const math::Matrix<double> & dstate_leftdProp_vars,
                                 const doubleType & step_size,
                                 const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c9 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c9 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f9, this->df9dProp_vars, needSTM);

            this->stage9_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage9wSTM(const math::Matrix<doubleType> & state_left,
                                     const math::Matrix<double> & dstate_leftdProp_vars,
                                     const doubleType & step_size,
                                     const bool & needSTM)
        {
            this->stage9(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage9_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage9_w_control(const math::Matrix<doubleType> & state_left,
                                           const math::Matrix<double> & dstate_leftdProp_vars,
                                           const math::Matrix<doubleType> & control,
                                           const doubleType & step_size,
                                           const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c9 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c9 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f9, this->df9dProp_vars, control, needSTM);

            this->stage9_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage9wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                               const math::Matrix<double> & dstate_leftdProp_vars,
                                               const math::Matrix<doubleType> & control,
                                               const doubleType & step_size,
                                               const bool & needSTM)
        {
            this->stage9_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage9_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage9_state_update(const math::Matrix<doubleType> & state_left,
                                              const doubleType & step_size,
                                              const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                             this->a10_1*step_size*this->f1(k) +
                             this->a10_4*step_size*this->f4(k) +
                             this->a10_5*step_size*this->f5(k) +
                             this->a10_6*step_size*this->f6(k) +
                             this->a10_7*step_size*this->f7(k) +
                             this->a10_8*step_size*this->f8(k) +
                             this->a10_9*step_size*this->f9(k);
            }
        }

        void RungeKutta8::stage9_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
            const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->a10_1*(step_size*this->df1dProp_vars(k, 0))
                                         + this->a10_4*(step_size*this->df4dProp_vars(k, 0))
                                         + this->a10_5*(step_size*this->df5dProp_vars(k, 0))
                                         + this->a10_6*(step_size*this->df6dProp_vars(k, 0))
                                         + this->a10_7*(step_size*this->df7dProp_vars(k, 0))
                                         + this->a10_8*(step_size*this->df8dProp_vars(k, 0))
                                         + this->a10_9*(step_size*this->df9dProp_vars(k, 0))) _GETVALUE;

                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->a10_1*(dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))
                                         + this->a10_4*(dstep_sizedProp_var*this->f4(k) + step_size * this->df4dProp_vars(k, 1))
                                         + this->a10_5*(dstep_sizedProp_var*this->f5(k) + step_size * this->df5dProp_vars(k, 1))
                                         + this->a10_6*(dstep_sizedProp_var*this->f6(k) + step_size * this->df6dProp_vars(k, 1))
                                         + this->a10_7*(dstep_sizedProp_var*this->f7(k) + step_size * this->df7dProp_vars(k, 1))
                                         + this->a10_8*(dstep_sizedProp_var*this->f8(k) + step_size * this->df8dProp_vars(k, 1))
                                         + this->a10_9*(dstep_sizedProp_var*this->f9(k) + step_size * this->df9dProp_vars(k, 1))) _GETVALUE;
            }
        }

        //////////////
        // stage 10
        //////////////
        void RungeKutta8::stage10(const math::Matrix<doubleType> & state_left,
                                  const math::Matrix<double> & dstate_leftdProp_vars,
                                  const doubleType & step_size,
                                  const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c10 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c10 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f10, this->df10dProp_vars, needSTM);

            this->stage10_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage10wSTM(const math::Matrix<doubleType> & state_left,
                                      const math::Matrix<double> & dstate_leftdProp_vars,
                                      const doubleType & step_size,
                                      const bool & needSTM)
        {
            this->stage10(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage10_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage10_w_control(const math::Matrix<doubleType> & state_left,
                                            const math::Matrix<double> & dstate_leftdProp_vars,
                                            const math::Matrix<doubleType> & control,
                                            const doubleType & step_size,
                                            const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c10 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c10 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f10, this->df10dProp_vars, control, needSTM);

            this->stage10_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage10wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                                const math::Matrix<double> & dstate_leftdProp_vars,
                                                const math::Matrix<doubleType> & control,
                                                const doubleType & step_size,
                                                const bool & needSTM)
        {
            this->stage10_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage10_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage10_state_update(const math::Matrix<doubleType> & state_left,
                                               const doubleType & step_size,
                                               const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                           this->a11_1 *  step_size*this->f1(k) +
                           this->a11_4 *  step_size*this->f4(k) +
                           this->a11_5 *  step_size*this->f5(k) +
                           this->a11_6 *  step_size*this->f6(k) +
                           this->a11_7 *  step_size*this->f7(k) +
                           this->a11_8 *  step_size*this->f8(k) +
                           this->a11_9 *  step_size*this->f9(k) +
                           this->a11_10 * step_size*this->f10(k);
            }
        }

        void RungeKutta8::stage10_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                                                 const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->a11_1*(step_size*this->df1dProp_vars(k, 0))
                                         + this->a11_4*(step_size*this->df4dProp_vars(k, 0))
                                         + this->a11_5*(step_size*this->df5dProp_vars(k, 0))
                                         + this->a11_6*(step_size*this->df6dProp_vars(k, 0))
                                         + this->a11_7*(step_size*this->df7dProp_vars(k, 0))
                                         + this->a11_8*(step_size*this->df8dProp_vars(k, 0))
                                         + this->a11_9*(step_size*this->df9dProp_vars(k, 0))
                                         + this->a11_10*(step_size*this->df10dProp_vars(k, 0))) _GETVALUE;

                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->a11_1*(dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))
                                         + this->a11_4*(dstep_sizedProp_var*this->f4(k) + step_size * this->df4dProp_vars(k, 1))
                                         + this->a11_5*(dstep_sizedProp_var*this->f5(k) + step_size * this->df5dProp_vars(k, 1))
                                         + this->a11_6*(dstep_sizedProp_var*this->f6(k) + step_size * this->df6dProp_vars(k, 1))
                                         + this->a11_7*(dstep_sizedProp_var*this->f7(k) + step_size * this->df7dProp_vars(k, 1))
                                         + this->a11_8*(dstep_sizedProp_var*this->f8(k) + step_size * this->df8dProp_vars(k, 1))
                                         + this->a11_9*(dstep_sizedProp_var*this->f9(k) + step_size * this->df9dProp_vars(k, 1))
                                         + this->a11_10*(dstep_sizedProp_var*this->f10(k) + step_size * this->df10dProp_vars(k, 1))) _GETVALUE;
            }
        }


        //////////////
        // stage 11
        //////////////
        void RungeKutta8::stage11(const math::Matrix<doubleType> & state_left,
                                  const math::Matrix<double> & dstate_leftdProp_vars,
                                  const doubleType & step_size,
                                  const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c11 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c11 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f11, this->df11dProp_vars, needSTM);

            this->stage11_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage11wSTM(const math::Matrix<doubleType> & state_left,
                                      const math::Matrix<double> & dstate_leftdProp_vars,
                                      const doubleType & step_size,
                                      const bool & needSTM)
        {
            this->stage11(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage11_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage11_w_control(const math::Matrix<doubleType> & state_left,
                                            const math::Matrix<double> & dstate_leftdProp_vars,
                                            const math::Matrix<doubleType> & control,
                                            const doubleType & step_size,
                                            const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c11 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c11 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f11, this->df11dProp_vars, control, needSTM);

            this->stage11_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage11wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                                const math::Matrix<double> & dstate_leftdProp_vars,
                                                const math::Matrix<doubleType> & control,
                                                const doubleType & step_size,
                                                const bool & needSTM)
        {
            this->stage11_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage11_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage11_state_update(const math::Matrix<doubleType> & state_left,
                                               const doubleType & step_size,
                                               const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                           this->a12_1*step_size*this->f1(k) +
                           this->a12_4*step_size*this->f4(k) +
                           this->a12_5*step_size*this->f5(k) +
                           this->a12_6*step_size*this->f6(k) +
                           this->a12_7*step_size*this->f7(k) +
                           this->a12_8*step_size*this->f8(k) +
                           this->a12_9*step_size*this->f9(k) +
                           this->a12_10*step_size*this->f10(k) +
                           this->a12_11*step_size*this->f11(k);
            }
        }

        void RungeKutta8::stage11_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
            const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->a12_1*(step_size*this->df1dProp_vars(k, 0))
                                         + this->a12_4*(step_size*this->df4dProp_vars(k, 0))
                                         + this->a12_5*(step_size*this->df5dProp_vars(k, 0))
                                         + this->a12_6*(step_size*this->df6dProp_vars(k, 0))
                                         + this->a12_7*(step_size*this->df7dProp_vars(k, 0))
                                         + this->a12_8*(step_size*this->df8dProp_vars(k, 0))
                                         + this->a12_9*(step_size*this->df9dProp_vars(k, 0))
                                         + this->a12_10*(step_size*this->df10dProp_vars(k, 0))
                                         + this->a12_11*(step_size*this->df11dProp_vars(k, 0))) _GETVALUE;

                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->a12_1*(dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))
                                         + this->a12_4*(dstep_sizedProp_var*this->f4(k) + step_size * this->df4dProp_vars(k, 1))
                                         + this->a12_5*(dstep_sizedProp_var*this->f5(k) + step_size * this->df5dProp_vars(k, 1))
                                         + this->a12_6*(dstep_sizedProp_var*this->f6(k) + step_size * this->df6dProp_vars(k, 1))
                                         + this->a12_7*(dstep_sizedProp_var*this->f7(k) + step_size * this->df7dProp_vars(k, 1))
                                         + this->a12_8*(dstep_sizedProp_var*this->f8(k) + step_size * this->df8dProp_vars(k, 1))
                                         + this->a12_9*(dstep_sizedProp_var*this->f9(k) + step_size * this->df9dProp_vars(k, 1))
                                         + this->a12_10*(dstep_sizedProp_var*this->f10(k) + step_size * this->df10dProp_vars(k, 1))
                                         + this->a12_11*(dstep_sizedProp_var*this->f11(k) + step_size * this->df11dProp_vars(k, 1))) _GETVALUE;
            }
        }
        

        //////////////
        // stage 12
        //////////////
        void RungeKutta8::stage12(const math::Matrix<doubleType> & state_left,
                                  const math::Matrix<double> & dstate_leftdProp_vars,
                                  const doubleType & step_size,
                                  const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c12 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c12 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f12, this->df12dProp_vars, needSTM);

            this->stage12_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage12wSTM(const math::Matrix<doubleType> & state_left,
                                      const math::Matrix<double> & dstate_leftdProp_vars,
                                      const doubleType & step_size,
                                      const bool & needSTM)
        {
            this->stage12(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage12_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage12_w_control(const math::Matrix<doubleType> & state_left,
                                            const math::Matrix<double> & dstate_leftdProp_vars,
                                            const math::Matrix<doubleType> & control,
                                            const doubleType & step_size,
                                            const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c12 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c12 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f12, this->df12dProp_vars, control, needSTM);

            this->stage12_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage12wSTM_w_control(const math::Matrix<doubleType> & state_left,
                                                const math::Matrix<double> & dstate_leftdProp_vars,
                                                const math::Matrix<doubleType> & control,
                                                const doubleType & step_size,
                                                const bool & needSTM)
        {
            this->stage12_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage12_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage12_state_update(const math::Matrix<doubleType> & state_left,
                                               const doubleType & step_size,
                                               const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                           this->a13_1*step_size*this->f1(k) +
                           this->a13_4*step_size*this->f4(k) +
                           this->a13_5*step_size*this->f5(k) +
                           this->a13_6*step_size*this->f6(k) +
                           this->a13_7*step_size*this->f7(k) +
                           this->a13_8*step_size*this->f8(k) +
                           this->a13_9*step_size*this->f9(k) +
                           this->a13_10*step_size*this->f10(k) +
                           this->a13_11*step_size*this->f11(k);
            }
        }

        void RungeKutta8::stage12_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
            const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->a13_1*(step_size*this->df1dProp_vars(k, 0))
                                         + this->a13_4*(step_size*this->df4dProp_vars(k, 0))
                                         + this->a13_5*(step_size*this->df5dProp_vars(k, 0))
                                         + this->a13_6*(step_size*this->df6dProp_vars(k, 0))
                                         + this->a13_7*(step_size*this->df7dProp_vars(k, 0))
                                         + this->a13_8*(step_size*this->df8dProp_vars(k, 0))
                                         + this->a13_9*(step_size*this->df9dProp_vars(k, 0))
                                         + this->a13_10*(step_size*this->df10dProp_vars(k, 0))
                                         + this->a13_11*(step_size*this->df11dProp_vars(k, 0))) _GETVALUE;

                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->a13_1*(dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))
                                         + this->a13_4*(dstep_sizedProp_var*this->f4(k) + step_size * this->df4dProp_vars(k, 1))
                                         + this->a13_5*(dstep_sizedProp_var*this->f5(k) + step_size * this->df5dProp_vars(k, 1))
                                         + this->a13_6*(dstep_sizedProp_var*this->f6(k) + step_size * this->df6dProp_vars(k, 1))
                                         + this->a13_7*(dstep_sizedProp_var*this->f7(k) + step_size * this->df7dProp_vars(k, 1))
                                         + this->a13_8*(dstep_sizedProp_var*this->f8(k) + step_size * this->df8dProp_vars(k, 1))
                                         + this->a13_9*(dstep_sizedProp_var*this->f9(k) + step_size * this->df9dProp_vars(k, 1))
                                         + this->a13_10*(dstep_sizedProp_var*this->f10(k) + step_size * this->df10dProp_vars(k, 1))
                                         + this->a13_11*(dstep_sizedProp_var*this->f11(k) + step_size * this->df11dProp_vars(k, 1))) _GETVALUE;
            }
        }

        ////////////////////////
        // stage 13 Seventh Order
        ////////////////////////
        void RungeKutta8::stage13seventhOrder(const math::Matrix<doubleType> & state_left,
                                              const math::Matrix<double> & dstate_leftdProp_vars,
                                              const math::Matrix<doubleType> & control,
                                              const doubleType & step_size,
                                              const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            // Seventh order does not require a new gradient information (i.e. an integrand call)

            // 7th order solution
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->x_seventh(k) = state_left(k) +
                                     this->b1lower*step_size*this->f1(k) +
                                     this->b6lower*step_size*this->f6(k) +
                                     this->b7lower*step_size*this->f7(k) +
                                     this->b8lower*step_size*this->f8(k) +
                                     this->b9lower*step_size*this->f9(k) +
                                     this->b10lower*step_size*this->f10(k) +
                                     this->b11lower*step_size*this->f11(k) +
                                     this->b12lower*step_size*this->f12(k);
            }
        }

        /*
        void RungeKutta8::stage13seventhOrderwSTM(const math::Matrix<doubleType> & state_left,
                                                  const math::Matrix<double> & dstate_leftdProp_vars,
                                                  const math::Matrix<doubleType> & control,
                                                  const doubleType & step_size,
                                                  const bool & needSTM)
        {

            this->stage13seventhOrder(state_left, dstate_leftdProp_vars, control, step_size, needSTM);

            for (size_t k = 0; k < num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->b1lower*(step_size*this->df1dProp_vars(k, 0))
                                         + this->b6lower*(step_size*this->df6dProp_vars(k, 0))
                                         + this->b7lower*(step_size*this->df7dProp_vars(k, 0))
                                         + this->b8lower*(step_size*this->df8dProp_vars(k, 0))
                                         + this->b9lower*(step_size*this->df9dProp_vars(k, 0))
                                         + this->b10lower*(step_size*this->df10dProp_vars(k, 0))
                                         + this->b11lower*(step_size*this->df11dProp_vars(k, 0))
                                         + this->b12lower*(step_size*this->df12dProp_vars(k, 0))) _GETVALUE;

                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->b1lower*(dstep_sizedProp_var*this->f1(k) + step_size*this->df1dProp_vars(k, 1))
                                         + this->b6lower*(dstep_sizedProp_var*this->f6(k) + step_size*this->df6dProp_vars(k, 1))
                                         + this->b7lower*(dstep_sizedProp_var*this->f7(k) + step_size*this->df7dProp_vars(k, 1))
                                         + this->b8lower*(dstep_sizedProp_var*this->f8(k) + step_size*this->df8dProp_vars(k, 1))
                                         + this->b9lower*(dstep_sizedProp_var*this->f9(k) + step_size*this->df9dProp_vars(k, 1))
                                         + this->b10lower*(dstep_sizedProp_var*this->f10(k) + step_size*this->df10dProp_vars(k, 1))
                                         + this->b11lower*(dstep_sizedProp_var*this->f11(k) + step_size*this->df11dProp_vars(k, 1))
                                         + this->b12lower*(dstep_sizedProp_var*this->f12(k) + step_size*this->df12dProp_vars(k, 1))) _GETVALUE;
            }
        }
        */

        ////////////////////////
        // stage 13 eighth Order
        ////////////////////////
        void RungeKutta8::stage13eighthOrder(const math::Matrix<doubleType> & state_left,
                                             const math::Matrix<double> & dstate_leftdProp_vars,
                                             const doubleType & step_size,
                                             const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c13 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c13 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f13, this->df13dProp_vars, needSTM);

            this->stage13eighthOrder_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage13eighthOrderwSTM(const math::Matrix<doubleType> & state_left,
                                                 const math::Matrix<double> & dstate_leftdProp_vars,
                                                 const doubleType & step_size,
                                                 const bool & needSTM)
        {
            this->stage13eighthOrder(state_left, dstate_leftdProp_vars, step_size, needSTM);
            this->stage13eighthOrder_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage13eighthOrder_w_control(const math::Matrix<doubleType> & state_left,
                                                       const math::Matrix<double> & dstate_leftdProp_vars,
                                                       const math::Matrix<doubleType> & control,
                                                       const doubleType & step_size,
                                                       const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;

            this->integrand->setCurrentIndependentVariable(*this->left_hand_independent_variable + this->c13 * step_size);
            this->integrand->setdCurrentIndVardPropVar(*this->dleft_hand_ind_var_dProp_var + this->c13 * dstep_sizedProp_var);
            this->integrand->evaluate(this->y, this->dydProp_vars, this->f13, this->df13dProp_vars, control, needSTM);

            this->stage13eighthOrder_state_update(state_left, step_size, num_states_to_integrate);
        }

        void RungeKutta8::stage13eighthOrderwSTM_w_control(const math::Matrix<doubleType> & state_left,
                                                           const math::Matrix<double> & dstate_leftdProp_vars,
                                                           const math::Matrix<doubleType> & control,
                                                           const doubleType & step_size,
                                                           const bool & needSTM)
        {
            this->stage13eighthOrder_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage13eighthOrder_propvar_update(dstate_leftdProp_vars, step_size);
        }

        void RungeKutta8::stage13eighthOrder_state_update(const math::Matrix<doubleType> & state_left,
                                                          const doubleType & step_size,
                                                          const size_t & num_states_to_integrate)
        {
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                this->y(k) = state_left(k) +
                             this->b1upper*step_size*this->f1(k) +
                             this->b6upper*step_size*this->f6(k) +
                             this->b7upper*step_size*this->f7(k) +
                             this->b8upper*step_size*this->f8(k) +
                             this->b9upper*step_size*this->f9(k) +
                             this->b10upper*step_size*this->f10(k) +
                             this->b11upper*step_size*this->f11(k) +
                             this->b12upper*step_size*this->f12(k) +
                             this->b13upper*step_size*this->f13(k);
            }
            this->x_eighth = this->y;
        }

        void RungeKutta8::stage13eighthOrder_propvar_update(const math::Matrix<double> & dstate_leftdProp_vars,
                                                            const doubleType & step_size)
        {
            for (size_t k = 0; k < this->num_prop_var_deriv_states; ++k)
            {
                this->dydProp_vars(k, 0) = (dstate_leftdProp_vars(k, 0)
                                         + this->b1upper*(step_size*this->df1dProp_vars(k, 0))
                                         + this->b6upper*(step_size*this->df6dProp_vars(k, 0))
                                         + this->b7upper*(step_size*this->df7dProp_vars(k, 0))
                                         + this->b8upper*(step_size*this->df8dProp_vars(k, 0))
                                         + this->b9upper*(step_size*this->df9dProp_vars(k, 0))
                                         + this->b10upper*(step_size*this->df10dProp_vars(k, 0))
                                         + this->b11upper*(step_size*this->df11dProp_vars(k, 0))
                                         + this->b12upper*(step_size*this->df12dProp_vars(k, 0))
                                         + this->b13upper*(step_size*this->df13dProp_vars(k, 0))) _GETVALUE;

                this->dydProp_vars(k, 1) = (dstate_leftdProp_vars(k, 1)
                                         + this->b1upper*(dstep_sizedProp_var*this->f1(k) + step_size * this->df1dProp_vars(k, 1))
                                         + this->b6upper*(dstep_sizedProp_var*this->f6(k) + step_size * this->df6dProp_vars(k, 1))
                                         + this->b7upper*(dstep_sizedProp_var*this->f7(k) + step_size * this->df7dProp_vars(k, 1))
                                         + this->b8upper*(dstep_sizedProp_var*this->f8(k) + step_size * this->df8dProp_vars(k, 1))
                                         + this->b9upper*(dstep_sizedProp_var*this->f9(k) + step_size * this->df9dProp_vars(k, 1))
                                         + this->b10upper*(dstep_sizedProp_var*this->f10(k) + step_size * this->df10dProp_vars(k, 1))
                                         + this->b11upper*(dstep_sizedProp_var*this->f11(k) + step_size * this->df11dProp_vars(k, 1))
                                         + this->b12upper*(dstep_sizedProp_var*this->f12(k) + step_size * this->df12dProp_vars(k, 1))
                                         + this->b13upper*(dstep_sizedProp_var*this->f13(k) + step_size * this->df13dProp_vars(k, 1))) _GETVALUE;
            }
            this->dx_eighthdProp_vars = this->dydProp_vars;
        }

        void RungeKutta8::step(const math::Matrix<doubleType> & state_left,
                               const math::Matrix<double> & dstate_leftdProp_vars,
                               math::Matrix<doubleType> & state_right,
                               math::Matrix<double> & dstate_rightdProp_vars,
                               const doubleType & step_size,
                               const double & dstep_sizedProp_var,
                               const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;
            this->dstep_sizedProp_var = dstep_sizedProp_var;
            this->integrand->setdCurrentIndVardPropVarPrevious(*this->dleft_hand_ind_var_dProp_var_previous);
            this->num_prop_var_deriv_states = this->dydProp_vars.get_n();

            if (needSTM)
            {
                this->stage1wSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage2wSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage3wSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage4wSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage5wSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage6wSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage7wSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage8wSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage9wSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage10wSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage11wSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage12wSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage13eighthOrderwSTM(state_left, dstate_leftdProp_vars, step_size, needSTM);

            }
            else
            {
                this->stage1(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage2(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage3(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage4(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage5(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage6(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage7(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage8(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage9(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage10(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage11(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage12(state_left, dstate_leftdProp_vars, step_size, needSTM);
                this->stage13eighthOrder(state_left, dstate_leftdProp_vars, step_size, needSTM);
            }

            // pass the final answers to the outbound containers
            dstate_rightdProp_vars = this->dx_eighthdProp_vars;
            state_right = this->x_eighth;
        }

        void RungeKutta8::step_w_control(const math::Matrix<doubleType> & state_left,
                                         const math::Matrix<double> & dstate_leftdProp_vars,
                                         math::Matrix<doubleType> & state_right,
                                         math::Matrix<double> & dstate_rightdProp_vars,
                                         const math::Matrix <doubleType> & control,
                                         const doubleType & step_size,
                                         const double & dstep_sizedProp_var,
                                         const bool & needSTM)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;
            this->dstep_sizedProp_var = dstep_sizedProp_var;
            this->integrand->setdCurrentIndVardPropVarPrevious(*this->dleft_hand_ind_var_dProp_var_previous);
            this->num_prop_var_deriv_states = this->dydProp_vars.get_n();

            if (needSTM)
            {
                this->stage1wSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage2wSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage3wSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage4wSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage5wSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage6wSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage7wSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage8wSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage9wSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage10wSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage11wSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage12wSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage13eighthOrderwSTM_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);

            }
            else
            {
                this->stage1_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage2_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage3_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage4_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage5_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage6_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage7_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage8_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage9_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage10_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage11_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage12_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
                this->stage13eighthOrder_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            }

            // pass the final answers to the outbound containers
            dstate_rightdProp_vars = this->dx_eighthdProp_vars;
            state_right = this->x_eighth;
        }

        void RungeKutta8::computeError(doubleType & error, math::Matrix<double> & error_scaling_factors)
        {
            error = 0.0;
            size_t & num_states_to_integrate = *this->num_states_to_integrate;
            doubleType current_error = 0.0;
            for (size_t k = 0; k < num_states_to_integrate; ++k)
            {
                current_error = fabs(x_eighth(k) - x_seventh(k)) * error_scaling_factors(k);
                if (current_error > error)
                {
                    error = current_error;
                }
            }
        }

        void RungeKutta8::errorControlledStep(const math::Matrix<doubleType> & state_left,
            const math::Matrix<double> & dstate_leftdProp_vars,
            math::Matrix<doubleType> & state_right,
            math::Matrix<double> & dstate_rightdProp_vars,
            const math::Matrix <doubleType> & control,
            const doubleType & step_size,
            const double & dstep_sizedProp_var,
            const bool & needSTM,
            doubleType & error,
            math::Matrix<double> & error_scaling_factors)
        {
            size_t & num_states_to_integrate = *this->num_states_to_integrate;
            this->dstep_sizedProp_var = dstep_sizedProp_var;
            this->integrand->setdCurrentIndVardPropVarPrevious(*this->dleft_hand_ind_var_dProp_var_previous);
            this->num_prop_var_deriv_states = this->dydProp_vars.get_n();

            // We are not computing analytic propagation variable partials for
            // adaptive step RK87 at this time

            //if (needSTM)
            //{
            //    this->stage1wSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage2wSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage3wSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage4wSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage5wSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage6wSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage7wSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage8wSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage9wSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage10wSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage11wSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage12wSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage13seventhOrder(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    //this->stage13seventhOrderwSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //    this->stage13eighthOrderwSTM(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //
            //}
            //else
            //{
            this->stage1_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage2_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage3_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage4_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage5_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage6_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage7_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage8_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage9_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage10_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage11_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage12_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage13seventhOrder(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            this->stage13eighthOrder_w_control(state_left, dstate_leftdProp_vars, control, step_size, needSTM);
            //}

            // compute the relative error between the 7th and 8th order solutions for this RK step
            this->computeError(error, error_scaling_factors);

            // pass the final answers to the outbound containers
            dstate_rightdProp_vars = this->dx_eighthdProp_vars;
            state_right = this->x_eighth;
        }

    } // end integration namespace
} // end EMTG namespace