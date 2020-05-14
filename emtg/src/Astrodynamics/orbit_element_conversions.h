// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2014 - 2016 United States Government as represented by the
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

//header file for orbit element conversions
//Jacob Englander 12-29-2015
//templated so that derivatives may be checked using algorithmic differentiation



#pragma once
#include "EMTG_math.h"
#include "EMTG_Matrix.h"

#include<vector>

namespace EMTG
{
    namespace Astrodynamics
    {
        // template<class T> void inertial2COE(const std::vector<T>& stateCartesianVector,//6x1
            // const double& mu,
            // std::vector<T>& stateCOEVector,//6x1
            // const bool& GenerateDerivatives,
            // math::Matrix<T>& COE_derivatives)//6x6
        // {
            // static math::Matrix<T> stateCartesian(6, 1);
            // static math::Matrix<T> stateCOE(6, 1);

            // stateCartesian.assign_all(stateCartesianVector);
            // stateCOE.assign_all(stateCOEVector);
            // inertial2COE(stateCartesian, mu, stateCOE, GenerateDerivatives, COE_derivatives);
        // }
        
        template<class T> void inertial2COE(const math::Matrix<T>& stateCartesian,//6x1
                                            const double& mu, 
                                            math::Matrix<T>& stateCOE,//6x1
                                            const bool& GenerateDerivatives,
                                            math::Matrix<T>& COE_derivatives)//6x6
        {
            T r, v, h, n, rdotv, s;
            math::Matrix<T> evec(3, 1), nvec(3, 1), hvec(3, 1);
            math::Matrix<T> ihat(3, 1, { 1, 0, 0 });
            math::Matrix<T> jhat(3, 1, { 0, 1, 0 });
            math::Matrix<T> khat(3, 1, { 0, 0, 1 });

            math::Matrix<T> R = stateCartesian.getSubMatrix1D(0, 2);
            math::Matrix<T> V = stateCartesian.getSubMatrix1D(3, 5);
            r = R.norm();
            v = V.norm();

            //SMA
            stateCOE(0) = r / (2 - r*v*v / mu);

            //eccentricity vector
            rdotv = R.dot(V);
            s = (v*v - mu / r);

            evec(0) = 1 / mu * (s*stateCartesian(0) - rdotv*stateCartesian(3));
            evec(1) = 1 / mu * (s*stateCartesian(1) - rdotv*stateCartesian(4));
            evec(2) = 1 / mu * (s*stateCartesian(2) - rdotv*stateCartesian(5));

            //ECC
            T ECC = evec.norm();
            stateCOE(1) = ECC;

            //angular momentum vector and scalar
            hvec = R.cross(V);
            h = hvec.norm();

            //INC
            stateCOE(2) = acos(hvec(2) / h);

            //nodal vector
            nvec = khat.cross(hvec) / h;

            n = nvec.norm();

            //RAAN

            if (nvec(1) >= 0)
                stateCOE(3) = acos(nvec(0) / n);
            else
                stateCOE(3) = 2 * EMTG::math::PI - acos(nvec(0) / n);

            if (n == 0)
                stateCOE(3) = 0;

            //AOP
            T ndote = nvec.dot(evec);
            if (evec(2) >= 0)
                stateCOE(4) = acos(ndote / (n * ECC));
            else
                stateCOE(4) = 2 * EMTG::math::PI - acos(ndote / (n * ECC));

            if (n == 0) //if no inclination, then eccentricity vector points to the periapse
                stateCOE(4) = atan2(evec(1), evec(0));

            //TA
            T edotR = evec.dot(R);
            if (rdotv >= 0)
                stateCOE(5) = acos(edotR / (r * ECC));
            else
                stateCOE(5) = 2 * EMTG::math::PI - acos(edotR / (r * ECC));

            //derivatives
            if (GenerateDerivatives)
            {
                T R = r _GETVALUE;
                T V = v _GETVALUE;
                T x = stateCartesian(0)_GETVALUE;
                T y = stateCartesian(1)_GETVALUE;
                T z = stateCartesian(2)_GETVALUE;
                T xdot = stateCartesian(3)_GETVALUE;
                T ydot = stateCartesian(4)_GETVALUE;
                T zdot = stateCartesian(5)_GETVALUE;

                static math::Matrix<T> dr_drstate(3, 1);
                static math::Matrix<T> dv_dvstate(3, 1);
                dr_drstate(0) = x / R;
                dr_drstate(1) = y / R;
                dr_drstate(2) = z / R;
                dv_dvstate(0) = xdot / V;
                dv_dvstate(1) = ydot / V;
                dv_dvstate(2) = zdot / V;

                //derivatives of SMA
                T dSMAdr = 2 * mu * mu / (R*V*V - 2 * mu) / (R*V*V - 2 * mu);
                COE_derivatives(0, 0) = dSMAdr * dr_drstate(0);
                COE_derivatives(1, 0) = dSMAdr * dr_drstate(1);
                COE_derivatives(2, 0) = dSMAdr * dr_drstate(2);
                T dSMAdv = 2 * mu * R*R*V / (R*V*V - 2 * mu) / (R*V*V - 2 * mu);
                COE_derivatives(3, 0) = dSMAdv * dv_dvstate(0);
                COE_derivatives(4, 0) = dSMAdv * dv_dvstate(1);
                COE_derivatives(5, 0) = dSMAdv * dv_dvstate(2);

                //derivatives of ECC
                //first need derivatives of rdotv and s
                static math::Matrix<T> drdotv_dstate(6, 1);
                drdotv_dstate(0) = xdot;
                drdotv_dstate(1) = ydot;
                drdotv_dstate(2) = zdot;
                drdotv_dstate(3) = x;
                drdotv_dstate(4) = y;
                drdotv_dstate(5) = z;
                T dsdv = 2 * V;
                T dsdr = mu / (R * R);
                static math::Matrix<T> ds_dstate(6, 1);
                ds_dstate(0) = dsdr * dr_drstate(0);
                ds_dstate(1) = dsdr * dr_drstate(1);
                ds_dstate(2) = dsdr * dr_drstate(2);
                ds_dstate(3) = dsdv * dv_dvstate(0);
                ds_dstate(4) = dsdv * dv_dvstate(1);
                ds_dstate(5) = dsdv * dv_dvstate(2);

                //derivatives of ECC vector
                static math::Matrix<T> devec_dstate(6, 3);
                //with respect to position
                for (size_t i = 0; i < 3; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                    {
                        devec_dstate(i, j) = 1 / mu * (ds_dstate(i) * stateCartesian(j) + (i == j ? s : 0) - drdotv_dstate(i) * stateCartesian(j + 3))_GETVALUE;
                    }
                }
                //with respect to velocity
                for (size_t i = 0; i < 3; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                    {
                        devec_dstate(i + 3, j) = 1 / mu * (ds_dstate(i + 3) * stateCartesian(j) - drdotv_dstate(i + 3) * stateCartesian(j + 3) - (i == j ? rdotv : 0))_GETVALUE;
                    }
                }

#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                //check of eccentricity vector derivatives, do be removed
                static math::Matrix<T> devec_dstate_algorithmic(6, 3);
                for (size_t i = 0; i < 6; ++i)
                    for (size_t j = 0; j < 3; ++j)
                        devec_dstate_algorithmic(i, j) = evec[j] _GETDERIVATIVE;

                (devec_dstate - devec_dstate_algorithmic).element_divide(devec_dstate_algorithmic).print_to_file("devec_dstate_relative_error.txt");
#endif

                //derivatives of ECC scalar
                for (size_t i = 0; i < 6; ++i)
                {
                    COE_derivatives(i, 1) = 0.0;
                    for (size_t j = 0; j < 3; ++j)
                        COE_derivatives(i, 1) += (evec(j) / stateCOE(1))_GETVALUE * devec_dstate(i, j);
                }

                //derivatives of INC
                //first need the derivatives of the h-vector
                static math::Matrix<T> dhvec_dstate(6, 3);
                dhvec_dstate(0, 0) = 0.0;
                dhvec_dstate(0, 1) = -zdot;
                dhvec_dstate(0, 2) = ydot;
                dhvec_dstate(1, 0) = zdot;
                dhvec_dstate(1, 1) = 0.0;
                dhvec_dstate(1, 2) = -xdot;
                dhvec_dstate(2, 0) = -ydot;
                dhvec_dstate(2, 1) = xdot;
                dhvec_dstate(2, 2) = 0.0;
                dhvec_dstate(3, 0) = 0.0;
                dhvec_dstate(3, 1) = z;
                dhvec_dstate(3, 2) = -y;
                dhvec_dstate(4, 0) = -z;
                dhvec_dstate(4, 1) = 0.0;
                dhvec_dstate(4, 2) = x;
                dhvec_dstate(5, 0) = y;
                dhvec_dstate(5, 1) = -x;
                dhvec_dstate(5, 2) = 0.0;

#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                //quick check of h-vector derivatives, to be removed later
                math::Matrix<T> h_vector_derivatives_algorthmic(6, 3);
                for (size_t i = 0; i < 6; ++i)
                    for (size_t j = 0; j < 3; ++j)
                        h_vector_derivatives_algorthmic(i, j) = hvec[j] _GETDERIVATIVE;

                (dhvec_dstate - h_vector_derivatives_algorthmic).element_divide(h_vector_derivatives_algorthmic).print_to_file("dhvec_dstate_relative_error.txt");
#endif

                //derivatives of h scalar
                static math::Matrix<T> dhscalar_dstate(6, 1);
                for (size_t i = 0; i < 6; ++i)
                {
                    dhscalar_dstate(i) = 0.0;
                    for (size_t j = 0; j < 3; ++j)
                        dhscalar_dstate(i) += (hvec(j) / h)_GETVALUE * dhvec_dstate(i, j);
                }

#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                //quick check of h-scalar derivatives, to be removed later
                static math::Matrix<T> dhscalar_dstate_algorithmic(6, 1);
                for (size_t i = 0; i < 6; ++i)
                    dhscalar_dstate_algorithmic(i) = h _GETDERIVATIVE;
                (dhscalar_dstate - dhscalar_dstate_algorithmic).element_divide(dhscalar_dstate_algorithmic).print_to_file("dhscalar_dstate_relative_error.txt");
#endif

                //derivatives of INC
                T INCterm = (hvec(2) / h)_GETVALUE;
                T diffINCterm = -1.0 / sqrt(1 - INCterm*INCterm) / (h*h)_GETVALUE;
                for (size_t i = 0; i < 6; ++i)
                {
                    COE_derivatives(i, 2) = diffINCterm * (dhvec_dstate(i, 2) * h - dhscalar_dstate(i) * hvec(2))_GETVALUE;
                }

                //derivatives of RAAN
                //we need derivatives of the nodal vector
                static math::Matrix<T> dnvec_dstate(6, 3, 0.0);
                static math::Matrix<T> dkcrossh_dstate(6, 3, 0.0);
                dkcrossh_dstate(0, 0) = zdot;
                dkcrossh_dstate(1, 1) = zdot;
                dkcrossh_dstate(2, 0) = -xdot;
                dkcrossh_dstate(2, 1) = -ydot;
                dkcrossh_dstate(3, 0) = -z;
                dkcrossh_dstate(4, 1) = -z;
                dkcrossh_dstate(5, 0) = x;
                dkcrossh_dstate(5, 1) = y;
                for (size_t i = 0; i < 6; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                    {
                        dnvec_dstate(i, j) = ((dkcrossh_dstate(i, j) * h - dhscalar_dstate(i) * nvec(j) * h) / (h*h))_GETVALUE;
                    }
                }

#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                //quick check of n-vector derivatives, to be removed later
                math::Matrix<T> n_vector_derivatives_algorthmic(6, 3);
                for (size_t i = 0; i < 6; ++i)
                    for (size_t j = 0; j < 3; ++j)
                        n_vector_derivatives_algorthmic(i, j) = nvec[j] _GETDERIVATIVE;

                (dnvec_dstate - n_vector_derivatives_algorthmic).element_divide(n_vector_derivatives_algorthmic).print_to_file("dnvec_dstate_relative_error.txt");
#endif

                //derivatives of n scalar
                static math::Matrix<T> dnscalar_dstate(6, 1);
                for (size_t i = 0; i < 6; ++i)
                {
                    dnscalar_dstate(i) = 0.0;
                    for (size_t j = 0; j < 3; ++j)
                        dnscalar_dstate(i) += (nvec(j) / n)_GETVALUE * dnvec_dstate(i, j);
                }

#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                //quick check of n-scalar derivatives, to be removed later
                static math::Matrix<T> dnscalar_dstate_algorithmic(6, 1);
                for (size_t i = 0; i < 6; ++i)
                    dnscalar_dstate_algorithmic(i) = n _GETDERIVATIVE;
                (dnscalar_dstate - dnscalar_dstate_algorithmic).element_divide(dnscalar_dstate_algorithmic).print_to_file("dnscalar_dstate_relative_error.txt");
#endif

                //now derivatives of RAAN
                T RAANterm = (nvec(0) / n)_GETVALUE;
                T diffRAANterm = -1.0 / sqrt(1 - RAANterm*RAANterm) / (n*n)_GETVALUE;
                if (nvec(1) >= 0)
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        COE_derivatives(i, 3) = diffRAANterm * (dnvec_dstate(i, 0) * n - dnscalar_dstate(i) * nvec(0))_GETVALUE;
                    }
                }
                else if (n == 0) //pathological crash case
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        COE_derivatives(i, 3) = 0.0;
                    }
                }
                else
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        COE_derivatives(i, 3) = -diffRAANterm * (dnvec_dstate(i, 0) * n - dnscalar_dstate(i) * nvec(0))_GETVALUE;
                    }
                }

                //derivatives of AOP
                T AOPterm = (ndote / n / ECC)_GETVALUE;
                T diffAOPterm = -1.0 / sqrt(1 - AOPterm*AOPterm) / (n*n*ECC*ECC)_GETVALUE;
                if (evec(2) >= 0)
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        T dndote_dstate = 0.0;
                        for (size_t j = 0; j < 3; ++j)
                        {
                            dndote_dstate += (dnvec_dstate(i, j) * evec(j) + devec_dstate(i, j) * nvec(j))_GETVALUE;
                        }
                        COE_derivatives(i, 4) = diffAOPterm * (dndote_dstate * n * ECC - (dnscalar_dstate(i)*ECC + COE_derivatives(i, 1) * n) * ndote)_GETVALUE;
                    }
                }
                else if (n == 0)//if no inclination, then eccentricity vector points to the periapse
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        COE_derivatives(i, 4) = ((-evec(1) * devec_dstate(i, 0) + evec(0) * devec_dstate(i, 1)) / (evec(0) * evec(0) + evec(1) * evec(1)))_GETVALUE;
                    }
                }
                else
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        T dndote_dstate = 0.0;
                        for (size_t j = 0; j < 3; ++j)
                        {
                            dndote_dstate += (dnvec_dstate(i, j) * evec(j) + devec_dstate(i, j) * nvec(j))_GETVALUE;
                        }
                        COE_derivatives(i, 4) = -diffAOPterm * (dndote_dstate * n * ECC - (dnscalar_dstate(i)*ECC + COE_derivatives(i, 1) * n) * ndote)_GETVALUE;
                    }
                }

                //derivatives of true anomaly
                T TAterm = (edotR / r / ECC)_GETVALUE;
                T diffTAterm = -1.0 / sqrt(1 - TAterm*TAterm) / (r*r*ECC*ECC)_GETVALUE;
                if (rdotv >= 0)
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        T dedotstate_dstate = 0.0;
                        for (size_t j = 0; j < 3; ++j)
                        {
                            dedotstate_dstate += (devec_dstate(i, j) * stateCartesian(j) + (i == j ? 1 : 0) * evec(j))_GETVALUE;
                        }
                        COE_derivatives(i, 5) = diffTAterm * (dedotstate_dstate * r * ECC - (COE_derivatives(i, 1)*r + (i < 3 ? dr_drstate(i) : 0.0) * ECC) * edotR)_GETVALUE;
                    }
                }
                else
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        T dedotstate_dstate = 0.0;
                        for (size_t j = 0; j < 3; ++j)
                        {
                            dedotstate_dstate += (devec_dstate(i, j) * stateCartesian(j) + (i == j ? 1 : 0) * evec(j))_GETVALUE;
                        }
                        COE_derivatives(i, 5) = -diffTAterm * (dedotstate_dstate * r * ECC - (COE_derivatives(i, 1)*r + (i < 3 ? dr_drstate(i) : 0.0) * ECC) * edotR)_GETVALUE;
                    }
                }
                    
            }//end code to generate derivatives of COE with respect to inertial
        }//end intertial2COE

        template<class T> void inertial2COE(const math::Matrix<T>& stateCartesian,//6x1
                                            const double& mu, 
                                            math::Matrix<T>& stateCOE)
        {
            static math::Matrix<T> dummy_derivatives(1, 1, 0.0);
            
            inertial2COE(stateCartesian,//6x1
                         mu, 
                         stateCOE,//6x1
                         false,
                         dummy_derivatives);//6x6
        }
        

        template<class T> void COE2inertial(const math::Matrix<T>& E_COE,//6x1
                                                    const double& mu,
                                                    math::Matrix<T>& state,//6x1
                                                    const bool& GenerateDerivatives = false,
                                                    math::Matrix<T>& inertial_derivatives = math::Matrix<T>(1, 1, 0.0))//6x6
        {
            T SMA = E_COE(0);
            T ECC = E_COE(1);
            T INC = E_COE(2);
            T RAAN = E_COE(3);
            T AOP = E_COE(4);
            T TA = E_COE(5);
            T THETA = AOP + TA;

            T cosTA = cos(TA);
            T cosINC = cos(INC);
            T sinINC = sin(INC);
            T cosRAAN = cos(RAAN);
            T sinRAAN = sin(RAAN);
            T cosAOP = cos(AOP);
            T sinAOP = sin(AOP);
            T cosTHETA = cos(THETA);
            T sinTHETA = sin(THETA);

            T p = SMA * (1 - ECC*ECC);

            if (p < 1.0e-30)
            {
                throw std::runtime_error("Error converting parabolic orbit to Cartesian coordinates. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            T r = p / (1 + ECC * cosTA);

            T h = sqrt(mu * SMA * (1 - ECC*ECC));

            //position
            state(0) = r * (cosRAAN*cosTHETA - sinRAAN*sinTHETA*cosINC);
            state(1) = r * (sinRAAN*cosTHETA + cosRAAN*sinTHETA*cosINC);
            state(2) = r * sinTHETA * sinINC;

            //velocity
            state(3) = -mu / h * (cosRAAN * (sinTHETA + ECC*sinAOP) + sinRAAN * (cosTHETA + ECC*cosAOP) * cosINC);
            state(4) = -mu / h * (sinRAAN * (sinTHETA + ECC*sinAOP) - cosRAAN * (cosTHETA + ECC*cosAOP) * cosINC);
            state(5) = mu / h * (cosTHETA + ECC*cosAOP) * sinINC;

            if (GenerateDerivatives)
            {
                //first we need derivatives of r and h with respect to their components
                T dr_dSMA = (-(ECC * ECC - 1) / (ECC*cosTA + 1))_GETVALUE;
                T dr_dECC = (-(SMA*(cosTA * ECC * ECC + 2 * ECC + cosTA)) / (ECC*cosTA + 1) / (ECC*cosTA + 1))_GETVALUE;
                T dr_dTA = (-(ECC*SMA*sin(TA)*(ECC*ECC - 1)) / (ECC*cosTA + 1) / (ECC*cosTA + 1))_GETVALUE;
                T dh_dSMA = (-(mu*(ECC*ECC - 1)) / (2 * sqrt(-SMA*mu*(ECC*ECC - 1))))_GETVALUE;
                T dh_dECC = (-(ECC*SMA*mu) / sqrt(-SMA*mu*(ECC *ECC - 1)))_GETVALUE;

                //derivatives of x
                inertial_derivatives(0, 0) = dr_dSMA * (cosRAAN*cosTHETA - sinRAAN*sinTHETA*cosINC)_GETVALUE;
                inertial_derivatives(1, 0) = dr_dECC * (cosRAAN*cosTHETA - sinRAAN*sinTHETA*cosINC)_GETVALUE;
                inertial_derivatives(2, 0) = (r * sinRAAN*sinTHETA*sinINC)_GETVALUE;
                inertial_derivatives(3, 0) = (r * (-sinRAAN*cosTHETA - cosRAAN*sinTHETA*cosINC))_GETVALUE;
                inertial_derivatives(4, 0) = (r * (-cosRAAN*sinTHETA - sinRAAN*cosTHETA*cosINC))_GETVALUE;
                inertial_derivatives(5, 0) = dr_dTA * (cosRAAN*cosTHETA - sinRAAN*sinTHETA*cosINC)_GETVALUE + (r * (-cosRAAN*sinTHETA - sinRAAN*cosTHETA*cosINC))_GETVALUE;

                //derivatives of y
                inertial_derivatives(0, 1) = dr_dSMA * (sinRAAN*cosTHETA + cosRAAN*sinTHETA*cosINC)_GETVALUE;
                inertial_derivatives(1, 1) = dr_dECC * (sinRAAN*cosTHETA + cosRAAN*sinTHETA*cosINC)_GETVALUE;
                inertial_derivatives(2, 1) = (r * -cosRAAN*sinTHETA*sinINC)_GETVALUE;
                inertial_derivatives(3, 1) = (r * (cosRAAN*cosTHETA - sinRAAN*sinTHETA*cosINC))_GETVALUE;
                inertial_derivatives(4, 1) = (r * (-sinRAAN*sinTHETA + cosRAAN*cosTHETA*cosINC))_GETVALUE;
                inertial_derivatives(5, 1) = dr_dTA * (sinRAAN*cosTHETA + cosRAAN*sinTHETA*cosINC)_GETVALUE + (r * (-sinRAAN*sinTHETA + cosRAAN*cosTHETA*cosINC))_GETVALUE;

                //derivatives of z
                inertial_derivatives(0, 2) = dr_dSMA * (sinTHETA * sinINC)_GETVALUE;
                inertial_derivatives(1, 2) = dr_dECC * (sinTHETA * sinINC)_GETVALUE;
                inertial_derivatives(2, 2) = (r * sinTHETA * cosINC)_GETVALUE;
                inertial_derivatives(3, 2) = 0.0;
                inertial_derivatives(4, 2) = (r * cosTHETA * sinINC)_GETVALUE;
                inertial_derivatives(5, 2) = dr_dTA * (sinTHETA * sinINC)_GETVALUE + (r * (cosTHETA * sinINC))_GETVALUE;

                //derivatives of xdot
                inertial_derivatives(0, 3) = (-dh_dSMA / (h*h) * -mu * (cosRAAN * (sinTHETA + ECC*sinAOP) + sinRAAN * (cosTHETA + ECC*cosAOP) * cosINC))_GETVALUE;
                inertial_derivatives(1, 3) = (-dh_dECC / (h*h) * -mu * (cosRAAN * (sinTHETA + ECC*sinAOP) + sinRAAN * (cosTHETA + ECC*cosAOP) * cosINC) 
                    - mu / h * (cosRAAN * sinAOP + sinRAAN * cosAOP * cosINC))_GETVALUE;
                inertial_derivatives(2, 3) = (-mu / h * (sinRAAN * (cosTHETA + ECC*cosAOP) * -sinINC))_GETVALUE;
                inertial_derivatives(3, 3) = (-mu / h * (-sinRAAN * (sinTHETA + ECC*sinAOP) + cosRAAN * (cosTHETA + ECC*cosAOP) * cosINC))_GETVALUE;
                inertial_derivatives(4, 3) = (-mu / h * (cosRAAN * (cosTHETA + ECC*cosAOP) + sinRAAN * (-sinTHETA + ECC*-sinAOP) * cosINC))_GETVALUE;
                inertial_derivatives(5, 3) = (-mu / h * (cosRAAN * (cosTHETA) + sinRAAN * (-sinTHETA) * cosINC))_GETVALUE;

                //derivatives of ydot
                inertial_derivatives(0, 4) = (-dh_dSMA / (h*h) *-mu * (sinRAAN * (sinTHETA + ECC*sinAOP) - cosRAAN * (cosTHETA + ECC*cosAOP) * cosINC))_GETVALUE;
                inertial_derivatives(1, 4) = (-dh_dECC / (h*h) *-mu * (sinRAAN * (sinTHETA + ECC*sinAOP) - cosRAAN * (cosTHETA + ECC*cosAOP) * cosINC)
                    - mu / h * (sinRAAN * sinAOP - cosRAAN * cosAOP * cosINC))_GETVALUE;
                inertial_derivatives(2, 4) = (-mu / h * (- cosRAAN * (cosTHETA + ECC*cosAOP) * -sinINC))_GETVALUE;
                inertial_derivatives(3, 4) = (-mu / h * (cosRAAN * (sinTHETA + ECC*sinAOP) + sinRAAN * (cosTHETA + ECC*cosAOP) * cosINC))_GETVALUE;
                inertial_derivatives(4, 4) = (-mu / h * (sinRAAN * (cosTHETA + ECC*cosAOP) - cosRAAN * (-sinTHETA + ECC*-sinAOP) * cosINC))_GETVALUE;
                inertial_derivatives(5, 4) = (-mu / h * (sinRAAN * (cosTHETA) - cosRAAN * (-sinTHETA) * cosINC))_GETVALUE;

                //derivatives of zdot
                inertial_derivatives(0, 5) = (-dh_dSMA / (h*h) * mu  * (cosTHETA + ECC*cosAOP) * sinINC)_GETVALUE;
                inertial_derivatives(1, 5) = (-dh_dECC / (h*h) * mu  * (cosTHETA + ECC*cosAOP) * sinINC + mu / h * cosAOP * sinINC)_GETVALUE;
                inertial_derivatives(2, 5) = (mu / h * (cosTHETA + ECC*cosAOP) * cosINC)_GETVALUE;
                inertial_derivatives(3, 5) = 0.0;
                inertial_derivatives(4, 5) = (mu / h * (-sinTHETA - ECC*sinAOP) * sinINC)_GETVALUE;
                inertial_derivatives(5, 5) = (mu / h * (-sinTHETA) * sinINC)_GETVALUE;
            }//end derivative code for COE2inertial
        }//end COE2inertial

        //template<class T> void COE2inertial(const std::vector<T>& E_COE,//6x1
        //    const double& mu,
        //    std::vector<T>& state)//6x1)
        //{
        //    static math::Matrix<double> COE_derivatives;
        //    COE2inertial(E_COE, mu, state, false, COE_derivatives);
        //}

        template<class T> void IncomingAsymptote2COE(const std::vector<T>& IncomingAsymptote,//6x1
            const double& mu,
            std::vector<T>& E_COE,//6x1
            const bool& GenerateDerivatives,
            math::Matrix<double>& COE_derivatives)//6x6
        {
            T RP = IncomingAsymptote[0];
            T C3 = IncomingAsymptote[1];
            T RHA = IncomingAsymptote[2];
            T DHA = IncomingAsymptote[3];
            T BVA = IncomingAsymptote[4];
            T TA = IncomingAsymptote[5];

            T cosRHA = cos(RHA);
            T sinRHA = sin(RHA);
            T cosDHA = cos(DHA);
            T sinDHA = sin(DHA);
            T cosBVA = cos(BVA);
            T sinBVA = sin(BVA);
            T cosTA = cos(TA);
            T sinTA = sin(TA);

            if (C3 < 1e-7)
            {
                std::cout << "Warning, IncomingAsymptote orbit is elliptic so using apsides vector for asymptote, C3 < 1e-7" << std::endl;
            }
            
            T SMA = -mu / C3;
            T ECC = 1 - RP / SMA;

            if (ECC < 1e-7)
            {
                throw std::runtime_error("Illegal conversion from IncomingAsymptote to Keplerian elements, ECC < 1e-7. IncomingAsymptote is undefined on a circular orbit. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            static math::Matrix<T> sVhat(3, 1);
            sVhat(0) = cosDHA*cosRHA;
            sVhat(1) = cosDHA*sinRHA;
            sVhat(2) = sinDHA;
            T khati[] = { 0.0,0.0,1.0 };
            static math::Matrix<T> khat(3, 1, khati);

            if (1.0 - sVhat.dot(khat) < 1.0e-7)
            {
                throw std::runtime_error("Illegal conversion from IncomingAsymptote to Keplerian elements, 1.0 - vhat.dot(khat) < 1e-7. IncomingAsymptote is undefined when asymptote is aligned with khat. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            static math::Matrix<T> eaVec = khat.cross(sVhat);
            static math::Matrix<T> eaVhat = eaVec.unitize();
            static math::Matrix<T> noVhat = sVhat.cross(eaVhat);
            T AMI = math::TwoPI / 4 - BVA; // angular azimuth at infinity 
            T cosAMI = cos(AMI);
            T sinAMI = sin(AMI);
            static math::Matrix<T> hVhat = eaVhat*sinAMI + noVhat*cosAMI;
            static math::Matrix<T> nodeVec = khat.cross(hVhat);
            T nMag = nodeVec.norm();
            static math::Matrix<T> eccVhat;
            T TAmax;
            T sinTAmax;
            T cosTAmax;
            static math::Matrix<T> oVhat;

            if (C3 <= -1.0e-7)
            {
                eccVhat = -sVhat;
            }
            else if (C3 >= 1.0e-7)
            {
                TAmax = acos(-1 / ECC);
                oVhat = hVhat.cross(sVhat);
                sinTAmax = sin(TAmax);
                cosTAmax = cos(TAmax);
                eccVhat = oVhat * sinTAmax + sVhat * cosTAmax;
            }

            T khatdothVhat = khat.dot(hVhat);
            T INC = acos(khatdothVhat);
            T RAAN;
            T AOP;

            if (ECC >= 1E-11 && INC >= 1E-11)  // CASE 1: Non-circular, Inclined Orbit
            {
                if (nMag == 0.0)
                {
                    throw std::runtime_error("Cannot convert from Incoming asymptote elements to Cartesian elements - line-of-nodes vector is a zero vector.");
                }
                RAAN = acos(nodeVec(0) / nMag);
                if (nodeVec(1) < 0)
                    RAAN = math::TwoPI - RAAN;

                AOP = acos((nodeVec.dot(eccVhat)) / nMag);
                if (eccVhat(2) < 0)
                    AOP = math::TwoPI - AOP;
            }
            else if (ECC >= 1E-11 && INC < 1E-7)  // CASE 2: Non-circular, Equatorial Orbit
            {
                RAAN = 0;
                AOP = acos(eccVhat(1));
                if (eccVhat(1) < 0)
                    AOP = math::TwoPI - AOP;

            }
            else if (ECC > 1E-11 && INC >= math::PI - 1E-7)  // CASE 3: Non-circular, Equatorial Retrograde Orbit
            {
                RAAN = 0;
                AOP = -acos(eccVhat(0));
                if (eccVhat(1) < 0)
                    AOP = math::TwoPI - AOP;
            }

            //assign values
            E_COE[0] = SMA;
            E_COE[1] = ECC;
            E_COE[2] = INC;
            E_COE[3] = RAAN;
            E_COE[4] = AOP;
            E_COE[5] = TA;

            //compute derivatives of the IncomingAsymptote to COE conversion
            if (GenerateDerivatives)
            {
                //SMA
                COE_derivatives(0, 0) = 0.0;
                COE_derivatives(1, 0) = mu / (C3 * C3)_GETVALUE;
                COE_derivatives(2, 0) = 0.0;
                COE_derivatives(3, 0) = 0.0;
                COE_derivatives(4, 0) = 0.0;
                COE_derivatives(5, 0) = 0.0;

                //ECC
                COE_derivatives(0, 1) = (C3 / mu)_GETVALUE;
                COE_derivatives(1, 1) = (RP / mu)_GETVALUE;
                COE_derivatives(2, 1) = 0.0;
                COE_derivatives(3, 1) = 0.0;
                COE_derivatives(4, 1) = 0.0;
                COE_derivatives(5, 1) = 0.0;

                //INC
                //first we need derivatives of eaVec
                static math::Matrix<double> deaVec_dIncomingAsymptote(6, 3, 0.0);
                //no dependencies of eaVec on RP and C3
                deaVec_dIncomingAsymptote(2, 0) = (-cosDHA * cosRHA)_GETVALUE;
                deaVec_dIncomingAsymptote(2, 1) = (-cosDHA * sinRHA)_GETVALUE;
                deaVec_dIncomingAsymptote(3, 0) = (sinDHA * sinRHA)_GETVALUE;
                deaVec_dIncomingAsymptote(3, 1) = (-cosRHA * sinDHA)_GETVALUE;
                //no dependencies of eaVec on BVA or TA

                //then derivatives of eaVhat
                T ea = eaVec.norm();
                double deascalar_dRHA = 0;
                double deascalar_dDHA = -math::sgn(cosDHA) * sinDHA _GETVALUE;
                
                static math::Matrix<double> deaVhat_dIncomingAsymptote(6, 3, 0.0);
                deaVhat_dIncomingAsymptote(2, 0) = (1 / (ea * ea) * (deaVec_dIncomingAsymptote(2, 0) * ea - deascalar_dRHA * eaVec(0)))_GETVALUE;
                deaVhat_dIncomingAsymptote(2, 1) = (1 / (ea * ea) * (deaVec_dIncomingAsymptote(2, 1) * ea - deascalar_dRHA * eaVec(1)))_GETVALUE;
                deaVhat_dIncomingAsymptote(3, 0) = (1 / (ea * ea) * (deaVec_dIncomingAsymptote(3, 0) * ea - deascalar_dDHA * eaVec(0)))_GETVALUE;
                deaVhat_dIncomingAsymptote(3, 1) = (1 / (ea * ea) * (deaVec_dIncomingAsymptote(3, 1) * ea - deascalar_dDHA * eaVec(1)))_GETVALUE;

                //and derivatives of noVhat
                static math::Matrix<T> noVec = sVhat.cross(eaVec);
                static math::Matrix<double> dnoVec_dRHA(3, 1, 0.0);
                static math::Matrix<double> dnoVec_dDHA(3, 1, 0.0);
                dnoVec_dRHA(0) = (cosDHA * sinDHA * sinRHA)_GETVALUE;
                dnoVec_dRHA(1) = (-cosDHA * cosRHA * sinDHA)_GETVALUE;
                dnoVec_dDHA(0) = (cosRHA * sinDHA * sinDHA - cosDHA*cosDHA * cosRHA)_GETVALUE;
                dnoVec_dDHA(1) = (sinDHA * sinDHA * sinRHA - cosDHA*cosDHA * sinRHA)_GETVALUE;
                dnoVec_dDHA(2) = (-2 * cosDHA * cosRHA * cosRHA * sinDHA - 2 * cosDHA * sinDHA * sinRHA * sinRHA)_GETVALUE;
                static math::Matrix<double> dnoVhat_dIncomingAsymptote(6, 3, 0.0);

                dnoVhat_dIncomingAsymptote(2, 0) = (1 / (ea*ea) * (dnoVec_dRHA(0) * ea - deascalar_dRHA * noVec(0)))_GETVALUE;
                dnoVhat_dIncomingAsymptote(2, 1) = (1 / (ea*ea) * (dnoVec_dRHA(1) * ea - deascalar_dRHA * noVec(1)))_GETVALUE;
                dnoVhat_dIncomingAsymptote(2, 2) = (1 / (ea*ea) * (dnoVec_dRHA(2) * ea - deascalar_dRHA * noVec(2)))_GETVALUE;
                dnoVhat_dIncomingAsymptote(3, 0) = (1 / (ea*ea) * (dnoVec_dDHA(0) * ea - deascalar_dDHA * noVec(0)))_GETVALUE;
                dnoVhat_dIncomingAsymptote(3, 1) = (1 / (ea*ea) * (dnoVec_dDHA(1) * ea - deascalar_dDHA * noVec(1)))_GETVALUE;
                dnoVhat_dIncomingAsymptote(3, 2) = (1 / (ea*ea) * (dnoVec_dDHA(2) * ea - deascalar_dDHA * noVec(2)))_GETVALUE;

                //check on noVhat derivatives
#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                math::Matrix<double> dnoVhat_dIncomingAsymptote_algorithmic(6, 3, 0.0);
                for (size_t ii = 0; ii < 6; ++ii)
                {
                    size_t i = ii + 12;
                    for (size_t j = 0; j < 3; ++j)
                        dnoVhat_dIncomingAsymptote_algorithmic(ii, j) = noVhat(j) _GETDERIVATIVE;
                }
                (dnoVhat_dIncomingAsymptote - dnoVhat_dIncomingAsymptote_algorithmic).element_divide(dnoVhat_dIncomingAsymptote_algorithmic).print_to_file("dnoVhat_dIncomingAsymptote_relative_error.txt");
                dnoVhat_dIncomingAsymptote.print_to_file("dnoVhat_dIncomingAsymptote_analytical.txt");
                dnoVhat_dIncomingAsymptote_algorithmic.print_to_file("dnoVhat_dIncomingAsymptote_algorithmic.txt");
#endif
                //derivatives of AMI
                double dAMI_dBVA = -1.0;

                //derivatives of hVhat
                static math::Matrix<double> dhVhat_dIncomingAsymptote(6, 3, 0.0);
                for (size_t j = 0; j < 3; ++j)
                {
                    dhVhat_dIncomingAsymptote(2, j) = (deaVhat_dIncomingAsymptote(2, j) * sinAMI + dnoVhat_dIncomingAsymptote(2, j) * cosAMI)_GETVALUE;
                    dhVhat_dIncomingAsymptote(3, j) = (deaVhat_dIncomingAsymptote(3, j) * sinAMI + dnoVhat_dIncomingAsymptote(3, j) * cosAMI)_GETVALUE;
                    dhVhat_dIncomingAsymptote(4, j) = ((eaVhat(j) * cosAMI - noVhat(j) * sinAMI) * dAMI_dBVA)_GETVALUE;
                }
#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                math::Matrix<double> dhVhat_dIncomingAsymptote_algorithmic(6, 3, 0.0);
                for (size_t ii = 0; ii < 6; ++ii)
                {
                    size_t i = ii + 12;
                    for (size_t j = 0; j < 3; ++j)
                        dhVhat_dIncomingAsymptote_algorithmic(ii, j) = hVhat(j) _GETDERIVATIVE;
                }
                (dhVhat_dIncomingAsymptote - dhVhat_dIncomingAsymptote_algorithmic).element_divide(dhVhat_dIncomingAsymptote_algorithmic).print_to_file("dhVhat_dIncomingAsymptote_relative_error.txt");
                dhVhat_dIncomingAsymptote.print_to_file("dhVhat_dIncomingAsymptote_analytical.txt");
                dhVhat_dIncomingAsymptote_algorithmic.print_to_file("dhVhat_dIncomingAsymptote_algorithmic.txt");
#endif                

                //now we finally have the ingredients to get derivatives of INC
                double INCmultiplier = -1.0 / sqrt(1.0 - khatdothVhat*khatdothVhat)_GETVALUE;
                COE_derivatives(0, 2) = 0.0;
                COE_derivatives(1, 2) = 0.0;
                COE_derivatives(2, 2) = INCmultiplier * dhVhat_dIncomingAsymptote(2, 2);
                COE_derivatives(3, 2) = INCmultiplier * dhVhat_dIncomingAsymptote(3, 2);
                COE_derivatives(4, 2) = INCmultiplier * dhVhat_dIncomingAsymptote(4, 2);
                COE_derivatives(5, 2) = 0.0;

                //RAAN
                //first we need derivatives of nodeVec
                static math::Matrix<double> dnodeVec_dIncomingAsymptote(6, 3, 0.0);
                for (size_t i = 0; i < 6; ++i)
                {
                    dnodeVec_dIncomingAsymptote(i, 0) = -dhVhat_dIncomingAsymptote(i, 1);
                    dnodeVec_dIncomingAsymptote(i, 1) = dhVhat_dIncomingAsymptote(i, 0);
                    dnodeVec_dIncomingAsymptote(i, 2) = 0.0;
                }
#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                math::Matrix<double> dnodeVec_dIncomingAsymptote_algorithmic(6, 3, 0.0);
                for (size_t ii = 0; ii < 6; ++ii)
                {
                    size_t i = ii + 12;
                    for (size_t j = 0; j < 3; ++j)
                        dnodeVec_dIncomingAsymptote_algorithmic(ii, j) = nodeVec(j) _GETDERIVATIVE;
                }
                (dnodeVec_dIncomingAsymptote - dnodeVec_dIncomingAsymptote_algorithmic).element_divide(dnodeVec_dIncomingAsymptote_algorithmic).print_to_file("dnodeVec_dIncomingAsymptote_relative_error.txt");
                dnodeVec_dIncomingAsymptote.print_to_file("dnodeVec_dIncomingAsymptote_analytical.txt");
                dnodeVec_dIncomingAsymptote_algorithmic.print_to_file("dnodeVec_dIncomingAsymptote_algorithmic.txt");
#endif
                static math::Matrix<double> dnMag_dIncomingAsymptote(6, 1, 0.0);
                for (size_t i = 0; i < 6; ++i)
                {
                    dnMag_dIncomingAsymptote(i) = (1.0 / nMag * (nodeVec(0) * dnodeVec_dIncomingAsymptote(i, 0) + nodeVec(1) * dnodeVec_dIncomingAsymptote(i, 1)))_GETVALUE;
                }

                //now we can compute derivatives of RAAN
                if (ECC >= 1E-11 && INC >= 1E-11)  // CASE 1: Non-circular, Inclined Orbit
                {
                    double RAANterm = -(nodeVec(1) < 0.0 ? -1 : 1) / sqrt(1.0 - (nodeVec(0) / nMag)*(nodeVec(0) / nMag))_GETVALUE;
                    COE_derivatives(0, 3) = 0.0;
                    COE_derivatives(1, 3) = 0.0;
                    COE_derivatives(2, 3) = (RAANterm / (nMag * nMag) * (dnodeVec_dIncomingAsymptote(2, 0) * nMag - dnMag_dIncomingAsymptote(2) * nodeVec(0)))_GETVALUE;
                    COE_derivatives(3, 3) = (RAANterm / (nMag * nMag) * (dnodeVec_dIncomingAsymptote(3, 0) * nMag - dnMag_dIncomingAsymptote(3) * nodeVec(0)))_GETVALUE;
                    COE_derivatives(4, 3) = (RAANterm / (nMag * nMag) * (dnodeVec_dIncomingAsymptote(4, 0) * nMag - dnMag_dIncomingAsymptote(4) * nodeVec(0)))_GETVALUE;
                    COE_derivatives(5, 3) = 0.0;
                }
                else
                {
                    COE_derivatives(0, 3) = 0.0;
                    COE_derivatives(1, 3) = 0.0;
                    COE_derivatives(2, 3) = 0.0;
                    COE_derivatives(3, 3) = 0.0;
                    COE_derivatives(4, 3) = 0.0;
                    COE_derivatives(5, 3) = 0.0;
                }

                //AOP
                //do get AOP derivatives we need derivatives of eccVhat
                static math::Matrix<double> dsVhat_dIncomingAsymptote(6, 3, 0.0);
                dsVhat_dIncomingAsymptote(2, 0) = (-cosDHA*sinRHA)_GETVALUE;
                dsVhat_dIncomingAsymptote(2, 1) = (cosDHA*cosRHA)_GETVALUE;
                dsVhat_dIncomingAsymptote(3, 0) = (-cosRHA*sinDHA)_GETVALUE;
                dsVhat_dIncomingAsymptote(3, 1) = (-sinDHA*sinRHA)_GETVALUE;
                dsVhat_dIncomingAsymptote(3, 2) = (cosDHA)_GETVALUE;
                static math::Matrix<double> deccVhat_dIncomingAsymptote(6, 3, 0.0);
                if (C3 <= -1.0e-7)
                {
                    deccVhat_dIncomingAsymptote(2, 0) = dsVhat_dIncomingAsymptote(2, 0);
                    deccVhat_dIncomingAsymptote(2, 1) = dsVhat_dIncomingAsymptote(2, 1);
                    deccVhat_dIncomingAsymptote(3, 0) = dsVhat_dIncomingAsymptote(3, 0);
                    deccVhat_dIncomingAsymptote(3, 1) = dsVhat_dIncomingAsymptote(3, 1);
                    deccVhat_dIncomingAsymptote(3, 2) = dsVhat_dIncomingAsymptote(3, 2);
                }
                else if (C3 >= 1.0e-7)
                {
                    //need derivatives of oVhat
                    static math::Matrix<double> doVhat_dIncomingAsymptote(6, 3);
                    for (size_t i = 0; i < 6; ++i)
                    {
                        doVhat_dIncomingAsymptote(i, 0) = (dhVhat_dIncomingAsymptote(i, 1) * sVhat(2) + hVhat(1) * dsVhat_dIncomingAsymptote(i, 2) - dhVhat_dIncomingAsymptote(i, 2) * sVhat(1) - hVhat(2) * dsVhat_dIncomingAsymptote(i, 1))_GETVALUE;
                        doVhat_dIncomingAsymptote(i, 1) = (dhVhat_dIncomingAsymptote(i, 2) * sVhat(0) + hVhat(2) * dsVhat_dIncomingAsymptote(i, 0) - dhVhat_dIncomingAsymptote(i, 0) * sVhat(2) - hVhat(0) * dsVhat_dIncomingAsymptote(i, 2))_GETVALUE;
                        doVhat_dIncomingAsymptote(i, 2) = (dhVhat_dIncomingAsymptote(i, 0) * sVhat(1) + hVhat(0) * dsVhat_dIncomingAsymptote(i, 1) - dhVhat_dIncomingAsymptote(i, 1) * sVhat(0) - hVhat(1) * dsVhat_dIncomingAsymptote(i, 0))_GETVALUE;
                    }

                    //need derivatives of TAmax
                    static math::Matrix<double> dTAmax_dIncomingAsymptote(6, 1, 0.0);
                    dTAmax_dIncomingAsymptote(0) = (-C3 / (ECC*ECC * mu*sqrt(1 - 1 / (ECC*ECC))))_GETVALUE;
                    dTAmax_dIncomingAsymptote(1) = (-RP / (ECC*ECC * mu*sqrt(1 - 1 / (ECC*ECC))))_GETVALUE;

                    //now we can compute derivatives of eccVhat
                    for (size_t i = 0; i < 6; ++i)
                    {
                        for (size_t j = 0; j < 3; ++j)
                        {
                            deccVhat_dIncomingAsymptote(i, j) = (doVhat_dIncomingAsymptote(i, j) * sinTAmax
                                + oVhat(j) * cosTAmax * dTAmax_dIncomingAsymptote(i)
                                + dsVhat_dIncomingAsymptote(i, j) * cosTAmax
                                - sVhat(j) * sinTAmax * dTAmax_dIncomingAsymptote(i))_GETVALUE;
                        }
                    }
                }

#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                math::Matrix<double> deccVhat_dIncomingAsymptote_algorithmic(6, 3, 0.0);
                for (size_t ii = 0; ii < 6; ++ii)
                {
                    size_t i = ii + 12;
                    for (size_t j = 0; j < 3; ++j)
                        deccVhat_dIncomingAsymptote_algorithmic(ii, j) = eccVhat(j) _GETDERIVATIVE;
                }
                (deccVhat_dIncomingAsymptote - deccVhat_dIncomingAsymptote_algorithmic).element_divide(deccVhat_dIncomingAsymptote_algorithmic).print_to_file("deccVhat_dIncomingAsymptote_relative_error.txt");
                deccVhat_dIncomingAsymptote.print_to_file("deccVhat_dIncomingAsymptote_analytical.txt");
                deccVhat_dIncomingAsymptote_algorithmic.print_to_file("deccVhat_dIncomingAsymptote_algorithmic.txt");
#endif    
                //and finally the derivatives of AOP
                if (ECC >= 1E-11 && INC >= 1E-11)  // CASE 1: Non-circular, Inclined Orbit
                {
                    T nodeVecdoteccVhat = nodeVec.dot(eccVhat);
                    static math::Matrix<double> dnodeVecdoteccVhat_dIncomingAsymptote(6, 1, 0.0);
                    for (size_t i = 0; i < 6; ++i)
                    {
                        dnodeVecdoteccVhat_dIncomingAsymptote(i) = (nodeVec(0) * deccVhat_dIncomingAsymptote(i, 0) + eccVhat(0) * dnodeVec_dIncomingAsymptote(i, 0)
                            + nodeVec(1) * deccVhat_dIncomingAsymptote(i, 1) + eccVhat(1) * dnodeVec_dIncomingAsymptote(i, 1)
                            + nodeVec(2) * deccVhat_dIncomingAsymptote(i, 2) + eccVhat(2) * dnodeVec_dIncomingAsymptote(i, 2))_GETVALUE;
                    }
                    double AOPterm = -(eccVhat(2) < 0.0 ? -1 : 1) / sqrt(1.0 - (nodeVecdoteccVhat / nMag)*(nodeVecdoteccVhat / nMag))_GETVALUE;
                    COE_derivatives(0, 4) = (AOPterm / (nMag * nMag) * (dnodeVecdoteccVhat_dIncomingAsymptote(0) * nMag - dnMag_dIncomingAsymptote(0) * nodeVecdoteccVhat))_GETVALUE;
                    COE_derivatives(1, 4) = (AOPterm / (nMag * nMag) * (dnodeVecdoteccVhat_dIncomingAsymptote(1) * nMag - dnMag_dIncomingAsymptote(1) * nodeVecdoteccVhat))_GETVALUE;
                    COE_derivatives(2, 4) = (AOPterm / (nMag * nMag) * (dnodeVecdoteccVhat_dIncomingAsymptote(2) * nMag - dnMag_dIncomingAsymptote(2) * nodeVecdoteccVhat))_GETVALUE;
                    COE_derivatives(3, 4) = (AOPterm / (nMag * nMag) * (dnodeVecdoteccVhat_dIncomingAsymptote(3) * nMag - dnMag_dIncomingAsymptote(3) * nodeVecdoteccVhat))_GETVALUE;
                    COE_derivatives(4, 4) = (AOPterm / (nMag * nMag) * (dnodeVecdoteccVhat_dIncomingAsymptote(4) * nMag - dnMag_dIncomingAsymptote(4) * nodeVecdoteccVhat))_GETVALUE;
                    COE_derivatives(5, 4) = 0.0;
                }
                else if (ECC >= 1E-11 && INC < 1E-7)  // CASE 2: Non-circular, Equatorial Orbit
                {
                    double AOPterm = (eccVhat(0) < 0.0 ? -1 : 1) / sqrt(1.0 - (eccVhat(0) * eccVhat(0)))_GETVALUE;
                    COE_derivatives(0, 4) = (AOPterm * deccVhat_dIncomingAsymptote(0, 0));
                    COE_derivatives(1, 4) = (AOPterm * deccVhat_dIncomingAsymptote(1, 0));
                    COE_derivatives(2, 4) = (AOPterm * deccVhat_dIncomingAsymptote(2, 0));
                    COE_derivatives(3, 4) = (AOPterm * deccVhat_dIncomingAsymptote(3, 0));
                    COE_derivatives(4, 4) = (AOPterm * deccVhat_dIncomingAsymptote(4, 0));
                    COE_derivatives(5, 4) = 0.0;

                }
                else if (ECC > 1E-11 && INC >= math::PI - 1E-7)  // CASE 3: Non-circular, Equatorial Retrograde Orbit
                {
                    double AOPterm = -(eccVhat(1) < 0.0 ? -1 : 1) / sqrt(1.0 - (eccVhat(1) * eccVhat(1)))_GETVALUE;
                    COE_derivatives(0, 4) = (AOPterm * deccVhat_dIncomingAsymptote(0, 1));
                    COE_derivatives(1, 4) = (AOPterm * deccVhat_dIncomingAsymptote(1, 1));
                    COE_derivatives(2, 4) = (AOPterm * deccVhat_dIncomingAsymptote(2, 1));
                    COE_derivatives(3, 4) = (AOPterm * deccVhat_dIncomingAsymptote(3, 1));
                    COE_derivatives(4, 4) = (AOPterm * deccVhat_dIncomingAsymptote(4, 1));
                    COE_derivatives(5, 4) = 0.0;
                }

                //TA
                COE_derivatives(0, 5) = 0.0;
                COE_derivatives(1, 5) = 0.0;
                COE_derivatives(2, 5) = 0.0;
                COE_derivatives(3, 5) = 0.0;
                COE_derivatives(4, 5) = 0.0;
                COE_derivatives(5, 5) = 1.0;
            }//end derivatives of the IncomingAsymptote to COE conversion
        }//end IncomingAsymptote2COE

        template<class T> void Inertial2IncomingAsymptote(const std::vector<T>& Inertial_state,//6x1
            const double& mu,
            std::vector<T>& IncomingAsymptote_state,//6x1
            const bool& GenerateDerivatives,
            math::Matrix<double>& IncomingAsymptote_derivatives)//6x6
        {
            // Calculate the magnitude of the position vector, right ascension, and declination
            static math::Matrix<T> R(3, 1, 0.0);
            static math::Matrix<T> V(3, 1, 0.0);
            for (size_t k = 0; k < 3; ++k)
            {
                R(k) = Inertial_state[k];
                V(k) = Inertial_state[k + 3];
            }
            
            T r = R.norm();
            T v = V.norm();
            T rdotv = R.dot(V);

            static math::Matrix<T> hVec = R.cross(V);
            T hMag = hVec.norm();
            T s = (v*v - mu / r);
            static math::Matrix<T> eccVec = (R * s - V * rdotv) / mu;
            T ECC = eccVec.norm();
            T C3 = v*v - 2.0 * mu / r;

            if (fabs(C3) < 1E-7)
            {
                throw std::runtime_error("Error, IncomingAsymptote elements are not defined for nearly parabolic orbits (C3 < 1.0e-7). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            if (v < 1E-7)
            {
                throw std::runtime_error("Error, IncomingAsymptote elements are not defined for near-zero velocity orbits (v < 1.0e-7). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            if (ECC <= 1E-7)
            {
                throw std::runtime_error("Error, IncomingAsymptote elements are not defined for near-zero eccentricity orbits (ECC < 1.0e-7). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            T SMA = -mu / C3;
            T RP = SMA * (1 - ECC);
            
            static math::Matrix<T> sVHat(3, 1, 0.0);
            T fac1;
            static math::Matrix<T> hcrosse(3, 1, 0.0);

            // sVHat is the asymptote vector of position 
            if (C3 > 1E-7)
            {
                fac1 = 1 / (1 + C3*hMag*hMag / mu / mu);
                hcrosse = hVec.cross(eccVec);
                sVHat = (-hcrosse * sqrt(C3) / mu - eccVec) * fac1;
            }
            else if (C3 < -1E-7)
            {
                std::cout << "Warning: Orbit is elliptic so using Apsides vector for asymptote." << std::endl;
                
                sVHat = -eccVec / ECC;
            }


            T khati[] = { 0.0, 0.0, 1.0 };
            static math::Matrix<T> khat(3, 1, khati);
            if (acos(sVHat.dot(khat)) < 1E-7) // 1e-7, 1e-11, criteria??
            {
                throw std::runtime_error("Error, IncomingAsymptote vector is aligned with z-direction. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            static math::Matrix<T> eaVec = khat.cross(sVHat);
            static math::Matrix<T>  eaVhat = eaVec / eaVec.norm();
            static math::Matrix<T>  noVhat = sVHat.cross(eaVhat);
            static math::Matrix<T>  bVec = hVec.cross(sVHat);
            T sinBVA = bVec.dot(eaVhat) / hMag;
            T cosBVA = bVec.dot(-noVhat) / hMag;
            T BVA = atan2(sinBVA, cosBVA);

            if (BVA < 0.0)
                BVA += math::TwoPI;

            T DHA = asin(sVHat(2));
            T RHA = atan2(sVHat(1), sVHat(0));
            if (RHA < 0.0)
                RHA += math::TwoPI;

            T TA = acos((eccVec.dot(R)) / (ECC*r));
            if (rdotv < 0.0)
                TA = math::TwoPI - TA;

            IncomingAsymptote_state[0] = RP;
            IncomingAsymptote_state[1] = C3;
            IncomingAsymptote_state[2] = RHA;
            IncomingAsymptote_state[3] = DHA;
            IncomingAsymptote_state[4] = BVA;
            IncomingAsymptote_state[5] = TA;


            if (GenerateDerivatives)
            {
                //derivatives of SMA and ECC, which are necessary to do everything else

                static math::Matrix<double> dr_dstate(6, 1, 0.0);
                static math::Matrix<double> dv_dstate(6, 1, 0.0);
                dr_dstate(0) = (Inertial_state[0] / r)_GETVALUE;
                dr_dstate(1) = (Inertial_state[1] / r)_GETVALUE;
                dr_dstate(2) = (Inertial_state[2] / r)_GETVALUE;
                dv_dstate(3) = (Inertial_state[3] / v)_GETVALUE;
                dv_dstate(4) = (Inertial_state[4] / v)_GETVALUE;
                dv_dstate(5) = (Inertial_state[5] / v)_GETVALUE;

                //derivatives of SMA
                static math::Matrix<double> dSMAdstate(6, 1, 0.0);
                double dSMAdr = (2 * mu * mu / (r*v*v - 2 * mu) / (r*v*v - 2 * mu))_GETVALUE;
                double dSMAdv = (2 * mu * r*r*v / (r*v*v - 2 * mu) / (r*v*v - 2 * mu))_GETVALUE;
                for (size_t k = 0; k < 3; ++k)
                {
                    dSMAdstate(k) = dSMAdr * dr_dstate(k);
                    dSMAdstate(k + 3) = dSMAdv * dv_dstate(k + 3);
                }

                //derivatives of ECC
                //first need derivatives of rdotv and s
                static math::Matrix<double> drdotv_dstate(6, 1);
                drdotv_dstate(0) = Inertial_state[3]_GETVALUE;
                drdotv_dstate(1) = Inertial_state[4]_GETVALUE;
                drdotv_dstate(2) = Inertial_state[5]_GETVALUE;
                drdotv_dstate(3) = Inertial_state[0]_GETVALUE;
                drdotv_dstate(4) = Inertial_state[1]_GETVALUE;
                drdotv_dstate(5) = Inertial_state[2]_GETVALUE;
                double dsdv = 2 * v _GETVALUE;
                double dsdr = mu / (r * r)_GETVALUE;
                static math::Matrix<double> ds_dstate(6, 1);
                ds_dstate(0) = dsdr * dr_dstate(0);
                ds_dstate(1) = dsdr * dr_dstate(1);
                ds_dstate(2) = dsdr * dr_dstate(2);
                ds_dstate(3) = dsdv * dv_dstate(3);
                ds_dstate(4) = dsdv * dv_dstate(4);
                ds_dstate(5) = dsdv * dv_dstate(5);

                //derivatives of ECC vector
                static math::Matrix<double> deccVec_dstate(6, 3);
                //with respect to position
                for (size_t i = 0; i < 3; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                    {
                        deccVec_dstate(i, j) = 1 / mu * (ds_dstate(i) * Inertial_state[j]_GETVALUE + (i == j ? s _GETVALUE : 0) - drdotv_dstate(i) * Inertial_state[j + 3]_GETVALUE);
                    }
                }
                //with respect to velocity
                for (size_t i = 0; i < 3; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                    {
                        deccVec_dstate(i + 3, j) = 1 / mu * (ds_dstate(i + 3) * Inertial_state[j]_GETVALUE - drdotv_dstate(i + 3) * (Inertial_state[j + 3] - (i == j ? rdotv : 0))_GETVALUE);
                    }
                }

                //derivatives of ECC scalar
                static math::Matrix<double> dECC_dstate(6, 1, 0.0);
                for (size_t i = 0; i < 6; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                        dECC_dstate(i) += (eccVec(j) / ECC)_GETVALUE * deccVec_dstate(i, j);
                }


                //derivatives of RP
                double dRP_dSMA = 1.0 - ECC _GETVALUE;
                double dRP_dECC = -SMA _GETVALUE;
                for (size_t i = 0; i < 6; ++i)
                    IncomingAsymptote_derivatives(i, 0) = dRP_dSMA * dSMAdstate(i) + dRP_dECC * dECC_dstate(i);

                //derivatives of C3
                double dC3dv = 2.0 * v _GETVALUE;
                double dC3dr = 2.0 * mu / (r * r)_GETVALUE;
                static math::Matrix<double> dC3_dstate(6, 1, 0.0);

                for (size_t i = 0; i < 6; ++i)
                {
                    dC3_dstate(i) = dC3dr * dr_dstate(i) + dC3dv * dv_dstate(i);
                    IncomingAsymptote_derivatives(i, 1) = dC3_dstate(i);
                }

                //next we need derivatives of sVHat and therefore of hVec
                static math::Matrix<double> dhVec_dstate(6, 3, 0.0);
                dhVec_dstate(0, 1) = -Inertial_state[5]_GETVALUE;
                dhVec_dstate(0, 2) = Inertial_state[4]_GETVALUE;
                dhVec_dstate(1, 0) = Inertial_state[5]_GETVALUE;
                dhVec_dstate(1, 2) = -Inertial_state[3]_GETVALUE;
                dhVec_dstate(2, 0) = -Inertial_state[4]_GETVALUE;
                dhVec_dstate(2, 1) = Inertial_state[3]_GETVALUE;
                dhVec_dstate(3, 1) = Inertial_state[2]_GETVALUE;
                dhVec_dstate(3, 2) = -Inertial_state[1]_GETVALUE;
                dhVec_dstate(4, 0) = -Inertial_state[2]_GETVALUE;
                dhVec_dstate(4, 2) = Inertial_state[0]_GETVALUE;
                dhVec_dstate(5, 0) = Inertial_state[1]_GETVALUE;
                dhVec_dstate(5, 1) = -Inertial_state[0]_GETVALUE;

#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                //quick check of sVHat derivatives, to be removed later
                math::Matrix<double> dhVec_dstate_algorthmic(6, 3);
                for (size_t i = 0; i < 6; ++i)
                    for (size_t j = 0; j < 3; ++j)
                        dhVec_dstate_algorthmic(i, j) = hVec(j).getDerivative(i + 18);

                (dhVec_dstate - dhVec_dstate_algorthmic).element_divide(dhVec_dstate_algorthmic).print_to_file("dhVec_dstate_relative_error.txt");
                dhVec_dstate.print_to_file("dhVec_dstate_analytical.txt");
                dhVec_dstate_algorthmic.print_to_file("dhVec_dstate_algorthmic.txt");
#endif

                //derivatives of hMag
                static math::Matrix<double> dhMag_dstate(6, 1);
                for (size_t i = 0; i < 6; ++i)
                {
                    dhMag_dstate(i) = 0.0;
                    for (size_t j = 0; j < 3; ++j)
                        dhMag_dstate(i) += (hVec(j) / hMag)_GETVALUE * dhVec_dstate(i, j);
                }

                //derivatives of sVhat
                static math::Matrix<double> dsVHat_dstate(6, 3, 0.0);
                if (C3 > 1E-7)
                {
                    //derivatives of fac1
                    static math::Matrix<T> dfac1_dstate(6, 1, 0.0);
                    double dfac1_dC3 = (-hMag * hMag / (mu * mu * ((C3 * hMag * hMag) / (mu * mu) + 1) * ((C3 * hMag * hMag) / (mu * mu) + 1)))_GETVALUE;
                    double dfac1_dhMag = (-(2 * C3*hMag) / (mu * mu * ((C3 * hMag *hMag) / (mu * mu) + 1) * ((C3 * hMag *hMag) / (mu * mu) + 1)))_GETVALUE;
                    for (size_t i = 0; i < 6; ++i)
                        dfac1_dstate(i) = dfac1_dC3 * dC3_dstate(i);
                    
                    //first need derivatives of hVec.cross(eccVec)
                    static math::Matrix<double> dhcrosse_dhVec(3, 3, 0.0);
                    static math::Matrix<double> dhcrosse_deccVec(3, 3, 0.0);
                    static math::Matrix<double> dhcrosse_dstate(6, 3, 0.0);

                    dhcrosse_dhVec(0, 0) = 0.0;
                    dhcrosse_dhVec(0, 1) = -eccVec(2)_GETVALUE;
                    dhcrosse_dhVec(0, 2) = eccVec(1)_GETVALUE;
                    dhcrosse_dhVec(1, 0) = eccVec(2)_GETVALUE;
                    dhcrosse_dhVec(1, 1) = 0.0;
                    dhcrosse_dhVec(1, 2) = -eccVec(0)_GETVALUE;
                    dhcrosse_dhVec(2, 0) = -eccVec(1)_GETVALUE;
                    dhcrosse_dhVec(2, 1) = eccVec(0)_GETVALUE;
                    dhcrosse_dhVec(2, 2) = 0.0;

                    dhcrosse_deccVec(0, 0) = 0.0;
                    dhcrosse_deccVec(0, 1) = hVec(2)_GETVALUE;
                    dhcrosse_deccVec(0, 2) = -hVec(1)_GETVALUE;
                    dhcrosse_deccVec(1, 0) = -hVec(2)_GETVALUE;
                    dhcrosse_deccVec(1, 1) = 0.0;
                    dhcrosse_deccVec(1, 2) = hVec(0)_GETVALUE;
                    dhcrosse_deccVec(2, 0) = hVec(1)_GETVALUE;
                    dhcrosse_deccVec(2, 1) = -hVec(0)_GETVALUE;
                    dhcrosse_deccVec(2, 2) = 0.0;


                    dhcrosse_dstate = dhcrosse_dhVec * dhVec_dstate + dhcrosse_deccVec * deccVec_dstate;
                    /*for (size_t i = 0; i < 6; ++i)
                    {
                        for (size_t j = 0; j < 3; ++j)
                        {
                            dhcrosse_dstate(i, j) = 0.0;
                            for (size_t k = 0; k < 3; ++k)
                            {
                                dhcrosse_dstate(i, k) += dhcrosse_dhVec(k, j) * dhVec_dstate(i, j) + dhcrosse_deccVec(k, j) * deccVec_dstate(i, j);
                            }
                        }
                    }*/
#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                    //quick check of sVHat derivatives, to be removed later
                    math::Matrix<double> dhcrosse_dstate_algorthmic(6, 3);
                    for (size_t i = 0; i < 6; ++i)
                        for (size_t j = 0; j < 3; ++j)
                            dhcrosse_dstate_algorthmic(i, j) = hcrosse(j).getDerivative(i + 18);

                    (dhcrosse_dstate - dhcrosse_dstate_algorthmic).element_divide(dhcrosse_dstate_algorthmic).print_to_file("dhcrosse_dstate_relative_error.txt");
                    dhcrosse_dstate.print_to_file("dhcrosse_dstate_analytical.txt");
                    dhcrosse_dstate_algorthmic.print_to_file("dhcrosse_dstate_algorthmic.txt");
#endif

                    //sVHat = (-hVec.cross(eccVec) * sqrt(C3) / mu - eccVec) * fac1;
                    //this one will stink
                    static math::Matrix<double> dsVHat_dstate(6, 3, 0.0);
                    for (size_t i = 0; i < 6; ++i)
                    {
                        for (size_t j = 0; j < 3; ++j)
                        {
                            dsVHat_dstate(i, j) = ( fac1 * (-sqrt(C3) / mu * dhcrosse_dstate(i, j)
                                                            - 1.0 / (2.0 * mu * sqrt(C3)) * hcrosse(j) * dC3_dstate(i)
                                                            - deccVec_dstate(i, j))
                                                    + dfac1_dstate(i) * (-hcrosse(j) * sqrt(C3) / mu - eccVec(j)))_GETVALUE;
                        }
                    }

#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                    //quick check of sVHat derivatives, to be removed later
                    math::Matrix<double> sVHat_derivatives_algorthmic(6, 3);
                    for (size_t i = 0; i < 6; ++i)
                        for (size_t j = 0; j < 3; ++j)
                            sVHat_derivatives_algorthmic(i, j) = sVHat(j).getDerivative(i + 18);

                    (dsVHat_dstate - sVHat_derivatives_algorthmic).element_divide(sVHat_derivatives_algorthmic).print_to_file("dsVHat_dstate_relative_error.txt");
                    dsVHat_dstate.print_to_file("dsVHat_dstate_analytical.txt");
                    sVHat_derivatives_algorthmic.print_to_file("dsVHat_dstate_algorithmic.txt");
#endif
                }
                else if (C3 < -1E-7)
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        for (size_t j = 0; j < 6; ++j)
                        {
                            dsVHat_dstate(i, j) = (-(eccVec(j) * dECC_dstate(i) - ECC * deccVec_dstate(i, j)) / (ECC * ECC))_GETVALUE;
                        }
                    }
                }

                //derivatives of RHA

                //derivatives of DHA

                //derivatives of BVA

                //derivatives of TA
                //T TA = acos((eccVec.dot(R)) / (ECC*r));
                //if (rdotv < 0.0)
                //    TA = math::TwoPI - TA;
                double TAterm = (eccVec.dot(R) / (ECC*r)) _GETVALUE;
                double diffTAterm = -1.0 / sqrt(1.0 - TAterm * TAterm);
                int TAsign = (rdotv _GETVALUE < 0.0 ? -1 : 1);

            }//end derivatives for Inertial2IncomingAsymptote
        }//end Inertial2IncomingAsymptote

        template<class T> void IncomingAsymptote2Inertial(const std::vector<T>& IncomingAsymptote_state,//6x1
            const double& mu,
            std::vector<T>& Inertial_state,//6x1
            const bool& GenerateDerivatives,
            math::Matrix<T>& Inertial_state_derivatives)//6x6
        {
            //first convert from IncomingAsymptote to COE
            static std::vector<T> COE_state(6, 0.0);
            static math::Matrix<T> COE_derivatives(6, 6, 0.0);
            IncomingAsymptote2COE(IncomingAsymptote_state, mu, COE_state, GenerateDerivatives, COE_derivatives);

            //then convert from COE to Inertial
            COE2inertial(COE_state, mu, Inertial_state, GenerateDerivatives, Inertial_state_derivatives);
            if (GenerateDerivatives)
            {
                Inertial_state_derivatives = COE_derivatives * Inertial_state_derivatives;
            }//end IncomingAsymptote2Inertial derivatives
        }//end IncomingAsymptote2Inertial

        template<class T> void COE2IncomingAsymptote(const std::vector<T>& COE_state,//6x1
            const double& mu,
            std::vector<T>& IncomingAsymptote_state,//6x1
            const bool& GenerateDerivatives,
            math::Matrix<T>& IncomingAsymptote_derivatives)//6x6
        {
            //first convert from COE to inertial
            static std::vector<T> Inertial_state(6, 0.0);
            static math::Matrix<T> Inertial_derivatives(6, 6, 0.0);
            COE2inertial(COE_state, mu, Inertial_state, GenerateDerivatives, Inertial_derivatives);

            //then convert from COE to Inertial
            Inertial2IncomingAsymptote(Inertial_state, mu, IncomingAsymptote_state, GenerateDerivatives, IncomingAsymptote_derivatives);
            if (GenerateDerivatives)
            {
                IncomingAsymptote_derivatives = Inertial_derivatives * IncomingAsymptote_derivatives;
            }//end COE2IncomingAsymptote derivatives
        }//end COE2IncomingAsymptote

    }//end namespace Astrodynamics
}//end namespace EMTG