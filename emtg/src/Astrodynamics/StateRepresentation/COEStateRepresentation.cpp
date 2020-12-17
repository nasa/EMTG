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

#include "COEStateRepresentation.h"
#include "EMTG_math.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        COEStateRepresentation::COEStateRepresentation(const double& mu) : StateRepresentationBase::StateRepresentationBase(mu)
        {
            this->name = "COE";
            this->stateNames = std::vector<std::string>({ "SMA", "ECC", "INC", "RAAN", "AOP", "TA" });
        };

        math::Matrix<doubleType> COEStateRepresentation::convertFromRepresentationToCartesian(const math::Matrix<doubleType>& StateVectorThisRepresentationIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorThisRepresentation(StateVectorThisRepresentationIn);

            //Step 2: convert the state           
            const doubleType& SMA = this->StateVectorThisRepresentation(0);
            const doubleType& ECC = this->StateVectorThisRepresentation(1);
            const doubleType& INC = this->StateVectorThisRepresentation(2);
            const doubleType& RAAN = this->StateVectorThisRepresentation(3);
            const doubleType& AOP = this->StateVectorThisRepresentation(4);
            const doubleType& TA = this->StateVectorThisRepresentation(5);
            doubleType THETA = AOP + TA;

            doubleType cosTA = cos(TA);
            doubleType cosINC = cos(INC);
            doubleType sinINC = sin(INC);
            doubleType cosRAAN = cos(RAAN);
            doubleType sinRAAN = sin(RAAN);
            doubleType cosAOP = cos(AOP);
            doubleType sinAOP = sin(AOP);
            doubleType cosTHETA = cos(THETA);
            doubleType sinTHETA = sin(THETA);

            doubleType p = SMA * (1 - ECC * ECC);

            if (p < 1.0e-30)
            {
                throw std::runtime_error("Error converting parabolic orbit to Cartesian coordinates. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            doubleType r = p / (1 + ECC * cosTA);

            doubleType h = sqrt(mu * SMA * (1 - ECC * ECC));

            this->StateVectorCartesian(0) = r * (cosRAAN*cosTHETA - sinRAAN * sinTHETA*cosINC);
            this->StateVectorCartesian(1) = r * (sinRAAN*cosTHETA + cosRAAN * sinTHETA*cosINC);
            this->StateVectorCartesian(2) = r * sinTHETA * sinINC;
            this->StateVectorCartesian(3) = -mu / h * (cosRAAN * (sinTHETA + ECC * sinAOP) + sinRAAN * (cosTHETA + ECC * cosAOP) * cosINC);
            this->StateVectorCartesian(4) = -mu / h * (sinRAAN * (sinTHETA + ECC * sinAOP) - cosRAAN * (cosTHETA + ECC * cosAOP) * cosINC);
            this->StateVectorCartesian(5) = mu / h * (cosTHETA + ECC * cosAOP) * sinINC;


            //Step 3: construct the partial derivatives
            if (needG)
            {
                this->RepresentationToCartesianTransitionMatrix.assign_zeros();

                //first we need derivatives of r and h with respect to their components
                doubleType dr_dSMA = (-(ECC * ECC - 1) / (ECC*cosTA + 1))_GETVALUE;
                doubleType dr_dECC = (-(SMA*(cosTA * ECC * ECC + 2 * ECC + cosTA)) / (ECC*cosTA + 1) / (ECC*cosTA + 1))_GETVALUE;
                doubleType dr_dTA = (-(ECC*SMA*sin(TA)*(ECC*ECC - 1)) / (ECC*cosTA + 1) / (ECC*cosTA + 1))_GETVALUE;
                doubleType dh_dSMA = (-(mu*(ECC*ECC - 1)) / (2 * sqrt(-SMA * mu*(ECC*ECC - 1))))_GETVALUE;
                doubleType dh_dECC = (-(ECC*SMA*mu) / sqrt(-SMA * mu*(ECC *ECC - 1)))_GETVALUE;

                //transition matrix has NewThing for rows and OldThing for columns
                //derivatives of x
                this->RepresentationToCartesianTransitionMatrix(0, 0) = dr_dSMA * (cosRAAN*cosTHETA - sinRAAN * sinTHETA*cosINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(0, 1) = dr_dECC * (cosRAAN*cosTHETA - sinRAAN * sinTHETA*cosINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(0, 2) = (r * sinRAAN*sinTHETA*sinINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(0, 3) = (r * (-sinRAAN * cosTHETA - cosRAAN * sinTHETA*cosINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(0, 4) = (r * (-cosRAAN * sinTHETA - sinRAAN * cosTHETA*cosINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(0, 5) = dr_dTA * (cosRAAN*cosTHETA - sinRAAN * sinTHETA*cosINC)_GETVALUE + (r * (-cosRAAN * sinTHETA - sinRAAN * cosTHETA*cosINC))_GETVALUE;

                //derivatives of y
                this->RepresentationToCartesianTransitionMatrix(1, 0) = dr_dSMA * (sinRAAN*cosTHETA + cosRAAN * sinTHETA*cosINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(1, 1) = dr_dECC * (sinRAAN*cosTHETA + cosRAAN * sinTHETA*cosINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(1, 2) = (r * -cosRAAN * sinTHETA*sinINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(1, 3) = (r * (cosRAAN*cosTHETA - sinRAAN * sinTHETA*cosINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(1, 4) = (r * (-sinRAAN * sinTHETA + cosRAAN * cosTHETA*cosINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(1, 5) = dr_dTA * (sinRAAN*cosTHETA + cosRAAN * sinTHETA*cosINC)_GETVALUE + (r * (-sinRAAN * sinTHETA + cosRAAN * cosTHETA*cosINC))_GETVALUE;

                //derivatives of z
                this->RepresentationToCartesianTransitionMatrix(2, 0) = dr_dSMA * (sinTHETA * sinINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(2, 1) = dr_dECC * (sinTHETA * sinINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(2, 2) = (r * sinTHETA * cosINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(2, 3) = 0.0;
                this->RepresentationToCartesianTransitionMatrix(2, 4) = (r * cosTHETA * sinINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(2, 5) = dr_dTA * (sinTHETA * sinINC)_GETVALUE + (r * (cosTHETA * sinINC))_GETVALUE;

                //derivatives of xdot
                this->RepresentationToCartesianTransitionMatrix(3, 0) = (-dh_dSMA / (h*h) * -mu * (cosRAAN * (sinTHETA + ECC * sinAOP) + sinRAAN * (cosTHETA + ECC * cosAOP) * cosINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(3, 1) = (-dh_dECC / (h*h) * -mu * (cosRAAN * (sinTHETA + ECC * sinAOP) + sinRAAN * (cosTHETA + ECC * cosAOP) * cosINC)
                    - mu / h * (cosRAAN * sinAOP + sinRAAN * cosAOP * cosINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(3, 2) = (-mu / h * (sinRAAN * (cosTHETA + ECC * cosAOP) * -sinINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(3, 3) = (-mu / h * (-sinRAAN * (sinTHETA + ECC * sinAOP) + cosRAAN * (cosTHETA + ECC * cosAOP) * cosINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(3, 4) = (-mu / h * (cosRAAN * (cosTHETA + ECC * cosAOP) + sinRAAN * (-sinTHETA + ECC * -sinAOP) * cosINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(3, 5) = (-mu / h * (cosRAAN * (cosTHETA)+sinRAAN * (-sinTHETA) * cosINC))_GETVALUE;

                //derivatives of ydot
                this->RepresentationToCartesianTransitionMatrix(4, 0) = (-dh_dSMA / (h*h) *-mu * (sinRAAN * (sinTHETA + ECC * sinAOP) - cosRAAN * (cosTHETA + ECC * cosAOP) * cosINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(4, 1) = (-dh_dECC / (h*h) *-mu * (sinRAAN * (sinTHETA + ECC * sinAOP) - cosRAAN * (cosTHETA + ECC * cosAOP) * cosINC)
                    - mu / h * (sinRAAN * sinAOP - cosRAAN * cosAOP * cosINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(4, 2) = (-mu / h * (-cosRAAN * (cosTHETA + ECC * cosAOP) * -sinINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(4, 3) = (-mu / h * (cosRAAN * (sinTHETA + ECC * sinAOP) + sinRAAN * (cosTHETA + ECC * cosAOP) * cosINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(4, 4) = (-mu / h * (sinRAAN * (cosTHETA + ECC * cosAOP) - cosRAAN * (-sinTHETA + ECC * -sinAOP) * cosINC))_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(4, 5) = (-mu / h * (sinRAAN * (cosTHETA)-cosRAAN * (-sinTHETA) * cosINC))_GETVALUE;

                //derivatives of zdot
                this->RepresentationToCartesianTransitionMatrix(5, 0) = (-dh_dSMA / (h*h) * mu  * (cosTHETA + ECC * cosAOP) * sinINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(5, 1) = (-dh_dECC / (h*h) * mu  * (cosTHETA + ECC * cosAOP) * sinINC + mu / h * cosAOP * sinINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(5, 2) = (mu / h * (cosTHETA + ECC * cosAOP) * cosINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(5, 3) = 0.0;
                this->RepresentationToCartesianTransitionMatrix(5, 4) = (mu / h * (-sinTHETA - ECC * sinAOP) * sinINC)_GETVALUE;
                this->RepresentationToCartesianTransitionMatrix(5, 5) = (mu / h * (-sinTHETA) * sinINC)_GETVALUE;
            }

            //Step 4: return the converted vector
            return this->StateVectorCartesian;
        }//end convertFromRepresentationToCartesian()

        math::Matrix<doubleType> COEStateRepresentation::convertFromCartesianToRepresentation(const math::Matrix<doubleType>& StateVectorCartesianIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorCartesian(StateVectorCartesianIn);

            //Step 2: convert the state           
            //this is trivial for Cartesian
            const doubleType& x = this->StateVectorCartesian(0);
            const doubleType& y = this->StateVectorCartesian(1);
            const doubleType& z = this->StateVectorCartesian(2);
            const doubleType& vx = this->StateVectorCartesian(3);
            const doubleType& vy = this->StateVectorCartesian(4);
            const doubleType& vz = this->StateVectorCartesian(5);

            doubleType r, v, h, n, rdotv, s;
            math::Matrix<doubleType> evec(3, 1), nvec(3, 1), hvec(3, 1);
            math::Matrix<doubleType> ihat(3, 1, { 1, 0, 0 });
            math::Matrix<doubleType> jhat(3, 1, { 0, 1, 0 });
            math::Matrix<doubleType> khat(3, 1, { 0, 0, 1 });

            math::Matrix<doubleType> R = this->StateVectorCartesian.getSubMatrix1D(0, 2);
            math::Matrix<doubleType> V = this->StateVectorCartesian.getSubMatrix1D(3, 5);
            r = R.norm();
            v = V.norm();

            //SMA
            this->StateVectorThisRepresentation(0) = r / (2 - r * v*v / mu);

            //eccentricity vector
            rdotv = R.dot(V);
            s = (v*v - mu / r);

            evec(0) = 1 / mu * (s*this->StateVectorCartesian(0) - rdotv * this->StateVectorCartesian(3));
            evec(1) = 1 / mu * (s*this->StateVectorCartesian(1) - rdotv * this->StateVectorCartesian(4));
            evec(2) = 1 / mu * (s*this->StateVectorCartesian(2) - rdotv * this->StateVectorCartesian(5));

            //ECC
            doubleType ECC = evec.norm();
            this->StateVectorThisRepresentation(1) = ECC;

            //angular momentum vector and scalar
            hvec = R.cross(V);
            h = hvec.norm();

            //INC
            this->StateVectorThisRepresentation(2) = math::safe_acos(hvec(2) / h);

            //nodal vector
            nvec = khat.cross(hvec) / h;

            n = nvec.norm();

            //RAAN

            if (nvec(1) >= 0)
                this->StateVectorThisRepresentation(3) = math::safe_acos(nvec(0) / n);
            else
                this->StateVectorThisRepresentation(3) = 2 * EMTG::math::PI - math::safe_acos(nvec(0) / n);

            if (n == 0)
                this->StateVectorThisRepresentation(3) = 0;

            //AOP
            doubleType ndote = nvec.dot(evec);
            if (evec(2) >= 0)
                this->StateVectorThisRepresentation(4) = math::safe_acos(ndote / (n * ECC));
            else
                this->StateVectorThisRepresentation(4) = 2 * EMTG::math::PI - math::safe_acos(ndote / (n * ECC));

            if (n == 0) //if no inclination, then eccentricity vector points to the periapse
                this->StateVectorThisRepresentation(4) = atan2(evec(1), evec(0));

            //TA
            doubleType edotR = evec.dot(R);
            if (rdotv >= 0)
                this->StateVectorThisRepresentation(5) = math::safe_acos(edotR / (r * ECC));
            else
                this->StateVectorThisRepresentation(5) = 2 * EMTG::math::PI - math::safe_acos(edotR / (r * ECC));

            //Step 3: construct the partial derivatives
            if (needG)
            {
                this->CartesianToRepresentationTransitionMatrix.assign_zeros();

                doubleType R = r _GETVALUE;
                doubleType V = v _GETVALUE;
                doubleType x = this->StateVectorCartesian(0)_GETVALUE;
                doubleType y = this->StateVectorCartesian(1)_GETVALUE;
                doubleType z = this->StateVectorCartesian(2)_GETVALUE;
                doubleType xdot = this->StateVectorCartesian(3)_GETVALUE;
                doubleType ydot = this->StateVectorCartesian(4)_GETVALUE;
                doubleType zdot = this->StateVectorCartesian(5)_GETVALUE;

                math::Matrix<doubleType> dr_drstate(3, 1);
                math::Matrix<doubleType> dv_dvstate(3, 1);
                dr_drstate(0) = x / R;
                dr_drstate(1) = y / R;
                dr_drstate(2) = z / R;
                dv_dvstate(0) = xdot / V;
                dv_dvstate(1) = ydot / V;
                dv_dvstate(2) = zdot / V;

                //transition matrix has NewThing for rows and OldThing for columns
                //derivatives of SMA
                doubleType dSMAdr = 2 * mu * mu / (R*V*V - 2 * mu) / (R*V*V - 2 * mu);
                this->CartesianToRepresentationTransitionMatrix(0, 0) = dSMAdr * dr_drstate(0);
                this->CartesianToRepresentationTransitionMatrix(0, 1) = dSMAdr * dr_drstate(1);
                this->CartesianToRepresentationTransitionMatrix(0, 2) = dSMAdr * dr_drstate(2);
                doubleType dSMAdv = 2 * mu * R*R*V / (R*V*V - 2 * mu) / (R*V*V - 2 * mu);
                this->CartesianToRepresentationTransitionMatrix(0, 3) = dSMAdv * dv_dvstate(0);
                this->CartesianToRepresentationTransitionMatrix(0, 4) = dSMAdv * dv_dvstate(1);
                this->CartesianToRepresentationTransitionMatrix(0, 5) = dSMAdv * dv_dvstate(2);

                //derivatives of ECC
                //first need derivatives of rdotv and s
                math::Matrix<doubleType> drdotv_dstate(6, 1);
                drdotv_dstate(0) = xdot;
                drdotv_dstate(1) = ydot;
                drdotv_dstate(2) = zdot;
                drdotv_dstate(3) = x;
                drdotv_dstate(4) = y;
                drdotv_dstate(5) = z;
                doubleType dsdv = 2 * V;
                doubleType dsdr = mu / (R * R);
                math::Matrix<doubleType> ds_dstate(6, 1);
                ds_dstate(0) = dsdr * dr_drstate(0);
                ds_dstate(1) = dsdr * dr_drstate(1);
                ds_dstate(2) = dsdr * dr_drstate(2);
                ds_dstate(3) = dsdv * dv_dvstate(0);
                ds_dstate(4) = dsdv * dv_dvstate(1);
                ds_dstate(5) = dsdv * dv_dvstate(2);

                //derivatives of ECC vector
                math::Matrix<doubleType> devec_dstate(6, 3);
                //with respect to position
                for (size_t i = 0; i < 3; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                    {
                        devec_dstate(i, j) = 1 / mu * (ds_dstate(i) * this->StateVectorCartesian(j) + (i == j ? s : 0) - drdotv_dstate(i) * this->StateVectorCartesian(j + 3))_GETVALUE;
                    }
                }
                //with respect to velocity
                for (size_t i = 0; i < 3; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                    {
                        devec_dstate(i + 3, j) = 1 / mu * (ds_dstate(i + 3) * this->StateVectorCartesian(j) - drdotv_dstate(i + 3) * this->StateVectorCartesian(j + 3) - (i == j ? rdotv : 0))_GETVALUE;
                    }
                }

#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                //check of eccentricity vector derivatives, do be removed
                math::Matrix<doubleType> devec_dstate_algorithmic(6, 3);
                for (size_t i = 0; i < 6; ++i)
                    for (size_t j = 0; j < 3; ++j)
                        devec_dstate_algorithmic(i, j) = evec[j] _GETDERIVATIVE;

                (devec_dstate - devec_dstate_algorithmic).element_divide(devec_dstate_algorithmic).print_to_file("devec_dstate_relative_error.txt");
#endif

                //derivatives of ECC scalar
                for (size_t i = 0; i < 6; ++i)
                {
                    this->CartesianToRepresentationTransitionMatrix(1, i) = 0.0;
                    for (size_t j = 0; j < 3; ++j)
                        this->CartesianToRepresentationTransitionMatrix(1, i) += (evec(j) / this->StateVectorThisRepresentation(1))_GETVALUE * devec_dstate(i, j);
                }

                //derivatives of INC
                //first need the derivatives of the h-vector
                math::Matrix<doubleType> dhvec_dstate(6, 3);
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
                math::Matrix<doubleType> h_vector_derivatives_algorthmic(6, 3);
                for (size_t i = 0; i < 6; ++i)
                    for (size_t j = 0; j < 3; ++j)
                        h_vector_derivatives_algorthmic(i, j) = hvec[j] _GETDERIVATIVE;

                (dhvec_dstate - h_vector_derivatives_algorthmic).element_divide(h_vector_derivatives_algorthmic).print_to_file("dhvec_dstate_relative_error.txt");
#endif

                //derivatives of h scalar
                math::Matrix<doubleType> dhscalar_dstate(6, 1);
                for (size_t i = 0; i < 6; ++i)
                {
                    dhscalar_dstate(i) = 0.0;
                    for (size_t j = 0; j < 3; ++j)
                        dhscalar_dstate(i) += (hvec(j) / h)_GETVALUE * dhvec_dstate(i, j);
                }

#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                //quick check of h-scalar derivatives, to be removed later
                math::Matrix<doubleType> dhscalar_dstate_algorithmic(6, 1);
                for (size_t i = 0; i < 6; ++i)
                    dhscalar_dstate_algorithmic(i) = h _GETDERIVATIVE;
                (dhscalar_dstate - dhscalar_dstate_algorithmic).element_divide(dhscalar_dstate_algorithmic).print_to_file("dhscalar_dstate_relative_error.txt");
#endif

                //derivatives of INC
                doubleType INCterm = (hvec(2) / h)_GETVALUE;
                doubleType diffINCterm = -1.0 / sqrt(1 - INCterm * INCterm) / (h*h)_GETVALUE;
                for (size_t i = 0; i < 6; ++i)
                {
                    this->CartesianToRepresentationTransitionMatrix(2, i) = diffINCterm * (dhvec_dstate(i, 2) * h - dhscalar_dstate(i) * hvec(2))_GETVALUE;
                }

                //derivatives of RAAN
                //we need derivatives of the nodal vector
                math::Matrix<doubleType> dnvec_dstate(6, 3, 0.0);
                math::Matrix<doubleType> dkcrossh_dstate(6, 3, 0.0);
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
                math::Matrix<doubleType> n_vector_derivatives_algorthmic(6, 3);
                for (size_t i = 0; i < 6; ++i)
                    for (size_t j = 0; j < 3; ++j)
                        n_vector_derivatives_algorthmic(i, j) = nvec[j] _GETDERIVATIVE;

                (dnvec_dstate - n_vector_derivatives_algorthmic).element_divide(n_vector_derivatives_algorthmic).print_to_file("dnvec_dstate_relative_error.txt");
#endif

                //derivatives of n scalar
                math::Matrix<doubleType> dnscalar_dstate(6, 1);
                for (size_t i = 0; i < 6; ++i)
                {
                    dnscalar_dstate(i) = 0.0;
                    for (size_t j = 0; j < 3; ++j)
                        dnscalar_dstate(i) += (nvec(j) / n)_GETVALUE * dnvec_dstate(i, j);
                }

#ifdef ORBIT_ELEMENT_COMPONENT_DERIVATIVE_CHECK
                //quick check of n-scalar derivatives, to be removed later
                math::Matrix<doubleType> dnscalar_dstate_algorithmic(6, 1);
                for (size_t i = 0; i < 6; ++i)
                    dnscalar_dstate_algorithmic(i) = n _GETDERIVATIVE;
                (dnscalar_dstate - dnscalar_dstate_algorithmic).element_divide(dnscalar_dstate_algorithmic).print_to_file("dnscalar_dstate_relative_error.txt");
#endif

                //now derivatives of RAAN
                doubleType RAANterm = (nvec(0) / n)_GETVALUE;
                doubleType diffRAANterm = -1.0 / sqrt(1 - RAANterm * RAANterm) / (n*n)_GETVALUE;
                if (nvec(1) >= 0)
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        this->CartesianToRepresentationTransitionMatrix(3, i) = diffRAANterm * (dnvec_dstate(i, 0) * n - dnscalar_dstate(i) * nvec(0))_GETVALUE;
                    }
                }
                else if (n == 0) //pathological crash case
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        this->CartesianToRepresentationTransitionMatrix(3, i) = 0.0;
                    }
                }
                else
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        this->CartesianToRepresentationTransitionMatrix(3, i) = -diffRAANterm * (dnvec_dstate(i, 0) * n - dnscalar_dstate(i) * nvec(0))_GETVALUE;
                    }
                }

                //derivatives of AOP
                doubleType AOPterm = (ndote / n / ECC)_GETVALUE;
                doubleType diffAOPterm = -1.0 / sqrt(1 - AOPterm * AOPterm) / (n*n*ECC*ECC)_GETVALUE;
                if (evec(2) >= 0)
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        doubleType dndote_dstate = 0.0;
                        for (size_t j = 0; j < 3; ++j)
                        {
                            dndote_dstate += (dnvec_dstate(i, j) * evec(j) + devec_dstate(i, j) * nvec(j))_GETVALUE;
                        }
                        this->CartesianToRepresentationTransitionMatrix(4, i) = diffAOPterm * (dndote_dstate * n * ECC - (dnscalar_dstate(i)*ECC + this->CartesianToRepresentationTransitionMatrix(1, i) * n) * ndote)_GETVALUE;
                    }
                }
                else if (n == 0)//if no inclination, then eccentricity vector points to the periapse
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        this->CartesianToRepresentationTransitionMatrix(4, i) = ((-evec(1) * devec_dstate(i, 0) + evec(0) * devec_dstate(i, 1)) / (evec(0) * evec(0) + evec(1) * evec(1)))_GETVALUE;
                    }
                }
                else
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        doubleType dndote_dstate = 0.0;
                        for (size_t j = 0; j < 3; ++j)
                        {
                            dndote_dstate += (dnvec_dstate(i, j) * evec(j) + devec_dstate(i, j) * nvec(j))_GETVALUE;
                        }
                        this->CartesianToRepresentationTransitionMatrix(4, i) = -diffAOPterm * (dndote_dstate * n * ECC - (dnscalar_dstate(i)*ECC + this->CartesianToRepresentationTransitionMatrix(1, i) * n) * ndote)_GETVALUE;
                    }
                }

                //derivatives of true anomaly
                doubleType TAterm = math::absclip((edotR / r / ECC)_GETVALUE, 1.0 - math::SMALL);
                doubleType diffTAterm = -1.0 / sqrt(1 - TAterm * TAterm) / (r*r*ECC*ECC)_GETVALUE;
                if (rdotv >= 0)
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        doubleType dedotstate_dstate = 0.0;
                        for (size_t j = 0; j < 3; ++j)
                        {
                            dedotstate_dstate += (devec_dstate(i, j) * this->StateVectorCartesian(j) + (i == j ? 1 : 0) * evec(j))_GETVALUE;
                        }
                        this->CartesianToRepresentationTransitionMatrix(5, i) = diffTAterm * (dedotstate_dstate * r * ECC - (this->CartesianToRepresentationTransitionMatrix(1, i)*r + (i < 3 ? dr_drstate(i) : 0.0) * ECC) * edotR)_GETVALUE;
                    }
                }
                else
                {
                    for (size_t i = 0; i < 6; ++i)
                    {
                        doubleType dedotstate_dstate = 0.0;
                        for (size_t j = 0; j < 3; ++j)
                        {
                            dedotstate_dstate += (devec_dstate(i, j) * this->StateVectorCartesian(j) + (i == j ? 1 : 0) * evec(j))_GETVALUE;
                        }
                        this->CartesianToRepresentationTransitionMatrix(5, i) = -diffTAterm * (dedotstate_dstate * r * ECC - (this->CartesianToRepresentationTransitionMatrix(1, i)*r + (i < 3 ? dr_drstate(i) : 0.0) * ECC) * edotR)_GETVALUE;
                    }
                }
            }

            //Step 4: return the converted vector
            return this->StateVectorThisRepresentation;
        }//end convertFromCartesianToRepresentation()
    }//end namespace StateRepresentation
}//end namespace EMTG