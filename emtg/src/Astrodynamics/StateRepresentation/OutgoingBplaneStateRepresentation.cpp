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

#include "OutgoingBplaneStateRepresentation.h"
#include "EMTG_math.h"

#include <stdexcept>

//OutgoingBplane state representation
//we expect that people will ask what our b-plane vector is. It is current coordinate system zhat.
//We have to use a constant vector because otherwise the transformation from Bplane to cartesian is self-referential.
//So we picked current coordinate system zhat because it requires the least work. So there.

namespace EMTG
{
    namespace Astrodynamics
    {
        OutgoingBplaneStateRepresentation::OutgoingBplaneStateRepresentation(const double& mu) : StateRepresentationBase::StateRepresentationBase(mu)
        {
            this->name = "OutgoingBplane";
            this->stateNames = std::vector<std::string>({ "VINF", "RHA", "DHA", "BRADIUS", "BTHETA", "TA" });
        };

        math::Matrix<doubleType> OutgoingBplaneStateRepresentation::convertFromRepresentationToCartesian(const math::Matrix<doubleType>& StateVectorThisRepresentationIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorThisRepresentation(StateVectorThisRepresentationIn);

            //Step 2: convert the state           
            const doubleType& VINF = this->StateVectorThisRepresentation(0);
            const doubleType& RHA = this->StateVectorThisRepresentation(1);
            const doubleType& DHA = this->StateVectorThisRepresentation(2);
            const doubleType& BRADIUS = this->StateVectorThisRepresentation(3);
            const doubleType& BTHETA = this->StateVectorThisRepresentation(4);
            const doubleType& TA = this->StateVectorThisRepresentation(5);

            doubleType cosTA = cos(TA);
            doubleType sinTA = sin(TA);
            doubleType cosRHA = cos(RHA);
            doubleType sinRHA = sin(RHA);
            doubleType cosDHA = cos(DHA);
            doubleType sinDHA = sin(DHA);
            doubleType cosBTHETA = cos(BTHETA);
            doubleType sinBTHETA = sin(BTHETA);

            //working from Noble's math spec in the docs/bplane_parameterization folder
            //components first

            //eq 70
            doubleType ECC = sqrt(1.0 + VINF*VINF*VINF*VINF * BRADIUS*BRADIUS / mu / mu);

            if (ECC < 1.0 + math::SMALL)
            {
                throw std::runtime_error("Error converting from IncomingBplane to Cartesian. IncomingBplaneStateRepresentation can only work if the state has an ECC of > 1.0. ECC is currently " + std::to_string(ECC _GETVALUE) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            //eq 71
            doubleType h = VINF * BRADIUS;

            //eq 72
            math::Matrix<doubleType> Vinfinity(3, 1, {VINF * cosRHA * cosDHA,
                                                      VINF * sinRHA * cosDHA,
                                                      VINF * sinDHA});
            
            //eq 73
            this->BdotR = BRADIUS * sinBTHETA;

            //eq 74
            this->BdotT = BRADIUS * cosBTHETA;

            //eq 75
            math::Matrix<doubleType> Shat = Vinfinity.unitize();

            //for the b-plane reference vector, we will use z-hat
            math::Matrix<doubleType> phi(3, 1, std::vector<doubleType>({ 0, 0, 1.0 }));

            //eq 76
            math::Matrix<doubleType> That = Shat.unitcross(phi);

            //eq 77
            math::Matrix<doubleType> Rhat = Shat.unitcross(That);

            //eq 78
            math::Matrix<doubleType> Bvector = Rhat * this->BdotR + That * this->BdotT;

            //eq 79
            math::Matrix<doubleType> Hhat = Bvector.unitcross(Shat);

            //eq 80
            math::Matrix<doubleType> Hvec = Hhat * h;

            //eq 155
            doubleType TA_inf_out = math::safe_acos(-1.0 / ECC);

            //eq 156 - note we don't have Bhat yet
            math::Matrix<doubleType> Bhat = Bvector.unitize();
            math::Matrix<doubleType> Ehat = (Shat * cos(TA_inf_out) + Bhat * sin(TA_inf_out)).unitize();

            //eq 83
            math::Matrix<doubleType> Evec = Ehat * ECC;


            //position and velocity
            //eq 68 - need to define Phat
            math::Matrix<doubleType> Phat = Hvec.unitcross(Evec);
            math::Matrix<doubleType> pos = (Ehat * cosTA + Phat * sinTA) * h * h / mu / (1.0 + ECC * cosTA);

            //eq 69
            math::Matrix<doubleType> vel = (Ehat * sinTA - Phat * (ECC + cosTA)) * -mu / h;

            //now put it all in one place
            this->StateVectorCartesian(0) = pos(0);
            this->StateVectorCartesian(1) = pos(1);
            this->StateVectorCartesian(2) = pos(2);
            this->StateVectorCartesian(3) = vel(0);
            this->StateVectorCartesian(4) = vel(1);
            this->StateVectorCartesian(5) = vel(2);          


            //Step 3: construct the partial derivatives
            if (needG)
            {
                this->RepresentationToCartesianTransitionMatrix.assign_zeros();

                //once again, following Noble's math spec
                //x := b-plane state representation of state
                math::Matrix<doubleType> x = this->StateVectorThisRepresentation;

                //need a 3x3 identity
                math::Matrix<doubleType> I3(3, math::identity);

                //eq 84
                math::Matrix<doubleType> dBRADIUS_dx(1, 6, std::vector<doubleType>({ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 }));

                //eq 85
                math::Matrix<doubleType> dBTHETA_dx(1, 6, std::vector<doubleType>({ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 }));

                //eq 86
                math::Matrix<doubleType> dTA_dx(1, 6, std::vector<doubleType>({ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 }));

                //eq 88
                math::Matrix<doubleType> dECC_dx(1, 6, std::vector<doubleType>({ 4 * VINF*VINF*VINF*BRADIUS*BRADIUS / mu / mu, 0.0, 0.0, 2 * VINF*VINF*VINF*VINF*BRADIUS / mu / mu, 0.0, 0.0 }));
                dECC_dx *= 1.0 / (2.0 * ECC);

                //a disgression to get dEhat_dx
                //eq 157
                math::Matrix<doubleType> dTA_inf_out_dx = -dECC_dx * 1.0 / ECC / sqrt(ECC * ECC - 1.0);

                doubleType cosTA_inf_out = cos(TA_inf_out);

                doubleType sinTA_inf_out = sin(TA_inf_out);

                //eq 159
                math::Matrix<doubleType> dcosTA_inf_out_dx = dTA_inf_out_dx  * -sinTA_inf_out;

                //eq 160
                math::Matrix<doubleType> dsinTA_inf_out_dx = dTA_inf_out_dx * cosTA_inf_out;

                //eq 123
                math::Matrix<doubleType> dBhat_dB = (I3 - Bvector * Bvector.transpose() / BRADIUS / BRADIUS) / BRADIUS;

                //derivatives of Shat
                //eq 93
                math::Matrix<doubleType> dShat_dx(3, 6, 0.0);
                dShat_dx(0, 1) = -cosDHA * sinRHA;
                dShat_dx(1, 1) = cosDHA * cosRHA;
                dShat_dx(0, 2) = -sinDHA * cosRHA;
                dShat_dx(1, 2) = -sinDHA * sinRHA;
                dShat_dx(2, 2) = cosDHA;

                //derivatives of That
                //eq 97
                math::Matrix<doubleType> Tvec = Shat.cross(phi);
                doubleType T = Tvec.norm();

                //eq 98
                math::Matrix<doubleType> dTvec_dx = -math::Matrix<doubleType>(phi, math::skewsymmetric) * dShat_dx; // + Shat.skew(dphi_dx) but phi is a constant so this term is zero

                //eq 99
                math::Matrix<doubleType> Xi2 = (Tvec.transpose() * dTvec_dx * -1.0 / T / T / T).transpose();

                //eq 100
                math::Matrix<doubleType> dThat_dx = dTvec_dx / T + Tvec * Xi2.transpose();

                //derivatives of Rhat
                //eq 101
                math::Matrix<doubleType> Rvec = Shat.cross(That);
                doubleType R = Rvec.norm();

                //eq 102
                math::Matrix<doubleType> dRvec_dx = -math::Matrix<doubleType>(That, math::skewsymmetric) * dShat_dx
                                                  +  math::Matrix<doubleType>(Shat, math::skewsymmetric) * dThat_dx;

                //eq 103
                Xi2 = (Rvec.transpose() * dRvec_dx * -1.0 / R / R / R).transpose();

                //eq 104
                math::Matrix<doubleType> dRhat_dx = dRvec_dx / R + Rvec * Xi2.transpose();

                //derivatives of B vector
                //eq 94
                math::Matrix<doubleType> dsinBTHETA_dxb(1, 6, std::vector<doubleType>({ 0.0, 0.0, 0.0, 0.0, cosBTHETA, 0.0 }));

                //eq 95
                math::Matrix<doubleType> dcosBTHETA_dxb(1, 6, std::vector<doubleType>({ 0.0, 0.0, 0.0, 0.0, -sinBTHETA, 0.0 }));

                //eq 96
                math::Matrix<doubleType> dBvector_dx = Rhat * dBRADIUS_dx * sinBTHETA
                                                     + Rhat * dsinBTHETA_dxb * BRADIUS
                                                     + dRhat_dx * BRADIUS * sinBTHETA
                                                     + That * dBRADIUS_dx * cosBTHETA
                                                     + That * BRADIUS * dcosBTHETA_dxb
                                                     + dThat_dx * BRADIUS * cosBTHETA;

                //eq 124
                math::Matrix<doubleType> dBhat_dx = dBhat_dB * dBvector_dx;
                
                //eq 158
                math::Matrix<doubleType> dEhat_dx = dShat_dx * cosTA_inf_out
                                                  + Shat * dcosTA_inf_out_dx
                                                  + dBhat_dx * sinTA_inf_out
                                                  + Bhat * dsinTA_inf_out_dx;

                ///and now back to the derivatives of eccentricity...
                //eq 87
                math::Matrix<doubleType> dEvec_dx = Ehat * dECC_dx + dEhat_dx * ECC;

                //derivatives of angular momentum magnitude
                //eq 90
                math::Matrix<doubleType> dh_dx(1, 6, { BRADIUS, 0.0, 0.0, VINF, 0.0, 0.0 });

                //angular momentum vector
                //eq 91
                math::Matrix<doubleType> Gamma = Bvector.cross(Shat);
                doubleType gamma = Gamma.norm();

                //eq 92
                math::Matrix<doubleType> Sskew = math::Matrix<doubleType>(Shat, math::skewsymmetric);
                math::Matrix<doubleType> Bskew = math::Matrix<doubleType>(Bvector, math::skewsymmetric);
                math::Matrix<doubleType> dHhat_dx = (-Gamma * Gamma.transpose() / gamma / gamma / gamma + I3 / gamma)
                    * (-Sskew * dBvector_dx + Bskew * dShat_dx);

                //eq 89
                math::Matrix<doubleType> dHvec_dx = Hhat * dh_dx + dHhat_dx * h;

                //derivatives of position wrt angular momentum
                //eq 130
                math::Matrix<doubleType> Xi1 = Ehat * Hvec.transpose() * 2 * cosTA;
                
                //eq 131 - first we need Pvec and P
                math::Matrix<doubleType> Pvec = Hvec.cross(Evec);
                doubleType P = Pvec.norm();
                math::Matrix<doubleType> Eskew(Evec, math::skewsymmetric);
                Xi2 = (Pvec * Hvec.transpose() * 2.0 + (-I3 + Pvec * Pvec.transpose() / P / P) * Eskew * h * h ) * sinTA / P;

                //eq 132
                math::Matrix<doubleType> dpos_dHvec = (Xi1 + Xi2) / (mu * (1 + ECC * cosTA));

                //derivatives of position vector wrt eccentricity vector
                //eq 133
                Xi1 = (Ehat * cosTA + Phat * sinTA) * (Ehat.transpose() * (-cosTA / (1 + ECC * cosTA) / (1 + ECC * cosTA)));

                //eq 134
                Xi2 = ((I3 - Evec * Evec.transpose() / ECC / ECC) * cosTA / ECC + (I3 - Pvec * Pvec.transpose() / P / P) * math::Matrix<doubleType>(Hvec, math::skewsymmetric) * sinTA / P)
                    / (1 + ECC * cosTA);

                //eq 135
                math::Matrix<doubleType> dpos_dEvec = (Xi1 + Xi2) * h * h / mu;

                //derivatives of position vector wrt true anomaly
                //eq 136
                Xi1 = (Ehat * cosTA + Phat * sinTA) * ECC * sinTA / (1 + ECC * cosTA) / (1 + ECC * cosTA);

                //eq 137
                Xi2 = (-Ehat * sinTA + Phat * cosTA) / (1 + ECC * cosTA);

                //eq 138
                math::Matrix<doubleType> dpos_dTA = (Xi1 + Xi2) * h * h / mu;

                //.....and the position derivatives themselves!
                //eq 139
                math::Matrix<doubleType> dpos_dx = dpos_dHvec * dHvec_dx + dpos_dEvec * dEvec_dx + dpos_dTA * dTA_dx;

                //derivatives of velocity vector wrt angular momentum vector
                //eq 140
                Xi1 = (Ehat * sinTA - Phat * (ECC + cosTA)) * Hvec.transpose() * -1.0 / h / h / h;

                //eq 141
                Xi2 = (-Eskew + Pvec * Pvec.transpose() * Eskew / P / P) * -(ECC + cosTA) / h / P;

                //eq 142
                math::Matrix<doubleType> dvel_dHvec = -(Xi1 + Xi2) * mu;

                //derivatives of velocity wrt eccentricity vector
                //eq 143
                Xi1 = (I3 - Evec * Evec.transpose() / ECC / ECC) * sinTA / ECC;

                //eq 144
                math::Matrix<doubleType> Xi21 = Phat * Ehat.transpose();

                //eq 145
                math::Matrix<doubleType> Hskew = math::Matrix<doubleType>(Hvec, math::skewsymmetric);
                math::Matrix<doubleType> Xi22 = (Hskew - Pvec * Pvec.transpose() * Hskew / P / P) * (ECC + cosTA) / P;

                //eq 146
                Xi2 = -(Xi21 + Xi22);

                //eq 147
                math::Matrix<doubleType> dvel_dEvec = -(Xi1 + Xi2) * mu / h;

                //derivatives of velocity wrt true anomaly
                //eq 148
                Xi1 = Ehat * cosTA;
                
                //eq 149
                Xi2 = Phat * sinTA;

                //eq 150
                math::Matrix<doubleType> dvel_dTA = -(Xi1 + Xi2) * mu / h;

                //...and the velocity derivative itself!
                //eq 139
                math::Matrix<doubleType> dvel_dx = dvel_dHvec * dHvec_dx + dvel_dEvec * dEvec_dx + dvel_dTA * dTA_dx;

                //finally, put it all in the big matrix
                for (size_t i = 0; i < 3; ++i)
                {
                    for (size_t j = 0; j < 6; ++j)
                    {
                        this->RepresentationToCartesianTransitionMatrix(i, j) = dpos_dx(i, j);
                        this->RepresentationToCartesianTransitionMatrix(i + 3, j) = dvel_dx(i, j);
                    }
                }
            }

            //Step 4: return the converted vector
            return this->StateVectorCartesian;
        }//end convertFromRepresentationToCartesian()

        math::Matrix<doubleType> OutgoingBplaneStateRepresentation::convertFromCartesianToRepresentation(const math::Matrix<doubleType>& StateVectorCartesianIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorCartesian(StateVectorCartesianIn);

            //need a 3x3 identity
            math::Matrix<doubleType> I3(3, math::identity);

            //Step 2: convert the state           
            //this is trivial for Cartesian
            const doubleType& x = this->StateVectorCartesian(0);
            const doubleType& y = this->StateVectorCartesian(1);
            const doubleType& z = this->StateVectorCartesian(2);
            const doubleType& vx = this->StateVectorCartesian(3);
            const doubleType& vy = this->StateVectorCartesian(4);
            const doubleType& vz = this->StateVectorCartesian(5);

            math::Matrix<doubleType> rvec = this->StateVectorCartesian.getSubMatrix1D(0, 2);
            math::Matrix<doubleType> vvec = this->StateVectorCartesian.getSubMatrix1D(3, 5);
            doubleType r = rvec.norm();
            doubleType v = vvec.norm();

            //eq 3
            math::Matrix<doubleType> Evec = (rvec * (v * v - mu / r) - vvec * (rvec.transpose() * vvec)) / mu;
            math::Matrix<doubleType> Ehat = Evec.unitize();
            doubleType ECC = Evec.norm();

            if (ECC < 1.0 + math::SMALL)
            {
                throw std::runtime_error("Error converting from Cartesian to IncomingBplane. IncomingBplaneStateRepresentation can only work if the state has an ECC of > 1.0. ECC is currently " + std::to_string(ECC  _GETVALUE) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            //eq 4
            math::Matrix<doubleType> Hvec = rvec.cross(vvec);
            doubleType h = Hvec.norm();

            //eq 5
            math::Matrix<doubleType> Pvec = Hvec.cross(Evec);
            math::Matrix<doubleType> Phat = Pvec.unitize();

            //eq 151
            math::Matrix<doubleType> Shat = -Ehat / ECC + Phat * sqrt(1.0 - 1.0 / ECC / ECC);

            //eq 7
            //first we need the reference vector phi
            math::Matrix<doubleType> phi(3, 1, { 0.0, 0.0, 1.0 });
            math::Matrix<doubleType> Tvec = Shat.cross(phi);
            math::Matrix<doubleType> That = Tvec.unitize();

            //eq 8
            math::Matrix<doubleType> Rvec = Shat.cross(Tvec);
            math::Matrix<doubleType> Rhat = Rvec.unitize();

            //eq 152:
            math::Matrix<doubleType> Bhat = Phat / ECC + Ehat * sqrt(1.0 - 1.0 / ECC / ECC);;

            //eq 11
            doubleType BRADIUS = h * h / mu / sqrt(ECC * ECC - 1.0);
            
            //eq 12
            math::Matrix<doubleType> Bvec = Bhat * BRADIUS;
            this->BdotT = (Bvec.transpose() * That)(0); //this is a matrix equation but a scalar comes out

            //eq 13
            this->BdotR = (Bvec.transpose() * Rhat)(0); //this is a matrix equation but a scalar comes out

            //eq 14:
            doubleType BTHETA = atan2(this->BdotR, this->BdotT);

            //eq 15:
            doubleType VINF = sqrt(v * v - 2.0 * mu / r);

            //eq 16:
            doubleType rp = mu * (ECC - 1.0) / VINF * VINF;

            //eq 17
            doubleType RHA = atan2(Shat(1), Shat(0));

            //eq 18
            doubleType DHA = math::safe_asin(Shat(2));

            //eq 19
            doubleType TA = atan2((Evec.cross(rvec)).norm(), Evec.dot(rvec));

            doubleType rdotv = rvec.dot(vvec);

            if (rdotv < 0.0)
            {
                TA = math::TwoPI - TA;
            }

            //put it in the state vector

            this->StateVectorThisRepresentation(0) = VINF;
            this->StateVectorThisRepresentation(1) = RHA;
            this->StateVectorThisRepresentation(2) = DHA;
            this->StateVectorThisRepresentation(3) = BRADIUS;
            this->StateVectorThisRepresentation(4) = BTHETA;
            this->StateVectorThisRepresentation(5) = TA;

            

            //Step 3: construct the partial derivatives
            if (needG)
            {
                this->CartesianToRepresentationTransitionMatrix.assign_zeros();

                //once again following Noble's math spec

                //eq 21 - rvec wrt x
                math::Matrix<doubleType> drvec_dx(3, 6, 0.0);
                drvec_dx(0, 0) = 1.0;
                drvec_dx(1, 1) = 1.0;
                drvec_dx(2, 2) = 1.0;

                //eq 22 - vvec wrt x
                math::Matrix<doubleType> dvvec_dx(3, 6, 0.0);
                dvvec_dx(0, 3) = 1.0;
                dvvec_dx(1, 4) = 1.0;
                dvvec_dx(2, 5) = 1.0;

                //derivatives of eccentricity vector wrt x
                //eq 62
                math::Matrix<doubleType> Sigma1 = rvec.transpose() * dvvec_dx + vvec.transpose() * drvec_dx;

                //eq 63
                math::Matrix<doubleType> dr_drvec = rvec.transpose() / r;

                //eq 64
                math::Matrix<doubleType> dv_dvvec = vvec.transpose() / v;

                //eq 65
                math::Matrix<doubleType> Xi1 = rvec * (dv_dvvec * dvvec_dx * 2.0 * v + dr_drvec * drvec_dx * mu / r / r) + drvec_dx * (v*v - mu / r);

                //eq 66
                math::Matrix<doubleType> Xi2 = vvec * Sigma1 + dvvec_dx * rvec.dot(vvec);

                //eq 67
                math::Matrix<doubleType> dEvec_dx = (Xi1 - Xi2) / mu;

                //derivatives of the angular momentum vector
                //eq 60
                math::Matrix<doubleType> Rskew(rvec, math::skewsymmetric);
                math::Matrix<doubleType> Vskew(vvec, math::skewsymmetric);
                math::Matrix<doubleType> dHvec_dx = -Vskew.horz_cat(-Rskew);

                //derivatives of the P vector
                //eq 61
                math::Matrix<doubleType> Eskew(Evec, math::skewsymmetric);
                math::Matrix<doubleType> Hskew(Hvec, math::skewsymmetric);
                math::Matrix<doubleType> dPvec_dx = -Eskew * dHvec_dx + Hskew * dEvec_dx;

                //derivatives of Shat
                //eq 52
                doubleType P = Pvec.norm();
                math::Matrix<doubleType> dPhat_dx = (I3 - Pvec * Pvec.transpose() / P / P) * dPvec_dx / P;

                //eq 53
                doubleType sigma1 = sqrt(1.0 - 1.0 / ECC / ECC);

                //eq 54
                math::Matrix<doubleType> dsigma1_dx = Ehat.transpose() * dEvec_dx / ECC / ECC / ECC / sigma1;

                //eq 55
                Xi1 = -Ehat * Ehat.transpose() * dEvec_dx / ECC / ECC;

                //eq 37, which we need out of order
                math::Matrix<doubleType> dEhat_dx = (I3 - Evec * Evec.transpose() / ECC / ECC) * dEvec_dx / ECC;

                //eq 56
                Xi2 = dEhat_dx / ECC;

                //eq 57
                math::Matrix <doubleType> Xi3 = Phat * dsigma1_dx;

                //eq 58
                math::Matrix<doubleType> Xi4 = dPhat_dx * sigma1;

                //eq 59
                math::Matrix<doubleType> dShat_dx = -Xi1 - Xi2 + Xi3 + Xi4;


                //RHA wrt x
                //eq 23
                math::Matrix<doubleType> dRHA_dShat(1, 3, { -Shat(1) / (Shat(0)*Shat(0) + Shat(1)*Shat(1)), Shat(0) / (Shat(0)*Shat(0) + Shat(1)*Shat(1)), 0.0 });

                //eq 24
                math::Matrix<doubleType> dRHA_dx = dRHA_dShat * dShat_dx;

                //DHA wrt x
                //eq 25
                math::Matrix<doubleType> dDHA_dShat(1, 3, { 0.0, 0.0, 1.0 / sqrt(1.0 - Shat(2) * Shat(2)) });

                //eq 26
                math::Matrix<doubleType> dDHA_dx = dDHA_dShat * dShat_dx;

                //derivatives of VINF
                //eq 27
                math::Matrix<doubleType> dVINF_dx = (vvec.transpose() * dvvec_dx + rvec.transpose() * drvec_dx * mu / r / r / r) / sqrt(v * v - 2.0 * mu / r);

                //derivatives of BRADIUS
                //eq 32
                doubleType e2m1 = ECC * ECC - 1.0; //"e squared minus 1"
                Xi1 = -Evec.transpose() * dEvec_dx / sqrt(e2m1 * e2m1 * e2m1);

                //eq 33
                Xi2 = Hvec.transpose() * dHvec_dx * 2.0;

                //eq 34
                math::Matrix<doubleType> dBRADIUS_dx = (Xi1 * h * h + Xi2 / sqrt(e2m1)) / mu;

                //derivatives of Bhat
                //eq 36, 37, 38, and 39 are duplicates of things we've already computed
                
                //eq 40
                math::Matrix<doubleType> Sigma2 = -Ehat.transpose() * dEvec_dx / ECC / ECC;

                //eq 41
                Xi1 = dEhat_dx * sigma1;

                //eq 42
                Xi2 = Ehat * dsigma1_dx;

                //eq 43
                Xi3 = -dPhat_dx / ECC;

                //eq 44
                Xi4 = -Phat * Sigma2;

                //eq 45
                math::Matrix<doubleType> dBhat_dx = Xi1 + Xi2 - Xi3 - Xi4;

                //derivatives of B vector
                //eq 35
                math::Matrix<doubleType> dBvec_dx = Bhat * dBRADIUS_dx + dBhat_dx * BRADIUS;

                //derivatives of That
                //eq 49 is already done
                //eq 50
                math::Matrix<doubleType> Sskew(Shat, math::skewsymmetric);
                math::Matrix<doubleType> phiskew(phi, math::skewsymmetric);
                math::Matrix<doubleType> dTvec_dx = -phiskew * dShat_dx;// + Sskew * dphi_dx, but dphi_dx is zero

                //eq 51
                doubleType T = Tvec.norm();
                math::Matrix<doubleType> dThat_dx = (I3 - Tvec * Tvec.transpose() / T / T) * dTvec_dx / T;

                //derivatives of Rhat
                //eq 46 is already done
                //eq 47
                math::Matrix<doubleType> Tskew(Tvec, math::skewsymmetric);
                math::Matrix<doubleType> dRvec_dx = -Tskew * dShat_dx + Sskew * dTvec_dx;

                //eq 48
                doubleType R = Rvec.norm();
                math::Matrix<doubleType> dRhat_dx = (I3 - Rvec * Rvec.transpose() / R / R) * dRvec_dx / R;

                //derivatives of BdotT
                //eq 30
                this->dBdotT_dx_cartesian = That.transpose() * dBvec_dx + Bvec.transpose() * dThat_dx;

                //derivatives of BdotR
                //eq 31
                this->dBdotR_dx_cartesian = Rhat.transpose() * dBvec_dx + Bvec.transpose() * dRhat_dx;

                //derivatives of rp
                //eq 28
                math::Matrix<doubleType> drp_dx = (Ehat.transpose() * dEvec_dx / VINF / VINF - dVINF_dx * 2.0 * (ECC - 1.0) / VINF / VINF / VINF) / mu;


                //derivatives of BTHETA
                //eq 29
                math::Matrix<doubleType> dBTHETA_dx = this->dBdotR_dx_cartesian * BdotT / (BdotR * BdotR + BdotT * BdotT)
                                                    - this->dBdotT_dx_cartesian * BdotR / (BdotR * BdotR + BdotT * BdotT);

                //derivatives of TA - borrowed from the COE derivatives
                //derivatives of true anomaly
                //TA = atan2(B, A)
                //dTA_dx = datan2(B,A)_dA * dA_dx + datan2(B,A)_dB * dB_dx

                doubleType A = Evec.dot(rvec);

                doubleType B = (Evec.cross(rvec)).norm();

                math::Matrix<doubleType> EcrossR = Evec.cross(rvec);

                doubleType dTA_dA = -B / (A * A + B * B);

                doubleType dTA_dB = A / (A * A + B * B);

                math::Matrix<doubleType> dA_dx = Evec.transpose() * drvec_dx + rvec.transpose() * dEvec_dx;

                math::Matrix<doubleType> dEcrossR_dEvec = -Rskew;

                math::Matrix<doubleType> dEcrossR_drvec = Eskew;

                math::Matrix<doubleType> dEcrossR_dx = dEcrossR_dEvec * dEvec_dx + dEcrossR_drvec * drvec_dx;

                math::Matrix<doubleType> dB_dx = EcrossR.transpose() * dEcrossR_dx / EcrossR.norm();

                math::Matrix<doubleType> dTA_dx = dA_dx * dTA_dA + dB_dx * dTA_dB;
                
                if (rdotv < 0.0)
                {
                    dTA_dx *= -1.0;
                }

                //now put it all in the transformation matrix
                for (size_t i : {0, 1, 2, 3, 4, 5})
                {
                    this->CartesianToRepresentationTransitionMatrix(0, i) = dVINF_dx(i);
                    this->CartesianToRepresentationTransitionMatrix(1, i) = dRHA_dx(i);
                    this->CartesianToRepresentationTransitionMatrix(2, i) = dDHA_dx(i);
                    this->CartesianToRepresentationTransitionMatrix(3, i) = dBRADIUS_dx(i);
                    this->CartesianToRepresentationTransitionMatrix(4, i) = dBTHETA_dx(i);
                    this->CartesianToRepresentationTransitionMatrix(5, i) = dTA_dx(i);
                }
            }//end derivatives

            //Step 4: return the converted vector
            return this->StateVectorThisRepresentation;
        }//end convertFromCartesianToRepresentation()
    }//end namespace StateRepresentation
}//end namespace EMTG