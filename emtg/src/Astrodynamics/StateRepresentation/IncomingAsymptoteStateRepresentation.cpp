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

#pragma once

#include "IncomingAsymptoteStateRepresentation.h"
#include "COEStateRepresentation.h"
#include "EMTG_math.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        IncomingAsymptote::IncomingAsymptote(const double& mu) : StateRepresentationBase::StateRepresentationBase(mu)
        {
            this->name = "IncomingAsymptote";
            this->stateNames = std::vector<std::string>({ "RP", "C3", "INC", "RHA", "DHA", "TA" });

            this->myCOEStateRepresentation = COEStateRepresentation(mu);
        };

        math::Matrix<doubleType> IncomingAsymptote::convertFromRepresentationToCartesian(const math::Matrix<doubleType>& StateVectorThisRepresentationIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorThisRepresentation(StateVectorThisRepresentationIn);

            //Step 2: convert the state
            //Step 2.1: one cannot convert directly from IncomingAsymptote to Cartesian. Instead we convert to COE first and then from there to Cartesian
            math::Matrix<doubleType> stateCOE;

            doubleType RP =  StateVectorThisRepresentationIn(0);
            doubleType C3 =  StateVectorThisRepresentationIn(1);
            doubleType RHA = StateVectorThisRepresentationIn(2);
            doubleType DHA = StateVectorThisRepresentationIn(3);
            doubleType BVA = StateVectorThisRepresentationIn(4);
            doubleType TA =  StateVectorThisRepresentationIn(5);

            doubleType cosRHA = cos(RHA);
            doubleType sinRHA = sin(RHA);
            doubleType cosDHA = cos(DHA);
            doubleType sinDHA = sin(DHA);
            doubleType cosBVA = cos(BVA);
            doubleType sinBVA = sin(BVA);
            doubleType cosTA = cos(TA);
            doubleType sinTA = sin(TA);

            if (C3 < 1e-7)
            {
                std::cout << "Warning, IncomingAsymptote orbit is elliptic so using apsides vector for asymptote, C3 < 1e-7" << std::endl;
            }

            doubleType SMA = -mu / C3;
            doubleType ECC = 1 - RP / SMA;

            if (ECC < 1e-7)
            {
                throw std::runtime_error("Illegal conversion from IncomingAsymptote to Keplerian elements, ECC < 1e-7. IncomingAsymptote is undefined on a circular orbit. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            math::Matrix<doubleType> sVhat(3, 1);
            sVhat(0) = cosDHA * cosRHA;
            sVhat(1) = cosDHA * sinRHA;
            sVhat(2) = sinDHA;
            doubleType khati[] = { 0.0,0.0,1.0 };
            math::Matrix<doubleType> khat(3, 1, khati);

            if (1.0 - sVhat.dot(khat) < 1.0e-7)
            {
                throw std::runtime_error("Illegal conversion from IncomingAsymptote to Keplerian elements, 1.0 - vhat.dot(khat) < 1e-7. IncomingAsymptote is undefined when asymptote is aligned with khat. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            math::Matrix<doubleType> eaVec = khat.cross(sVhat);
            math::Matrix<doubleType> eaVhat = eaVec.unitize();
            math::Matrix<doubleType> noVhat = sVhat.cross(eaVhat);
            doubleType AMI = math::TwoPI / 4 - BVA; // angular azimuth at infinity 
            doubleType cosAMI = cos(AMI);
            doubleType sinAMI = sin(AMI);
            math::Matrix<doubleType> hVhat = eaVhat * sinAMI + noVhat * cosAMI;
            math::Matrix<doubleType> nodeVec = khat.cross(hVhat);
            doubleType nMag = nodeVec.norm();
            math::Matrix<doubleType> eccVhat;
            doubleType TAmax;
            doubleType sinTAmax;
            doubleType cosTAmax;
            math::Matrix<doubleType> oVhat;

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

            doubleType khatdothVhat = khat.dot(hVhat);
            doubleType INC = acos(khatdothVhat);
            doubleType RAAN;
            doubleType AOP;

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
            stateCOE(0) = SMA;
            stateCOE(1) = ECC;
            stateCOE(2) = INC;
            stateCOE(3) = RAAN;
            stateCOE(4) = AOP;
            stateCOE(5) = TA;

            //Step 2.2: convert from COE to Cartesian
            this->StateVectorCartesian = this->myCOEStateRepresentation.convertFromRepresentationToCartesian(stateCOE);


            //Step 3: construct the partial derivatives
            if (needG)
            {
                this->RepresentationToCartesianTransitionMatrix.assign_zeros();

                //Step 3.1: once again, we have to convert to COE first and then from COE to cartesian
                math::Matrix<doubleType> TransformationFromIncomingAsymptoteToCOE(6, 6, math::identity);

                //SMA
                TransformationFromIncomingAsymptoteToCOE(0, 0) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(1, 0) = mu / (C3 * C3)_GETVALUE;
                TransformationFromIncomingAsymptoteToCOE(2, 0) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(3, 0) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(4, 0) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(5, 0) = 0.0;

                //ECC
                TransformationFromIncomingAsymptoteToCOE(0, 1) = (C3 / mu)_GETVALUE;
                TransformationFromIncomingAsymptoteToCOE(1, 1) = (RP / mu)_GETVALUE;
                TransformationFromIncomingAsymptoteToCOE(2, 1) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(3, 1) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(4, 1) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(5, 1) = 0.0;

                //INC
                //first we need derivatives of eaVec
                math::Matrix<double> deaVec_dIncomingAsymptote(6, 3, 0.0);
                //no dependencies of eaVec on RP and C3
                deaVec_dIncomingAsymptote(2, 0) = (-cosDHA * cosRHA)_GETVALUE;
                deaVec_dIncomingAsymptote(2, 1) = (-cosDHA * sinRHA)_GETVALUE;
                deaVec_dIncomingAsymptote(3, 0) = (sinDHA * sinRHA)_GETVALUE;
                deaVec_dIncomingAsymptote(3, 1) = (-cosRHA * sinDHA)_GETVALUE;
                //no dependencies of eaVec on BVA or TA

                //then derivatives of eaVhat
                doubleType ea = eaVec.norm();
                double deascalar_dRHA = 0;
                double deascalar_dDHA = -math::sgn(cosDHA) * sinDHA _GETVALUE;

                math::Matrix<double> deaVhat_dIncomingAsymptote(6, 3, 0.0);
                deaVhat_dIncomingAsymptote(2, 0) = (1 / (ea * ea) * (deaVec_dIncomingAsymptote(2, 0) * ea - deascalar_dRHA * eaVec(0)))_GETVALUE;
                deaVhat_dIncomingAsymptote(2, 1) = (1 / (ea * ea) * (deaVec_dIncomingAsymptote(2, 1) * ea - deascalar_dRHA * eaVec(1)))_GETVALUE;
                deaVhat_dIncomingAsymptote(3, 0) = (1 / (ea * ea) * (deaVec_dIncomingAsymptote(3, 0) * ea - deascalar_dDHA * eaVec(0)))_GETVALUE;
                deaVhat_dIncomingAsymptote(3, 1) = (1 / (ea * ea) * (deaVec_dIncomingAsymptote(3, 1) * ea - deascalar_dDHA * eaVec(1)))_GETVALUE;

                //and derivatives of noVhat
                math::Matrix<doubleType> noVec = sVhat.cross(eaVec);
                math::Matrix<double> dnoVec_dRHA(3, 1, 0.0);
                math::Matrix<double> dnoVec_dDHA(3, 1, 0.0);
                dnoVec_dRHA(0) = (cosDHA * sinDHA * sinRHA)_GETVALUE;
                dnoVec_dRHA(1) = (-cosDHA * cosRHA * sinDHA)_GETVALUE;
                dnoVec_dDHA(0) = (cosRHA * sinDHA * sinDHA - cosDHA * cosDHA * cosRHA)_GETVALUE;
                dnoVec_dDHA(1) = (sinDHA * sinDHA * sinRHA - cosDHA * cosDHA * sinRHA)_GETVALUE;
                dnoVec_dDHA(2) = (-2 * cosDHA * cosRHA * cosRHA * sinDHA - 2 * cosDHA * sinDHA * sinRHA * sinRHA)_GETVALUE;
                math::Matrix<double> dnoVhat_dIncomingAsymptote(6, 3, 0.0);

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
                math::Matrix<double> dhVhat_dIncomingAsymptote(6, 3, 0.0);
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
                double INCmultiplier = -1.0 / sqrt(1.0 - khatdothVhat * khatdothVhat)_GETVALUE;
                TransformationFromIncomingAsymptoteToCOE(0, 2) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(1, 2) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(2, 2) = INCmultiplier * dhVhat_dIncomingAsymptote(2, 2);
                TransformationFromIncomingAsymptoteToCOE(3, 2) = INCmultiplier * dhVhat_dIncomingAsymptote(3, 2);
                TransformationFromIncomingAsymptoteToCOE(4, 2) = INCmultiplier * dhVhat_dIncomingAsymptote(4, 2);
                TransformationFromIncomingAsymptoteToCOE(5, 2) = 0.0;

                //RAAN
                //first we need derivatives of nodeVec
                math::Matrix<double> dnodeVec_dIncomingAsymptote(6, 3, 0.0);
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
                math::Matrix<double> dnMag_dIncomingAsymptote(6, 1, 0.0);
                for (size_t i = 0; i < 6; ++i)
                {
                    dnMag_dIncomingAsymptote(i) = (1.0 / nMag * (nodeVec(0) * dnodeVec_dIncomingAsymptote(i, 0) + nodeVec(1) * dnodeVec_dIncomingAsymptote(i, 1)))_GETVALUE;
                }

                //now we can compute derivatives of RAAN
                if (ECC >= 1E-11 && INC >= 1E-11)  // CASE 1: Non-circular, Inclined Orbit
                {
                    double RAANterm = -(nodeVec(1) < 0.0 ? -1 : 1) / sqrt(1.0 - (nodeVec(0) / nMag)*(nodeVec(0) / nMag))_GETVALUE;
                    TransformationFromIncomingAsymptoteToCOE(0, 3) = 0.0;
                    TransformationFromIncomingAsymptoteToCOE(1, 3) = 0.0;
                    TransformationFromIncomingAsymptoteToCOE(2, 3) = (RAANterm / (nMag * nMag) * (dnodeVec_dIncomingAsymptote(2, 0) * nMag - dnMag_dIncomingAsymptote(2) * nodeVec(0)))_GETVALUE;
                    TransformationFromIncomingAsymptoteToCOE(3, 3) = (RAANterm / (nMag * nMag) * (dnodeVec_dIncomingAsymptote(3, 0) * nMag - dnMag_dIncomingAsymptote(3) * nodeVec(0)))_GETVALUE;
                    TransformationFromIncomingAsymptoteToCOE(4, 3) = (RAANterm / (nMag * nMag) * (dnodeVec_dIncomingAsymptote(4, 0) * nMag - dnMag_dIncomingAsymptote(4) * nodeVec(0)))_GETVALUE;
                    TransformationFromIncomingAsymptoteToCOE(5, 3) = 0.0;
                }
                else
                {
                    TransformationFromIncomingAsymptoteToCOE(0, 3) = 0.0;
                    TransformationFromIncomingAsymptoteToCOE(1, 3) = 0.0;
                    TransformationFromIncomingAsymptoteToCOE(2, 3) = 0.0;
                    TransformationFromIncomingAsymptoteToCOE(3, 3) = 0.0;
                    TransformationFromIncomingAsymptoteToCOE(4, 3) = 0.0;
                    TransformationFromIncomingAsymptoteToCOE(5, 3) = 0.0;
                }

                //AOP
                //do get AOP derivatives we need derivatives of eccVhat
                math::Matrix<double> dsVhat_dIncomingAsymptote(6, 3, 0.0);
                dsVhat_dIncomingAsymptote(2, 0) = (-cosDHA * sinRHA)_GETVALUE;
                dsVhat_dIncomingAsymptote(2, 1) = (cosDHA*cosRHA)_GETVALUE;
                dsVhat_dIncomingAsymptote(3, 0) = (-cosRHA * sinDHA)_GETVALUE;
                dsVhat_dIncomingAsymptote(3, 1) = (-sinDHA * sinRHA)_GETVALUE;
                dsVhat_dIncomingAsymptote(3, 2) = (cosDHA)_GETVALUE;
                math::Matrix<double> deccVhat_dIncomingAsymptote(6, 3, 0.0);
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
                    math::Matrix<double> doVhat_dIncomingAsymptote(6, 3);
                    for (size_t i = 0; i < 6; ++i)
                    {
                        doVhat_dIncomingAsymptote(i, 0) = (dhVhat_dIncomingAsymptote(i, 1) * sVhat(2) + hVhat(1) * dsVhat_dIncomingAsymptote(i, 2) - dhVhat_dIncomingAsymptote(i, 2) * sVhat(1) - hVhat(2) * dsVhat_dIncomingAsymptote(i, 1))_GETVALUE;
                        doVhat_dIncomingAsymptote(i, 1) = (dhVhat_dIncomingAsymptote(i, 2) * sVhat(0) + hVhat(2) * dsVhat_dIncomingAsymptote(i, 0) - dhVhat_dIncomingAsymptote(i, 0) * sVhat(2) - hVhat(0) * dsVhat_dIncomingAsymptote(i, 2))_GETVALUE;
                        doVhat_dIncomingAsymptote(i, 2) = (dhVhat_dIncomingAsymptote(i, 0) * sVhat(1) + hVhat(0) * dsVhat_dIncomingAsymptote(i, 1) - dhVhat_dIncomingAsymptote(i, 1) * sVhat(0) - hVhat(1) * dsVhat_dIncomingAsymptote(i, 0))_GETVALUE;
                    }

                    //need derivatives of TAmax
                    math::Matrix<double> dTAmax_dIncomingAsymptote(6, 1, 0.0);
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
                    doubleType nodeVecdoteccVhat = nodeVec.dot(eccVhat);
                    math::Matrix<double> dnodeVecdoteccVhat_dIncomingAsymptote(6, 1, 0.0);
                    for (size_t i = 0; i < 6; ++i)
                    {
                        dnodeVecdoteccVhat_dIncomingAsymptote(i) = (nodeVec(0) * deccVhat_dIncomingAsymptote(i, 0) + eccVhat(0) * dnodeVec_dIncomingAsymptote(i, 0)
                            + nodeVec(1) * deccVhat_dIncomingAsymptote(i, 1) + eccVhat(1) * dnodeVec_dIncomingAsymptote(i, 1)
                            + nodeVec(2) * deccVhat_dIncomingAsymptote(i, 2) + eccVhat(2) * dnodeVec_dIncomingAsymptote(i, 2))_GETVALUE;
                    }
                    double AOPterm = -(eccVhat(2) < 0.0 ? -1 : 1) / sqrt(1.0 - (nodeVecdoteccVhat / nMag)*(nodeVecdoteccVhat / nMag))_GETVALUE;
                    TransformationFromIncomingAsymptoteToCOE(0, 4) = (AOPterm / (nMag * nMag) * (dnodeVecdoteccVhat_dIncomingAsymptote(0) * nMag - dnMag_dIncomingAsymptote(0) * nodeVecdoteccVhat))_GETVALUE;
                    TransformationFromIncomingAsymptoteToCOE(1, 4) = (AOPterm / (nMag * nMag) * (dnodeVecdoteccVhat_dIncomingAsymptote(1) * nMag - dnMag_dIncomingAsymptote(1) * nodeVecdoteccVhat))_GETVALUE;
                    TransformationFromIncomingAsymptoteToCOE(2, 4) = (AOPterm / (nMag * nMag) * (dnodeVecdoteccVhat_dIncomingAsymptote(2) * nMag - dnMag_dIncomingAsymptote(2) * nodeVecdoteccVhat))_GETVALUE;
                    TransformationFromIncomingAsymptoteToCOE(3, 4) = (AOPterm / (nMag * nMag) * (dnodeVecdoteccVhat_dIncomingAsymptote(3) * nMag - dnMag_dIncomingAsymptote(3) * nodeVecdoteccVhat))_GETVALUE;
                    TransformationFromIncomingAsymptoteToCOE(4, 4) = (AOPterm / (nMag * nMag) * (dnodeVecdoteccVhat_dIncomingAsymptote(4) * nMag - dnMag_dIncomingAsymptote(4) * nodeVecdoteccVhat))_GETVALUE;
                    TransformationFromIncomingAsymptoteToCOE(5, 4) = 0.0;
                }
                else if (ECC >= 1E-11 && INC < 1E-7)  // CASE 2: Non-circular, Equatorial Orbit
                {
                    double AOPterm = (eccVhat(0) < 0.0 ? -1 : 1) / sqrt(1.0 - (eccVhat(0) * eccVhat(0)))_GETVALUE;
                    TransformationFromIncomingAsymptoteToCOE(0, 4) = (AOPterm * deccVhat_dIncomingAsymptote(0, 0));
                    TransformationFromIncomingAsymptoteToCOE(1, 4) = (AOPterm * deccVhat_dIncomingAsymptote(1, 0));
                    TransformationFromIncomingAsymptoteToCOE(2, 4) = (AOPterm * deccVhat_dIncomingAsymptote(2, 0));
                    TransformationFromIncomingAsymptoteToCOE(3, 4) = (AOPterm * deccVhat_dIncomingAsymptote(3, 0));
                    TransformationFromIncomingAsymptoteToCOE(4, 4) = (AOPterm * deccVhat_dIncomingAsymptote(4, 0));
                    TransformationFromIncomingAsymptoteToCOE(5, 4) = 0.0;

                }
                else if (ECC > 1E-11 && INC >= math::PI - 1E-7)  // CASE 3: Non-circular, Equatorial Retrograde Orbit
                {
                    double AOPterm = -(eccVhat(1) < 0.0 ? -1 : 1) / sqrt(1.0 - (eccVhat(1) * eccVhat(1)))_GETVALUE;
                    TransformationFromIncomingAsymptoteToCOE(0, 4) = (AOPterm * deccVhat_dIncomingAsymptote(0, 1));
                    TransformationFromIncomingAsymptoteToCOE(1, 4) = (AOPterm * deccVhat_dIncomingAsymptote(1, 1));
                    TransformationFromIncomingAsymptoteToCOE(2, 4) = (AOPterm * deccVhat_dIncomingAsymptote(2, 1));
                    TransformationFromIncomingAsymptoteToCOE(3, 4) = (AOPterm * deccVhat_dIncomingAsymptote(3, 1));
                    TransformationFromIncomingAsymptoteToCOE(4, 4) = (AOPterm * deccVhat_dIncomingAsymptote(4, 1));
                    TransformationFromIncomingAsymptoteToCOE(5, 4) = 0.0;
                }

                //TA
                TransformationFromIncomingAsymptoteToCOE(0, 5) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(1, 5) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(2, 5) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(3, 5) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(4, 5) = 0.0;
                TransformationFromIncomingAsymptoteToCOE(5, 5) = 1.0;

                //Step 3.2: now convert from COE to cartesian
                this->RepresentationToCartesianTransitionMatrix = this->myCOEStateRepresentation.getRepresentationToCartesianTransitionMatrix() * TransformationFromIncomingAsymptoteToCOE;
            }

            //Step 4: return the converted vector
            return this->StateVectorCartesian;
        }//end convertFromRepresentationToCartesian()

        math::Matrix<doubleType> IncomingAsymptote::convertFromCartesianToRepresentation(const math::Matrix<doubleType>& StateVectorCartesianIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorCartesian(StateVectorCartesianIn);

            //Step 2: convert the state
            // Calculate the magnitude of the position vector, right ascension, and declination
            math::Matrix<doubleType> R(3, 1, 0.0);
            math::Matrix<doubleType> V(3, 1, 0.0);
            for (size_t k = 0; k < 3; ++k)
            {
                R(k) = StateVectorCartesianIn(k);
                V(k) = StateVectorCartesianIn(k + 3);
            }

            doubleType r = R.norm();
            doubleType v = V.norm();
            doubleType rdotv = R.dot(V);

            math::Matrix<doubleType> hVec = R.cross(V);
            doubleType hMag = hVec.norm();
            doubleType s = (v*v - mu / r);
            math::Matrix<doubleType> eccVec = (R * s - V * rdotv) / mu;
            doubleType ECC = eccVec.norm();
            doubleType C3 = v * v - 2.0 * mu / r;

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

            doubleType SMA = -mu / C3;
            doubleType RP = SMA * (1 - ECC);

            math::Matrix<doubleType> sVHat(3, 1, 0.0);
            doubleType fac1;
            math::Matrix<doubleType> hcrosse(3, 1, 0.0);

            // sVHat is the asymptote vector of position 
            if (C3 > 1E-7)
            {
                fac1 = 1 / (1 + C3 * hMag*hMag / mu / mu);
                hcrosse = hVec.cross(eccVec);
                sVHat = (-hcrosse * sqrt(C3) / mu - eccVec) * fac1;
            }
            else if (C3 < -1E-7)
            {
                std::cout << "Warning: Orbit is elliptic so using Apsides vector for asymptote." << std::endl;

                sVHat = -eccVec / ECC;
            }


            doubleType khati[] = { 0.0, 0.0, 1.0 };
            math::Matrix<doubleType> khat(3, 1, khati);
            if (acos(sVHat.dot(khat)) < 1E-7) // 1e-7, 1e-11, criteria??
            {
                throw std::runtime_error("Error, IncomingAsymptote vector is aligned with z-direction. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            math::Matrix<doubleType> eaVec = khat.cross(sVHat);
            math::Matrix<doubleType>  eaVhat = eaVec / eaVec.norm();
            math::Matrix<doubleType>  noVhat = sVHat.cross(eaVhat);
            math::Matrix<doubleType>  bVec = hVec.cross(sVHat);
            doubleType sinBVA = bVec.dot(eaVhat) / hMag;
            doubleType cosBVA = bVec.dot(-noVhat) / hMag;
            doubleType BVA = atan2(sinBVA, cosBVA);

            if (BVA < 0.0)
                BVA += math::TwoPI;

            doubleType DHA = asin(sVHat(2));
            doubleType RHA = atan2(sVHat(1), sVHat(0));
            if (RHA < 0.0)
                RHA += math::TwoPI;

            doubleType TA = acos((eccVec.dot(R)) / (ECC*r));
            if (rdotv < 0.0)
                TA = math::TwoPI - TA;

           this->StateVectorThisRepresentation(0) = RP;
           this->StateVectorThisRepresentation(1) = C3;
           this->StateVectorThisRepresentation(2) = RHA;
           this->StateVectorThisRepresentation(3) = DHA;
           this->StateVectorThisRepresentation(4) = BVA;
           this->StateVectorThisRepresentation(5) = TA;

            //Step 3: construct the partial derivatives
            if (needG)
            {
                this->CartesianToRepresentationTransitionMatrix.assign_zeros();

                //derivatives of SMA and ECC, which are necessary to do everything else

                math::Matrix<double> dr_dstate(6, 1, 0.0);
                math::Matrix<double> dv_dstate(6, 1, 0.0);
                dr_dstate(0) = (this->StateVectorCartesian(0) / r)_GETVALUE;
                dr_dstate(1) = (this->StateVectorCartesian(1) / r)_GETVALUE;
                dr_dstate(2) = (this->StateVectorCartesian(2) / r)_GETVALUE;
                dv_dstate(3) = (this->StateVectorCartesian(3) / v)_GETVALUE;
                dv_dstate(4) = (this->StateVectorCartesian(4) / v)_GETVALUE;
                dv_dstate(5) = (this->StateVectorCartesian(5) / v)_GETVALUE;

                //derivatives of SMA
                math::Matrix<double> dSMAdstate(6, 1, 0.0);
                double dSMAdr = (2 * mu * mu / (r*v*v - 2 * mu) / (r*v*v - 2 * mu))_GETVALUE;
                double dSMAdv = (2 * mu * r*r*v / (r*v*v - 2 * mu) / (r*v*v - 2 * mu))_GETVALUE;
                for (size_t k = 0; k < 3; ++k)
                {
                    dSMAdstate(k) = dSMAdr * dr_dstate(k);
                    dSMAdstate(k + 3) = dSMAdv * dv_dstate(k + 3);
                }

                //derivatives of ECC
                //first need derivatives of rdotv and s
                math::Matrix<double> drdotv_dstate(6, 1);
                drdotv_dstate(0) = this->StateVectorCartesian(3)_GETVALUE;
                drdotv_dstate(1) = this->StateVectorCartesian(4)_GETVALUE;
                drdotv_dstate(2) = this->StateVectorCartesian(5)_GETVALUE;
                drdotv_dstate(3) = this->StateVectorCartesian(0)_GETVALUE;
                drdotv_dstate(4) = this->StateVectorCartesian(1)_GETVALUE;
                drdotv_dstate(5) = this->StateVectorCartesian(2)_GETVALUE;
                double dsdv = 2 * v _GETVALUE;
                double dsdr = mu / (r * r)_GETVALUE;
                math::Matrix<double> ds_dstate(6, 1);
                ds_dstate(0) = dsdr * dr_dstate(0);
                ds_dstate(1) = dsdr * dr_dstate(1);
                ds_dstate(2) = dsdr * dr_dstate(2);
                ds_dstate(3) = dsdv * dv_dstate(3);
                ds_dstate(4) = dsdv * dv_dstate(4);
                ds_dstate(5) = dsdv * dv_dstate(5);

                //derivatives of ECC vector
                math::Matrix<double> deccVec_dstate(6, 3);
                //with respect to position
                for (size_t i = 0; i < 3; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                    {
                        deccVec_dstate(i, j) = 1 / mu * (ds_dstate(i) * this->StateVectorCartesian(j)_GETVALUE + (i == j ? s _GETVALUE : 0) - drdotv_dstate(i) * this->StateVectorCartesian(j + 3)_GETVALUE);
                    }
                }
                //with respect to velocity
                for (size_t i = 0; i < 3; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                    {
                        deccVec_dstate(i + 3, j) = 1 / mu * (ds_dstate(i + 3) * this->StateVectorCartesian(j)_GETVALUE - drdotv_dstate(i + 3) * (this->StateVectorCartesian(j + 3) - (i == j ? rdotv : 0))_GETVALUE);
                    }
                }

                //derivatives of ECC scalar
                math::Matrix<double> dECC_dstate(6, 1, 0.0);
                for (size_t i = 0; i < 6; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                        dECC_dstate(i) += (eccVec(j) / ECC)_GETVALUE * deccVec_dstate(i, j);
                }


                //derivatives of RP
                double dRP_dSMA = 1.0 - ECC _GETVALUE;
                double dRP_dECC = -SMA _GETVALUE;
                for (size_t i = 0; i < 6; ++i)
                {
                    this->CartesianToRepresentationTransitionMatrix(i, 0) = dRP_dSMA * dSMAdstate(i) + dRP_dECC * dECC_dstate(i);
                }

                //derivatives of C3
                double dC3dv = 2.0 * v _GETVALUE;
                double dC3dr = 2.0 * mu / (r * r)_GETVALUE;
                math::Matrix<double> dC3_dstate(6, 1, 0.0);

                for (size_t i = 0; i < 6; ++i)
                {
                    dC3_dstate(i) = dC3dr * dr_dstate(i) + dC3dv * dv_dstate(i);
                    this->CartesianToRepresentationTransitionMatrix(i, 1) = dC3_dstate(i);
                }

                //next we need derivatives of sVHat and therefore of hVec
                math::Matrix<double> dhVec_dstate(6, 3, 0.0);
                dhVec_dstate(0, 1) = -this->StateVectorCartesian(5)_GETVALUE;
                dhVec_dstate(0, 2) =  this->StateVectorCartesian(4)_GETVALUE;
                dhVec_dstate(1, 0) =  this->StateVectorCartesian(5)_GETVALUE;
                dhVec_dstate(1, 2) = -this->StateVectorCartesian(3)_GETVALUE;
                dhVec_dstate(2, 0) = -this->StateVectorCartesian(4)_GETVALUE;
                dhVec_dstate(2, 1) =  this->StateVectorCartesian(3)_GETVALUE;
                dhVec_dstate(3, 1) =  this->StateVectorCartesian(2)_GETVALUE;
                dhVec_dstate(3, 2) = -this->StateVectorCartesian(1)_GETVALUE;
                dhVec_dstate(4, 0) = -this->StateVectorCartesian(2)_GETVALUE;
                dhVec_dstate(4, 2) =  this->StateVectorCartesian(0)_GETVALUE;
                dhVec_dstate(5, 0) =  this->StateVectorCartesian(1)_GETVALUE;
                dhVec_dstate(5, 1) = -this->StateVectorCartesian(0)_GETVALUE;

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
                math::Matrix<double> dhMag_dstate(6, 1);
                for (size_t i = 0; i < 6; ++i)
                {
                    dhMag_dstate(i) = 0.0;
                    for (size_t j = 0; j < 3; ++j)
                        dhMag_dstate(i) += (hVec(j) / hMag)_GETVALUE * dhVec_dstate(i, j);
                }

                //derivatives of sVhat
                math::Matrix<double> dsVHat_dstate(6, 3, 0.0);
                if (C3 > 1E-7)
                {
                    //derivatives of fac1
                    math::Matrix<doubleType> dfac1_dstate(6, 1, 0.0);
                    double dfac1_dC3 = (-hMag * hMag / (mu * mu * ((C3 * hMag * hMag) / (mu * mu) + 1) * ((C3 * hMag * hMag) / (mu * mu) + 1)))_GETVALUE;
                    double dfac1_dhMag = (-(2 * C3*hMag) / (mu * mu * ((C3 * hMag *hMag) / (mu * mu) + 1) * ((C3 * hMag *hMag) / (mu * mu) + 1)))_GETVALUE;
                    for (size_t i = 0; i < 6; ++i)
                        dfac1_dstate(i) = dfac1_dC3 * dC3_dstate(i);

                    //first need derivatives of hVec.cross(eccVec)
                    math::Matrix<double> dhcrosse_dhVec(3, 3, 0.0);
                    math::Matrix<double> dhcrosse_deccVec(3, 3, 0.0);
                    math::Matrix<double> dhcrosse_dstate(6, 3, 0.0);

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
                    math::Matrix<double> dsVHat_dstate(6, 3, 0.0);
                    for (size_t i = 0; i < 6; ++i)
                    {
                        for (size_t j = 0; j < 3; ++j)
                        {
                            dsVHat_dstate(i, j) = (fac1 * (-sqrt(C3) / mu * dhcrosse_dstate(i, j)
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
            }

            //Step 4: return the converted vector
            return this->StateVectorThisRepresentation;
        }//end convertFromCartesianToRepresentation()
    }//end namespace StateRepresentation
}//end namespace EMTG