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

#include "MEEStateRepresentation.h"
#include "EMTG_math.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        MEEStateRepresentation::MEEStateRepresentation(const double& mu) : StateRepresentationBase::StateRepresentationBase(mu)
        {
            this->name = "MEE";
            this->stateNames = std::vector<std::string>({ "P", "F", "G", "H", "K", "L" });
            this->myCOEStateRepresentation.setmu(mu);
        };

        math::Matrix<doubleType> MEEStateRepresentation::convertFromRepresentationToCartesian(const math::Matrix<doubleType>& StateVectorThisRepresentationIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorThisRepresentation(StateVectorThisRepresentationIn);

            //Step 2: convert the state
            //Step 2.1: convert from MEE to COE
            const doubleType& P = this->StateVectorThisRepresentation(0);
            const doubleType& F = this->StateVectorThisRepresentation(1);
            const doubleType& G = this->StateVectorThisRepresentation(2);
            const doubleType& H = this->StateVectorThisRepresentation(3);
            const doubleType& K = this->StateVectorThisRepresentation(4);
            const doubleType& L = this->StateVectorThisRepresentation(5);

            doubleType SMA = P / (1.0 - F * F - G * G);
            doubleType ECC = sqrt(F * F + G * G);

            doubleType A = sqrt(H * H + K * K);
            doubleType B = atan2(G / ECC, F / ECC);

            doubleType RAAN = atan2(K / A, H / A);

            bool prograde = true; // TODO: how do you check an MEE state for whether it is prograde or retrograde?

            doubleType INC, AOP;
            if (prograde)
            {
                INC = 2.0 * atan(A);
                AOP = B - RAAN;
            }
            else
            {
                INC = math::PI - 2.0 * atan(A);
                AOP = B + RAAN;
            }

            doubleType TA = L - B;

            math::Matrix<doubleType> StateVectorCOE(6, 1, { SMA, ECC, INC, RAAN, AOP, TA });

            //Step 2.2: convert from COE to Cartesian
            this->StateVectorCartesian = this->myCOEStateRepresentation.convertFromRepresentationToCartesian(StateVectorCOE, needG);


            //Step 3: construct the partial derivatives
            if (needG)
            {
                this->RepresentationToCartesianTransitionMatrix.assign_zeros();

                //Step 3.1: compute the derivatives from MEE to COE
                doubleType dSMA_dP = 1.0 / (1.0 - F * F - G * G);
                doubleType dSMA_dF = 2.0 * F * P / (1.0 - F * F - G * G) / (1.0 - F * F - G * G);
                doubleType dSMA_dG = 2.0 * G * P / (1.0 - F * F - G * G) / (1.0 - F * F - G * G);

                doubleType dECC_dF = F / ECC;
                doubleType dECC_dG = G / ECC;

                doubleType dRAAN_dH = -K / A / A;
                doubleType dRAAN_dK = H / A / A;

                doubleType dB_dF = -G / ECC / ECC;
                doubleType dB_dG = F / ECC / ECC;

                doubleType dAOP_dF = dB_dF;
                doubleType dAOP_dG = dB_dG;

                doubleType dINC_dH, dINC_dK, dAOP_dH, dAOP_dK;
                if (prograde)
                {
                    dINC_dH = 2.0 * H / A / (1.0 + A * A);
                    dINC_dK = 2.0 * K / A / (1.0 + A * A);
                    dAOP_dH = -dRAAN_dH;
                    dAOP_dK = -dRAAN_dK;

                }
                else
                {
                    dINC_dH = -2.0 * H / A / (1.0 + A * A);
                    dINC_dK = -2.0 * K / A / (1.0 + A * A);
                    dAOP_dH = dRAAN_dH;
                    dAOP_dK = dRAAN_dK;
                }

                doubleType dTA_dF = -dB_dF;
                doubleType dTA_dG = -dB_dG;
                doubleType dTA_dL = 1.0;

                math::Matrix<doubleType> TransformationFromMEEtoCOE(6, 6, 0.0);

                TransformationFromMEEtoCOE(0, 0) = dSMA_dP;
                TransformationFromMEEtoCOE(0, 1) = dSMA_dF;
                TransformationFromMEEtoCOE(0, 2) = dSMA_dG;

                TransformationFromMEEtoCOE(1, 1) = dECC_dF;
                TransformationFromMEEtoCOE(1, 2) = dECC_dG;

                TransformationFromMEEtoCOE(2, 3) = dINC_dH;
                TransformationFromMEEtoCOE(2, 4) = dINC_dK;
                
                TransformationFromMEEtoCOE(3, 3) = dRAAN_dH;
                TransformationFromMEEtoCOE(3, 4) = dRAAN_dK;

                TransformationFromMEEtoCOE(4, 1) = dAOP_dF;
                TransformationFromMEEtoCOE(4, 2) = dAOP_dG;
                TransformationFromMEEtoCOE(4, 3) = dAOP_dH;
                TransformationFromMEEtoCOE(4, 4) = dAOP_dK;

                TransformationFromMEEtoCOE(5, 1) = dTA_dF;
                TransformationFromMEEtoCOE(5, 2) = dTA_dG;
                TransformationFromMEEtoCOE(5, 5) = 1.0;


                //Step 3.2: assemble the full transformation matrix
                this->RepresentationToCartesianTransitionMatrix = this->myCOEStateRepresentation.getRepresentationToCartesianTransitionMatrix() * TransformationFromMEEtoCOE;
            }

            //Step 4: return the converted vector
            return this->StateVectorCartesian;
        }//end convertFromRepresentationToCartesian()

        math::Matrix<doubleType> MEEStateRepresentation::convertFromCartesianToRepresentation(const math::Matrix<doubleType>& StateVectorCartesianIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorCartesian(StateVectorCartesianIn);

            //Step 2: convert the state           
            //Step 2.1: convert from Cartesian to COE
            math::Matrix<doubleType> stateCOE = this->myCOEStateRepresentation.convertFromCartesianToRepresentation(this->StateVectorCartesian, needG);

            //Step 2.2: convert from COE to MEE
            const doubleType& SMA = stateCOE(0);
            const doubleType& ECC = stateCOE(1);
            const doubleType& INC = stateCOE(2);
            const doubleType& RAAN = stateCOE(3);
            const doubleType& AOP = stateCOE(4);
            const doubleType& TA = stateCOE(5);

            doubleType tanINCover2 = tan(INC / 2.0);
            doubleType cosINCover2 = cos(INC / 2.0);
            doubleType cosRAANplusAOP = cos(RAAN + AOP);
            doubleType sinRAANplusAOP = sin(RAAN + AOP);
            doubleType cosRAAN = cos(RAAN);
            doubleType sinRAAN = sin(RAAN);
            
            doubleType P = SMA * (1 - ECC * ECC);
            doubleType F = ECC * cosRAANplusAOP;
            doubleType G = ECC * sinRAANplusAOP;
            doubleType H = tanINCover2 * cosRAAN;
            doubleType K = tanINCover2 * sinRAAN;
            doubleType L = RAAN + AOP + TA;

            this->StateVectorThisRepresentation(0) = P;
            this->StateVectorThisRepresentation(1) = F;
            this->StateVectorThisRepresentation(2) = G;
            this->StateVectorThisRepresentation(3) = H;
            this->StateVectorThisRepresentation(4) = K;
            this->StateVectorThisRepresentation(5) = L;


            //Step 3: construct the partial derivatives
            if (needG)
            {
                this->CartesianToRepresentationTransitionMatrix.assign_zeros();

                //Step 3.1: compute the derivatives from COE to MEE
                doubleType dP_dSMA = 1.0 - ECC * ECC;
                doubleType dP_dECC = -2.0 * SMA * ECC;
                
                doubleType dF_dECC = cosRAANplusAOP;
                doubleType dF_dRAAN = -ECC * sinRAANplusAOP;
                doubleType dF_dAOP = -ECC * sinRAANplusAOP;

                doubleType dG_dECC = sinRAANplusAOP;
                doubleType dG_dRAAN = ECC * cosRAANplusAOP;
                doubleType dG_dAOP = ECC * cosRAANplusAOP;

                doubleType dH_dINC = cosRAAN / 2.0 / cosINCover2 / cosINCover2;
                doubleType dH_dRAAN = -tanINCover2 * sinRAAN;

                doubleType dK_dINC = sinRAAN / 2.0 / cosINCover2 / cosINCover2;
                doubleType dK_dRAAN = tanINCover2 * cosRAAN;

                doubleType dL_dRAAN = 1.0;
                doubleType dL_dAOP = 1.0;
                doubleType dL_dTA = 1.0;

                math::Matrix<doubleType> TransformationFromCOEtoMEE(6, 6, 0.0);
                
                TransformationFromCOEtoMEE(0, 0) = dP_dSMA;
                TransformationFromCOEtoMEE(0, 1) = dP_dECC;

                TransformationFromCOEtoMEE(1, 1) = dF_dECC;
                TransformationFromCOEtoMEE(1, 3) = dF_dRAAN;
                TransformationFromCOEtoMEE(1, 4) = dF_dAOP;

                TransformationFromCOEtoMEE(2, 1) = dG_dECC;
                TransformationFromCOEtoMEE(2, 3) = dG_dRAAN;
                TransformationFromCOEtoMEE(2, 4) = dG_dAOP;

                TransformationFromCOEtoMEE(3, 2) = dH_dINC;
                TransformationFromCOEtoMEE(3, 3) = dH_dRAAN;

                TransformationFromCOEtoMEE(4, 2) = dK_dINC;
                TransformationFromCOEtoMEE(4, 3) = dK_dRAAN;

                TransformationFromCOEtoMEE(5, 3) = dL_dRAAN;
                TransformationFromCOEtoMEE(5, 4) = dL_dAOP;
                TransformationFromCOEtoMEE(5, 5) = dL_dTA;

                //Step 3.2: assemble the full transformation matrix
                this->CartesianToRepresentationTransitionMatrix = TransformationFromCOEtoMEE * this->myCOEStateRepresentation.getCartesianToRepresentationTransitionMatrix();
            }//end derivatives

            //Step 4: return the converted vector
            return this->StateVectorThisRepresentation;
        }//end convertFromCartesianToRepresentation()
    }//end namespace StateRepresentation
}//end namespace EMTG