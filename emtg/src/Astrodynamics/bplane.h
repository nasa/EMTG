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

//header file for b-plane class
//for use with EMTGv9
//Jacob Englander 1-10-2013


#pragma once

#include "doubleType.h"

#include "EMTG_Matrix.h"

namespace EMTG 
{ 
    namespace Astrodynamics
    {
        class bplane
        {
        public:
            //constructor
            bplane();
            bplane(const double& mu);

            //destructor
            virtual ~bplane();

            //methods
            inline doubleType getBdotR() const { return this->BdotR; }
            inline doubleType getBdotT() const { return this->BdotT; }
            inline doubleType getBradius() const { return this->Bradius; }
            inline doubleType getBtheta() const { return this->Btheta; }

            //functions to define the b-plane coordinate system
            void define_bplane(math::Matrix<doubleType>& V_infinity_in, math::Matrix<doubleType>& bodyState);
            void define_bplane(math::Matrix<doubleType> periapseState);
            void define_bplane_without_bend_angle(math::Matrix<doubleType> periapseState, math::Matrix<doubleType>& V_infinity_in);

            //function to compute the b-plane Bradius, Btheta, BdotR, BdotT, and periapse distance from an outbound V-infinity vector, for flybys
            void compute_bplane_coordinates_from_Vinfinity_out(const double& mu, math::Matrix<doubleType>& V_infinity_out, doubleType& Bradius, doubleType& Btheta, doubleType& BdotR, doubleType& BdotT, doubleType& rp);

            //function to compute the outbound V-infinity vector and periapse distance from Bradius and Btheta
            void compute_Vinfinity_out_from_Bradius_Btheta(const double& mu, const doubleType& Bradius, const doubleType& Btheta, math::Matrix<doubleType>& V_infinity_out, doubleType& rp);

            //function to compute the outbound V-infinity vector and periapse distance from BdotR and BdotT
            void compute_Vinfinity_out_from_BdotR_BdotT(const double& mu, const doubleType& BdotR, const doubleType& BdotT, math::Matrix<doubleType>& V_infinity_out, doubleType& rp);

            //function to compute Bradius and Btheta from periapse position vector
            void compute_Bradius_Btheta_from_periapse_position(const double& mu, math::Matrix<doubleType>& Vinfinity_in, math::Matrix<doubleType>& Rp, doubleType& Bradius, doubleType& Btheta);

            //function to compute BdotR and BdotT from periapse position vector
            void compute_BdotR_BdotT_from_periapse_position(const double& mu, math::Matrix<doubleType>& Vinfinity_in, math::Matrix<doubleType>& Rp, doubleType& BdotR, doubleType& BdotT);
            
            //function to convert from Bradius, Btheta to BdotR, BdotT
            void convert_polar_to_cartesian(const doubleType& Bradius, const doubleType& Btheta, doubleType& BdotR, doubleType& BdotT);

            //function to convert from BdotR, BdotT to Bradius, Btheta
            void convert_cartesian_to_polar(const doubleType& BdotR, const doubleType& BdotT, doubleType& Bradius, doubleType& Btheta);


        private:
            //fields
            math::Matrix<doubleType> S_hat, h_hat, k_hat, B_hat, T_hat, R_hat, B;
            doubleType Bmag, h;
            doubleType V_infinity_mag;
            doubleType BdotR, BdotT, Bradius, Btheta;
            double mu;
            double rbody;

        };
    }//close namespace Astrodynamics
}//close namespace EMTG