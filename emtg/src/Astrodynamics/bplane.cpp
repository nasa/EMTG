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

//source file for b-plane class
//for use with EMTGv9
//Jacob Englander 1-10-2013

#include "EMTG_Matrix.h"
#include "EMTG_math.h"
#include "bplane.h"

namespace EMTG 
{ 
    namespace Astrodynamics 
    {
        //default constructor
        bplane::bplane() {}

        //constructor with mu
        bplane::bplane(const double& mu) : bplane()
        {
            this->mu = mu;
        }

        //destructor doesn't currently do anything
        bplane::~bplane() {}


        //function to define the b-plane coordinate system. B-plane is defined by the vectors T_hat and R_hat
        //S_hat is the incoming velocity asymptote unit vector
        //from http://ccar.colorado.edu/asen5519/imd/documents/BPlaneHandout.pdf
        void bplane::define_bplane(math::Matrix<doubleType>& V_infinity_in,
            math::Matrix<doubleType>& bodyState)
        {
            //create Rbody and Vbody
            math::Matrix<doubleType> R_body = bodyState.getSubMatrix1D(0, 2);
            math::Matrix<doubleType> V_body = bodyState.getSubMatrix1D(3, 5);

            //compute S_hat
            this->V_infinity_mag = V_infinity_in.norm();
            this->S_hat = V_infinity_in.unitize();

            //compute the body orbit normal vector - this is the STK/GMAT "reference vector" from which the b-plane is constructed
            //note: reference frame is ICRF
            this->k_hat = R_body.unitcross(V_body);

            //T_hat and R_hat
            this->T_hat = S_hat.unitcross(k_hat);
            this->R_hat = S_hat.unitcross(T_hat);
        }//end define_bplane() with v-infinity in and body state
        
        void bplane::define_bplane(math::Matrix<doubleType> periapseState)
        {
            math::Matrix<doubleType> R = periapseState.getSubMatrix1D(0, 2);
            math::Matrix<doubleType> V = periapseState.getSubMatrix1D(3, 5);

            //Step 1: compute eccentricity vector
            doubleType rdotv = R.dot(V);
            doubleType r = R.norm();
            doubleType v = V.norm();
            doubleType s = (v*v - this->mu / r);

            math::Matrix<doubleType> e_vec = (R * s - V * rdotv) / mu;

            doubleType e = e_vec.norm();
            math::Matrix<doubleType> e_hat = e_vec.unitize();
            
            //Step 2: compute angular momentum
            math::Matrix<doubleType> h_vec = R.cross(V);
            this->h = h_vec.norm();
            this->h_hat = h_vec.unitize();

            //Step 3: compute S_hat, the vector in the direction of v-infinity
            this->S_hat = e_hat / e + h_hat.unitcross(e_hat) * sqrt(1.0 - 1.0 / (e * e));

            //Step 4: compute Bradius
            this->Bradius = this->h * this->h / (this->mu * sqrt(e * e - 1.0));

            //Step 5: compute B vector
            this->B = S_hat.unitcross(h_hat) * this->Bradius;

            //Step 6: compute T_hat and R_hat
            this->T_hat = math::Matrix<doubleType>(3, 1, std::vector<doubleType>({ this->S_hat(1), -this->S_hat(0), 0.0 })).unitize();
            this->R_hat = S_hat.unitcross(T_hat);

            //Step 7: compute BdotT and BdotR
            this->BdotT = this->B.dot(this->T_hat);
            this->BdotR = this->B.dot(this->R_hat);

            //Step 8: compute Btheta
            this->Btheta = atan2(BdotR, BdotT);
        }//end define_bplane() with periapse state

        void bplane::define_bplane_without_bend_angle(math::Matrix<doubleType> periapseState, 
            math::Matrix<doubleType>& V_infinity_in)
        {
            math::Matrix<doubleType> R = periapseState.getSubMatrix1D(0, 2);
            math::Matrix<doubleType> V = periapseState.getSubMatrix1D(3, 5);
            
            //Step 1: compute angular momentum
            math::Matrix<doubleType> h_vec = R.cross(V);
            this->h = h_vec.norm();
            this->h_hat = h_vec.unitize();

            //Step 2: compute S_hat, the vector in the direction of v-infinity
            this->S_hat = V_infinity_in.unitize();

            //Step 3: compute Bradius
            this->Bradius = S_hat.cross(R.cross(S_hat)).norm();

            //Step 4: compute B vector
            this->B = S_hat.unitcross(h_hat) * this->Bradius;

            //Step 5: compute T_hat and R_hat
            this->T_hat = math::Matrix<doubleType>(3, 1, std::vector<doubleType>({ this->S_hat(1), -this->S_hat(0), 0.0 }));
            this->R_hat = S_hat.unitcross(T_hat);

            //Step 6: compute BdotT and BdotR
            this->BdotT = this->B.dot(this->T_hat);
            this->BdotR = this->B.dot(this->R_hat);

            //Step 8: compute Btheta
            this->Btheta = atan2(BdotR, BdotT);
        }//end define_bplane() with periapse state

        //function to convert from Bradius, Btheta to BdotR, BdotT
        void bplane::convert_polar_to_cartesian(const doubleType& Bradius, const doubleType& Btheta, doubleType& BdotR, doubleType& BdotT)
        {
            BdotT = Bradius * cos(Btheta);
            BdotR = Bradius * sin(Btheta);
        }

        //function to convert from BdotR, BdotT to Bradius, Btheta
        void bplane::convert_cartesian_to_polar(const doubleType& BdotR, const doubleType& BdotT, doubleType& Bradius, doubleType& Btheta)
        {
            Bradius = sqrt(BdotR*BdotR + BdotT*BdotT);
            Btheta = atan2(BdotR, BdotT);
        }

        //function to compute the b-plane Bradius, Btheta, BdotR, BdotT, and periapse distance from an outbound V-infinity vector, for flybys
        //from http://ccar.colorado.edu/asen5519/imd/documents/BPlaneHandout.pdf
        void bplane::compute_bplane_coordinates_from_Vinfinity_out(const double& mu, math::Matrix<doubleType>& V_infinity_out, doubleType& Bradius, doubleType& Btheta, doubleType& BdotR, doubleType& BdotT, doubleType& rp)
        {
            //compute the square of the V_infinity magnitude
            doubleType C3 = V_infinity_mag*V_infinity_mag;
        
            //compute h_hat, the unit vector in the direction of spacecraft orbit normal
            //note that both the incoming and outgoing V_infinity will be in the plane of the orbit, so the angular momentum unit vector is the normalized 
            //cross product of the two V_infinity vectors
            h_hat = S_hat.unitcross(V_infinity_out);
        
            //compute unit vector pointing toward b-plane crossing
            B_hat = S_hat.unitcross(h_hat);

            //compute phi
            doubleType phi = math::safe_acos(S_hat.dot(V_infinity_out) / V_infinity_mag);

            //compute rp
            rp = mu / C3 * ( 1 / cos( (math::PI - phi) / 2) - 1);

            //compute Bmag and B
            doubleType A = 1 + C3 * rp / mu;

            Bradius = mu / C3 * sqrt(A*A - 1);
            B = B_hat * Bradius;

            //compute Btheta
            Btheta = B_hat.dot(R_hat) >= 0 ?
                math::safe_acos(T_hat.dot(B_hat))
                :
                2 * math::PI - math::safe_acos(T_hat.dot(B_hat));

            //compute BdotR, BdotT
            BdotR = B.dot(R_hat);
            BdotT = B.dot(T_hat);
        }

        //function to compute the outbound V-infinity vector and periapse distance from Bradius and Btheta
        void bplane::compute_Vinfinity_out_from_Bradius_Btheta(const double& mu, const doubleType& Bradius, const doubleType& Btheta, math::Matrix<doubleType>& V_infinity_out, doubleType& rp)
        {
            //first convert from (Bradius, Btheta) to (BdotR, BdotT)
            doubleType BdotR, BdotT;
            this->convert_polar_to_cartesian(Bradius, Btheta, BdotR, BdotT);

            this->compute_Vinfinity_out_from_BdotR_BdotT(mu, BdotR, BdotT, V_infinity_out, rp);
        }

        //function to compute the outbound V-infinity vector and periapse distance from BdotR and BdotT
        void bplane::compute_Vinfinity_out_from_BdotR_BdotT(const double& mu, const doubleType& BdotR, const doubleType& BdotT, math::Matrix<doubleType>& V_infinity_out, doubleType& rp)
        {
            throw std::runtime_error("bplane::compute_Vinfinity_out_from_BdotR_BdotT is not yet implemented!");

            //basically, compute B_hat and negate the component of V_infinity_in in the direction of Bhat
            //i.e. (V_infinity_in dot B_hat) = (V_infinity_out dot B_hat)
        }

        //function to compute Bradius and Btheta from periapse position vector
        void bplane::compute_Bradius_Btheta_from_periapse_position(const double& mu, math::Matrix<doubleType>& Vinfinity_in, math::Matrix<doubleType>& Rp, doubleType& Bradius, doubleType& Btheta)
        {
            //first compute periapse distance
            doubleType rp = Rp.norm();

            //compute C3
            doubleType C3 = V_infinity_mag*V_infinity_mag;
        
            //compute velocity magnitude at periapse
            doubleType vp = sqrt(C3 + 2*mu/rp);

            //compute Bradius
            Bradius = vp * rp / V_infinity_mag;

            //compute h_hat
            h_hat = Rp.unitcross(Vinfinity_in);

            //compute B
            B_hat = S_hat.unitcross(h_hat);
            B = B_hat * Bradius;

            //compute Btheta
            Btheta = B_hat.dot(R_hat) >= 0 ? 
                math::safe_acos(T_hat.dot(B_hat))
                :
                2 * math::PI - math::safe_acos(T_hat.dot(B_hat));
        }

        //function to compute BdotR and BdotT from periapse position vector
        void bplane::compute_BdotR_BdotT_from_periapse_position(const double& mu, math::Matrix<doubleType>& Vinfinity_in, math::Matrix<doubleType>& Rp, doubleType& BdotR, doubleType& BdotT)
        {
            doubleType Bradius, Btheta;
            this->compute_Bradius_Btheta_from_periapse_position(mu, Vinfinity_in, Rp, Bradius, Btheta);

            this->convert_polar_to_cartesian(Bradius, Btheta, BdotR, BdotT);
        }
    }//close namespace Astrodynamics
}//close namespace EMTG