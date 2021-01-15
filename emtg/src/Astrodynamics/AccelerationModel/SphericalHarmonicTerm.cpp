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

// J2 oblateness term

#include "body.h"
#include "CentralBodyGravityTerm.h"
#include "SphericalHarmonicTerm.h"
#include "doubleType.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        //constructors
        SphericalHarmonicTerm::SphericalHarmonicTerm(SpacecraftAccelerationModel* acceleration_model_in, body* my_body_in, CentralBodyGravityTerm* parent_gravity_term_in, const size_t& degree_in, const size_t& order_in) :
            SpacecraftAccelerationModelTerm::SpacecraftAccelerationModelTerm(acceleration_model_in),
            my_body(my_body_in),
            parent_gravity_term(parent_gravity_term_in),
            degree(degree_in),
            order(order_in)
        {
            this->r_body2sc_BCF.resize(3, 1, 0.0);
            this->spherical_harmonic_term_acceleration_BCF.resize(3, 1, 0.0);
        }

        SphericalHarmonicTerm::~SphericalHarmonicTerm() {}

        //methods
        void SphericalHarmonicTerm::computeAccelerationTerm()
        {
            // J2 calculations must be performed in the central-body fixed equatorial frame
            static EMTG::math::Matrix <doubleType> temp_matrix_in(3, 1, 0.0);
            static EMTG::math::Matrix <doubleType> temp_matrix_out(3, 1, 0.0);

            // rotate position ICRF->BCF
            for (size_t k = 0; k < 3; ++k)
            {
                temp_matrix_in(k, 0) = this->parent_gravity_term->r_body2sc(k);
            }
            this->acceleration_model->my_universe->LocalFrame.rotate_frame_to_frame(EMTG::ReferenceFrame::ICRF,
                temp_matrix_in,
                EMTG::ReferenceFrame::TrueOfDate_BCF,
                temp_matrix_out,
                this->acceleration_model->current_epoch);
            for (size_t k = 0; k < 3; ++k)
            {
                this->r_body2sc_BCF(k) = temp_matrix_out(k, 0);
            }

            this->r_body2sc_BCF_norm = sqrt(this->r_body2sc_BCF(0) * this->r_body2sc_BCF(0) +
                this->r_body2sc_BCF(1) * this->r_body2sc_BCF(1) +
                this->r_body2sc_BCF(2) * this->r_body2sc_BCF(2));

            doubleType spacecraft_distance_from_central_body_BCF2 = this->r_body2sc_BCF_norm * this->r_body2sc_BCF_norm;
            doubleType spacecraft_distance_from_central_body_BCF5 = this->r_body2sc_BCF_norm * spacecraft_distance_from_central_body_BCF2 * spacecraft_distance_from_central_body_BCF2;

            doubleType z2 = this->r_body2sc_BCF(2) * this->r_body2sc_BCF(2);
            doubleType zcoeff = 5.0 * z2 / (spacecraft_distance_from_central_body_BCF2);
            doubleType J2coeff = 3.0 * this->acceleration_model->my_universe->central_body_J2 * this->my_body->mu * this->my_body->J2_ref_radius * this->my_body->J2_ref_radius / (2.0 * spacecraft_distance_from_central_body_BCF5);

            // compute the J2 acceleration in BCF coordinates
            this->spherical_harmonic_term_acceleration_BCF(0) = -J2coeff * (1.0 - zcoeff) * this->r_body2sc_BCF(0);
            this->spherical_harmonic_term_acceleration_BCF(1) = -J2coeff * (1.0 - zcoeff) * this->r_body2sc_BCF(1);
            this->spherical_harmonic_term_acceleration_BCF(2) = -J2coeff * (3.0 - zcoeff) * this->r_body2sc_BCF(2);

            // rotate the J2 acceleration vector BCF->ICRF
            for (size_t k = 0; k < 3; ++k)
            {
                temp_matrix_in(k, 0) = this->spherical_harmonic_term_acceleration_BCF(k);
            }
            this->acceleration_model->my_universe->LocalFrame.rotate_frame_to_frame(EMTG::ReferenceFrame::TrueOfDate_BCF,
                temp_matrix_in,
                EMTG::ReferenceFrame::ICRF,
                temp_matrix_out,
                this->acceleration_model->current_epoch);

            // dump the spherical harmonic contribution into the parent gravity term
            for (size_t k = 0; k < 3; ++k)
            {
                this->term_acceleration(k) = temp_matrix_out(k, 0);
                // for now, we are adding the spherical harmonic influence into the parent GravityTerm's term_acceleration
                this->parent_gravity_term->term_acceleration(k) += temp_matrix_out(k, 0);
                // add the J2 acceleration into the overall acceleration vector
                this->acceleration_model->acceleration(k) += temp_matrix_out(k, 0);
            }

        }// end computeAccelerationTerm()

        void SphericalHarmonicTerm::computeAccelerationTerm(const bool& generate_derivatives)
        {
            this->computeAccelerationTerm();

            doubleType xy = this->r_body2sc_BCF(0) * this->r_body2sc_BCF(1);
            doubleType xz = this->r_body2sc_BCF(0) * this->r_body2sc_BCF(2);
            doubleType yz = this->r_body2sc_BCF(1) * this->r_body2sc_BCF(2);
            doubleType x2 = this->r_body2sc_BCF(0) * this->r_body2sc_BCF(0);
            doubleType y2 = this->r_body2sc_BCF(1) * this->r_body2sc_BCF(1);
            doubleType x4 = x2 * x2;
            doubleType y4 = y2 * y2;
            doubleType z2 = this->r_body2sc_BCF(2) * this->r_body2sc_BCF(2);
            doubleType z4 = z2 * z2;
            doubleType x2y2 = x2 * y2;
            doubleType x2z2 = x2 * z2;
            doubleType y2z2 = y2 * z2;
            doubleType spacecraft_distance_from_central_body_BCF2 = this->r_body2sc_BCF_norm * this->r_body2sc_BCF_norm;
            doubleType spacecraft_distance_from_central_body_BCF5 = this->r_body2sc_BCF_norm * spacecraft_distance_from_central_body_BCF2 * spacecraft_distance_from_central_body_BCF2;
            doubleType spacecraft_distance_from_central_body_BCF9 = spacecraft_distance_from_central_body_BCF5 * spacecraft_distance_from_central_body_BCF2 * spacecraft_distance_from_central_body_BCF2;
            doubleType J2_deriv_coeff = -this->acceleration_model->my_universe->central_body_J2 * this->my_body->mu * this->my_body->J2_ref_radius * this->my_body->J2_ref_radius / (2.0 * spacecraft_distance_from_central_body_BCF9);
            doubleType J2_deriv_coeff3 = J2_deriv_coeff * J2_deriv_coeff * J2_deriv_coeff;


            // grab the ICRF->BCF and BCF->ICRF frame rotation matrices
            static EMTG::math::Matrix <doubleType> R_from_ICRF_to_BCF;
            static EMTG::math::Matrix <doubleType> R_from_BCF_to_ICRF;
            static EMTG::math::Matrix <doubleType> accel_position_Jacobian_BCF(3, 3, 0.0);
            static EMTG::math::Matrix <doubleType> accel_position_Jacobian_ICRF(3, 3, 0.0);
            // we need to compute the matrices because we need their explicit time partials, and to get those we need to pass the generate_derivatives flag
            this->acceleration_model->my_universe->LocalFrame.construct_rotation_matrices(this->acceleration_model->current_epoch, generate_derivatives);

            R_from_ICRF_to_BCF = this->acceleration_model->my_universe->LocalFrame.get_R_from_ICRF_to_BCF();
            R_from_BCF_to_ICRF = this->acceleration_model->my_universe->LocalFrame.get_R_from_BCF_to_ICRF();

            accel_position_Jacobian_BCF(0, 0) = -3.0 * J2_deriv_coeff * (4.0 * x4 + 3.0 * x2y2 - 27.0 * x2z2 - y4 + 3.0 * y2z2 + 4.0 * z4);
            accel_position_Jacobian_BCF(0, 1) = -15.0 * J2_deriv_coeff * xy * (x2 + y2 - 6.0 * z2);
            accel_position_Jacobian_BCF(0, 2) = -15.0 * J2_deriv_coeff * xz * (3.0 * x2 + 3.0 * y2 - 4.0 * z2);
            accel_position_Jacobian_BCF(1, 0) = -15.0 * J2_deriv_coeff * xy * (x2 + y2 - 6.0 * z2);
            accel_position_Jacobian_BCF(1, 1) = -3.0 * J2_deriv_coeff * (-1.0 * x4 + 3.0 * x2y2 + 3.0 * x2z2 + 4.0 * y4 - 27.0 * y2z2 + 4.0 * z4);
            accel_position_Jacobian_BCF(1, 2) = -15.0 * J2_deriv_coeff * yz * (3.0 * x2 + 3.0 * y2 - 4.0 * z2);
            accel_position_Jacobian_BCF(2, 0) = -15.0 * J2_deriv_coeff * xz * (3.0 * x2 + 3.0 * y2 - 4.0 * z2);
            accel_position_Jacobian_BCF(2, 1) = -15.0 * J2_deriv_coeff * yz * (3.0 * x2 + 3.0 * y2 - 4.0 * z2);
            accel_position_Jacobian_BCF(2, 2) = 3.0 * J2_deriv_coeff * (3.0 * x4 + 6.0 * x2y2 - 24.0 * x2z2 + 3.0 * y4 - 24.0 * y2z2 + 8.0 * z4);

            accel_position_Jacobian_ICRF = (R_from_BCF_to_ICRF * accel_position_Jacobian_BCF * R_from_ICRF_to_BCF);

            //doubleType dax_BCFdx_ICRF = dax_BCFdx_BCF * R_from_ICRF_to_BCF(0, 0) + dax_BCFdy_BCF * R_from_ICRF_to_BCF(1, 0) + dax_BCFdz_BCF * R_from_ICRF_to_BCF(2, 0);
            //doubleType dax_BCFdy_ICRF = dax_BCFdx_BCF * R_from_ICRF_to_BCF(0, 1) + dax_BCFdy_BCF * R_from_ICRF_to_BCF(1, 1) + dax_BCFdz_BCF * R_from_ICRF_to_BCF(2, 1);
            //doubleType dax_BCFdz_ICRF = dax_BCFdx_BCF * R_from_ICRF_to_BCF(0, 2) + dax_BCFdy_BCF * R_from_ICRF_to_BCF(1, 2) + dax_BCFdz_BCF * R_from_ICRF_to_BCF(2, 2);
            //doubleType day_BCFdx_ICRF = day_BCFdx_BCF * R_from_ICRF_to_BCF(0, 0) + day_BCFdy_BCF * R_from_ICRF_to_BCF(1, 0) + day_BCFdz_BCF * R_from_ICRF_to_BCF(2, 0);
            //doubleType day_BCFdy_ICRF = day_BCFdx_BCF * R_from_ICRF_to_BCF(0, 1) + day_BCFdy_BCF * R_from_ICRF_to_BCF(1, 1) + day_BCFdz_BCF * R_from_ICRF_to_BCF(2, 1);
            //doubleType day_BCFdz_ICRF = day_BCFdx_BCF * R_from_ICRF_to_BCF(0, 2) + day_BCFdy_BCF * R_from_ICRF_to_BCF(1, 2) + day_BCFdz_BCF * R_from_ICRF_to_BCF(2, 2);
            //doubleType daz_BCFdx_ICRF = daz_BCFdx_BCF * R_from_ICRF_to_BCF(0, 0) + daz_BCFdy_BCF * R_from_ICRF_to_BCF(1, 0) + daz_BCFdz_BCF * R_from_ICRF_to_BCF(2, 0);
            //doubleType daz_BCFdy_ICRF = daz_BCFdx_BCF * R_from_ICRF_to_BCF(0, 1) + daz_BCFdy_BCF * R_from_ICRF_to_BCF(1, 1) + daz_BCFdz_BCF * R_from_ICRF_to_BCF(2, 1);
            //doubleType daz_BCFdz_ICRF = daz_BCFdx_BCF * R_from_ICRF_to_BCF(0, 2) + daz_BCFdy_BCF * R_from_ICRF_to_BCF(1, 2) + daz_BCFdz_BCF * R_from_ICRF_to_BCF(2, 2);
            //
            //doubleType dax_ICRFdx_ICRF = R_from_BCF_to_ICRF(0, 0) * dax_BCFdx_ICRF + R_from_BCF_to_ICRF(0, 1) * day_BCFdx_ICRF + R_from_BCF_to_ICRF(0, 2) * daz_BCFdx_ICRF;
            //doubleType dax_ICRFdy_ICRF = R_from_BCF_to_ICRF(0, 0) * dax_BCFdy_ICRF + R_from_BCF_to_ICRF(0, 1) * day_BCFdy_ICRF + R_from_BCF_to_ICRF(0, 2) * daz_BCFdy_ICRF;
            //doubleType dax_ICRFdz_ICRF = R_from_BCF_to_ICRF(0, 0) * dax_BCFdz_ICRF + R_from_BCF_to_ICRF(0, 1) * day_BCFdz_ICRF + R_from_BCF_to_ICRF(0, 2) * daz_BCFdz_ICRF;
            //doubleType day_ICRFdx_ICRF = R_from_BCF_to_ICRF(1, 0) * dax_BCFdx_ICRF + R_from_BCF_to_ICRF(1, 1) * day_BCFdx_ICRF + R_from_BCF_to_ICRF(1, 2) * daz_BCFdx_ICRF;
            //doubleType day_ICRFdy_ICRF = R_from_BCF_to_ICRF(1, 0) * dax_BCFdy_ICRF + R_from_BCF_to_ICRF(1, 1) * day_BCFdy_ICRF + R_from_BCF_to_ICRF(1, 2) * daz_BCFdy_ICRF;
            //doubleType day_ICRFdz_ICRF = R_from_BCF_to_ICRF(1, 0) * dax_BCFdz_ICRF + R_from_BCF_to_ICRF(1, 1) * day_BCFdz_ICRF + R_from_BCF_to_ICRF(1, 2) * daz_BCFdz_ICRF;
            //doubleType daz_ICRFdx_ICRF = R_from_BCF_to_ICRF(2, 0) * dax_BCFdx_ICRF + R_from_BCF_to_ICRF(2, 1) * day_BCFdx_ICRF + R_from_BCF_to_ICRF(2, 2) * daz_BCFdx_ICRF;
            //doubleType daz_ICRFdy_ICRF = R_from_BCF_to_ICRF(2, 0) * dax_BCFdy_ICRF + R_from_BCF_to_ICRF(2, 1) * day_BCFdy_ICRF + R_from_BCF_to_ICRF(2, 2) * daz_BCFdy_ICRF;
            //doubleType daz_ICRFdz_ICRF = R_from_BCF_to_ICRF(2, 0) * dax_BCFdz_ICRF + R_from_BCF_to_ICRF(2, 1) * day_BCFdz_ICRF + R_from_BCF_to_ICRF(2, 2) * daz_BCFdz_ICRF;

            // A21 dadr
            this->acceleration_model->fx(3, 0) += (accel_position_Jacobian_ICRF(0, 0)) _GETVALUE;
            this->acceleration_model->fx(3, 1) += (accel_position_Jacobian_ICRF(0, 1)) _GETVALUE;
            this->acceleration_model->fx(3, 2) += (accel_position_Jacobian_ICRF(0, 2)) _GETVALUE;
            this->acceleration_model->fx(4, 0) += (accel_position_Jacobian_ICRF(1, 0)) _GETVALUE;
            this->acceleration_model->fx(4, 1) += (accel_position_Jacobian_ICRF(1, 1)) _GETVALUE;
            this->acceleration_model->fx(4, 2) += (accel_position_Jacobian_ICRF(1, 2)) _GETVALUE;
            this->acceleration_model->fx(5, 0) += (accel_position_Jacobian_ICRF(2, 0)) _GETVALUE;
            this->acceleration_model->fx(5, 1) += (accel_position_Jacobian_ICRF(2, 1)) _GETVALUE;
            this->acceleration_model->fx(5, 2) += (accel_position_Jacobian_ICRF(2, 2)) _GETVALUE;


            // dadt
            // grab the explicit time partial derivatives of the rotation matrices
            static EMTG::math::Matrix <doubleType> d_spherical_accel_ICRF_dt(3, 3, 0.0);
            static EMTG::math::Matrix <doubleType> dR_from_ICRF_to_BCF_dt;
            static EMTG::math::Matrix <doubleType> dR_from_BCF_to_ICRF_dt;
            dR_from_ICRF_to_BCF_dt = this->acceleration_model->my_universe->LocalFrame.get_dR_from_ICRF_to_BCF_dt();
            dR_from_BCF_to_ICRF_dt = this->acceleration_model->my_universe->LocalFrame.get_dR_from_BCF_to_ICRF_dt();

            static EMTG::math::Matrix <doubleType> dR_from_BCF_to_;

            // accel_ICRF = R_BCF_2_ICRF * accel_BCF
            // (daccel_BCF / dt) = (daccel_BCF / dr_BCF) * (dr_BCF / dt) 
            // r_BCF = R_ICRF_2_BCF * r_ICRF
            // (dr_BCF / dt) = (dR_ICRF_2_BCF / dt) * r_ICRF + (R_ICRF_2_BCF * dr_ICRF / dt)

            // daccel_ICRF / dt = (dR_BCF_2_ICRF / dt) * accel_BCF + (daccel_BCF / dt)
            //                  = (dR_BCF_2_ICRF / dt) * accel_BCF + (daccel_BCF / dr_BCF) * (dr_BCF / dt)
            //                  = (dR_BCF_2_ICRF / dt) * accel_BCF + (daccel_BCF / dr_BCF) * (dR_ICRF_2_BCF / dt) * r_ICRF + (R_ICRF_2_BCF * dr_ICRF / dt)
            //                  = (dR_BCF_2_ICRF / dt) * accel_BCF + (daccel_BCF / dr_BCF) * (dR_ICRF_2_BCF / dt) * r_ICRF + (R_ICRF_2_BCF * 0)
            //                  = (dR_BCF_2_ICRF / dt) * accel_BCF + (daccel_BCF / dr_BCF) * (dR_ICRF_2_BCF / dt) * r_ICRF
            d_spherical_accel_ICRF_dt = (dR_from_BCF_to_ICRF_dt * this->spherical_harmonic_term_acceleration_BCF
                + R_from_BCF_to_ICRF * accel_position_Jacobian_BCF * dR_from_ICRF_to_BCF_dt * this->parent_gravity_term->r_body2sc);

            this->acceleration_model->fx(3, 7) += (d_spherical_accel_ICRF_dt(0)) _GETVALUE;
            this->acceleration_model->fx(4, 7) += (d_spherical_accel_ICRF_dt(1)) _GETVALUE;
            this->acceleration_model->fx(5, 7) += (d_spherical_accel_ICRF_dt(2)) _GETVALUE;



            // dump them into this->acceleration_model->da_cb2scdPropVars(0, propVar)

        }// end computeAccelerationTerm(bool)

        void SphericalHarmonicTerm::populateInstrumentationFile(std::ofstream& acceleration_model_file)
        {
            acceleration_model_file << "," << sqrt(this->term_acceleration(0) * this->term_acceleration(0)
                + this->term_acceleration(1) * this->term_acceleration(1)
                + this->term_acceleration(2) * this->term_acceleration(2));
            for (size_t k = 0; k < this->term_acceleration.get_n(); ++k)
            {
                acceleration_model_file << "," << this->term_acceleration(k);
            }
        }
    }//close namespace Astrodynamics
}//close namespace EMTG