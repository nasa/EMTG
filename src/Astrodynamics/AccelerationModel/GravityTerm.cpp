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

//central force term

#include "SpacecraftAccelerationModel.h"
#include "body.h"
#include "doubleType.h"
#include "GravityTerm.h"
#include "missionoptions.h"
#include "Spacecraft.h"
#include "universe.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        //constructors
        GravityTerm::GravityTerm(SpacecraftAccelerationModel * acceleration_model_in, body * body_in) :
            SpacecraftAccelerationModelTerm::SpacecraftAccelerationModelTerm(acceleration_model_in), my_body(body_in)
        {
            this->r_body2sc.resize(3, 1, 0.0);
            this->v_body2sc.resize(3, 1, 0.0);
            this->r_cb2body.resize(3, 1, 0.0);
            this->v_cb2body.resize(3, 1, 0.0);
            this->a_body_on_cb.resize(3, 1, 0.0);

            // derivative containers
            this->da_cb2scdr_cb2body.resize(3, 3, 0.0);
            this->dr_cb2bodydcurrent_epoch.resize(3, 1, 0.0);
        }

        void GravityTerm::computeScBodyCBtriangle(const bool & generate_derivatives)
        {
            // determine where the gravitating body is relative to the central body 
            //std::cout << this->my_body->name << std::setprecision(16) << " " << "epoch: " << this->acceleration_model->current_epoch << " ";
            std::vector<doubleType> temp_state(12, 0.0);
            this->my_body->locate_body(this->acceleration_model->current_epoch,
                temp_state.data(),
                generate_derivatives,
                *this->acceleration_model->my_options);

            for (size_t k = 0; k < 3; ++k)
            {
                this->r_cb2body(k) = temp_state[k];
                this->v_cb2body(k) = temp_state[k + 3];
                this->dr_cb2bodydcurrent_epoch(k) = (temp_state[k + 6]) _GETVALUE;
            }


            this->r_cb2body_norm = this->r_cb2body.norm();
            this->v_cb2body_norm = this->v_cb2body.norm();


            // determine where the spacecraft is relative to the gravitating body and its velocity w.r.t. that body
            for (size_t k = 0; k < 3; ++k)
            {
                this->r_body2sc(k) = this->acceleration_model->r_cb2sc(k) - this->r_cb2body(k);
                this->v_body2sc(k) = this->acceleration_model->v_cb2sc(k) - this->v_cb2body(k);
            }
            this->r_body2sc_norm = this->r_body2sc.norm();
            this->v_body2sc_norm = this->v_body2sc.norm();

        }

        void GravityTerm::computeFrameDragAcceleration()
        {
            for (size_t k = 0; k < 3; ++k)
            {
                this->a_body_on_cb(k) = -(this->my_body->mu) / (this->r_cb2body_norm * this->r_cb2body_norm) * this->r_cb2body(k) / this->r_cb2body_norm;
            }
        }

        void GravityTerm::computePointMassGravityAcceleration(const bool & generate_derivatives)
        {
            // If we are the Sun, then the main SpacecraftAccelerationModel code will have 
            // already taken care of the ephemeris call and vector math for us
            if (this->my_body->name == "Sun")
            {
                this->r_cb2body = this->acceleration_model->r_cb2sun;
                this->v_cb2body = this->acceleration_model->v_cb2sun;
                this->r_cb2body_norm = this->r_cb2body.norm();
                this->v_cb2body_norm = this->v_cb2body.norm();
                this->dr_cb2bodydcurrent_epoch = this->acceleration_model->dr_cb2sundcurrent_epoch;
                this->r_body2sc = this->acceleration_model->r_sun2sc;
                this->r_body2sc_norm = this->acceleration_model->r_sun2sc_norm;
                this->v_body2sc = this->acceleration_model->v_sun2sc;
                this->v_body2sc_norm = this->acceleration_model->v_sun2sc_norm;
            }
            else 
            {
                this->computeScBodyCBtriangle(generate_derivatives);
            }
                         
            //if (this->r_body2sc_norm > this->my_body->radius)
            //{
                this->a_body_on_sc_norm = -(this->my_body->mu) / (this->r_body2sc_norm * this->r_body2sc_norm);

                this->computeFrameDragAcceleration();

                //********************************************************
                //
                // Acceleration vector modification: Perturbing gravitational body
                //
                //********************************************************
                for (size_t k = 0; k < 3; ++k)
                {
                    this->term_acceleration(k) = this->a_body_on_sc_norm * this->r_body2sc(k) / this->r_body2sc_norm + this->a_body_on_cb(k);
                    this->acceleration_model->acceleration(k) += this->term_acceleration(k);
                }
            //}
        }

        //methods
        void GravityTerm::computeAccelerationTerm()
        {
            bool generate_derivatives = false;

            this->computePointMassGravityAcceleration(generate_derivatives);

        }//end computeAccelerationTerm()

        void GravityTerm::computeAccelerationTerm(const bool & generate_derivatives)
        {
            this->computePointMassGravityAcceleration(generate_derivatives);

            this->computePointMassGravityDerivatives();

        } // end computeAccelerationTerm()

        void GravityTerm::computePointMassGravityDerivatives()
        {
            if (this->r_body2sc_norm > this->my_body->radius)
            {
                double mu_body = this->my_body->mu;
                doubleType r_body2sc_norm3 = this->r_body2sc_norm * this->r_body2sc_norm * this->r_body2sc_norm;
                doubleType r_body2sc_norm5 = r_body2sc_norm3 * r_body2sc_norm * this->r_body2sc_norm;
                doubleType three_mu_over_r_body2sc_norm5 = 3.0 * mu_body / r_body2sc_norm5;
                doubleType mu_over_r_body2sc_norm3 = mu_body / r_body2sc_norm3;


                //****************************************
                //
                //State Propagation (fx) matrix calculation
                //
                //****************************************

                this->acceleration_model->fx(3, 0) += (three_mu_over_r_body2sc_norm5 * this->r_body2sc(0) * this->r_body2sc(0) - mu_over_r_body2sc_norm3) _GETVALUE;
                this->acceleration_model->fx(3, 1) += (three_mu_over_r_body2sc_norm5 * this->r_body2sc(0) * this->r_body2sc(1)) _GETVALUE;
                this->acceleration_model->fx(3, 2) += (three_mu_over_r_body2sc_norm5 * this->r_body2sc(0) * this->r_body2sc(2)) _GETVALUE;
                this->acceleration_model->fx(4, 0) += (three_mu_over_r_body2sc_norm5 * this->r_body2sc(1) * this->r_body2sc(0)) _GETVALUE;
                this->acceleration_model->fx(4, 1) += (three_mu_over_r_body2sc_norm5 * this->r_body2sc(1) * this->r_body2sc(1) - mu_over_r_body2sc_norm3) _GETVALUE;
                this->acceleration_model->fx(4, 2) += (three_mu_over_r_body2sc_norm5 * this->r_body2sc(1) * this->r_body2sc(2)) _GETVALUE;
                this->acceleration_model->fx(5, 0) += (three_mu_over_r_body2sc_norm5 * this->r_body2sc(2) * this->r_body2sc(0)) _GETVALUE;
                this->acceleration_model->fx(5, 1) += (three_mu_over_r_body2sc_norm5 * this->r_body2sc(2) * this->r_body2sc(1)) _GETVALUE;
                this->acceleration_model->fx(5, 2) += (three_mu_over_r_body2sc_norm5 * this->r_body2sc(2) * this->r_body2sc(2) - mu_over_r_body2sc_norm3) _GETVALUE;               

                this->computePointMassGravityTimeDerivatives();
            }
        }

        void GravityTerm::computePointMassGravityTimeDerivatives()
        {
            double mu_body = this->my_body->mu;
            doubleType r_body2sc_norm3 = this->r_body2sc_norm * this->r_body2sc_norm * this->r_body2sc_norm;
            doubleType r_body2sc_norm5 = r_body2sc_norm3 * r_body2sc_norm * this->r_body2sc_norm;
            doubleType three_mu_over_r_body2sc_norm5 = 3.0 * mu_body / r_body2sc_norm5;
            doubleType mu_over_r_body2sc_norm3 = mu_body / r_body2sc_norm3;
            doubleType r_cb2body_norm3 = this->r_cb2body_norm * this->r_cb2body_norm * this->r_cb2body_norm;
            doubleType mu_over_r_cb2body_norm3 = mu_body / r_cb2body_norm3;
            doubleType r_cb2body_norm5 = r_cb2body_norm3 * this->r_cb2body_norm * this->r_cb2body_norm;
            doubleType three_mu3B_over_r_cb2body_norm5 = 3.0 * mu_body / r_cb2body_norm5;

            // daxdr3B
            this->da_cb2scdr_cb2body(0, 0) = (-three_mu_over_r_body2sc_norm5 * this->r_body2sc(0) * this->r_body2sc(0) + mu_over_r_body2sc_norm3 + three_mu3B_over_r_cb2body_norm5 * this->r_cb2body(0) * this->r_cb2body(0) - mu_over_r_cb2body_norm3) _GETVALUE;
            this->da_cb2scdr_cb2body(0, 1) = (-three_mu_over_r_body2sc_norm5 * this->r_body2sc(0) * this->r_body2sc(1)                           + three_mu3B_over_r_cb2body_norm5 * this->r_cb2body(0) * this->r_cb2body(1)) _GETVALUE;
            this->da_cb2scdr_cb2body(0, 2) = (-three_mu_over_r_body2sc_norm5 * this->r_body2sc(0) * this->r_body2sc(2)                           + three_mu3B_over_r_cb2body_norm5 * this->r_cb2body(0) * this->r_cb2body(2)) _GETVALUE;

            // daydr3B
            this->da_cb2scdr_cb2body(1, 0) = (-three_mu_over_r_body2sc_norm5 * this->r_body2sc(1) * this->r_body2sc(0)                           + three_mu3B_over_r_cb2body_norm5 * this->r_cb2body(1) * this->r_cb2body(0)) _GETVALUE;
            this->da_cb2scdr_cb2body(1, 1) = (-three_mu_over_r_body2sc_norm5 * this->r_body2sc(1) * this->r_body2sc(1) + mu_over_r_body2sc_norm3 + three_mu3B_over_r_cb2body_norm5 * this->r_cb2body(1) * this->r_cb2body(1) - mu_over_r_cb2body_norm3) _GETVALUE;
            this->da_cb2scdr_cb2body(1, 2) = (-three_mu_over_r_body2sc_norm5 * this->r_body2sc(1) * this->r_body2sc(2)                           + three_mu3B_over_r_cb2body_norm5 * this->r_cb2body(1) * this->r_cb2body(2)) _GETVALUE;

            // dazdr3B
            this->da_cb2scdr_cb2body(2, 0) = (-three_mu_over_r_body2sc_norm5 * this->r_body2sc(2) * this->r_body2sc(0)                           + three_mu3B_over_r_cb2body_norm5 * this->r_cb2body(2) * this->r_cb2body(0)) _GETVALUE;
            this->da_cb2scdr_cb2body(2, 1) = (-three_mu_over_r_body2sc_norm5 * this->r_body2sc(2) * this->r_body2sc(1)                           + three_mu3B_over_r_cb2body_norm5 * this->r_cb2body(2) * this->r_cb2body(1)) _GETVALUE;
            this->da_cb2scdr_cb2body(2, 2) = (-three_mu_over_r_body2sc_norm5 * this->r_body2sc(2) * this->r_body2sc(2) + mu_over_r_body2sc_norm3 + three_mu3B_over_r_cb2body_norm5 * this->r_cb2body(2) * this->r_cb2body(2) - mu_over_r_cb2body_norm3) _GETVALUE;

            this->acceleration_model->fx(3, 7) += this->da_cb2scdr_cb2body(0, 0) * this->dr_cb2bodydcurrent_epoch(0)
                                               + this->da_cb2scdr_cb2body(0, 1) * this->dr_cb2bodydcurrent_epoch(1)
                                               + this->da_cb2scdr_cb2body(0, 2) * this->dr_cb2bodydcurrent_epoch(2);

            this->acceleration_model->fx(4, 7) += this->da_cb2scdr_cb2body(1, 0) * this->dr_cb2bodydcurrent_epoch(0)
                                               + this->da_cb2scdr_cb2body(1, 1) * this->dr_cb2bodydcurrent_epoch(1)
                                               + this->da_cb2scdr_cb2body(1, 2) * this->dr_cb2bodydcurrent_epoch(2);

            this->acceleration_model->fx(5, 7) += this->da_cb2scdr_cb2body(2, 0) * this->dr_cb2bodydcurrent_epoch(0)
                                               + this->da_cb2scdr_cb2body(2, 1) * this->dr_cb2bodydcurrent_epoch(1)
                                               + this->da_cb2scdr_cb2body(2, 2) * this->dr_cb2bodydcurrent_epoch(2);

        }

        void GravityTerm::populateInstrumentationFile(std::ofstream & acceleration_model_file)
        {
            acceleration_model_file << "," << this->r_cb2body_norm;
            for (size_t k = 0; k < this->term_acceleration.get_n(); ++k)
            {
                acceleration_model_file << "," << this->r_cb2body(k);
            }
            acceleration_model_file << "," << this->v_cb2body_norm;
            for (size_t k = 0; k < this->v_cb2body.get_n(); ++k)
            {
                acceleration_model_file << "," << this->v_cb2body(k);
            }
            acceleration_model_file << "," << this->r_body2sc_norm;
            for (size_t k = 0; k < this->r_body2sc.get_n(); ++k)
            {
                acceleration_model_file << "," << this->r_body2sc(k);
            }
            acceleration_model_file << "," << this->v_body2sc_norm;
            for (size_t k = 0; k < this->v_body2sc.get_n(); ++k)
            {
                acceleration_model_file << "," << this->v_body2sc(k);
            }
            acceleration_model_file << "," << sqrt(this->term_acceleration(0)*this->term_acceleration(0)
                                                 + this->term_acceleration(1)*this->term_acceleration(1)
                                                 + this->term_acceleration(2)*this->term_acceleration(2));
            for (size_t k = 0; k < this->term_acceleration.get_n(); ++k)
            {
                acceleration_model_file << "," << this->term_acceleration(k);
            }
        }

    } // close namespace Astrodynamics
} // close namespace EMTG