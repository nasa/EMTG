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
#include "SpacecraftAccelerationModelTerm.h"
#include "doubleType.h"
#include "SolarRadiationPressureTerm.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        // constructors
        SolarRadiationPressureTerm::SolarRadiationPressureTerm(SpacecraftAccelerationModel * acceleration_model_in) :
            SpacecraftAccelerationModelTerm::SpacecraftAccelerationModelTerm(acceleration_model_in)
        {
            this->sc_area = this->acceleration_model->my_options->spacecraft_area / (1000.0 * 1000.0); 
            this->K = this->acceleration_model->my_options->solar_percentage; // percentage of Sun seen by spacecraft [0, 1]
            this->Cr = this->acceleration_model->my_options->coefficient_of_reflectivity; // [0, 2]
            this->solar_flux = this->acceleration_model->my_options->solar_flux; // W/m^2 = kg/s^3 --> equal to 3.846e26 Watts / (4 * pi * r_km^2)
            this->speed_of_light_vac = this->acceleration_model->my_options->speed_of_light_vac / 1000.0; // speed of light in a vacuum km/s

            this->SRP_coeff = this->Cr * this->sc_area * this->K * this->solar_flux / this->speed_of_light_vac 
                              * (this->acceleration_model->my_options->AU * this->acceleration_model->my_options->AU);
        }

        SolarRadiationPressureTerm::~SolarRadiationPressureTerm() {}

        // methods
        void SolarRadiationPressureTerm::computeAccelerationTerm()
        {

            doubleType spacecraft_distance_from_sun = this->acceleration_model->r_sun2sc_norm;
            this->SRP_force_norm = this->SRP_coeff / (spacecraft_distance_from_sun * spacecraft_distance_from_sun); // normalized force

            doubleType SRP_force[3];
            for (size_t k = 0; k < 3; ++k) 
            {
                SRP_force[k] = SRP_force_norm * this->acceleration_model->r_sun2sc(k) / spacecraft_distance_from_sun;
                this->term_acceleration(k) = SRP_force[k] / this->acceleration_model->spacecraft_mass;
                this->acceleration_model->acceleration(k) += this->term_acceleration(k);
            }

        }// end computeAccelerationTerm()

        void SolarRadiationPressureTerm::computeAccelerationTerm(const bool & generate_derivatives)
        {
            
            this->computeAccelerationTerm();

            doubleType spacecraft_mass = this->acceleration_model->spacecraft_mass;
            doubleType spacecraft_distance_from_sun = this->acceleration_model->r_sun2sc_norm;
            doubleType one_over_r_sun_sc3 = 1.0 / (spacecraft_distance_from_sun * spacecraft_distance_from_sun * spacecraft_distance_from_sun);
            doubleType one_over_r_sun_sc5 = one_over_r_sun_sc3 * 1.0 / (spacecraft_distance_from_sun * spacecraft_distance_from_sun);
            doubleType SRPcoeff_over_msc = this->SRP_coeff / spacecraft_mass;
            
            static math::Matrix<double> thing(3, 3, 0.0);
            thing(0, 0) = (SRPcoeff_over_msc * (-3.0 * this->acceleration_model->r_sun2sc(0) * this->acceleration_model->r_sun2sc(0) * one_over_r_sun_sc5 + one_over_r_sun_sc3)) _GETVALUE;
            thing(0, 1) = (SRPcoeff_over_msc * (-3.0 * this->acceleration_model->r_sun2sc(0) * this->acceleration_model->r_sun2sc(1) * one_over_r_sun_sc5)) _GETVALUE;
            thing(0, 2) = (SRPcoeff_over_msc * (-3.0 * this->acceleration_model->r_sun2sc(0) * this->acceleration_model->r_sun2sc(2) * one_over_r_sun_sc5)) _GETVALUE;
            thing(1, 0) = (SRPcoeff_over_msc * (-3.0 * this->acceleration_model->r_sun2sc(1) * this->acceleration_model->r_sun2sc(0) * one_over_r_sun_sc5)) _GETVALUE;
            thing(1, 1) = (SRPcoeff_over_msc * (-3.0 * this->acceleration_model->r_sun2sc(1) * this->acceleration_model->r_sun2sc(1) * one_over_r_sun_sc5 + one_over_r_sun_sc3)) _GETVALUE;
            thing(1, 2) = (SRPcoeff_over_msc * (-3.0 * this->acceleration_model->r_sun2sc(1) * this->acceleration_model->r_sun2sc(2) * one_over_r_sun_sc5)) _GETVALUE;
            thing(2, 0) = (SRPcoeff_over_msc * (-3.0 * this->acceleration_model->r_sun2sc(2) * this->acceleration_model->r_sun2sc(0) * one_over_r_sun_sc5)) _GETVALUE;
            thing(2, 1) = (SRPcoeff_over_msc * (-3.0 * this->acceleration_model->r_sun2sc(2) * this->acceleration_model->r_sun2sc(1) * one_over_r_sun_sc5)) _GETVALUE;
            thing(2, 2) = (SRPcoeff_over_msc * (-3.0 * this->acceleration_model->r_sun2sc(2) * this->acceleration_model->r_sun2sc(2) * one_over_r_sun_sc5 + one_over_r_sun_sc3)) _GETVALUE;

            // A21 dadr
            this->acceleration_model->fx(3, 0) += thing(0, 0);
            this->acceleration_model->fx(3, 1) += thing(0, 1);
            this->acceleration_model->fx(3, 2) += thing(0, 2);
            this->acceleration_model->fx(4, 0) += thing(1, 0);
            this->acceleration_model->fx(4, 1) += thing(1, 1);
            this->acceleration_model->fx(4, 2) += thing(1, 2);
            this->acceleration_model->fx(5, 0) += thing(2, 0);
            this->acceleration_model->fx(5, 1) += thing(2, 1);
            this->acceleration_model->fx(5, 2) += thing(2, 2);

            // A23 dadm
            this->acceleration_model->fx(3, 6) += (-this->SRP_force_norm / (spacecraft_mass * spacecraft_mass) * this->acceleration_model->r_sun2sc(0) / spacecraft_distance_from_sun) _GETVALUE;
            this->acceleration_model->fx(4, 6) += (-this->SRP_force_norm / (spacecraft_mass * spacecraft_mass) * this->acceleration_model->r_sun2sc(1) / spacecraft_distance_from_sun) _GETVALUE;
            this->acceleration_model->fx(5, 6) += (-this->SRP_force_norm / (spacecraft_mass * spacecraft_mass) * this->acceleration_model->r_sun2sc(2) / spacecraft_distance_from_sun) _GETVALUE;


            // Epoch and TOF derivative terms
            // TODO: refactor these
            static math::Matrix<double> da_cb2scdr_cb2sun_SRP_term(3, 3, 0.0);
            da_cb2scdr_cb2sun_SRP_term = -thing;

            this->acceleration_model->fx(3, 7) += da_cb2scdr_cb2sun_SRP_term(0, 0) * this->acceleration_model->dr_cb2sundcurrent_epoch(0)
                                               +  da_cb2scdr_cb2sun_SRP_term(0, 1) * this->acceleration_model->dr_cb2sundcurrent_epoch(1)
                                               +  da_cb2scdr_cb2sun_SRP_term(0, 2) * this->acceleration_model->dr_cb2sundcurrent_epoch(2);
                                                                                           
            this->acceleration_model->fx(4, 7) += da_cb2scdr_cb2sun_SRP_term(1, 0) * this->acceleration_model->dr_cb2sundcurrent_epoch(0)
                                               +  da_cb2scdr_cb2sun_SRP_term(1, 1) * this->acceleration_model->dr_cb2sundcurrent_epoch(1)
                                               +  da_cb2scdr_cb2sun_SRP_term(1, 2) * this->acceleration_model->dr_cb2sundcurrent_epoch(2);
                                                                                           
            this->acceleration_model->fx(5, 7) += da_cb2scdr_cb2sun_SRP_term(2, 0) * this->acceleration_model->dr_cb2sundcurrent_epoch(0)
                                               +  da_cb2scdr_cb2sun_SRP_term(2, 1) * this->acceleration_model->dr_cb2sundcurrent_epoch(1)
                                               +  da_cb2scdr_cb2sun_SRP_term(2, 2) * this->acceleration_model->dr_cb2sundcurrent_epoch(2);

        } // end computeAccelerationTerm(bool)

        void SolarRadiationPressureTerm::populateInstrumentationFile(std::ofstream & acceleration_model_file)
        {
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