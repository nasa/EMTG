// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2024 United States Government as represented by the
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
#include "CentralBodyGravityTerm.h"
#include "doubleType.h"
#include "GravityTerm.h"
#include "missionoptions.h"
#include "Spacecraft.h"
#include "SphericalHarmonicTerm.h"
#include "universe.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        // constructors
        CentralBodyGravityTerm::CentralBodyGravityTerm(SpacecraftAccelerationModel * acceleration_model_in, CentralBody * central_body_in) :
            GravityTerm::GravityTerm(acceleration_model_in, central_body_in)
        {
            /*
            if (this->acceleration_model->my_options->perturb_J2 && std::abs(this->acceleration_model->my_universe->central_body_J2) > 0.0)
            {
                this->gravitational_harmonic_terms.push_back(new SphericalHarmonicTerm(acceleration_model_in, this));
            }
            */
            if (this->acceleration_model->my_journey_options->perturb_central_body_gravity_harmonics)
            {
                this->spherical_harmonic_term = std::make_shared<SphericalHarmonicTerm>(acceleration_model_in, this);
                this->spherical_harmonic_term->setCentralBodyGM(this->my_body->harmonic_gravity_field->getCentralBodyGM());
                this->spherical_harmonic_term->setReferenceRadius(this->my_body->harmonic_gravity_field->getReferenceRadius());
            }
            else if (this->acceleration_model->my_options->perturb_J2)
            {
                this->spherical_harmonic_term = std::make_shared<SphericalHarmonicTerm>(acceleration_model_in, this);
                this->spherical_harmonic_term->setCentralBodyGM(this->my_body->mu);
                this->spherical_harmonic_term->setReferenceRadius(this->acceleration_model->my_universe->central_body_J2_reference_radius);
            }
        }

        CentralBodyGravityTerm::~CentralBodyGravityTerm()
        {
        }

        void CentralBodyGravityTerm::computeFrameDragAcceleration()
        {
            // Since we are the central body, then we do not want to compute the acceleration of
            // ourselves due to our own gravity. Things would...stop working.
            for (size_t k = 0; k < 3; ++k)
            {
                this->a_body_on_cb(k) = 0.0;
            }
        }

        //methods
        void CentralBodyGravityTerm::computeAccelerationTerm()
        {
            bool generate_derivatives = false;

            this->computePointMassGravityAcceleration(generate_derivatives);
            if (this->acceleration_model->my_journey_options->perturb_central_body_gravity_harmonics)
            {
                this->spherical_harmonic_term->computeAccelerationTerm();
            }
            else if (this->acceleration_model->my_options->perturb_J2)
            {
                this->spherical_harmonic_term->computeJ2Acceleration();
            }

        }//end computeAccelerationTerm()

        void CentralBodyGravityTerm::computeAccelerationTerm(const bool & generate_derivatives)
        {
            this->computePointMassGravityAcceleration(generate_derivatives);

            this->computePointMassGravityDerivatives();
            if (this->acceleration_model->my_journey_options->perturb_central_body_gravity_harmonics)
            {
                this->spherical_harmonic_term->computeAccelerationTerm(generate_derivatives);
            }
            else if (this->acceleration_model->my_options->perturb_J2)
            {
                this->spherical_harmonic_term->computeJ2Acceleration(generate_derivatives);
            }
        } // end computeAccelerationTerm()

        void CentralBodyGravityTerm::computePointMassGravityTimeDerivatives()
        {
            // DO NOTHING
            // For the case of the central body, these calculations are redundant, so
            // we zero out these terms in order to prevent double-counting them
            //this->acceleration_model->da_cb2scdPropVars.assign_zeros();
        }

        void CentralBodyGravityTerm::populateInstrumentationFile(std::ofstream & acceleration_model_file)
        {
            
            acceleration_model_file << "," << sqrt(this->term_acceleration(0)*this->term_acceleration(0)
                                                 + this->term_acceleration(1)*this->term_acceleration(1)
                                                 + this->term_acceleration(2)*this->term_acceleration(2));
            for (size_t k = 0; k < this->term_acceleration.get_n(); ++k)
            {
                acceleration_model_file << "," << this->term_acceleration(k);
            }

            // now request that the harmonic term report its acceleration instrumentation
            this->spherical_harmonic_term->populateInstrumentationFile(acceleration_model_file);
        }

        //methods
        
    }//close namespace Astrodynamics
}//close namespace EMTG