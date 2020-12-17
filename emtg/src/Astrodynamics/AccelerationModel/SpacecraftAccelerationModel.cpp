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

#include "SpacecraftAccelerationModel.h"
#include "body.h"
#include "CentralBodyGravityTerm.h"
#include "doubleType.h"
#include "GravityTerm.h"
#include "journeyoptions.h"
#include "missionoptions.h"
#include "Spacecraft.h"
#include "universe.h"

namespace EMTG {
    namespace Astrodynamics {

        SpacecraftAccelerationModel::SpacecraftAccelerationModel(missionoptions * options_in,
                                                                 JourneyOptions * journey_options_in,
                                                                 universe * universe_in, 
                                                                 std::vector<std::string>* Xdescriptions,
                                                                 HardwareModels::Spacecraft * spacecraft_in,
                                                                 const size_t & STM_size_in,                                                                 
                                                                 const bool& needCentralBody) :
                                                                 AccelerationModel(options_in, 
                                                                                   journey_options_in, 
                                                                                   universe_in, 
                                                                                   Xdescriptions),
                                                                 my_spacecraft(spacecraft_in),
                                                                 thrust_term(this),
                                                                 STM_size(STM_size_in)
        {
            this->constructorInitialize(needCentralBody);
        }

        SpacecraftAccelerationModel::~SpacecraftAccelerationModel() {}

        void SpacecraftAccelerationModel::constructorInitialize(const bool& needCentralBody)
        {
            this->isHeliocentric = my_universe->central_body_name == "Sun";
            this->acceleration.resize(3, 1, 0.0);
            this->r_cb2sc.resize(3, 1, 0.0);
            this->v_cb2sc.resize(3, 1, 0.0);
            this->r_sun2sc.resize(3, 1, 0.0);
            this->v_sun2sc.resize(3, 1, 0.0);
            this->r_sun2cb.resize(3, 1, 0.0);
            this->v_sun2cb.resize(3, 1, 0.0);
            this->r_cb2sun.resize(3, 1, 0.0);
            this->v_cb2sun.resize(3, 1, 0.0);
            this->control.resize(4, 1, 0.0);

            // size the derivative containers
            this->fx.resize(this->STM_size, this->STM_size, 0.0);
            this->dr_sun2sc_normdr_cb2sc.resize(3, 1, 0.0);
            this->dr_sun2sc_normdr_cb2sun.resize(3, 1, 0.0);
            this->dr_cb2sundcurrent_epoch.resize(3, 1, 0.0);

            // we will always have central body gravity
            if (needCentralBody)
            {
                this->natural_acceleration_terms.push_back(new CentralBodyGravityTerm(this, &this->my_universe->central_body));
            }

            // add a GravityTerm for each perturbing body that we want to consider
            if (this->my_options->perturb_thirdbody)
            {
                for (size_t k = 0; k < this->my_journey_options->perturbation_bodies.size(); ++k)
                {
                    size_t body_code_to_find = this->my_journey_options->perturbation_bodies[k];
                    std::vector<body>::iterator it = std::find_if(this->my_universe->bodies.begin(), this->my_universe->bodies.end(),
                        [body_code_to_find](const body & body_info) -> bool { return body_info.body_code == body_code_to_find; });

                    // if a perturbing body is requested in the options file, but is not present in the Universe file,
                    // then inform the user of this discrepancy
                    if (it != this->my_universe->bodies.end())
                    {
                        this->natural_acceleration_terms.push_back(new GravityTerm(this, &this->my_universe->bodies[it - this->my_universe->bodies.begin()]));
                    }
                    else
                    {
                        throw std::invalid_argument("Third-body perturber index number " + std::to_string(this->my_journey_options->perturbation_bodies[k])
                            + " does not appear in the universe file: " + this->my_options->universe_folder + "/" + this->my_journey_options->journey_central_body + ".emtg_universe");
                    }

                }
            }

            // add gravity terms into the actual gravity terms bin
            for (size_t k = 0; k < this->natural_acceleration_terms.size(); ++k)
            {
                this->gravity_terms.push_back(static_cast<GravityTerm*>(&this->natural_acceleration_terms[k]));
            }

            // add SRP if it was requested
            if (this->my_options->perturb_SRP)
            {
                this->SRP_term = new SolarRadiationPressureTerm(this);
                this->natural_acceleration_terms.push_back(this->SRP_term);
            }

			// add aerodynamic drag if it was requested
			if (this->my_journey_options->perturb_drag)
			{
				this->drag_term = new AerodynamicDragTerm(this);
				this->natural_acceleration_terms.push_back(this->drag_term);
			}
        }

        void SpacecraftAccelerationModel::initializeAndComputeGeometry(const math::Matrix<doubleType> & spacecraft_state,
                                                                       const bool & generate_derivatives)
        {
            // the total acceleration will be modified 
            // by each of the individual terms, so we need to start from zeros
            this->acceleration.assign_zeros();

            // If the control overload of this method is called, then it will set this
            this->control_norm = 0.0;

            // reset the derivative containers
            if (generate_derivatives)
            {
                this->fx.assign_zeros();
                // Top block 3x3 identity
                this->fx(0, 3) = 1.0;
                this->fx(1, 4) = 1.0;
                this->fx(2, 5) = 1.0;
            }


            for (size_t k = 0; k < 3; ++k)
            {
                this->r_cb2sc(k) = spacecraft_state(k);
                this->v_cb2sc(k) = spacecraft_state(k + 3);
            }
            this->r_cb2sc_norm = this->r_cb2sc.norm();
            this->v_cb2sc_norm = this->v_cb2sc.norm();
            this->spacecraft_mass = spacecraft_state(6);

            this->computeScCbSunTriangle(generate_derivatives);

            // TODO: compute the power state here, instead of in the thruster term
        }

        void SpacecraftAccelerationModel::computeAcceleration(const math::Matrix<doubleType> & spacecraft_state,
                                                              const bool & generate_derivatives)
        {
            this->initializeAndComputeGeometry(spacecraft_state, generate_derivatives);

            // compute the individual terms that comprise the acceleration model
            // also generate entries in the propagation matrix fx if requested
            if (generate_derivatives)
            {
                // NOTE: we must call the thrust term first in order to set all of the appropriate
                // quantities in the spacecraft object
                ThrustControlLaw control_law = this->thrust_term.getThrustControlLaw();
                if (control_law == ThrustControlLaw::Velocity ||
                    control_law == ThrustControlLaw::AntiVelocity)
                {
                    this->thrust_term.computeAccelerationTerm(generate_derivatives);
                }
                for (size_t k = 0; k < this->natural_acceleration_terms.size(); ++k)
                {
                    this->natural_acceleration_terms[k].computeAccelerationTerm(generate_derivatives);
                }
            }
            else
            {
                ThrustControlLaw control_law = this->thrust_term.getThrustControlLaw();
                if (control_law == ThrustControlLaw::Velocity ||
                    control_law == ThrustControlLaw::AntiVelocity)
                {
                    this->thrust_term.computeAccelerationTerm();
                }
                for (size_t k = 0; k < this->natural_acceleration_terms.size(); ++k)
                {
                    this->natural_acceleration_terms[k].computeAccelerationTerm();
                }
            }

            // check to see if any acceleration terms have NaN'd for some reason
            for (size_t k = 0; k < this->acceleration.get_n(); ++k)
            {
                if (this->acceleration(k) != this->acceleration(k))
                {
                    throw std::runtime_error("NaN alert in the SpacecraftAccelerationModel!!");
                }
            }

        }

        void SpacecraftAccelerationModel::computeAcceleration(const math::Matrix<doubleType> & spacecraft_state,
                                                              const math::Matrix<doubleType> & control_in,
                                                              const bool & generate_derivatives)
        {
            this->control = control_in;
            this->initializeAndComputeGeometry(spacecraft_state, generate_derivatives);
            if (generate_derivatives)
            {
                // NOTE: we must call the thrust term first in order to set all of the appropriate
                // quantities in the spacecraft object
                this->thrust_term.computeAccelerationTerm(generate_derivatives);
                for (size_t k = 0; k < this->natural_acceleration_terms.size(); ++k)
                {
                    this->natural_acceleration_terms[k].computeAccelerationTerm(generate_derivatives);
                }
            }
            else
            {
                this->thrust_term.computeAccelerationTerm();
                for (size_t k = 0; k < this->natural_acceleration_terms.size(); ++k)
                {
                    this->natural_acceleration_terms[k].computeAccelerationTerm();
                }
            }

            // check to see if any acceleration terms have NaN'd for some reason
            for (size_t k = 0; k < this->acceleration.get_n(); ++k)
            {
                if (this->acceleration(k) != this->acceleration(k))
                {
                    throw std::runtime_error("NaN alert in the SpacecraftAccelerationModel!!");
                }
            }
        }

        void SpacecraftAccelerationModel::computeScCbSunTriangle(const bool & generate_derivatives)
        {

            if (!this->isHeliocentric)
            {
                // get the position vector of the sun
                std::vector<doubleType> temp_state(12, 0.0);
                this->my_universe->locate_central_body(this->current_epoch, temp_state.data(), *this->my_options, generate_derivatives);

#ifdef AD_INSTRUMENTATION
                std::vector<size_t> timevars = this->current_epoch.getDerivativeIndicies();
#endif
                for (size_t k = 0; k < 3; ++k)
                {
                    this->r_sun2cb(k) = temp_state[k];
                    this->v_sun2cb(k) = temp_state[k + 3];

                    // locate_central_body gives us the state of the central body relative to the sun 
                    // we want the state of the sun relative to the central body though (so its negation)
                    this->r_cb2sun(k) = -this->r_sun2cb(k);
                    this->v_cb2sun(k) = -this->v_sun2cb(k);
                    this->dr_cb2sundcurrent_epoch(k) = -(temp_state[k + 6]) _GETVALUE;

                    // perform the vector subtraction to compute the position of the spacecraft relative to the sun
                    this->r_sun2sc(k) = this->r_cb2sc(k) - this->r_cb2sun(k);
                    this->v_sun2sc(k) = this->v_cb2sc(k) - this->v_cb2sun(k);
                }
                this->r_sun2sc_norm = this->r_sun2sc.norm();
                this->v_sun2sc_norm = this->v_sun2sc.norm();
            }
            else
            {
                // if we are dealing with a heliocentric trajectory, 
                // then the spacecraft's position relative to the sun is the same as its position relative to the central body
                this->r_sun2sc = this->r_cb2sc;
                for (size_t k = 0; k < 3; ++k)
                {
                    this->dr_cb2sundcurrent_epoch(k) = 0.0;
                }
                this->r_sun2sc_norm = this->r_cb2sc_norm;
            } // end if isHeliocentric

              // sensitivity of the s/c distance from the sun to changes in s/c position w.r.t. the central body
            doubleType one_over_r_sun2sc_norm = 1.0 / this->r_sun2sc_norm;
            this->dr_sun2sc_normdr_cb2sc(0) = (this->r_sun2sc(0) * one_over_r_sun2sc_norm) _GETVALUE;
            this->dr_sun2sc_normdr_cb2sc(1) = (this->r_sun2sc(1) * one_over_r_sun2sc_norm) _GETVALUE;
            this->dr_sun2sc_normdr_cb2sc(2) = (this->r_sun2sc(2) * one_over_r_sun2sc_norm) _GETVALUE;

            this->dr_sun2sc_normdr_cb2sun = -this->dr_sun2sc_normdr_cb2sc;
        }

        void SpacecraftAccelerationModel::populateInstrumentationFile(std::ofstream & acceleration_model_file, 
                                                                      const math::Matrix<doubleType> & spacecraft_state, 
                                                                      const doubleType & epoch)
        {
            this->setEpoch(epoch);
            // TODO: power state needs to be computed here if we move it out of the thruster term
            this->computeAcceleration(spacecraft_state, false);
            
            acceleration_model_file << "\n";
            acceleration_model_file << this->current_epoch_JD;
            acceleration_model_file << "," << sqrt(spacecraft_state(0)*spacecraft_state(0) + spacecraft_state(1)*spacecraft_state(1) + spacecraft_state(2)*spacecraft_state(2));
            acceleration_model_file << "," << spacecraft_state(0) << "," << spacecraft_state(1) << "," << spacecraft_state(2);
            acceleration_model_file << "," << sqrt(spacecraft_state(3)*spacecraft_state(3) + spacecraft_state(4)*spacecraft_state(4) + spacecraft_state(5)*spacecraft_state(5));
            acceleration_model_file << "," << spacecraft_state(3) << "," << spacecraft_state(4) << "," << spacecraft_state(5);
            acceleration_model_file << "," << this->spacecraft_mass;
            acceleration_model_file << "," << sqrt(this->acceleration(0)*this->acceleration(0) + this->acceleration(1)*this->acceleration(1) + this->acceleration(2)*this->acceleration(2));
            acceleration_model_file << "," << this->acceleration(0) << "," << this->acceleration(1) << "," << this->acceleration(2);
            if (this->my_options->perturb_SRP)
            {
                this->SRP_term->populateInstrumentationFile(acceleration_model_file);
            }
			if (this->my_journey_options->perturb_drag)
			{
				this->drag_term->populateInstrumentationFile(acceleration_model_file);
			}
            for (size_t k = 0; k < this->gravity_terms.size(); ++k)
            {
                this->gravity_terms[k]->populateInstrumentationFile(acceleration_model_file);
            }
        }

        void SpacecraftAccelerationModel::populateInstrumentationFile(std::ofstream & acceleration_model_file, 
                                                                      const math::Matrix<doubleType> & spacecraft_state, 
                                                                      const doubleType & epoch,
                                                                      const math::Matrix<doubleType> & control)
        {
            this->populateInstrumentationFile(acceleration_model_file, spacecraft_state, epoch);
            this->thrust_term.populateInstrumentationFile(acceleration_model_file);
        }
    } // end Astrodynamics namespace
} // end EMTG namespace