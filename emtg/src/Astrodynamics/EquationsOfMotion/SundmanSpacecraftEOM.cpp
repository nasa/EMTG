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

#include "doubleType.h"
#include "EMTG_Matrix.h"
#include "Integrand.h"
#include "SundmanSpacecraftEOM.h"

namespace EMTG {
    namespace Astrodynamics {

        // constructor
        SundmanSpacecraftEOM::SundmanSpacecraftEOM() : TimeDomainSpacecraftEOM()
        {
            this->acceleration.resize(3, 1, 0.0);
            this->velocity.resize(3, 1, 0.0);
            this->specific_angular_momentum.resize(3, 1, 0.0);
            this->dsundman_constantdState.resize(1, 6, 0.0);
			
			// default behaviour is to discretize by eccentric anomaly
			this->sundman_constant = 1.0;
			this->sundman_gamma = 2.0; // true anomaly
            //this->sundman_gamma = 1.0; // generalized anomaly
        }

        // destructor
		SundmanSpacecraftEOM::~SundmanSpacecraftEOM() {}

        void SundmanSpacecraftEOM::configureAccelerationModel(const math::Matrix<doubleType> & state)
        {
            this->spacecraft_acceleration_model->setEpoch(this->current_epoch);
        }

        void SundmanSpacecraftEOM::computeAcceleration(const math::Matrix<doubleType> & state,
                                                       const bool & STM_needed)
        {
            this->configureAccelerationModel(state);
            this->spacecraft_acceleration_model->computeAcceleration(state, STM_needed);
        }

        void SundmanSpacecraftEOM::computeAcceleration(const math::Matrix<doubleType> & state,
                                                       const math::Matrix<doubleType> & control,
                                                       const bool & STM_needed)
        {
            this->configureAccelerationModel(state);
            this->spacecraft_acceleration_model->computeAcceleration(state, control, STM_needed);
        }

        void SundmanSpacecraftEOM::evaluate(const math::Matrix<doubleType> & state,
                                            math::Matrix<doubleType> & state_dot,
                                            const bool & STM_needed)
        {
            // For s-domain EOM, the epoch is the integrated value in the state
            this->current_epoch = state(7);

			// call the time domain EOM first
			this->computeAcceleration(state, STM_needed);
			this->ballisticEOM(state, state_dot);
            this->velocity(0) = state_dot(0);
            this->velocity(1) = state_dot(1);
            this->velocity(2) = state_dot(2);
            this->mass_flow_rate = state_dot(6);

            //this->computeSundmanConstantTA(state);
            //this->computeSundmanConstantGA(state);
            this->computeSundmanConstantScaled(state);

			//  Sundman is just a multiplicative transformation of the time domain EOM
			state_dot(0) *= sundman_transform; // x
		    state_dot(1) *= sundman_transform; // y
		    state_dot(2) *= sundman_transform; // z
		    state_dot(3) *= sundman_transform; // vx
		    state_dot(4) *= sundman_transform; // vy
		    state_dot(5) *= sundman_transform; // vz
            state_dot(6) *= sundman_transform; // mass
            state_dot(7)  = sundman_transform; // time
			state_dot(8) *= sundman_transform; // fuel
			//state_dot(9) *= sundman_transform; // oxidizer
			

			if (STM_needed)
			{
				this->computeSundmanStatePropagationMatrix();
			}
        } // end Sundman EOM evaluate method

        void SundmanSpacecraftEOM::evaluate(const math::Matrix<doubleType> & state,
                                            math::Matrix<doubleType> & state_dot,
                                            const math::Matrix<doubleType> & control,
                                            const bool & STM_needed)
        {
            // For s-domain EOM, the epoch is the integrated value in the state
            this->current_epoch = state(7);

			// call the time domain EOM first
			this->computeAcceleration(state, control, STM_needed);
			this->ballisticEOM(state, state_dot);
			this->propulsionEOM(state_dot, control);
            this->velocity(0) = state_dot(0);
            this->velocity(1) = state_dot(1);
            this->velocity(2) = state_dot(2);
            this->mass_flow_rate = state_dot(6);

            //this->computeSundmanConstantTA(state);
            //this->computeSundmanConstantGA(state);
            this->computeSundmanConstantScaled(state);

			//  Sundman is just a multiplicative transformation of the time domain EOM
			state_dot(0) *= sundman_transform; // x
			state_dot(1) *= sundman_transform; // y
			state_dot(2) *= sundman_transform; // z
			state_dot(3) *= sundman_transform; // vx
			state_dot(4) *= sundman_transform; // vy
			state_dot(5) *= sundman_transform; // vz
			state_dot(6) *= sundman_transform; // mass
            state_dot(7)  = sundman_transform;  // time
			state_dot(8) *= sundman_transform; // fuel
			state_dot(9) *= sundman_transform; // oxidizer / electric propellant

			if (STM_needed)
			{
				this->computeSundmanStatePropagationMatrix();
			}

        } // end Sundman EOM with control evaluate method


        void SundmanSpacecraftEOM::computeSundmanConstantScaled(const math::Matrix<doubleType> & state)
        {
            // compute the specific Sundman constant for a scaled anomaly step distribution
            doubleType r = this->spacecraft_acceleration_model->getCB2SC();

            this->sundman_constant = 1.0 / this->spacecraft_acceleration_model->my_universe->LU;
            this->sundman_transform = this->sundman_constant * r;

            this->dsundman_constantdState(0) = (this->sundman_constant * state(0) / r) _GETVALUE;
            this->dsundman_constantdState(1) = (this->sundman_constant * state(1) / r) _GETVALUE;
            this->dsundman_constantdState(2) = (this->sundman_constant * state(2) / r) _GETVALUE;
            this->dsundman_constantdState(3) = 0.0;
            this->dsundman_constantdState(4) = 0.0;
            this->dsundman_constantdState(5) = 0.0;
        }

        void SundmanSpacecraftEOM::computeSundmanConstantGA(const math::Matrix<doubleType> & state)
        {
            // compute the specific Sundman constant for a generalized anomaly step distribution
            doubleType r = this->spacecraft_acceleration_model->getCB2SC();
            doubleType v = sqrt(state(3)*state(3) + state(4)*state(4) + state(5)*state(5));

            // semimajor axis
            double mu = this->spacecraft_acceleration_model->my_universe->central_body.mu;
            doubleType semimajor_axis = r / (2.0 - r * v * v / mu);

            this->sundman_constant = sqrt(fabs(semimajor_axis) / mu);
            this->sundman_transform = this->sundman_constant * r;

            doubleType thing = r*v*v - 2.0*mu;
            doubleType thing2 = thing * thing;
            
            this->dsundman_constantdState(0) = -((state(0) * (r*v*v - 3.0*mu)) / (sqrt(-r / thing) * thing2)) _GETVALUE;
            this->dsundman_constantdState(1) = -((state(1) * (r*v*v - 3.0*mu)) / (sqrt(-r / thing) * thing2)) _GETVALUE;
            this->dsundman_constantdState(2) = -((state(2) * (r*v*v - 3.0*mu)) / (sqrt(-r / thing) * thing2)) _GETVALUE;
            this->dsundman_constantdState(3) = (state(3)*r*r*r / (thing2 * this->sundman_constant)) _GETVALUE;
            this->dsundman_constantdState(4) = (state(4)*r*r*r / (thing2 * this->sundman_constant)) _GETVALUE;
            this->dsundman_constantdState(5) = (state(5)*r*r*r / (thing2 * this->sundman_constant)) _GETVALUE;
        }

        void SundmanSpacecraftEOM::computeSundmanConstantTA(const math::Matrix<doubleType> & state)
        {
            // compute the specific Sundman constant for a true anomaly step distribution

            // angular momentum = r cross v
            this->specific_angular_momentum(0) = state(1) * state(5) - state(2) * state(4); // y*zdot - z*ydot
            this->specific_angular_momentum(1) = state(2) * state(3) - state(0) * state(5); // z*xdot - x*zdot
            this->specific_angular_momentum(2) = state(0) * state(4) - state(1) * state(3); // x*ydot - y*xdot

            doubleType r = this->spacecraft_acceleration_model->getCB2SC();
            doubleType r2 = r * r;
            doubleType hnorm = specific_angular_momentum.norm();

            if (hnorm < 1.0e-7)
            {
                throw std::runtime_error("Sundman equations of motion is at a singularity due to angular momentum being < 1.0e-7.");
            }

            this->sundman_constant = 1.0 / hnorm;
            this->sundman_transform = this->sundman_constant * r2;
            doubleType coeff = -1.0 / (hnorm * hnorm * hnorm);

            this->dsundman_constantdState(0) = (coeff * (state(4) * this->specific_angular_momentum(2) - state(5) * this->specific_angular_momentum(1)) * r2
                                             + this->sundman_constant * 2.0 * state(0) ) _GETVALUE;
            this->dsundman_constantdState(1) = (coeff * (state(5) * this->specific_angular_momentum(0) - state(3) * this->specific_angular_momentum(2)) * r2
                                             + this->sundman_constant * 2.0 * state(1) ) _GETVALUE;
            this->dsundman_constantdState(2) = (coeff * (state(3) * this->specific_angular_momentum(1) - state(4) * this->specific_angular_momentum(0)) * r2
                                             + this->sundman_constant * 2.0 * state(2) ) _GETVALUE;
            this->dsundman_constantdState(3) = (coeff * (state(2) * this->specific_angular_momentum(1) - state(1) * this->specific_angular_momentum(2)) * r2 ) _GETVALUE;
            this->dsundman_constantdState(4) = (coeff * (state(0) * this->specific_angular_momentum(2) - state(2) * this->specific_angular_momentum(0)) * r2 ) _GETVALUE;
            this->dsundman_constantdState(5) = (coeff * (state(1) * this->specific_angular_momentum(0) - state(0) * this->specific_angular_momentum(1)) * r2 ) _GETVALUE;
        }

        void SundmanSpacecraftEOM::computeSundmanStatePropagationMatrix()
        {
            
            math::Matrix<double> thing1(this->STM.get_n(), this->STM.get_m(), 0.0);
            for (size_t i = 0; i < 3; ++i)
            {
                for (size_t j = 0; j < 6; ++j)
                {
                    thing1(i, j) = (this->velocity(i) * this->dsundman_constantdState(j)) _GETVALUE;
                    thing1(i + 3, j) = (this->acceleration(i) * this->dsundman_constantdState(j)) _GETVALUE;
                }
            }

            // mass row
            for (size_t j = 0; j < 6; ++j)
            {
                thing1(6, j) = (this->mass_flow_rate * this->dsundman_constantdState(j)) _GETVALUE;
            }
            
            // current epoch row
            for (size_t j = 0; j < 6; ++j)
            {
                thing1(7, j) = this->dsundman_constantdState(j);
            }

            // electric propellant / oxidizer row
            for (size_t j = 0; j < 6; ++j)
            {
                thing1(9, j) = (-this->mass_flow_rate * this->dsundman_constantdState(j)) _GETVALUE;
            }

            this->state_propagation_matrix = this->spacecraft_acceleration_model->getfx() * this->sundman_transform _GETVALUE;
            this->state_propagation_matrix += thing1;

        } // end Sundman variational equations method

    } // end Astrodynamics namespace
} // end EMTG namespace