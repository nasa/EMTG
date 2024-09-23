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
#include "TimeDomainSpacecraftEOM.h"

namespace EMTG {
    namespace Astrodynamics {

        // constructor
        TimeDomainSpacecraftEOM::TimeDomainSpacecraftEOM() : Integrand() 
        {
            this->acceleration.resize(3, 1, 0.0);
        }

        // destructor
        TimeDomainSpacecraftEOM::~TimeDomainSpacecraftEOM() {}


        void TimeDomainSpacecraftEOM::configureAccelerationModel(const math::Matrix<doubleType> & state)
        {
            this->spacecraft_acceleration_model->setEpoch(this->current_epoch);
        }

		void TimeDomainSpacecraftEOM::computeAcceleration(const math::Matrix<doubleType> & state, 
			                                              const bool & STM_needed)
		{
            this->configureAccelerationModel(state);
			this->spacecraft_acceleration_model->computeAcceleration(state, STM_needed);
		}

        void TimeDomainSpacecraftEOM::computeAcceleration(const math::Matrix<doubleType> & state, 
                                                          const math::Matrix<doubleType> & control,
			                                              const bool & STM_needed)
		{
            this->configureAccelerationModel(state);
			this->spacecraft_acceleration_model->computeAcceleration(state, control, STM_needed);
		}

        void TimeDomainSpacecraftEOM::evaluate(const math::Matrix<doubleType> & state,
                                               math::Matrix<doubleType> & state_dot,
                                               const bool & STM_needed)
        {
            // For t-domain EOM, the epoch IS the independent variable
            this->current_epoch = this->current_independent_variable;
            this->computeAcceleration(state, STM_needed);
            this->ballisticEOM(state, state_dot);

            ThrustControlLaw control_law = this->spacecraft_acceleration_model->getThrustControlLaw();
            if (control_law == ThrustControlLaw::Velocity ||
                control_law == ThrustControlLaw::AntiVelocity)
            {
                this->propulsionEOM(state_dot, this->spacecraft_acceleration_model->getControl());
            }

            state_dot(7) = 1.0; // epoch

			if (STM_needed)
			{
                this->state_propagation_matrix = this->spacecraft_acceleration_model->getfx();
				//this->computeVariationalEquations(state, state_dot);
			}
        } // end evaluate without control method

        void TimeDomainSpacecraftEOM::evaluate(const math::Matrix<doubleType> & state,
                                               math::Matrix<doubleType> & state_dot,
                                               const math::Matrix<doubleType> & control,
                                               const bool & STM_needed)
        {
            // For t-domain EOM, the epoch IS the independent variable
            this->current_epoch = this->current_independent_variable;
            this->computeAcceleration(state, control, STM_needed);
			this->ballisticEOM(state, state_dot);
			this->propulsionEOM(state_dot, control);

            state_dot(7) = 1.0;  // epoch

            if (STM_needed)
            {
                this->state_propagation_matrix = this->spacecraft_acceleration_model->getfx();
				//this->computeVariationalEquations(state, state_dot);
            }
        } // end evaluate with control method

        void TimeDomainSpacecraftEOM::ballisticEOM(const math::Matrix<doubleType> & state,
                                                   math::Matrix<doubleType> & state_dot)
        {
            // state_dot needs to be assigned to zero here, otherwise stale values
            // might persist if an instance of this EOM class is being re-used by multiple entities
            state_dot.assign_zeros();

            this->acceleration = this->spacecraft_acceleration_model->getAccelerationVec();
            state_dot(0) = state(3);        // x
            state_dot(1) = state(4);        // y
            state_dot(2) = state(5);        // z
            state_dot(3) = acceleration(0); // xddot
            state_dot(4) = acceleration(1); // yddot
            state_dot(5) = acceleration(2); // zddot
            state_dot(6) = 0.0;//mdot

            // Virtual chemical fuel equation
            // Even during a ballistic arc, we want to track ACS if applicable
            state_dot(8) = (this->spacecraft_acceleration_model->my_options->trackACS ?
                            this->spacecraft_acceleration_model->my_options->ACS_kg_per_day / 86400.0 : 0.0);
            //virtual electric propellant / oxidizer
            //state_dot(9) = 0.0;
        }

		void TimeDomainSpacecraftEOM::propulsionEOM(math::Matrix<doubleType> & state_dot, const math::Matrix<doubleType> & control)
		{
			// General propellant equation
			doubleType throttle = sqrt(control(0) * control(0)
				                     + control(1) * control(1)
				                     + control(2) * control(2));
			state_dot(6) = -throttle * this->spacecraft_acceleration_model->getThrusterMaxMassFlowRate();
			// TODO: verify that this is in the maxMassFlowRate 
			//       * this->spacecraft_acceleration_model->getDutyCycle();			

			// Virtual electric tank propellant / chemical oxidizer
			state_dot(9) = throttle * this->spacecraft_acceleration_model->getThrusterMaxMassFlowRate();
		}
    } // end Astrodynamics namespace
} // end EMTG namespace