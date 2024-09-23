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


#ifndef SPACECRAFT_ACCELERATION_MODEL_H
#define SPACECRAFT_ACCELERATION_MODEL_H

#include <iomanip>
#include <iostream>
#include <map>

#include "AccelerationModel.h"
#include "boost/ptr_container/ptr_vector.hpp"
#include "doubleType.h"
#include "EMTG_enums.h"
#include "GravityTerm.h"
#include "journeyoptions.h"
#include "missionoptions.h"
#include "SolarRadiationPressureTerm.h"
#include "Spacecraft.h"
#include "SphericalHarmonicTerm.h"
//#include "STM.h"
#include "ThrustTerm.h"
#include "AerodynamicDragTerm.h"
#include "universe.h"

namespace EMTG {
    namespace Astrodynamics {
        
        class SpacecraftAccelerationModel : public AccelerationModel
        {
            friend class GravityTerm;
            friend class CentralBodyGravityTerm;
            friend class SolarRadiationPressureTerm;
            friend class SphericalHarmonicTerm;
            friend class ThrustTerm;
			friend class AerodynamicDragTerm;

        public:

            // constructor
            SpacecraftAccelerationModel(missionoptions * options_in,
                                        JourneyOptions * journey_options_in,
                                        universe * universe_in,
                                        std::vector<std::string>* Xdescriptions,
                                        HardwareModels::Spacecraft * spacecraft_in,
                                        const size_t & STM_size_in,
                                        const bool& needCentralBody = true);

            // destructor
            virtual ~SpacecraftAccelerationModel();

            // set the launch epoch of the spaceraft
            inline void setLaunchEpoch(const doubleType & launch_epoch_in)  { this->launch_epoch = launch_epoch_in; }
            inline void setDutyCycle(const double & duty_cycle_in)  { this->duty_cycle = duty_cycle_in; }
            inline void setThrustControlLaw(const ThrustControlLaw & control_law_in) { this->thrust_term.setThrustControlLaw(control_law_in); };

            inline ThrustControlLaw getThrustControlLaw() const { return this->thrust_term.getThrustControlLaw(); };

            // returns the net acceleration vector of the spacecraft
            inline math::Matrix<doubleType> getAccelerationVec() const { return this->acceleration; }
            
            // returns the acceleration vector due to only gravitational sources
            
            inline math::Matrix<doubleType> getGravityAccelerationVec() 
            { 
                math::Matrix<doubleType> temp_grav_accel (3, 1, 0.0);
                for (GravityTerm* gravity_term : this->gravity_terms)
                {
                    temp_grav_accel += gravity_term->getTermAcceleration();
                }
                return temp_grav_accel; 
            }
            

            // returns a vector of 3-tuples, where each tuple contains a body name, body mu, acceleration of s/c due to body
            inline std::map<std::string, std::tuple<double, math::Matrix<doubleType>>> getGravityAccelSources() const
            {
                //std::vector<std::tuple<std::string, double, math::Matrix<doubleType>>> grav_sources;
                std::map<std::string, std::tuple<double, math::Matrix<doubleType>>> grav_sources;
                for (GravityTerm* gravity_term : this->gravity_terms)
                {
                    body * current_body = gravity_term->getBody();
                    grav_sources.emplace(current_body->name, std::make_tuple(current_body->mu, gravity_term->getTermAcceleration()));
                }
                return grav_sources;
            }
            

            // returns the acceleration vector due to only solar radiation pressure
            inline math::Matrix<doubleType> getSRPAccelerationVec() const 
            {
                if (this->my_options->perturb_SRP)
                {
                    return this->SRP_term->getTermAcceleration();
                }
                else
                {
                    math::Matrix<doubleType> zeros (3, 1, 0.0);
                    return zeros;
                }
            }

            // returns the acceleration vector due to the spacecraft propulstion system
            inline math::Matrix<doubleType> getThrustAccelerationVec() const { return this->thrust_term.getTermAcceleration(); }

            // returns the current maximum possible propellant mass flow rate for the propulsion system
            inline doubleType getThrusterMaxMassFlowRate() const { return this->thrust_term.getMaxMassFlowRate(); }

            // returns the norm of the control parameters
            inline doubleType getControlNorm() const { return this->control_norm; }

            inline math::Matrix<doubleType> getControl() const { return this->control; }

            // returns the current duty cycle of the propulsion system (modelled continuously)
            inline double getDutyCycle() const { return this->duty_cycle; }

            inline math::Matrix<double> getfx() const { return this->fx; }
            
            inline size_t getSTMsize() const { return this->STM_size; }
			inline doubleType getCB2SC() const { return this->r_cb2sc_norm; }
			inline doubleType getMaxMassFlowRate() const { return this->thrust_term.getMaxMassFlowRate(); }

            // compute the current net acceleration of the spacecraft
            virtual void computeAcceleration(const math::Matrix<doubleType> & state_in,
                                             const bool & generate_derivatives);

            // compute the current net acceleration of the spacecraft given the input throttle levels
            virtual void computeAcceleration(const math::Matrix<doubleType> & spacecraft_state,
                                             const math::Matrix<doubleType> & control,
                                             const bool & generate_derivatives);
            
            virtual void populateInstrumentationFile(std::ofstream & acceleration_model_file, 
                                                     const math::Matrix<doubleType> & spacecraft_state,
                                                     const doubleType & epoch);
            virtual void populateInstrumentationFile(std::ofstream & acceleration_model_file,
                                                     const math::Matrix<doubleType> & spacecraft_state,
                                                     const doubleType & epoch,
                                                     const math::Matrix<doubleType> & control);

        private:

            void constructorInitialize(const bool& needCentralBody = true);
            void initializeAndComputeGeometry(const math::Matrix<doubleType> & spacecraft_state,
                                              const bool & generate_derivatives);
            
        protected:

            void computeScCbSunTriangle(const bool & generate_derivatives);

            // fields

            // spacecraft object
            HardwareModels::Spacecraft * my_spacecraft;

            // flag indicating that the central body is the Sun
            bool isHeliocentric;

            // launch epoch in MJD
            doubleType launch_epoch;

            // thruster duty cycle
            double duty_cycle;

            // current total mass of the spacecraft
            doubleType spacecraft_mass;

            // vector containing all of the acceleration terms
            boost::ptr_vector< SpacecraftAccelerationModelTerm > natural_acceleration_terms;

            // acceleration components due to any gravitating bodies
            //boost::ptr_vector< GravityTerm > gravity_terms;
            std::vector< GravityTerm * > gravity_terms;

            // acceleration component due to solar radiation pressure
            SolarRadiationPressureTerm * SRP_term;

            // acceleration component due to the finite-burn thruster 
            ThrustTerm thrust_term;

			AerodynamicDragTerm * drag_term;

            // control input to the thruster, if applicable
            math::Matrix<doubleType> control;
            doubleType control_norm;

            // position vector of the spacecraft relative to the central body and its norm
            math::Matrix<doubleType> r_cb2sc;
            doubleType r_cb2sc_norm;

            // velocity vector of the spacecraft relative to the central body and its norm
            math::Matrix<doubleType> v_cb2sc;
            doubleType v_cb2sc_norm;

            // position vector of the spacecraft relative to the sun and its norm
            math::Matrix<doubleType> r_sun2sc;
            doubleType r_sun2sc_norm;

            // velocity vector of the spacecraft relative to the sun and its norm
            math::Matrix<doubleType> v_sun2sc;
            doubleType v_sun2sc_norm;

            // position vector of the central body relative to the sun
            math::Matrix<doubleType> r_sun2cb;

            // position vector of the central body relative to the sun
            math::Matrix<doubleType> v_sun2cb;

            // position vector of the sun relative to the central body
            math::Matrix<doubleType> r_cb2sun;

            // velocity vector of the sun relative to the central body
            math::Matrix<doubleType> v_cb2sun;

            /////////////////////////
            // derivative containers
            /////////////////////////

            // STM size
            size_t STM_size;

            // state propagation matrix 
            math::Matrix<double> fx;

            math::Matrix<double> dr_sun2sc_normdr_cb2sc;
            math::Matrix<double> dr_sun2sc_normdr_cb2sun;
            math::Matrix<double> dr_cb2sundcurrent_epoch;
        };
    } // end Astrodynamics namespace
} // end EMTG namespace

#endif // SPACECRAFT_ACCELERATION_MODEL_H