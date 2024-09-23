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

#ifndef THRUST_TERM_H
#define THRUST_TERM_H

#include "doubleType.h"
#include "SpacecraftAccelerationModel.h"
#include "SpacecraftAccelerationModelTerm.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        class ThrustTerm : public SpacecraftAccelerationModelTerm
        {
        public:
            //constructor
            ThrustTerm(SpacecraftAccelerationModel * acceleration_model_in);

            ~ThrustTerm();

            //methods

            inline void setThrustControlLaw(const ThrustControlLaw & control_law_in) { this->control_law = control_law_in; };
            inline ThrustControlLaw getThrustControlLaw() const { return this->control_law; };

            void processControlLogic();
            virtual void computeAccelerationTerm() override;
            virtual void computeAccelerationTerm(const bool & generate_derivatives) override;
            void populateInstrumentationFile(std::ofstream & acceleration_model_file) override;
            ThrustTerm* clone() const override { return new ThrustTerm(*this); } 

            inline doubleType getMaxMassFlowRate() const { return this->max_mass_flow_rate; }

        protected:
            //fields
            ThrustControlLaw control_law;

            // spacecraft state w.r.t. the central body
            math::Matrix<doubleType> r_cb2sc;
            doubleType r_cb2sc_norm;
            math::Matrix<doubleType> v_cb2sc;
            doubleType v_cb2sc_norm;            

            doubleType u_x;
            doubleType u_y;
            doubleType u_z;
            doubleType u_command;
            doubleType max_thrust;
            doubleType max_mass_flow_rate;
            doubleType power;
            doubleType number_of_active_engines;
            double dTdP;
            double dmdotdP;
            double dTdu_command;
            double dmdotdu_command;
            double dPdr_sun;
            double dPdt;

        };
    } // end Astrodynamics namespace
} // end EMTG namespace
#endif