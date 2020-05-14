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


#ifndef ACCELERATION_MODEL_H
#define ACCELERATION_MODEL_H

#include "boost/ptr_container/ptr_vector.hpp"
#include "doubleType.h"
#include "EMTG_Matrix.h"
#include "missionoptions.h"
#include "journeyoptions.h"
#include "universe.h"

namespace EMTG {
    namespace Astrodynamics {
        
        class AccelerationModel
        {

        public:

            // constructor
            AccelerationModel(missionoptions * options_in, 
                              JourneyOptions * journey_options_in, 
                              universe * universe_in,
                              std::vector<std::string>* Xdescriptions);

            // destructor
            ~AccelerationModel();

            inline void setEpoch(const doubleType & epoch_in) { this->current_epoch = epoch_in; 
                                                                this->current_epoch_JD = (this->current_epoch / 86400.0 + 2400000.5);}
            inline void setEpochJD(const doubleType & epoch_JD_in) { this->current_epoch_JD = epoch_JD_in;
                                                                     this->current_epoch = (this->current_epoch_JD - 2400000.5) * 86400.0;}
            inline void setdCurrentEpochdPropVar(const double & dcurrent_epochdProp_var_in) { this->dcurrent_epochdProp_var = dcurrent_epochdProp_var_in; }
            inline void setdCurrentEpochdPropVarPrevious(const double & dcurrent_epochdProp_var_previous_in) 
            { 
                this->dcurrent_epochdProp_var_previous = dcurrent_epochdProp_var_previous_in; 
            }

            inline void setdstate_dProp_vars(const math::Matrix<double> & dstate_dProp_vars_in)
            {
                this->dstate_dProp_vars = dstate_dProp_vars_in;
            }

            inline void setdstate_dotdProp_varsPtr(math::Matrix<double> & dstate_dotdProp_varsPtr_in) { this->dstate_dotdProp_varsPtr = &dstate_dotdProp_varsPtr_in; }

            inline math::Matrix<doubleType> getAcceleration() const { return this->acceleration; }

            std::vector<std::string>* getXescriptions() { return this->Xdescriptions; }

            virtual void computeAcceleration(const math::Matrix<doubleType> & state_in,
                                             const bool & generate_derivatives,
                                             const bool & generate_time_derivatives = true) = 0;

            virtual void computeAcceleration(const math::Matrix<doubleType> & spacecraft_state,
                                             const math::Matrix<doubleType> & control,
                                             const bool & generate_derivatives,
                                             const bool & generate_time_derivatives = true) = 0;

            virtual void populateInstrumentationFile(std::ofstream & acceleration_model_file, 
                                                     const math::Matrix<doubleType> & spacecraft_state,
                                                     const doubleType & epoch) = 0;
            virtual void populateInstrumentationFile(std::ofstream & acceleration_model_file, 
                                                     const math::Matrix<doubleType> & spacecraft_state,  
                                                     const doubleType & epoch,
                                                     const math::Matrix<doubleType> & control) = 0;
            
            // pointers to EMTG configuration objects
            missionoptions * my_options;
            JourneyOptions * my_journey_options;
            universe * my_universe;

        protected:
            
            math::Matrix<doubleType> acceleration;

            // current epoch on MJD
            doubleType current_epoch;

            // current epoch in JD
            doubleType current_epoch_JD;

            double dcurrent_epochdProp_var;
            double dcurrent_epochdProp_var_previous;
            math::Matrix<double> dstate_dProp_vars;

            // partial derivatives of the state propagation matrix 
            // w.r.t. the propagation variables
            math::Matrix<double> * dstate_dotdProp_varsPtr;

            //pointer to Xdescriptions, in case we need to manually insert derivatives due to black boxes like GSL/SPICE
            std::vector<std::string>* Xdescriptions;

        };
    } // end Astrodynamics namespace
} // end EMTG namespace

#endif // ACCELERATION_MODEL_H