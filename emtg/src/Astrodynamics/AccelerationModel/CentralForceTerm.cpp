// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2014 - 2017 United States Government as represented by the
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

#pragma once

#include "CentralForceTerm.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        //constructors
        CentralForceTerm::CentralForceTerm(missionoptions* myOptions,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& STMrows,
            const size_t& STMcolumns) :
            AccelerationModelTerm::AccelerationModelTerm(myOptions,
                myUniverse,
                mySpacecraft,
                journeyIndex,
                phaseIndex,
                STMrows,
                STMcolumns)
        {
            //central force gravity never affects the spacecraft mass
            this->mdot = 0.0;
        }

        //methods
        void CentralForceTerm::computeAccelerationTerm(const math::Matrix<doubleType>& spacecraft_state_relative_to_central_body,
            const math::Matrix <double> & dspacecraft_state_relative_to_central_bodydTOF,
            const doubleType& launch_epoch,
            const math::Matrix<doubleType>& control,
            const doubleType & epoch_step_left,
            std::vector <double> & depoch_left_segmentdTOF,
            const double & c,
            const doubleType & h,
            const double & dhdTOF,
            const bool& generate_derivatives)
        {
            //******************************acceleration vector
            doubleType spacecraft_distance_from_central_body = sqrt(spacecraft_state_relative_to_central_body(0) * spacecraft_state_relative_to_central_body(0)
                + spacecraft_state_relative_to_central_body(1) * spacecraft_state_relative_to_central_body(1) 
                + spacecraft_state_relative_to_central_body(2) * spacecraft_state_relative_to_central_body(2));
            
            doubleType& r = spacecraft_distance_from_central_body; //short name
            r = fabs(r) < EMTG::math::SMALL ? EMTG::math::sgn(r) * EMTG::math::SMALL : r;

            doubleType r3 = r*r*r;

            for (size_t i = 0; i < 3; ++i)
                acceleration_vector(i) = -1.0*spacecraft_state_relative_to_central_body(i) / r3;

            //******************************derivatives
            if (generate_derivatives) 
            {
                //***************************derivatives with respect to state
                //Top Row Identity
                //See Equation 45
                fx(0, 3) = 1.0;
                fx(1, 4) = 1.0;
                fx(2, 5) = 1.0;

                //A21 dadr (everything except for the third body terms, they are handled below)
                //Equation 46 (terms 1,2 and 5)
                double mu_CB = this->myUniverse->mu;
                doubleType spacecraft_distance_from_CB3 = spacecraft_distance_from_central_body * spacecraft_distance_from_central_body * spacecraft_distance_from_central_body;
                doubleType spacecraft_distance_from_CB5 = spacecraft_distance_from_CB3 * spacecraft_distance_from_central_body * spacecraft_distance_from_central_body;
                doubleType three_muCB_over_spacecraft_distance_from_CB5 = 3.0 * mu_CB / spacecraft_distance_from_CB5;
                doubleType muCB_over_spacecraft_distance_from_CB3 = 1.0 * mu_CB / spacecraft_distance_from_CB3;

                fx(3, 0) = three_muCB_over_spacecraft_distance_from_CB5 * spacecraft_state_relative_to_central_body(0) * spacecraft_state_relative_to_central_body(0) - muCB_over_spacecraft_distance_from_CB3;
                fx(3, 1) = three_muCB_over_spacecraft_distance_from_CB5 * spacecraft_state_relative_to_central_body(0) * spacecraft_state_relative_to_central_body(1);
                fx(3, 2) = three_muCB_over_spacecraft_distance_from_CB5 * spacecraft_state_relative_to_central_body(0) * spacecraft_state_relative_to_central_body(2);
                fx(4, 0) = three_muCB_over_spacecraft_distance_from_CB5 * spacecraft_state_relative_to_central_body(1) * spacecraft_state_relative_to_central_body(0);
                fx(4, 1) = three_muCB_over_spacecraft_distance_from_CB5 * spacecraft_state_relative_to_central_body(1) * spacecraft_state_relative_to_central_body(1) - muCB_over_spacecraft_distance_from_CB3;
                fx(4, 2) = three_muCB_over_spacecraft_distance_from_CB5 * spacecraft_state_relative_to_central_body(1) * spacecraft_state_relative_to_central_body(2);
                fx(5, 0) = three_muCB_over_spacecraft_distance_from_CB5 * spacecraft_state_relative_to_central_body(2) * spacecraft_state_relative_to_central_body(0);
                fx(5, 1) = three_muCB_over_spacecraft_distance_from_CB5 * spacecraft_state_relative_to_central_body(2) * spacecraft_state_relative_to_central_body(1);
                fx(5, 2) = three_muCB_over_spacecraft_distance_from_CB5 * spacecraft_state_relative_to_central_body(2) * spacecraft_state_relative_to_central_body(2) - muCB_over_spacecraft_distance_from_CB3;
            

                //***************************derivatives with respect to time

                //TODO!
            }// end derivative code
        }//end computeAccelerationTerm()
    }//close namespace Astrodynamics
}//close namespace EMTG