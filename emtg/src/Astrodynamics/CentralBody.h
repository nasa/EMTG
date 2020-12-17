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

//header file for EMTG body class

#pragma once

#include "body.h"
#include "doubleType.h"

#include <string>
#include <vector>
#include <cmath>
#include <fstream>

#include "missionoptions.h"
#include "frame.h"

#include "EMTG_math.h"
#include "orbit_element_conversions.h"

#include "SpiceUsr.h"

#ifdef SPLINE_EPHEM
#include "SplineEphem_universe.h"
#endif


namespace EMTG
{
    namespace Astrodynamics
    {
        class CentralBody : public body
        {
        public:
            //constructor
            CentralBody() {}; //default constructor
            CentralBody(const missionoptions& options
#ifdef SPLINE_EPHEM
                , SplineEphem::universe* SplineEphemUniverse
#endif
            );

            CentralBody(const int& ibody_code,
                        const std::string& iname,
                        const std::string& ishortname,
                        const int& ispice_ID,
                        const double& imininum_altitude,
                        const double& imass,
                        const double& iradius,
                        const double& iJ2,
                        const double& iJ2_ref_radius,
                        const double& iflattening_coefficient,
                        const double& iAbsoluteMagnitude,
                        const double& ialbedo,
                        const double& iepoch,
                        std::vector<double>& ireference_angles,
                        std::vector<double>& iclassical_orbit_elements,
                        const double& iuniverse_mu,
                        const int& icentral_body_SPICE_ID,
                        const std::string& icentral_body_name,
                        const double& icentral_body_radius,
                        const double& icentral_body_LU,
                        frame& central_body_frame,
                        const missionoptions& options
#ifdef SPLINE_EPHEM
                , SplineEphem::universe* SplineEphemUniverse
#endif
            );

            //destructor
            virtual ~CentralBody();


            //**************************************
            //methods

            //function to find the body state vector at epoch
            int locate_body(const doubleType& epoch,
                            doubleType* state,
                            const bool& need_deriv,
                            const missionoptions& options) const override;

        };//end class body

    }//close namespace Astrodynamics
}//close namespace EMTG