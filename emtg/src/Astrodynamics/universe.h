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

//header file for EMTG universe class

#pragma once

#include "doubleType.h"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "missionoptions.h"
#include "CentralBody.h"
#include "body.h"
#include "frame.h"
#include "EMTG_math.h"
#include "atmosphere.h"

#include "SpiceUsr.h"

#include "boost/ptr_container/ptr_vector.hpp"
#include "boost/algorithm/string.hpp"

#ifdef SPLINE_EPHEM
#include "SplineEphem_universe.h"
#endif

#ifdef AD_INSTRUMENTATION
#include <set>
#endif



namespace EMTG
{
    namespace Astrodynamics 
    {
        class universe
        {
        public:
            //constructor
            universe(); //default constructor
            universe(const size_t& j, 
                std::string universefile,
                const missionoptions& options
    #ifdef SPLINE_EPHEM
                , SplineEphem::universe* SplineEphemUniverse
    #endif
                );

            //destructor
            virtual ~universe();

            //**************************************
            //methods

            //get/set
            universe* get_nextUniverse() { return this->nextUniverse; }
            void set_nextUniverse(universe& inextUniverse) 
            {
                this->nextUniverse = &inextUniverse; 
            }

            //set scale factor pointer
            inline void set_X_scale_factors(std::vector<double>* X_scale_factors) 
            { 
                this->X_scale_factors = X_scale_factors;

                for (body& thisBody : this->bodies)
                    thisBody.set_X_scale_factors(X_scale_factors);
            }

            //function to find the central body state vector relative to the sun at epoch
            int locate_central_body(const doubleType& epoch, doubleType* state, const missionoptions& options, const bool& need_deriv) const;

            //function to print the universe
            void print_universe(std::string filename, const missionoptions& options) const;

            //method to return percent sun for a spacecraft orbiting the central body
            void get_percent_sun(math::Matrix<doubleType> spacecraft_state,
                math::Matrix<doubleType> sun_state, 
                const missionoptions& options,
                doubleType& percent_sun);

			// method to get central body reference angles
			void get_central_body_reference_angles(std::vector<double>& central_body_reference_angles_out);

            //**************************************
            //fields

            //the following fields are read in
            std::string central_body_name;
            int central_body_SPICE_ID;
            double central_body_radius;
            double central_body_J2;
            double central_body_J2_reference_radius;
            double mu;
            double central_body_flattening_coefficient;
            double LU;
            double r_SOI; //radius of the central body's sphere of influence
            double minimum_safe_distance; //minimum safe distance from the central body, in km
            
            //the following fields are computed internal to the class
            double TU;
            CentralBody central_body;
            std::vector<body> bodies;
            std::vector<size_t> flyby_menu; //vector containing the list of bodies, in the order that they appear in the Universe list, which can be used for flybys
            size_t size_of_flyby_menu; //this integer will be twice the number of flyby-capable bodies. If it is zero, then no flybys are possible in this universe.
            std::vector<size_t> perturbation_menu; //vector containing the list of bodies large enough to be considered for third-body perturbations
            frame LocalFrame; //local reference frame, for use in calculating rotation matrices
            math::Matrix<double> continuity_constraint_scale_factors;
            math::Matrix<double> COE_scale_factors;

			//boost::ptr_vector<EMTG::Astrodynamics::atmosphere> TheAtmosphere;
			std::shared_ptr<atmosphere> TheAtmosphere;



            //pointer to SplineEphem universe
    #ifdef SPLINE_EPHEM
            SplineEphem::universe* MySplineEphemUniverse;
    #endif
        private:
            // creates the CentralBody object, which is required by the acceleration model
            void createCentralBody(const missionoptions & options);
            std::vector<double> central_body_reference_angles;

        private:
            //function to load new data into the universe
            int load_universe_data(const size_t& j, std::string universefile, const missionoptions& options);

            //function to create the flyby menu - creates a list of bodies, by SPICE ID, which are flyby capable
            void create_flyby_and_perturbation_menus(const size_t& j, const missionoptions& options);

            //function to print the flyby menu
            void print_flyby_and_perturbation_menus(std::string filename, const missionoptions& options) const;

            universe* nextUniverse;
            std::vector<double>* X_scale_factors;

        };//end class universe

    }//close namespace Astrodynamics
}//close namespace EMTG