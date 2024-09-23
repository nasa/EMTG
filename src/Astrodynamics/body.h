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
        class body
        {
        public:
            //constructor
            body() {}; //default constructor
            body(const missionoptions& options
#ifdef SPLINE_EPHEM
                , SplineEphem::universe* SplineEphemUniverse
#endif
                );

            body(const int& ibody_code,
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
            virtual ~body();

            //comparator
            inline bool operator==(body OtherBody)
            {
                return this->spice_ID == OtherBody.spice_ID ? true : false;
            }

            //**************************************
            //methods
        
            //function to load new data into the body
            void load_body_data(const int& ibody_code,
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
                const missionoptions& options);



            //function to find the body state vector at epoch
            virtual int locate_body(const doubleType& epoch,
                                    doubleType* state,
                                    const bool& need_deriv,
                                    const missionoptions& options) const;

            //function to locate a point on the surface in BCI coordinates
            //assumes a spherical cow
            void locate_point_on_surface(const doubleType& ETepoch,
                                         const double& latitude,
                                         const double& longitude,
                                         math::Matrix<doubleType>& point_on_sphere_BCI);

            //function to determine the latitude and longitude of a point in BCI coordinates
            //also assumes a spherical cow
            void determine_latitude_longitude_for_point_in_BCI( const doubleType& ETepoch,
                                                                const math::Matrix<doubleType>& point_on_sphere_BCI,
                                                                doubleType& latitude,
                                                                doubleType& longitude);

            //function to determine the velocity vector of the surface (or atmosphere) from a BCI position vector
            //this assumes that the atmosphere is moving at the same angular rate as the body's spin
            //in BCI coordinates
            void determine_atmosphere_or_surface_velocity_vector(const doubleType& ETepoch,
                                                                const math::Matrix<doubleType>& position_BCI,
                                                                math::Matrix<doubleType>& atmosphere_or_surface_velocity);

            //function to print body to screen, for debug purposes
            void print_body_to_screen(std::string filename) const;

            //set scale factor pointer
            inline void set_X_scale_factors(std::vector<double>* X_scale_factors) { this->X_scale_factors = X_scale_factors; }

            //comparisons
            virtual bool operator== (const body& OtherBody) const;
            virtual bool operator!= (const body& OtherBody) const;

            double getEphemerisWindowOpen() const { return this->MySplineUniverse->getEphemerisWindowOpen(this->spice_ID, this->central_body_spice_ID); }
            double getEphemerisWindowClose() const { return this->MySplineUniverse->getEphemerisWindowClose(this->spice_ID, this->central_body_spice_ID); }

            //**************************************
            //fields

            //the following fields are read in
            int spice_ID;
            int body_code;
            int central_body_spice_ID;
            double central_body_radius;
            double central_body_LU;
            std::string central_body_name;
            std::string name;
            std::string short_name;
            double mass; //in kg
            double radius; //in km
            double flattening_coefficient;//flattening coefficient
            double J2;
            double J2_ref_radius;
            double AbsoluteMagnitude;
            double albedo;
            double reference_epoch; //in MJD
            double SMA; //in km
            double ECC;
            double INC; //in radians, convert this from input degrees
            double RAAN; //in radians, convert this from input degrees
            double AOP; //in radians, convert this from input degrees
            double MA; //in radians, convert this from input degrees
            double Period; //in days
            doubleType true_anomaly;
            double universe_mu; //gravitational constant for the local universe
            double minimum_safe_flyby_altitude; //in km, if this number is <= 0, then the object will not be a member of the flyby menu

            frame body_frame; //local reference frame, for use in calculating rotation matrices, can be evaluated at any epoch

            //the following fields are computed internal to the class
            std::string SPICE_frame;
            double mu; //gravitational constant
            double r_SOI; //radius of the sphere of influence/hill sphere
            double exclusion_radius; //radius of the sphere of exclusion inside which integrated propagators should turn off this body to avoid singularities

        protected:
            std::vector<double> reference_state;
            //pointer to spline ephem object
#ifdef SPLINE_EPHEM
            SplineEphem::universe *MySplineUniverse;
#endif
        private:
            void getCoverageWindow();      
            int body_ephemeris_source; //0: DE-405, 1: SPICE, 2: static
            double ephemeris_window_open;
            double ephemeris_window_close;
            double state_at_beginning_of_ephemeris[6];
            double state_at_end_of_ephemeris[6];
            std::vector<double>* X_scale_factors;

        };//end class body

    }//close namespace Astrodynamics
}//close namespace EMTG