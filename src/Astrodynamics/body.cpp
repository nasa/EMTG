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

#include "body.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        //constructor for default-Earth
        body::body(const missionoptions& options
#ifdef SPLINE_EPHEM
            , SplineEphem::universe* SplineEphemUniverse
#endif
        )
        {
#ifdef SPLINE_EPHEM
            this->MySplineUniverse = SplineEphemUniverse;
#endif
            std::vector<double> reference_angles;
            std::vector<double> orbit_elements;
            frame LocalFrame;
            reference_angles.push_back(0.0);
            reference_angles.push_back(-0.641);
            reference_angles.push_back(90.0);
            reference_angles.push_back(-0.557);
            reference_angles.push_back(190.147);
            reference_angles.push_back(360.9856235);
            orbit_elements.push_back(149653496.267);
            orbit_elements.push_back(0.0170423971756);
            orbit_elements.push_back(23.4390341482);
            orbit_elements.push_back(0.000185630044551);
            orbit_elements.push_back(101.741880245);
            orbit_elements.push_back(358.189140855);

            this->reference_state.resize(6, 0.0);

            this->load_body_data(3,
                "Earth",
                "E",
                399,
                300,
                3.986004418e+5,
                6378.136,
                0.0010826265,
                6378.136,
                1.0 / 298.257223563,
                10.0,
                0.367,
                51544.5 * 86400.0,
                reference_angles,
                orbit_elements,
                1.32712440018e+11,
                10,
                "Sun",
                4379000.0,
                149597870.691,
                LocalFrame,
                options);
        }

        //constructor takes in data and calls load_body_data()
        body::body(const int& ibody_code,
            const std::string& iname,
            const std::string& ishortname,
            const int& ispice_ID,
            const double& iminimum_altitude,
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
        )
        {
#ifdef SPLINE_EPHEM
            this->MySplineUniverse = SplineEphemUniverse;
#endif
            this->reference_state.resize(6, 0.0);

            this->load_body_data(ibody_code,
                iname,
                ishortname,
                ispice_ID,
                iminimum_altitude,
                imass,
                iradius,
                iJ2,
                iJ2_ref_radius,
                iflattening_coefficient,
                iAbsoluteMagnitude,
                ialbedo,
                iepoch,
                ireference_angles,
                iclassical_orbit_elements,
                iuniverse_mu,
                icentral_body_SPICE_ID,
                icentral_body_name,
                icentral_body_radius,
                icentral_body_LU,
                central_body_frame,
                options);
        }

        //destructor
        body::~body() {}


        //function to load new data into the body
        void body::load_body_data(const int& ibody_code,
            const std::string& iname,
            const std::string& ishortname,
            const int& ispice_ID,
            const double& imininum_altitude,
            const double& imu,
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
            frame& central_body_frame, const missionoptions& options)
        {
            //copy information from the inputs into the body
            this->name = iname;
            this->short_name = ishortname;
            this->universe_mu = iuniverse_mu;
            this->body_code = ibody_code;
            this->central_body_spice_ID = icentral_body_SPICE_ID;
            this->central_body_name = icentral_body_name;
            this->central_body_radius = icentral_body_radius;
            this->central_body_LU = icentral_body_LU;

            this->spice_ID = ispice_ID;
            this->minimum_safe_flyby_altitude = imininum_altitude;
            this->mu = imu;
            this->radius = iradius;
            this->reference_epoch = iepoch;

            this->J2 = iJ2;
            this->J2_ref_radius = iJ2_ref_radius;
            this->flattening_coefficient = iflattening_coefficient;
            this->AbsoluteMagnitude = iAbsoluteMagnitude;
            this->albedo = ialbedo;

            this->ephemeris_window_open = this->reference_epoch * 86400.0;
            this->ephemeris_window_close = 100000.0 * 86400.0;
            
            this->SMA = iclassical_orbit_elements[0];
            this->ECC = iclassical_orbit_elements[1];
            this->INC = iclassical_orbit_elements[2] * EMTG::math::deg2rad;
            this->RAAN = iclassical_orbit_elements[3] * EMTG::math::deg2rad;
            this->AOP = iclassical_orbit_elements[4] * EMTG::math::deg2rad;
            this->MA = iclassical_orbit_elements[5] * EMTG::math::deg2rad;
            this->Period = math::TwoPI * sqrt(SMA * SMA * SMA / this->universe_mu) / 86400.0;
            
            //determine which ephemeris to draw from
            if (options.ephemeris_source == 0)
            {
                this->body_ephemeris_source = 0; //use static ephemeris
            }
            else if (options.ephemeris_source == 1 || options.ephemeris_source == 2)
            {
                //first, check to see if the body exists in the currently loaded SPICE kernels
                double LT_dump;
                spkez_c(this->spice_ID, this->reference_epoch - (51544.5 * 86400.0), "J2000", "NONE", this->central_body_spice_ID, this->reference_state.data(), &LT_dump);

                if (!failed_c())
                {
                    //activate SPICE for this body

                    //replace the orbit elements with something drawn from SPICE
                    if (this->spice_ID != this->central_body_spice_ID)
                    {
                        std::vector<double> temp_elements(8);
                        oscelt_c(this->reference_state.data(), this->reference_epoch - (51544.5 * 86400.0), this->universe_mu, temp_elements.data());
                        this->SMA = temp_elements[0];
                        this->ECC = temp_elements[1];
                        this->INC = temp_elements[2];
                        this->RAAN = temp_elements[3];
                        this->AOP = temp_elements[4];
                        this->MA = temp_elements[5];
                        this->Period = math::TwoPI * sqrt(SMA * SMA * SMA / this->universe_mu) / 86400.0;
                    }

                    //activate SplineEphem if available
                    if (options.ephemeris_source == 2)
                    {
#ifdef SPLINE_EPHEM
                        //TODO some kind of check to see if the SplineEphem universe contains the information that we need
                        this->body_ephemeris_source = 2;
                        this->ephemeris_window_open = this->MySplineUniverse->getEphemerisWindowOpen(this->spice_ID, this->central_body_spice_ID);
                        this->ephemeris_window_close = this->MySplineUniverse->getEphemerisWindowClose(this->spice_ID, this->central_body_spice_ID);
#else
                        std::cout << "SplineEphem not available." << std::endl;
                        this->body_ephemeris_source = 1;
                        this->getCoverageWindow();
#endif
                    }
                    else
                    {
                        this->body_ephemeris_source = 1;

                        this->getCoverageWindow();
                    }
                }
                else
                {
                    std::cout << "Warning, body " << this->name << " is not defined in the kernel pool at reference epoch " << this->reference_epoch / 86400.0 << std::endl;
                    reset_c();
                    this->body_ephemeris_source = 0; //use static ephemeris
                }
            }



            this->body_frame.initialize(ireference_angles[0], ireference_angles[1], ireference_angles[2], ireference_angles[3], ireference_angles[4], ireference_angles[5]);

            //compute additional values
            this->mass = this->mu / options.G;
            this->r_SOI = this->SMA * (1.0 - this->ECC) * pow(this->mu / (3.0 * this->universe_mu), 0.333333333333333333333333);
            this->exclusion_radius = math::TwoPI * this->SMA * (1.0 - this->ECC) / 360.0;
        }

        //function to find the body state vector at epoch
        int body::locate_body(const doubleType& epoch, doubleType* state, const bool& need_deriv, const EMTG::missionoptions& options) const
        {
            
#ifdef AD_INSTRUMENTATION
            double statedouble[12];
            std::vector<size_t> timevars;
#endif

            switch (body_ephemeris_source)
            {
#ifdef SPLINE_EPHEM
            case 2: //SplineEphem

#ifdef AD_INSTRUMENTATION
                this->MySplineUniverse->getBody6StateAndDerivative(this->spice_ID, this->central_body_spice_ID, epoch _GETVALUE, statedouble);
                timevars = epoch.getDerivativeIndicies();
                for (size_t stateindex = 0; stateindex < 6; ++stateindex)
                {
                    state[stateindex].setValue(statedouble[stateindex]);
                    state[stateindex + 6].setValue(statedouble[stateindex + 6]);
                    //loop over derivative indices
                    for (size_t timeindex = 0; timeindex < timevars.size(); ++timeindex)
                    {
                        size_t timevar = timevars[timeindex];

                        double timeScale = epoch.getDerivative(timevar) / this->X_scale_factors->operator[](timevar);

                        state[stateindex].setDerivative(timevar, statedouble[stateindex + 6] * this->X_scale_factors->operator[](timevar) * timeScale);
                        //state[stateindex].setDerivative(timevar, statedouble[stateindex + 6]);
                    }
                }
#else
                if (need_deriv)
                    this->MySplineUniverse->getBody6StateAndDerivative(this->spice_ID, this->central_body_spice_ID, epoch, state);

                else
                    this->MySplineUniverse->getBody6State(this->spice_ID, this->central_body_spice_ID, epoch, state);
#endif
                break;
#endif
            case 1: //SPICE
                double LT_dump; 
                double statepert[6];

#ifdef AD_INSTRUMENTATION
                spkez_c(this->spice_ID, (epoch)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->central_body_spice_ID, statedouble, &LT_dump);
                spkez_c(this->spice_ID, (epoch)_GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", this->central_body_spice_ID, statepert, &LT_dump);
                timevars = epoch.getDerivativeIndicies();
                for (size_t stateindex = 0; stateindex < 6; ++stateindex)
                {
                    state[stateindex].setValue(statedouble[stateindex]);
                    state[stateindex + 6] = (statepert[stateindex] - state[stateindex]_GETVALUE) / 10.0;
                    //loop over derivative indices
                    for (size_t timeindex = 0; timeindex < timevars.size(); ++timeindex)
                    {
                        size_t timevar = timevars[timeindex];
                        state[stateindex].setDerivative(timevar, (statepert[stateindex] - state[stateindex]_GETVALUE) / (10.0));
                    }
                }
#else
                spkez_c(this->spice_ID, (epoch) _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->central_body_spice_ID, state, &LT_dump);

                if (need_deriv)
                {
                    spkez_c(this->spice_ID, (epoch) _GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", this->central_body_spice_ID, statepert, &LT_dump);
                    state[6] = (statepert[0] - state[0]) / (10.0);
                    state[7] = (statepert[1] - state[1]) / (10.0);
                    state[8] = (statepert[2] - state[2]) / (10.0);
                    state[9] = (statepert[3] - state[3]) / (10.0);
                    state[10] = (statepert[4] - state[4]) / (10.0);
                    state[11] = (statepert[5] - state[5]) / (10.0);
                }
#endif

                if (failed_c())//test for SPICE errors
                {
                    if (!options.quiet_NLP)
                        std::cout << "SPICE error detected" << std::endl;
                    reset_c();
                    throw std::runtime_error("SPICE error detected. SPICE failed to find body " + std::to_string(this->spice_ID)
                        + " with respect to " + std::to_string(this->central_body_spice_ID) + " on epoch " + std::to_string(epoch _GETVALUE / 86400.0) + ".");
                }
                break;
            }

            return 0;
        }


        //function to print body to screen (for debug purposes)
        void body::print_body_to_screen(std::string filename) const
        {
            std::ofstream outputfile(filename.c_str(), std::ios::app);
            outputfile << "Body name: " << this->name << std::endl;
            outputfile << "Short name: " << this->short_name << std::endl;
            outputfile << "Body position in menu: " << this->body_code << std::endl;
            outputfile << "SPICE ID: " << this->spice_ID << std::endl;
            outputfile << "Valid flyby target? " << (this->minimum_safe_flyby_altitude > 0.0 ? "True" : "False") << std::endl;
            if (this->minimum_safe_flyby_altitude > 0.0)
                outputfile << "Minimum safe flyby altitude (km) " << this->minimum_safe_flyby_altitude << std::endl;
            outputfile << "mu (km^3/s^2): " << this->mu << std::endl;
            outputfile << "Radius (km): " << this->radius << std::endl;
            outputfile << "Ephemeris source: " << this->body_ephemeris_source << std::endl;
            outputfile << "R_SOI: " << this->r_SOI << std::endl;
            outputfile << "Reference Epoch (MJD): " << this->reference_epoch << std::endl;
            outputfile << "SMA (km): " << this->SMA << std::endl;
            outputfile << "ECC: " << this->ECC << std::endl;
            outputfile << "INC (deg): " << this->INC * 180.0 / EMTG::math::PI << std::endl;
            outputfile << "RAAN (deg): " << this->RAAN * 180.0 / EMTG::math::PI << std::endl;
            outputfile << "AOP (deg): " << this->AOP * 180.0 / EMTG::math::PI << std::endl;
            outputfile << "MA (deg): " << this->MA * 180.0 / EMTG::math::PI << std::endl;
            outputfile << std::endl;

            outputfile.close();
        }

        //function to locate a point on the surface in BCI coordinates
        //assumes a spherical cow
        void body::locate_point_on_surface(const doubleType& ETepoch,
            const double& latitude,
            const double& longitude,
            math::Matrix<doubleType>& point_on_sphere_BCI)
        {
            //first locate the point in BCF coordinates
            math::Matrix<doubleType> point_on_sphere_BCF(3, 1, 0.0);
            point_on_sphere_BCF(0) = this->radius * cos(longitude) * cos(latitude);
            point_on_sphere_BCF(1) = this->radius * sin(longitude) * cos(latitude);
            point_on_sphere_BCF(2) = this->radius * sin(latitude);

            //then convert from BCF to BCI
            //first convert the point from BCI to BCF
            this->body_frame.rotate_frame_to_frame(ReferenceFrame::TrueOfDate_BCF, point_on_sphere_BCF, ReferenceFrame::TrueOfDate_BCI, point_on_sphere_BCI, ETepoch);
        }

        //function to determine the latitude and longitude of a point in BCI coordinates
        //also assumes a spherical cow
        void body::determine_latitude_longitude_for_point_in_BCI(const doubleType& ETepoch,
            const math::Matrix<doubleType>& point_on_sphere_BCI,
            doubleType& latitude,
            doubleType& longitude)
        {
            //first convert the point from BCI to BCF
            math::Matrix<doubleType> point_on_sphere_BCF;
            this->body_frame.rotate_frame_to_frame(ReferenceFrame::TrueOfDate_BCI, point_on_sphere_BCI, ReferenceFrame::TrueOfDate_BCF, point_on_sphere_BCF, ETepoch);

            //then compute the latitude and longitude of the point
            doubleType projection_to_equatorial_plane = sqrt(point_on_sphere_BCF(0) * point_on_sphere_BCF(0) + point_on_sphere_BCF(1) * point_on_sphere_BCF(1));
            latitude = atan2(point_on_sphere_BCF(2), projection_to_equatorial_plane);
            longitude = atan2(point_on_sphere_BCF(1), point_on_sphere_BCF(0));
        }

        //function to determine the velocity vector of the surface (or atmosphere) from a BCI position vector
        //this assumes that the atmosphere is moving at the same angular rate as the body's spin
        //in BCI coordinates
        void body::determine_atmosphere_or_surface_velocity_vector(const doubleType& ETepoch,
            const math::Matrix<doubleType>& position_BCI,
            math::Matrix<doubleType>& atmosphere_or_surface_velocity)
        {
            //the angular velocity of the body is given in radians per day
            //we just need to find the distance from the spin axis at this altitude and latitude
            doubleType distance_from_body_center = position_BCI.norm();
            doubleType latitude = 0.0, longitude = 0.0;
            this->determine_latitude_longitude_for_point_in_BCI(ETepoch,
                position_BCI,
                latitude,
                longitude);

            doubleType distance_from_spin_axis = distance_from_body_center * cos(latitude);

            //now we can find the magnitude of the atmosphere or surface velocity
            //note we have to convert from radians per day to radians per second
            doubleType velocity_magnitude = distance_from_spin_axis * this->body_frame.getWdot() / 86400.0;

            //the unit vector in the direction of velocity is the BCI position vector crossed with the spin axis vector
            //6-21-2016: changed the order of the cross product because we think the planet is rotating in the wrong direction
            atmosphere_or_surface_velocity = this->body_frame.get_zhat().unitcross(position_BCI) * velocity_magnitude;
        }

        //comparator
        bool body::operator== (const body& OtherBody) const
        {
            //compare three fields for accuracy
            if (this->name == OtherBody.name && this->spice_ID == OtherBody.spice_ID && this->mu == OtherBody.mu)
            {
                return true;
            }
            return false;
        }

        bool body::operator!= (const body& OtherBody) const
        {
            return !(*this == OtherBody);
        }

        //get coverage window
        void body::getCoverageWindow()
        {
            //use spkcov to get coverage window
            const size_t  FILSIZ = 256;
            const size_t  LNSIZE = 81;
            const size_t  MAXCOV = 100000;
            const size_t  WINSIZ = (2 * MAXCOV);
            const size_t  TIMLEN = 51;
            const size_t MAXOBJ = 10000;

            SPICEDOUBLE_CELL(cover, WINSIZ);
            SPICEINT_CELL(ids, MAXOBJ);

            SpiceBoolean found = false;
            bool FoundBody = false;

            SpiceChar file[FILSIZ];
            SpiceChar source[FILSIZ];
            SpiceChar type[LNSIZE];

            SpiceDouble b = 0;
            SpiceDouble e = 0;

            SpiceInt count;
            SpiceInt handle;
            SpiceInt i;
            SpiceInt idcode = this->spice_ID;
            SpiceInt niv;

            /*
            Find out how many kernels are loaded.Loop over the
            kernels : for each loaded SPK file, add its coverage
            for `idcode', if any, to the coverage window.
            */
            ktotal_c("SPK", &count);

            double EarliestEpoch = this->ephemeris_window_open;
            double LatestEpoch = this->ephemeris_window_open;
            double temp_latest = this->ephemeris_window_open;

            for (i = 0; i < count; ++i)
            {
                kdata_c(i, "SPK", FILSIZ, LNSIZE, FILSIZ,
                    file, type, source, &handle, &found);


                //let's see if the object we want is in this file
                scard_c(0, &ids);
                spkobj_c(file, &ids);

                for (size_t j = 0; j < card_c(&ids); ++j)
                {
                    if (SPICE_CELL_ELEM_I(&ids, j) == this->spice_ID)
                    {
                        FoundBody = true;
                        scard_c(0, &cover);
                        spkcov_c(file, idcode, &cover);

                        niv = wncard_c(&cover);
                        FoundBody = true;
                        for (size_t k = 0; k < niv; ++k)
                        {
                            //Get the endpoints of the ith interval.
                            wnfetd_c(&cover, k, &b, &e);

                            if (b + (51544.5 * 86400.0) < EarliestEpoch)
                                EarliestEpoch = b + (51544.5 * 86400.0);
                            if (e + (51544.5 * 86400.0) > LatestEpoch)
                                LatestEpoch = e + (51544.5 * 86400.0);
                        }
                    }
                }

                this->ephemeris_window_open = fmax(this->ephemeris_window_open, EarliestEpoch);
                temp_latest = fmax(temp_latest, LatestEpoch);
            }

            this->ephemeris_window_close = fmin(this->ephemeris_window_close, temp_latest);

            if (!found)
            {
                std::cout << "The SPICE kernel pool does not contain enough information to create an ephemeris for body " << this->spice_ID << " with respect to " << this->central_body_spice_ID << "." << std::endl;
            }

        }

    }//close namespace Astrodynamics
}//close namespace EMTG