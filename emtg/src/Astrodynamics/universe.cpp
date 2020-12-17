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
#include "CentralBody.h"
#include "universe.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        universe::universe()
        {
        
            this->central_body_J2 = 0.0;
            this->central_body_flattening_coefficient = 0.0;
            this->central_body_J2_reference_radius = 0.0;
        
        }

        //constructor to load a data file
        universe::universe(const size_t& j,
            std::string universe_file,
            const missionoptions& options
    #ifdef SPLINE_EPHEM
            , SplineEphem::universe* SplineEphemUniverse
    #endif
        )
        {

            this->central_body_J2 = 0.0;
            this->central_body_flattening_coefficient = 0.0;
            this->central_body_J2_reference_radius = 0.0;

    #ifdef SPLINE_EPHEM
            this->MySplineEphemUniverse = SplineEphemUniverse;
    #endif
            load_universe_data(j, universe_file, options);

            this->continuity_constraint_scale_factors = math::Matrix<double>(8, 1, std::vector<double>({ 1.0 / this->LU,
                                                                                                         1.0 / this->LU,
                                                                                                         1.0 / this->LU,
                                                                                                         1.0 / (this->LU / this->TU),
                                                                                                         1.0 / (this->LU / this->TU),
                                                                                                         1.0 / (this->LU / this->TU),
                                                                                                         1.0 / options.maximum_mass,
                                                                                                         1.0 / this->TU }));

            this->COE_scale_factors = math::Matrix<double>(8, 1, std::vector<double>({ 1.0 / this->LU,
                                                                                       1.0,
                                                                                       1.0,
                                                                                       1.0,
                                                                                       1.0,
                                                                                       1.0,
                                                                                       1.0 / options.maximum_mass,
                                                                                       1.0 / this->TU }));
        }

        //destructor
        universe::~universe() {}

        //function to load a data file
        int universe::load_universe_data(const size_t& j, std::string universefile, const missionoptions& options)
        {
            std::ifstream inputfile(universefile.c_str());
            int linenumber = 0;
            std::string choice;
            std::string peek;
            double value;
            char dump_buffer[1024];
            bool contains_J2_AbsoluteMagnitude_Albedo = false;
            bool J2_set = false;
            bool J2_ref_radius_set = false;

            if (!inputfile.is_open())
            {
                std::cout << "Failure to read " << universefile << std::endl;
                return 2;
            }

            while (!inputfile.eof())
            {
                peek = inputfile.peek();
                if (peek == "#" || peek == "\r" || peek == "\n") {
                    //comment or blank line, do not parse
                    inputfile.getline(dump_buffer, 1024);
                    if (boost::contains(dump_buffer, "J2"))
                    {
                        contains_J2_AbsoluteMagnitude_Albedo = true;
                    }
                    ++linenumber;
                }
                else
                {
                    inputfile >> choice;

                    if (choice == "central_body_name")
                        inputfile >> this->central_body_name;
                    else if (choice == "central_body_SPICE_ID")
                        inputfile >> this->central_body_SPICE_ID;
                    else if (choice == "mu")
                    {
                        inputfile >> this->mu;
                    }
                    else if (choice == "central_body_radius")
                    {
                        inputfile >> this->central_body_radius;
                    }
                    else if (choice == "central_body_J2")
                    {
                        J2_set = true;
                        inputfile >> this->central_body_J2;
                    }
                    else if (choice == "central_body_J2_reference_radius")
                    {
                        J2_ref_radius_set = true;
                        inputfile >> this->central_body_J2_reference_radius;
                    }
                    else if (choice == "central_body_flattening_coefficient")
                    {
                        inputfile >> this->central_body_flattening_coefficient;
                    }
                    else if (choice == "LU")
                    {
                        inputfile >> this->LU;

                        //compute TU, which can be overridden later
                        this->TU = sqrt(LU*LU*LU / mu);
                    }
                    else if (choice == "TU")
                    {
                        inputfile >> value;
                        if (value > 0.0)
                            this->TU = value;
                    }
                    else if (choice == "reference_angles")
                    {
                        this->central_body_reference_angles.resize(6);
                        inputfile >> this->central_body_reference_angles[0];
                        inputfile >> this->central_body_reference_angles[1];
                        inputfile >> this->central_body_reference_angles[2];
                        inputfile >> this->central_body_reference_angles[3];
                        inputfile >> this->central_body_reference_angles[4];
                        inputfile >> this->central_body_reference_angles[5];

                        if (central_body_name == "EARTH")
                            LocalFrame.initialize();
                        else
                            LocalFrame.initialize(this->central_body_reference_angles[0] * math::deg2rad, 
                                                  this->central_body_reference_angles[1] * math::deg2rad, 
                                                  this->central_body_reference_angles[2] * math::deg2rad, 
                                                  this->central_body_reference_angles[3] * math::deg2rad, 
                                                  this->central_body_reference_angles[4] * math::deg2rad, 
                                                  this->central_body_reference_angles[5] * math::deg2rad);
                    }
                    else if (choice == "r_SOI")
                    {
                        inputfile >> this->r_SOI;
                    }
                    else if (choice == "minimum_safe_distance")
                    {
                        inputfile >> this->minimum_safe_distance;
                    }
                    else if (choice == "convert_elements_from_central_body_frame")
                    {
                        inputfile >> value;
                        
                    }
                    else if (choice == "begin_body_list")
                    {
                        std::string tempname;
                        std::string tempshortname;
                        int temp_bodycode;
                        int temp_SPICEnum;
                        double temp_minimum_altitude;
                        double temp_mu;
                        double temp_radius;
                        double temp_epoch;
                        std::vector<double> temp_reference_angles(6);
                        std::vector<double> temp_elements(6);

                        double temp_J2 = 0.0010826265;
                        double temp_J2_ref_radius = 6378.0;
                        double temp_AbsoluteMagnitude = 10.0;
                        double temp_albedo = 0.367;

                        do
                        {
                            tempname.clear();

                            //if not yet reached end of list:
                            inputfile >> tempname;
                            inputfile >> tempshortname;
                            inputfile >> temp_bodycode;
                            if (tempname == "end_body_list")
                                break;
                            inputfile >> temp_SPICEnum;
                            inputfile >> temp_minimum_altitude;
                            inputfile >> temp_mu;
                            inputfile >> temp_radius;
                            if (contains_J2_AbsoluteMagnitude_Albedo)
                            {
                                inputfile >> temp_J2;
                                inputfile >> temp_AbsoluteMagnitude;
                                inputfile >> temp_albedo;
                            }
                            inputfile >> temp_epoch;
                            temp_epoch *= 86400.0;

                            for (int k = 0; k < 6; ++k)
                            {
                                inputfile >> temp_reference_angles[k];
                                temp_reference_angles[k] *= math::deg2rad;
                            }

                            for (int k = 0; k < 6; ++k)
                                inputfile >> temp_elements[k];

                            // TODO: right now we are not going to read a reference radius from 
                            // the universe file (so we don't break all universe parsing)
                            // since a Body doesn't actually ever use this, just give it the body radius for now
                            temp_J2_ref_radius = temp_radius;
                            bodies.push_back(body(temp_bodycode,
                                                  tempname,
                                                  tempshortname,
                                                  temp_SPICEnum,
                                                  temp_minimum_altitude,
                                                  temp_mu,
                                                  temp_radius,
                                                  temp_J2,
                                                  temp_J2_ref_radius,
                                                  0.0,
                                                  temp_AbsoluteMagnitude,
                                                  temp_albedo,
                                                  temp_epoch,
                                                  temp_reference_angles,
                                                  temp_elements,
                                                  this->mu,
                                                  this->central_body_SPICE_ID,
                                                  this->central_body_name,
                                                  this->central_body_radius,
                                                  this->LU,
                                                  this->LocalFrame,
                                                  options
    #ifdef SPLINE_EPHEM
                                               , this->MySplineEphemUniverse
    #endif
                                                ));

                        } while (true);
                        choice.clear();
                    }

                }
            }

            // we got to the end of the universe file and the user has specified a central body J2 value, but no reference radius....abort
            if (J2_set && !J2_ref_radius_set)
            {
                throw std::invalid_argument("The central_body_J2 parameter has been set in the Universe file: " + universefile + ", but the central_body_J2_reference_radius parameter has not.");
            }

            if (J2_ref_radius_set && !J2_set)
            {
                throw std::invalid_argument("The central_body_J2_reference_radius parameter has been set in the Universe file: " + universefile + ", but the central_body_J2 parameter has not.");
            }

            //create the flyby menu
            create_flyby_and_perturbation_menus(j, options);

            // create the CentralBody object
            this->createCentralBody(options);
            return 0;
        }

        void universe::createCentralBody(const missionoptions & options)
        {
            int temp_body_code = 0;
            std::string temp_name = this->central_body_name;
            std::string temp_short_name = "CB";
            int temp_SPICE_id = this->central_body_SPICE_ID;
            double temp_minimum_altitude = this->minimum_safe_distance - this->central_body_radius;
            double temp_mu = this->mu;
            double temp_radius = this->central_body_radius;
            double temp_J2 = this->central_body_J2;
            double temp_J2_reference_radius = this->central_body_J2_reference_radius;
            std::vector<double> temp_reference_angles = this->central_body_reference_angles;
            // TODO: I'm setting these to something trivial for now, but really, we should just
            // make these body fields private, so that central_body doesn't get them
            double temp_absolute_magnitude = -1.0; // not available from Universe file presently
            double temp_albedo = -1.0; // not available from Universe file presently
            double temp_epoch = -1.0; // not available from Universe file presently
            std::vector<double> temp_elements = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            // TODO: refactor several of these so that, for example universe.mu is universe.central_body.mu
            this->central_body = CentralBody(temp_body_code,
                                             temp_name,
                                             temp_short_name,
                                             temp_SPICE_id,
                                             temp_minimum_altitude,
                                             temp_mu,
                                             temp_radius,
                                             temp_J2,
                                             temp_J2_reference_radius,
                                             this->central_body_flattening_coefficient,
                                             temp_absolute_magnitude,
                                             temp_albedo,
                                             temp_epoch,
                                             temp_reference_angles,
                                             temp_elements,
                                             this->mu,
                                             this->central_body_SPICE_ID,
                                             this->central_body_name,
                                             this->central_body_radius,
                                             this->LU,
                                             this->LocalFrame,
                                             options
#ifdef SPLINE_EPHEM
                                           , this->MySplineEphemUniverse
#endif
                                             );

          
            this->central_body.exclusion_radius = this->central_body_radius;
            this->central_body.SPICE_frame = "Does not appear to be used presently...";

        }

        //function to locate the central body relative to the sun
        int universe::locate_central_body(const doubleType& epoch,
            doubleType* state,
            const missionoptions& options,
            const bool& need_deriv) const
        {
            if (!(boost::to_upper_copy(this->central_body_name) == "SUN"))
            {
                std::vector<size_t> timevars;

                switch (options.ephemeris_source)
                {
                case 1: //SPICE
                    double LT_dump;
                    double temp_state[6];
                    double statepert[6];

#ifdef AD_INSTRUMENTATION
                    double statedouble[12];
                    spkez_c(this->central_body_SPICE_ID, (epoch)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 10, statedouble, &LT_dump);
                    spkez_c(this->central_body_SPICE_ID, (epoch)_GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", 10, statepert, &LT_dump);
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
                    spkez_c(this->central_body_SPICE_ID, (epoch)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 10, state, &LT_dump);

                    if (need_deriv)
                    {
                        spkez_c(this->central_body_SPICE_ID, (epoch)_GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", 10, statepert, &LT_dump);
                        state[6] = (statepert[0] - state[0]) / (10.0);
                        state[7] = (statepert[1] - state[1]) / (10.0);
                        state[8] = (statepert[2] - state[2]) / (10.0);
                        state[9] = (statepert[3] - state[3]) / (10.0);
                        state[10] = (statepert[4] - state[4]) / (10.0);
                        state[11] = (statepert[5] - state[5]) / (10.0);
                    }
#endif



                    for (size_t k = 0; k < 6; ++k) {
                        temp_state[k] = (state[k]) _GETVALUE;
                    }
                    spkez_c(central_body_SPICE_ID, (epoch)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 10, temp_state, &LT_dump);

                    // REPOPULATE STATE VECTOR

                    if (need_deriv)
                    {
                        double statepert[6];
                        spkez_c(central_body_SPICE_ID, (epoch)_GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", 10, statepert, &LT_dump);
                        state[6] = (statepert[0] - state[0]) / (10.0);
                        state[7] = (statepert[1] - state[1]) / (10.0);
                        state[8] = (statepert[2] - state[2]) / (10.0);
                        state[9] = (statepert[3] - state[3]) / (10.0);
                        state[10] = (statepert[4] - state[4]) / (10.0);
                        state[11] = (statepert[5] - state[5]) / (10.0);
                    }

                    break;
#ifdef SPLINE_EPHEM
                case 2: //SplineEphem
                {
#ifdef AD_INSTRUMENTATION
                    this->MySplineEphemUniverse->getBody6StateAndDerivative(this->central_body.spice_ID, 10, epoch _GETVALUE, statedouble);
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
                        this->MySplineEphemUniverse->getBody6StateAndDerivative(this->central_body.spice_ID, 10, epoch, state);

                    else
                        this->MySplineEphemUniverse->getBody6State(this->central_body.spice_ID, 10, epoch, state);
#endif
                    break;
                }
#endif
                default:
                    std::cout << "Central body ephemeris source must be SPICE or SplineEphem" << std::endl;
                }//switch on ephemeris source
            }
            else
            {
                for (int k = 0; k < (need_deriv ? 12 : 6); ++k)
                    state[k] = 0.0;
            }

            return 0;
        }

        //function to create the flyby menu - creates a list of bodies, by SPICE ID, which are flyby capable
        void universe::create_flyby_and_perturbation_menus(const size_t& j, const missionoptions& options)
        {
            for (size_t k = 0; k < bodies.size(); ++k)
            {
                //flyby menu
                if (bodies[k].minimum_safe_flyby_altitude > 0.0)
                    flyby_menu.push_back(k); //store the position in the universe list for flyby-capable bodies

                                             //perturbation menu
                if (options.perturb_thirdbody)
                {
                    for (size_t b = 0; b < options.Journeys[j].perturbation_bodies.size(); ++b)
                    {
                        if (bodies[k].body_code == options.Journeys[j].perturbation_bodies[b])
                            perturbation_menu.push_back(k); //store the position in the universe list for perturbation bodies
                    }
                }
            }

            size_of_flyby_menu = flyby_menu.size() * 2;
        }

        //function to print the flyby menu
        void universe::print_flyby_and_perturbation_menus(std::string filename, const missionoptions& options) const
        {
            std::ofstream outputfile(filename.c_str(), std::ios::app);

            outputfile << std::endl;
            outputfile << "Flyby menu:" << std::endl;
            outputfile << "Position in flyby list | Position in Universe List | Name         | SPICE ID" << std::endl;
            outputfile << "----------------------------------------------------------------------------" << std::endl;
            for (size_t k = 0; k < size_of_flyby_menu / 2; ++k)
            {
                outputfile.width(29);
                outputfile << std::left << k + 1 << " | ";
                outputfile.width(25);
                outputfile << std::left << flyby_menu[k] + 1 << " | ";
                outputfile.width(11);
                outputfile << bodies[flyby_menu[k]].name;
                outputfile.width(3);
                outputfile << " | " << bodies[flyby_menu[k]].spice_ID << std::endl;
            }

            if (options.perturb_thirdbody)
            {
                outputfile << std::endl;
                outputfile << "Perturbation menu:" << std::endl;
                outputfile << "Position in perturbation list | Position in Universe List | Name         | SPICE ID" << std::endl;
                outputfile << "-----------------------------------------------------------------------------------" << std::endl;
                for (size_t k = 0; k < perturbation_menu.size(); ++k)
                {
                    outputfile.width(29);
                    outputfile << std::left << k + 1 << " | ";
                    outputfile.width(25);
                    outputfile << std::left << perturbation_menu[k] + 1 << " | ";
                    outputfile.width(11);
                    outputfile << bodies[perturbation_menu[k]].name;
                    outputfile.width(3);
                    outputfile << " | " << bodies[perturbation_menu[k]].spice_ID << std::endl;
                }
            }

            outputfile.close();
        }

        //function to print the universe to a file
        void universe::print_universe(std::string filename, const missionoptions& options) const
        {
            std::ofstream outputfile(filename.c_str(), std::ios::trunc);

            outputfile << "Central body name: " << this->central_body_name << std::endl;
            outputfile << "Central body SPICE ID: " << this->central_body_SPICE_ID << std::endl;
            outputfile << "mu (km^3/s^2) = " << this->mu << std::endl;
            outputfile << "LU (km) = " << this->LU << std::endl;
            outputfile << "TU (s) = " << this->TU << std::endl;
            outputfile << "alpha0 (degrees) = " << this->LocalFrame.get_alpha0() * 180.0 / math::PI << std::endl;
            outputfile << "alphadot (degrees/century) = " << this->LocalFrame.get_alphadot() * 180.0 / math::PI << std::endl;
            outputfile << "delta0 (degrees) = " << this->LocalFrame.get_delta0() * 180.0 / math::PI << std::endl;
            outputfile << "deltadot (degrees/century) = " << this->LocalFrame.get_deltadot() * 180.0 / math::PI << std::endl;
            outputfile << "Sphere of influence radius (km) = " << this->r_SOI << std::endl;
            outputfile << "Minimum safe distance (km) = " << this->minimum_safe_distance << std::endl;
            outputfile << "Orbit elements are in ICRF? (0: yes, 1: no, they are in central body frame)" << std::endl;
            outputfile << std::endl;
            outputfile << "Bodies:" << std::endl;
            outputfile << std::endl;

            outputfile.close();

            for (size_t k = 0; k < bodies.size(); ++k)
                bodies[k].print_body_to_screen(filename);

            print_flyby_and_perturbation_menus(filename, options);
        }

        void universe::get_percent_sun(math::Matrix<doubleType> spacecraft_position,
            math::Matrix<doubleType> sun_position,
            const missionoptions& options,
            doubleType& percent_sun)
        {
            //adapted from Matlab code by Steve Hughes, GSFC
            //outputs 1.0 for full shadow, 0.0 for umbra, (0, 1) for annular or penumbral eclipse

            //first of all, if the central body IS the sun then we just quit out
            if (this->central_body_SPICE_ID == 10)
            {
                percent_sun = 1.0;
                return;
            }

            //we need the radius of the Sun (hard coded) and also the radius of the current central body (known in the Universe)
            const static double SunRadius = 695990.0000;
            double BodyRadius = this->central_body_radius;

            //do some vector math
            doubleType norm_sun_position = sun_position.norm();
            doubleType norm_spacecraft_position = spacecraft_position.norm();
            math::Matrix<doubleType> sun_position_with_respect_to_spacecraft = sun_position - spacecraft_position;
            doubleType norm_sun_position_with_respect_to_spacecraft = sun_position_with_respect_to_spacecraft.norm();

            //  Calculate the apparent radius of Sun and the central body
            doubleType apparentSunRadius = EMTG::math::safe_asin(SunRadius / norm_sun_position_with_respect_to_spacecraft);
            doubleType apparentBodyRadius = EMTG::math::safe_asin(BodyRadius / norm_spacecraft_position);
            //  Calculate the apparent distance between Sun and Earth
            doubleType middle = -1.0 * (spacecraft_position / norm_spacecraft_position).dot(sun_position_with_respect_to_spacecraft / norm_sun_position_with_respect_to_spacecraft);
            doubleType apparentDistanceSunBody = EMTG::math::safe_acos(middle);


            //  Compute percent shadow based on the conditions
            if (apparentSunRadius + apparentBodyRadius <= apparentDistanceSunBody)
                percent_sun = 1.0;  // In full Sun
            else if (apparentDistanceSunBody <= apparentBodyRadius - apparentSunRadius)
                percent_sun = 0.0;  //  In Umbra
            else if (fabs(apparentSunRadius - apparentBodyRadius) < apparentDistanceSunBody
                && apparentDistanceSunBody < apparentSunRadius + apparentBodyRadius)
            {
                // Penumbra, not an annular eclipse
                doubleType c1 = (apparentDistanceSunBody * apparentDistanceSunBody + apparentSunRadius * apparentSunRadius - apparentBodyRadius * apparentBodyRadius)
                    / 2.0 / apparentDistanceSunBody;
                doubleType c2 = sqrt(apparentSunRadius * apparentSunRadius - c1 * c1);
                doubleType A = apparentSunRadius * apparentSunRadius * EMTG::math::safe_acos(c1 / apparentSunRadius)
                    + apparentBodyRadius * apparentBodyRadius * EMTG::math::safe_acos((apparentDistanceSunBody - c1) / apparentBodyRadius)
                    - apparentDistanceSunBody * c2;
                percent_sun = 1.0 - A / math::PI / (apparentSunRadius * apparentSunRadius);
            }
            else
                percent_sun = 1.0 - apparentBodyRadius * apparentBodyRadius / (apparentSunRadius * apparentSunRadius); //  Annular

        }

		// method to get central body reference angles
		void universe::get_central_body_reference_angles(std::vector<double>& central_body_reference_angles_out)
		{
			central_body_reference_angles_out = this->central_body_reference_angles;
		}

    }//close namespace Astrodynamics
}//close namespace EMTG