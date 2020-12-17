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

//master testbed
//Jacob Englander 6-20-2017

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>

#include "missionoptions.h"
#include "mission.h"
#include "chinchilla.h"
#include "EMTG_enums.h"
#include "EMTG_Matrix.h"

#include "LaunchVehicle.h"
#include "Spacecraft.h"
#include "LaunchVehicleOptionsFactory.h"
#include "SpacecraftOptionsFactory.h"

#include "universe.h"
#include "body.h"
#include "atmosphere.h"
#include "ExponentialAtmosphere.h"

#include "file_utilities.h"

#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include "boost/date_time.hpp"
#include "boost/date_time/local_time/local_date_time.hpp"
#include "boost/lexical_cast.hpp"

#include "SpiceUsr.h"


#ifdef SPLINE_EPHEM
#include "SplineEphem_universe.h"
#endif

#include "mission_testbed.h"
#include "MTMtestbed.h"
#include "propagatorTestbed.h"
#include "eventTestbed.h"
#include "mathTestbed.h"
#include "frame_testbed.h"
#include "MIRAGE_verification.h"
#include "StateRepresentation_testbed.h"

int main(int argc, char* argv[])
{
    try
    {
        //startup stuff
        std::cout << "program starting" << std::endl;

        //parse the options file
        std::string options_file_name;
        if (argc == 1)
            options_file_name = "default.emtgopt";
        else if (argc == 2)
            options_file_name.assign(argv[1]);

        std::cout << options_file_name << std::endl;

        EMTG::missionoptions options;
        try
        {
            options.parse_mission(options_file_name);
        }
        catch (std::exception error)
        {
            std::cout << error.what() << std::endl;
            std::cout << "Aborting program run." << std::endl;
            std::cin.ignore();
            return 0;
        }

        //configure the LaunchVehicleOptions and SpacecraftOptions objects


        //configure the LaunchVehicleOptions and SpacecraftOptions objects
        EMTG::HardwareModels::LaunchVehicleOptions myLaunchVehicleOptions = EMTG::HardwareModels::CreateLaunchVehicleOptions(options);
        EMTG::HardwareModels::SpacecraftOptions mySpacecraftOptions = EMTG::HardwareModels::CreateSpacecraftOptions(options);

        //load all ephemeris data if using SPICE
        std::vector<fs::path> SPICE_files_initial;
        std::vector<fs::path> SPICE_files_not_required;
        std::vector<fs::path> SPICE_files_required;
        std::vector<int> SPICE_bodies_required;
        std::string filestring;
        if (options.ephemeris_source >= 1)
        {
            //load all BSP files
            EMTG::file_utilities::get_all_files_with_extension(fs::path(options.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files_initial);

            for (size_t k = 0; k < SPICE_files_initial.size(); ++k)
            {
                filestring = options.universe_folder + "/ephemeris_files/" + SPICE_files_initial[k].string();
                furnsh_c(filestring.c_str());
                std::cout << filestring << std::endl;
            }

            //disable quit-on-SPICE-error so that we can see what happens if the leap second and/or frame kernels don't load properly
            erract_c((SpiceChar*)"SET", 100, (SpiceChar*)"RETURN");

            //SPICE reference frame kernel
            std::string leapsecondstring = options.universe_folder + "/ephemeris_files/" + options.SPICE_leap_seconds_kernel;
            std::string referenceframestring = options.universe_folder + "/ephemeris_files/" + options.SPICE_reference_frame_kernel;
            furnsh_c(leapsecondstring.c_str());
            furnsh_c(referenceframestring.c_str());

            //disable SPICE error printing. This is because we can, and will often, go off the edge of an ephemeris file.
            errprt_c((SpiceChar*)"SET", 100, (SpiceChar*)"NONE");

            SPICE_files_required = SPICE_files_initial;

            std::cout << "Completed loading SPICE kernels." << std::endl;
        }

        //RNG
        std::mt19937 RNG(0); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> UniformDouble(0.0, 1.0);

        //create a picture of a chinchilla
        draw_chinchilla();

        //if SplineEphem is enabled, create an empty SplineEphem universe
#ifdef SPLINE_EPHEM
        std::vector< std::tuple<int, int, int, double> > SplineUniverse_keyList;

        SplineEphem::universe SplineUniverse(SplineUniverse_keyList);
#endif

        //create a vector of universes for each journey
        std::vector<EMTG::Astrodynamics::universe > TheUniverse;
        options.TU = 0;
        for (int j = 0; j < options.number_of_journeys; ++j)
        {
#ifdef SPLINE_EPHEM
            TheUniverse.push_back(EMTG::Astrodynamics::universe(j, options.universe_folder + "//" + options.Journeys[j].journey_central_body + ".emtg_universe", options, &SplineUniverse));
#else
            TheUniverse.push_back(EMTG::Astrodynamics::universe<GSAD::adouble>(j, options.universe_folder + "//" + options.Journeys[j].journey_central_body + ".emtg_universe", options));
#endif
            std::stringstream universenamestream;

            universenamestream << options.Journeys[j].journey_central_body + "_Journey_" << j << ".universe_output";

            if (TheUniverse[j].TU > options.TU)
                options.TU = TheUniverse[j].TU;
        }

        for (int j = 0; j < options.number_of_journeys; ++j)
        {
            if (j > 0)
            {
                TheUniverse[j - 1].set_nextUniverse(TheUniverse[j]);
            }
        }

		//now that we have a Universe vector, we can set the atmosphere object for each universe
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			try
			{
				if (options.Journeys[j].perturb_drag ||
                    (options.Journeys[j].phase_type == EMTG::PhaseType::ProbeEntryPhase && (options.Journeys[j].perturb_drag_probe_AEI_to_end || options.Journeys[j].perturb_drag_probe_separation_to_AEI))
                    )
				{
					// need to choose the correct kind of atmosphere to create based on journey option
					if (options.Journeys[j].AtmosphericDensityModelKey == "Exponential")
					{
                        TheUniverse[j].TheAtmosphere = std::make_shared<EMTG::Astrodynamics::ExponentialAtmosphere>(j, options.Journeys[j].AtmosphericDensityModelDataFile, options);
					}
					else
					{
						std::cout << "Impermissible choice for AtmosphericDensityModelKey for Journey " << j << std::endl;
						throw std::exception();
					}
				}
				else
				{
					// we are not calculating drag, so just build a placeholder atmosphere that does nothing
					//TheUniverse[j].TheAtmosphere = new EMTG::Astrodynamics::ExponentialAtmosphere();
				}
			}
			catch (std::exception &myError)
			{
				std::cout << "Failure with configuring TheAtmosphere." << std::endl;
				throw;
			}
		}

        //now that we have a Universe vector, we can use it to populate the SplineEphem::universe
               //add every body that will we used in the mission to the SplineUniverse
#ifdef SPLINE_EPHEM
        SplineUniverse_keyList.clear();
        try
        {
            //double earliest_possible_epoch = options.launch_window_open_date + options.Journeys.front().wait_time_bounds[0];
            //double latest_possible_epoch = options.latestPossibleEpoch * 86400.0;

            size_t number_of_journeys_to_spline = std::min(options.number_of_journeys, options.stop_after_journey + 1);
            for (size_t j = 0; j < number_of_journeys_to_spline; ++j)
            {
                std::vector<int> body_index_array;

                //first boundary point
                if (options.Journeys[j].departure_class != EMTG::BoundaryClass::FreePoint
                    && options.Journeys[j].departure_elements_frame != EMTG::ReferenceFrame::ObjectReferenced)
                {
                    if (options.Journeys[j].destination_list[0] > 0)
                        body_index_array.push_back(options.Journeys[j].destination_list[0] - 1);
                }

                //last boundary point
                if (options.Journeys[j].arrival_class != EMTG::BoundaryClass::FreePoint
                    && options.Journeys[j].arrival_elements_frame != EMTG::ReferenceFrame::ObjectReferenced)
                {
                    if (options.Journeys[j].destination_list[1] > 0)
                        body_index_array.push_back(options.Journeys[j].destination_list[1] - 1);
                }

                //sequence
                for (int body : options.Journeys[j].sequence)
                    if (body > 0)
                        body_index_array.push_back(body - 1);

                //perturbation list
                if (options.perturb_thirdbody)
                {
                    for (size_t b = 0; b < TheUniverse[j].perturbation_menu.size(); ++b)
                        body_index_array.push_back(TheUniverse[j].perturbation_menu[b]);
                }

                //distance constraint list
                for (std::string& constraint : options.Journeys[j].PhaseDistanceConstraintDefinitions)
                {
                    std::vector<std::string> ConstraintDefinitionCell;
                    boost::split(ConstraintDefinitionCell,
                        constraint,
                        boost::is_any_of("_"),
                        boost::token_compress_on);

                    if (boost::to_lower_copy(ConstraintDefinitionCell[1]) != "cb")
                    {
                        int bodyIndex = std::stoi(ConstraintDefinitionCell[1]) - 1;

                        body_index_array.push_back(bodyIndex);
                    }
                }

                for (size_t b = 0; b < body_index_array.size(); ++b)
                {
                    //do we already have this body?
                    bool body_in_keylist = false;
                    for (size_t k = 0; k < SplineUniverse_keyList.size(); ++k)
                    {
                        if (std::get<0>(SplineUniverse_keyList[k]) == TheUniverse[j].bodies[body_index_array[b]].spice_ID
                            && std::get<1>(SplineUniverse_keyList[k]) == TheUniverse[j].central_body_SPICE_ID)
                        {
                            body_in_keylist = true;
                            break;
                        }
                    }

                    if (!body_in_keylist && body_index_array[b] >= 0)
                    {
                        SplineUniverse_keyList.push_back(std::make_tuple(
                            TheUniverse[j].bodies[body_index_array[b]].spice_ID,
                            TheUniverse[j].central_body_SPICE_ID,
                            options.SplineEphem_points_per_period,
                            TheUniverse[j].mu));
                    }
                }//end loop over bodies in the universe

                //is this universe's central body the sun? If not, let's add this body with respect to the sun. Let's add extra ephemeris points, too.
                if (!(TheUniverse[j].central_body_SPICE_ID == 10))
                {
                    bool body_in_keylist = false;
                    for (size_t k = 0; k < SplineUniverse_keyList.size(); ++k)
                    {
                        if (std::get<0>(SplineUniverse_keyList[k]) == TheUniverse[j].central_body_SPICE_ID
                            && std::get<1>(SplineUniverse_keyList[k]) == 10)
                        {
                            body_in_keylist = true;
                            break;
                        }
                    }

                    if (!body_in_keylist)
                    {
                        SplineUniverse_keyList.push_back(std::make_tuple(
                            TheUniverse[j].central_body_SPICE_ID,
                            10,
                            options.SplineEphem_non_central_body_sun_points_per_period,
                            1.32712440018e+11));
                    }
                }

                ////do we need to update the earliest or latest possible epoch?
                //if (options.Journeys[j].arrival_class == EMTG::BoundaryClass::FreePoint)
                //{
                //    earliest_possible_epoch = options.Journeys[j].arrival_elements_reference_epoch < earliest_possible_epoch ? options.Journeys[j].arrival_elements_reference_epoch : earliest_possible_epoch;
                //    latest_possible_epoch = options.Journeys[j].arrival_elements_reference_epoch > latest_possible_epoch ? options.Journeys[j].arrival_elements_reference_epoch : latest_possible_epoch;
                //}
                //if (options.Journeys[j].departure_class == EMTG::BoundaryClass::FreePoint && j == 0)
                //{
                //    earliest_possible_epoch = options.Journeys[j].departure_elements_reference_epoch < earliest_possible_epoch ? options.Journeys[j].departure_elements_reference_epoch : earliest_possible_epoch;
                //    latest_possible_epoch = options.Journeys[j].departure_elements_reference_epoch > latest_possible_epoch ? options.Journeys[j].departure_elements_reference_epoch : latest_possible_epoch;
                //}
            }

            double earliestPossibleEpoch = options.earliestPossibleEpoch * 86400.0;
            double latestPossibleEpoch = options.latestPossibleEpoch * 86400.0;

            if (options.SplineEphem_truncate_ephemeris_at_maximum_mission_epoch
                && latestPossibleEpoch < (options.launch_window_open_date + options.Journeys.front().wait_time_bounds[1] + options.total_flight_time_bounds[1]))
                latestPossibleEpoch = options.launch_window_open_date + options.Journeys.front().wait_time_bounds[1] + options.total_flight_time_bounds[1] * 86400.0;
            /*if (earliest_possible_epoch > options.earliestPossibleEpoch * 86400.0)
                earliest_possible_epoch = options.earliestPossibleEpoch * 86400.0;*/
            SplineUniverse.reinitialize(SplineUniverse_keyList,
                earliestPossibleEpoch - 10.0 * 86400.0,
                latestPossibleEpoch + 10.0 * 86400.0);
        }
        catch (std::exception &myError)
        {
            std::cout << "Failure while configuring SplineEphem." << std::endl;
            std::cout << myError.what() << std::endl;
            std::cout << "Submit this error message to the EMTG development team, along with your .emtgopt, .emtg_universe file(s), your hardware model files, any relevant ephemeris files, and which branch you are using. This information will allow us to properly help you." << std::endl;
#ifndef BACKGROUND_MODE //macro overrides if statement
            std::cout << "Press enter to close window." << std::endl;
            std::cin.ignore();
#endif
            throw;
        }
#endif

        EMTG::HardwareModels::Spacecraft mySpacecraft(mySpacecraftOptions);
        EMTG::HardwareModels::LaunchVehicle myLaunchVehicle(myLaunchVehicleOptions);

        //*****************************************************************

                //assemble the mission
        options.description.clear();

        for (size_t j = 0; j < options.number_of_journeys; ++j)
        {
            std::vector<int> phase_targets = options.Journeys[j].sequence;
            options.Journeys[j].sequence.clear();

            options.Journeys[j].sequence.push_back(options.Journeys[j].destination_list[0]);

            if (j > 0) //if not the first journey, insert an underscore
                options.description.append("_");
            options.description.append(TheUniverse[j].central_body_name + "(");

            if (options.Journeys[j].departure_class == EMTG::BoundaryClass::EphemerisPegged)
            {
                options.description.append(TheUniverse[j].bodies[options.Journeys[j].sequence[0] - 1].short_name);
            }
            else if (options.Journeys[j].departure_class == EMTG::BoundaryClass::EphemerisReferenced)
            {
                if (options.Journeys[j].sequence[0] == -1) //SOI
                {
                    options.description.append("s");
                }
                else if (options.Journeys[j].sequence[0] == 0) //central body
                {
                    options.description.append("c");
                }
                else //body in the universe
                {
                    options.description.append(TheUniverse[j].bodies[options.Journeys[j].sequence[0] - 1].short_name);
                }
            }
            else if (options.Journeys[j].departure_class == EMTG::BoundaryClass::Periapse)
            {
                options.description.append("c"); //always the central body
            }
            else //free point
            {
                options.description.append("f");
            }

            //first, how many phases are there in the journey?
            for (size_t p = 0; p < options.Journeys[j].number_of_phases - 1; ++p)
            {
                size_t bodyIndex = phase_targets[p];
                if (bodyIndex > 0 && bodyIndex < (TheUniverse[j].size_of_flyby_menu / 2) + 1) //this is a legitimate flyby
                {
                    if (bodyIndex - 1 > TheUniverse[j].flyby_menu.size())
                    {
                        throw std::invalid_argument("ERROR: Journey " + std::to_string(j) + " phase " + std::to_string(p) + " body index " + std::to_string(bodyIndex)
                            + " exceeds size of flyby menu. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }

                    //append the flyby
                    options.Journeys[j].sequence.push_back(phase_targets[p]);

                    //update the mission description
                    options.description.append(TheUniverse[j].bodies[TheUniverse[j].flyby_menu[options.Journeys[j].sequence.back() - 1]].short_name);
                }
            }

            options.Journeys[j].sequence.push_back(options.Journeys[j].destination_list[1]);


            if (options.Journeys[j].arrival_class == EMTG::BoundaryClass::EphemerisPegged)
            {
                options.description.append(TheUniverse[j].bodies[options.Journeys[j].sequence.back() - 1].short_name);
            }
            else if (options.Journeys[j].arrival_class == EMTG::BoundaryClass::EphemerisReferenced)
            {
                if (options.Journeys[j].sequence.back() == -1) //SOI
                {
                    options.description.append("s");
                }
                else if (options.Journeys[j].sequence.back() == 0) //central body
                {
                    options.description.append("c");
                }
                else //body in the universe
                {
                    options.description.append(TheUniverse[j].bodies[options.Journeys[j].sequence.back() - 1].short_name);
                }
            }
            else if (options.Journeys[j].arrival_class == EMTG::BoundaryClass::Periapse)
            {
                options.description.append("c"); //always the central body
            }
            else //free point
            {
                options.description.append("f");
            }
            options.description.append(")");

            options.Journeys[j].number_of_phases = options.Journeys[j].sequence.size() - 1;

        }
        //end startup stuff


        //***************************************************************************actually useful things

        //math testbed
        //mathTestbed(options, TheUniverse, mySpacecraft, myLaunchVehicle, RNG, UniformDouble);

        //event testbed
        //eventTestbed(options, TheUniverse, mySpacecraft, myLaunchVehicle, RNG, UniformDouble);

        //propagator testbed
        //propagatorTestbed(options, TheUniverse, mySpacecraft, RNG, UniformDouble);

        //MTM testbed
        //MTMtestbed(options, TheUniverse, mySpacecraft, RNG, UniformDouble);

        //MGAnDSMs subphases
        //MGAnDSMs_subphase_testbed(options, TheUniverse, mySpacecraft, RNG, UniformDouble);

        //phase testbed
        //phaseTestbed(options, TheUniverse, mySpacecraft, myLaunchVehicle, RNG, UniformDouble);

        //state representation testbed
        //StateRepresentation_testbed(options, TheUniverse, RNG, UniformDouble);

        //mission testbed
        missionTestbed(options, TheUniverse, mySpacecraft, myLaunchVehicle, RNG, UniformDouble);

        //frame testbed
        //frameTestbed(options, TheUniverse, RNG, UniformDouble);

        // MIRAGE verification of acceleration model accuracy
        //MIRAGE_verification(options, TheUniverse, mySpacecraft, myLaunchVehicle);
    }
    catch (std::exception exception)
    {
        std::cout << "EMTG testbed failed with error:" << std::endl;
        std::cout << exception.what() << std::endl;
        std::cout << "Submit this error message to the EMTG development team, along with your .emtgopt, .emtg_universe file(s), your hardware model files, any relevant ephemeris files, and which branch you are using. This information will allow us to properly help you." << std::endl;
        std::cin.ignore();
    }
}