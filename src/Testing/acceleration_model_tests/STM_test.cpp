// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2014 - 2016 United States Government as represented by the
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

// FBLT_derivatives_testbed.cpp : Defines the entry point for the console application.
//

#include <iomanip>
#include <string>
#include <sstream>
#include <chrono>

#include "boost/random/uniform_real.hpp"
#include "boost/random/mersenne_twister.hpp"


#include "LaunchVehicle.h"
#include "Spacecraft.h"

#include "LaunchVehicleOptionsFactory.h"
#include "SpacecraftOptionsFactory.h"
#include "doubleType.h"
#include "EMTG_Matrix.h"
#include "ExplicitRungeKutta.h"
#include "file_utilities.h"
#include "IntegratedAdaptiveStepPropagator.h"
#include "IntegratedFixedStepPropagator.h"
#include "mission.h"
#include "missionoptions.h"
#include "PropagatorFactory.h"
#include "SpacecraftAccelerationModel.h"
#include "SpiceUsr.h" 
#include "TimeDomainSpacecraftEOM.h"
#include "SundmanSpacecraftEOM.h"
#include "universe.h"
#include "atmosphere.h"
#include "ExponentialAtmosphere.h"


// needed for GSAD 4B
//size_t GSAD::adouble::reserve_size = 10;

// needed for GSAD A5
//GSAD::adouble GSAD::adouble::temp;
//std::vector<size_t>::size_type GSAD::adouble::point = 0;
//GSAD::adouble::derivative_t GSAD::adouble::pair_temp = std::make_pair(0, 0);
//size_t GSAD::adouble::reserve_size = 10;

int main(int argc, char* argv[])
{

    boost::mt19937 RNG;
    boost::uniform_real<> DoubleDistribution;
    DoubleDistribution = boost::uniform_real<>(0.0, 0.55);
    RNG.seed(time(NULL));

    //parse the options file
    std::string options_file_name;
    if (argc == 1)
        options_file_name = "default.emtgopt";
    else if (argc == 2)
        options_file_name.assign(argv[1]);

    std::cout << options_file_name << std::endl;

    EMTG::missionoptions options;
    options.parse_mission(options_file_name);

    //configure the LaunchVehicleOptions and SpacecraftOptions objects
    EMTG::HardwareModels::LaunchVehicleOptions myLaunchVehicleOptions = EMTG::HardwareModels::CreateLaunchVehicleOptions(options);
    EMTG::HardwareModels::SpacecraftOptions mySpacecraftOptions = EMTG::HardwareModels::CreateSpacecraftOptions(options);

    EMTG::HardwareModels::LaunchVehicle myLaunchVehicle(myLaunchVehicleOptions);
    EMTG::HardwareModels::Spacecraft mySpacecraft(mySpacecraftOptions);

    //Load appropriate SPICE kernels
    //load DE430, leap seconds kernel, frame kernel

    //load all ephemeris data if using SPICE
    std::vector<::boost::filesystem::path> SPICE_files_initial;
    std::vector<::boost::filesystem::path> SPICE_files_not_required;
    std::vector<::boost::filesystem::path> SPICE_files_required;
    std::vector<int> SPICE_bodies_required;
    std::string filestring;
    if (options.ephemeris_source >= 1)
    {
        //load all BSP files
        EMTG::file_utilities::get_all_files_with_extension(::boost::filesystem::path(options.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files_initial);

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

    //create a vector of Universes, one for each journey

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
        TheUniverse.push_back(EMTG::Astrodynamics::universe(j, options.universe_folder + "//" + options.Journeys[j].journey_central_body + ".emtg_universe", options));
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
            if (options.Journeys[j].destination_list[0] > 0)
                body_index_array.push_back(options.Journeys[j].destination_list[0] - 1);

            //last boundary point
            if (options.Journeys[j].destination_list[1] > 0)
                body_index_array.push_back(options.Journeys[j].destination_list[1] - 1);

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

	//now that we have a Universe vector, we can set the atmosphere object for each universe
	for (int j = 0; j < options.number_of_journeys; ++j)
	{
		try
		{
			if (options.Journeys[j].perturb_drag == 1)
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
        switch (options.Journeys[j].sequence[0])
        {
        case -1: //begin at SOI
        {
            options.description.append("s");
            break;
        }
        case 0: //begin at central body
        {
            options.description.append("c");
            break;
        }
        default:
            options.description.append(TheUniverse[j].bodies[options.Journeys[j].sequence[0] - 1].short_name);
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
                        + " exceeds size of flyby menu.");
                }

                //append the flyby
                options.Journeys[j].sequence.push_back(phase_targets[p]);

                //update the mission description
                options.description.append(TheUniverse[j].bodies[TheUniverse[j].flyby_menu[options.Journeys[j].sequence.back() - 1]].short_name);
            }
        }

        options.Journeys[j].sequence.push_back(options.Journeys[j].destination_list[1]);


        switch (options.Journeys[j].sequence.back())
        {
        case -1: //begin at SOI
        {
            options.description.append("s");
            break;
        }
        case 0: //begin at central body
        {
            options.description.append("c");
            break;
        }
        default:
            options.description.append(TheUniverse[j].bodies[options.Journeys[j].sequence.back() - 1].short_name);
        }
        options.description.append(")");

        options.Journeys[j].number_of_phases = options.Journeys[j].sequence.size() - 1;

    }

    EMTG::Mission myMission(options,
                            TheUniverse,
                            myLaunchVehicle,
                            mySpacecraft);

    
    //TheUniverse[0].set_X_scale_factors(&myMission.X_scale_factors);

    size_t state_vector_size = 10;
    size_t STM_size = 14;
    std::vector<std::string> Xdescriptions;
    EMTG::Astrodynamics::SpacecraftAccelerationModel test_acceleration_model(&myMission.options, 
                                                                             &myMission.options.Journeys[0], 
                                                                             &TheUniverse[0],
                                                                             &Xdescriptions,
                                                                             &mySpacecraft, 
                                                                             STM_size);
    
    doubleType test_launch_epoch = myMission.options.launch_window_open_date; // seconds past MJD J2000
    doubleType TOFprevious, TOF;
    if (options.mission_type == EMTG::PhaseType::SundmanCoastPhase)
    {
        TOFprevious = 0.0 * 86400.0;
        //TOF = 22750.0*EMTG::math::PI;
        TOF = 1000.0*EMTG::math::PI;
    }
    else
    {
        TOFprevious = 0.0 * 86400.0;
        TOF = 5.07152939736843 * 86400.0; // what it was before
		//TOF = 58.0; // for drag
        //TOF = 1.0 * 86400.0;

    }
	//test_launch_epoch -= 86400.*365.25;
    doubleType test_epoch;

#ifdef AD_INSTRUMENTATION
    size_t GSADindex = 0;
    test_launch_epoch.setDerivative(GSADindex++, 1.0);
    TOFprevious.setDerivative(GSADindex++, 1.0);
    TOF.setDerivative(GSADindex++, 1.0);
#endif

    test_epoch = test_launch_epoch + TOFprevious;

    //initialize state vector
    double LT_dump;
    std::vector<double> x(8);
    double shit = test_epoch _GETVALUE - (51544.5 * 86400.0);
    //spkez_c(399, shit, "J2000", "NONE", 10, x.data(), &LT_dump);
    spkez_c(504, test_epoch _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 599, x.data(), &LT_dump);
    //spkez_c(301, test_epoch _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 399, x.data(), &LT_dump);
    //spkez_c(499, test_epoch _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 4, x.data(), &LT_dump);
    //spkez_c(602, test_epoch _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 699, x.data(), &LT_dump);
    //x[0] = -1169636.86875674;
    //x[1] = -781158.68641361;
    //x[2] = -337861.83290009;
    //x[3] = 4.67341199;
    //x[4] = 3.04126260;
    //x[5] = 1.31500098;
    x[6] = myMission.options.maximum_mass;
    x[7] = test_epoch _GETVALUE;
    //x[0] = 4571.0;
    //x[1] = 3403.3;
    //x[2] = 3003.3;
    //x[3] = -0.1;
    //x[4] = 7.4;
    //x[5] = -0.5;
	//x[0] = 6571.0;
	//x[1] = 14.3;
	//x[2] = 100.3;
	//x[3] = -0.1;
	//x[4] = 7.4;
	//x[5] = -0.5;
    //x[6] = myMission.options.maximum_mass;

    //create perturbed initial state for TOF derivative use
    std::vector <double> x_pert(7, 0.0);
    std::vector <double> dx_dTOFprevious(7, 0.0);
    std::vector <double> dx_dTOFcurrent(7, 0.0);
    //spkez_c(399, shit + 10.0, "J2000", "NONE", 10, x_pert.data(), &LT_dump);
    spkez_c(504, test_epoch _GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", 599, x_pert.data(), &LT_dump);
    //spkez_c(301, test_epoch _GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", 399, x_pert.data(), &LT_dump);
    //spkez_c(602, test_epoch _GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", 699, x_pert.data(), &LT_dump);
    //dx_dTOFprevious[0] = (x_pert[0] - x[0]) / 10.0;
    //dx_dTOFprevious[1] = (x_pert[1] - x[1]) / 10.0;
    //dx_dTOFprevious[2] = (x_pert[2] - x[2]) / 10.0;
    //dx_dTOFprevious[3] = (x_pert[3] - x[3]) / 10.0;
    //dx_dTOFprevious[4] = (x_pert[4] - x[4]) / 10.0;
    //dx_dTOFprevious[5] = (x_pert[5] - x[5]) / 10.0;
    //dx_dTOFprevious[6] = 0.0;

    EMTG::math::Matrix<doubleType> test_state(state_vector_size, 1, 0.0);
    for (size_t k = 0; k < 7; ++k)
    {
        test_state(k) = x[k];
#ifdef AD_INSTRUMENTATION
        test_state(k).setDerivative(GSADindex++, 1.0);
#endif
        //test_state(k).setDerivative(0, dx_dTOFprevious[k]);
        //test_state(k).setDerivative(1, dx_dTOFprevious[k]);
        //test_state(k).setDerivative(2, dx_dTOFcurrent[k]);
    }

    // current epoch
    test_state(7) = test_epoch;
#ifdef AD_INSTRUMENTATION
    test_state(7).setDerivative(GSADindex++, 1.0);
#endif

// virtual chemical thruster propellant tank
    test_state(8) = 0.0;
#ifdef AD_INSTRUMENTATION
    test_state(8).setDerivative(GSADindex++, 1.0);
#endif

    // virtual electric thruster propellant tank 
    test_state(9) = 0.0;
#ifdef AD_INSTRUMENTATION
    test_state(9).setDerivative(GSADindex++, 1.0);
#endif

    // control decision variables
    EMTG::math::Matrix<doubleType> control(4, 1, 0.0);
    for (size_t i = 0; i < 3; ++i)
    {
        control(i) = DoubleDistribution(RNG);
        //control(i).setValue(0.0);
#ifdef AD_INSTRUMENTATION
        control(i).setDerivative(GSADindex++, 1.0);
#endif
    }

    test_acceleration_model.setDutyCycle(myMission.options.engine_duty_cycle);
    test_acceleration_model.setLaunchEpoch(test_launch_epoch);
    test_acceleration_model.setEpoch(test_epoch);
    //test_acceleration_model.setEpoch(test_state(7));
    //test_acceleration_model.computeAcceleration(test_state, control, true);
    //EMTG::math::Matrix<doubleType> acceleration = test_acceleration_model.getAcceleration();

    // create the spacecraft equations of motion
    EMTG::Integration::Integrand * thing; // create the integrand pointer
    EMTG::Astrodynamics::SundmanSpacecraftEOM SundmanEOM = EMTG::Astrodynamics::SundmanSpacecraftEOM();
    EMTG::Astrodynamics::TimeDomainSpacecraftEOM TimeEOM = EMTG::Astrodynamics::TimeDomainSpacecraftEOM();
    if (options.mission_type == EMTG::PhaseType::SundmanCoastPhase)
    {
        SundmanEOM.setSpacecraftAccelerationModel(&test_acceleration_model);
        thing = &SundmanEOM;
    }
    else
    {
        TimeEOM.setSpacecraftAccelerationModel(&test_acceleration_model);
        thing = &TimeEOM;
    }

    // create the integration scheme
    // {x, y, z, vx, vy, vz, m, epoch, tank1, tank2, ux, uy, uz}
    EMTG::Integration::ExplicitRungeKutta explicit_rk8(thing, EMTG::rkdp87, state_vector_size, STM_size);

    // create the propagator
    EMTG::math::Matrix <doubleType> state_left(state_vector_size, 1, 0.0);
    EMTG::math::Matrix <double> dstate_leftdTOF(10, 2, 0.0); // initial state partials w.r.t. previous and current TOF
    EMTG::math::Matrix <doubleType> state_right(state_vector_size, 1, 0.0);
    EMTG::math::Matrix <double> dstate_rightdTOF(10, 2, 0.0); // final state partials w.r.t. previous and current TOF
    EMTG::math::Matrix <double> STM(STM_size, STM_size, 0.0);
    double BoundaryTarget_dStepSizedPropVar = 1.0 / options.num_timesteps;

    // insert true state
    state_left = test_state;

    // seed derivatives for spacecraft state with respect to current and previous phase flight times
    for (size_t i = 0; i < 3; i++)
    {
        //dstate_leftdTOF(i, 0) = dx_dTOFprevious[i];
        //dstate_leftdTOF(i + 3, 0) = dx_dTOFprevious[i + 3];
        dstate_leftdTOF(i, 0) = 0.0;
        dstate_leftdTOF(i + 3, 0) = 0.0;
        dstate_leftdTOF(i, 1) = 0.0;
        dstate_leftdTOF(i + 3, 1) = 0.0;
    }


    EMTG::Astrodynamics::PropagatorBase * integratedPropagator = EMTG::Astrodynamics::CreatePropagator(&options,
        &TheUniverse[0],
        state_vector_size,
        STM_size,
        state_left,
        state_right,
        STM,
        dstate_leftdTOF,
        thing,
        &explicit_rk8,
        &BoundaryTarget_dStepSizedPropVar,
        options.Journeys[0].override_integration_step_size
        ? options.Journeys[0].integration_step_size
        : options.integration_time_step_size);

    integratedPropagator->setCurrentEpoch(test_state(7));
    integratedPropagator->setIndexOfEpochInStateVec(7); // where is epoch located in the state vector?
    if (options.mission_type == EMTG::PhaseType::SundmanCoastPhase)
    {
        integratedPropagator->setCurrentIndependentVariable(0.0);
    }
    else
    {
        integratedPropagator->setCurrentIndependentVariable(test_state(7));
    }

    //options.X_scale_factors[0] = 1.0;
    //options.X_scale_factors[1] = 1.0;
    //options.X_scale_factors[2] = 1.0;
    //integratedFSProp.propagate(TOF / options.num_timesteps, control, true);
    //integratedFSProp.propagate(TOF / options.num_timesteps, true);
    // Start a timer
    std::cout << std::setprecision(16);
    clock_t begin = clock();
    //integratedPropagator->propagate(TOF / options.num_timesteps, control, true);
    integratedPropagator->propagate(TOF / options.num_timesteps, true);
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Propagation time: " << std::setprecision(16) << elapsed_secs << " secs\n" << std::endl;

    dstate_rightdTOF = dstate_leftdTOF;    

    std::cout << "Numerically integrated STM:" << std::endl;
    for (size_t i = 0; i < STM_size; ++i)
    {
        for (size_t j = 0; j < STM_size; ++j)
        {
            std::cout << STM(i, j) << "     ";
        }
        std::cout << std::endl << std::endl;
    }

    std::cout << std::endl << std::endl;

#ifdef AD_INSTRUMENTATION
    std::cout << "AD computed STM:" << std::endl;
    EMTG::math::Matrix <double> AD_STM(STM_size, STM_size, 0.0);
    size_t derivIndex = 3;
    for (size_t i = 0; i < state_vector_size; ++i)
    {
        for (size_t derivIndex = 3; derivIndex < STM_size + 3; ++derivIndex)
        {
            AD_STM(i, derivIndex - 3) = state_right(i).getDerivative(derivIndex);
            std::cout << AD_STM(i, derivIndex - 3) << "     ";
        }
        // propagation variable derivative for this row
        AD_STM(i, 13) = state_right(i).getDerivative(2);
        std::cout << AD_STM(i, 13) << "     ";
        std::cout << std::endl << std::endl;
    }

    std::cout << std::endl << std::endl;

    std::cout << "Absolute STM error:" << std::endl;
    for (size_t i = 0; i < STM_size; ++i)
    {
        for (size_t j = 0; j < STM_size; ++j)
        {
            std::cout << AD_STM(i, j) - STM(i, j) << "     ";
        }
        std::cout << std::endl << std::endl;
    }

    std::cout << std::endl << std::endl;

    std::cout << "Relative STM error:" << std::endl;
    for (size_t i = 0; i < STM_size; ++i)
    {
        for (size_t j = 0; j < STM_size; ++j)
        {
            if (AD_STM(i, j) < 1.0e-20 && (AD_STM(i, j) - STM(i, j)) < 1.0e-20)
            {
                std::cout << 0.0 << "     ";
            }
            else
            {
                std::cout << (AD_STM(i, j) - STM(i, j)) / AD_STM(i, j) << "     ";
            }
        }
        std::cout << std::endl << std::endl;
    }

    std::cout << std::endl << std::endl;

    std::cout << "Independent variable derivatives" << std::endl;
    std::cout << "Format: analytical, GSAD, absolute error, relative error" << std::endl;
    std::cout << "dxdIndVarcurrent: "      << STM(0, 13) << " " << state_right(0).getDerivative(2) << " " << STM(0, 13) - state_right(0).getDerivative(2) << " " << (STM(0, 13) - state_right(0).getDerivative(2)) / state_right(0).getDerivative(2) << std::endl;
    std::cout << "dxdIndVarprevious: "     << STM(0, 7)  << " " << state_right(0).getDerivative(1) << " " << STM(0, 7)  - state_right(0).getDerivative(1) << " " << (STM(0, 7)  - state_right(0).getDerivative(1)) / state_right(0).getDerivative(1) << std::endl;
    std::cout << "dydIndVarcurrent: "      << STM(1, 13) << " " << state_right(1).getDerivative(2) << " " << STM(1, 13) - state_right(1).getDerivative(2) << " " << (STM(1, 13) - state_right(1).getDerivative(2)) / state_right(1).getDerivative(2) << std::endl;
    std::cout << "dydIndVarprevious: "     << STM(1, 7)  << " " << state_right(1).getDerivative(1) << " " << STM(1, 7)  - state_right(1).getDerivative(1) << " " << (STM(1, 7)  - state_right(1).getDerivative(1)) / state_right(1).getDerivative(1) << std::endl;
    std::cout << "dzdIndVarcurrent: "      << STM(2, 13) << " " << state_right(2).getDerivative(2) << " " << STM(2, 13) - state_right(2).getDerivative(2) << " " << (STM(2, 13) - state_right(2).getDerivative(2)) / state_right(2).getDerivative(2) << std::endl;
    std::cout << "dzdIndVarprevious: "     << STM(2, 7)  << " " << state_right(2).getDerivative(1) << " " << STM(2, 7)  - state_right(2).getDerivative(1) << " " << (STM(2, 7)  - state_right(2).getDerivative(1)) / state_right(2).getDerivative(1) << std::endl;
    std::cout << "dvxdIndVarcurrent: "     << STM(3, 13) << " " << state_right(3).getDerivative(2) << " " << STM(3, 13) - state_right(3).getDerivative(2) << " " << (STM(3, 13) - state_right(3).getDerivative(2)) / state_right(3).getDerivative(2) << std::endl;
    std::cout << "dvxdIndVarprevious: "    << STM(3, 7)  << " " << state_right(3).getDerivative(1) << " " << STM(3, 7)  - state_right(3).getDerivative(1) << " " << (STM(3, 7)  - state_right(3).getDerivative(1)) / state_right(3).getDerivative(1) << std::endl;
    std::cout << "dvydIndVarcurrent: "     << STM(4, 13) << " " << state_right(4).getDerivative(2) << " " << STM(4, 13) - state_right(4).getDerivative(2) << " " << (STM(4, 13) - state_right(4).getDerivative(2)) / state_right(4).getDerivative(2) << std::endl;
    std::cout << "dvydIndVarprevious: "    << STM(4, 7)  << " " << state_right(4).getDerivative(1) << " " << STM(4, 7)  - state_right(4).getDerivative(1) << " " << (STM(4, 7)  - state_right(4).getDerivative(1)) / state_right(4).getDerivative(1) << std::endl;
    std::cout << "dvzdIndVarcurrent: "     << STM(5, 13) << " " << state_right(5).getDerivative(2) << " " << STM(5, 13) - state_right(5).getDerivative(2) << " " << (STM(5, 13) - state_right(5).getDerivative(2)) / state_right(5).getDerivative(2) << std::endl;
    std::cout << "dvzdIndVarprevious: "    << STM(5, 7)  << " " << state_right(5).getDerivative(1) << " " << STM(5, 7)  - state_right(5).getDerivative(1) << " " << (STM(5, 7)  - state_right(5).getDerivative(1)) / state_right(5).getDerivative(1) << std::endl;
    std::cout << "dmdIndVarcurrent: "      << STM(6, 13) << " " << state_right(6).getDerivative(2) << " " << STM(6, 13) - state_right(6).getDerivative(2) << " " << (STM(6, 13) - state_right(6).getDerivative(2)) / state_right(6).getDerivative(2) << std::endl;
    std::cout << "dmdIndVarprevious: "     << STM(6, 7)  << " " << state_right(6).getDerivative(1) << " " << STM(6, 7)  - state_right(6).getDerivative(1) << " " << (STM(6, 7)  - state_right(6).getDerivative(1)) / state_right(6).getDerivative(1) << std::endl;
    std::cout << "depochdIndVarcurrent: "  << STM(7, 13) << " " << state_right(7).getDerivative(2) << " " << STM(7, 13) - state_right(7).getDerivative(2) << " " << (STM(7, 13) - state_right(7).getDerivative(2)) / state_right(7).getDerivative(2) << std::endl;
    std::cout << "depochdIndVarprevious: " << STM(7, 7)  << " " << state_right(7).getDerivative(1) << " " << STM(7, 7)  - state_right(7).getDerivative(1) << " " << (STM(7, 7)  - state_right(7).getDerivative(1)) / state_right(7).getDerivative(1) << std::endl;
    std::cout << "dchemdIndVarcurrent: "   << STM(8, 13) << " " << state_right(8).getDerivative(2) << " " << STM(8, 13) - state_right(8).getDerivative(2) << " " << (STM(8, 13) - state_right(8).getDerivative(2)) / state_right(8).getDerivative(2) << std::endl;
    std::cout << "dchemdIndVarprevious: "  << STM(8, 7)  << " " << state_right(8).getDerivative(1) << " " << STM(8, 7)  - state_right(8).getDerivative(1) << " " << (STM(8, 7)  - state_right(8).getDerivative(1)) / state_right(8).getDerivative(1) << std::endl;
    std::cout << "delecdIndVarcurrent: "   << STM(9, 13) << " " << state_right(9).getDerivative(2) << " " << STM(9, 13) - state_right(9).getDerivative(2) << " " << (STM(9, 13) - state_right(9).getDerivative(2)) / state_right(9).getDerivative(2) << std::endl;
    std::cout << "delecdIndVarprevious: "  << STM(9, 7)  << " " << state_right(9).getDerivative(1) << " " << STM(9, 7)  - state_right(9).getDerivative(1) << " " << (STM(9, 7)  - state_right(9).getDerivative(1)) / state_right(9).getDerivative(1) << std::endl;
#endif
    
    std::cout << "\nInitial State Information:\n" << std::endl;
    std::cout << "X [km]:    " << state_left(0) << std::endl;
    std::cout << "Y [km]:    " << state_left(1) << std::endl;
    std::cout << "Z [km]:    " << state_left(2) << std::endl;
    std::cout << "VX [km/s]: " << state_left(3) << std::endl;
    std::cout << "VY [km/s]: " << state_left(4) << std::endl;
    std::cout << "VZ [km/s]: " << state_left(5) << std::endl;
    std::cout << "Mass [kg]: " << state_left(6) << std::endl;
    std::cout << "Initial Epoch [JD TDB]: " << state_left(7) / 86400.0 << std::endl << std::endl;

    std::cout << "\nFinal State Information:\n" << std::endl;
    std::cout << "X [km]:    " << state_right(0) << std::endl;
    std::cout << "Y [km]:    " << state_right(1) << std::endl;
    std::cout << "Z [km]:    " << state_right(2) << std::endl;
    std::cout << "VX [km/s]: " << state_right(3) << std::endl;
    std::cout << "VY [km/s]: " << state_right(4) << std::endl;
    std::cout << "VZ [km/s]: " << state_right(5) << std::endl;
    std::cout << "Mass [kg]: " << state_right(6) << std::endl;
    std::cout << "Final Epoch [JD TDB]: " << state_right(7) / 86400.0 << std::endl << std::endl;

    std::cout << "Final State: " << state_right(0) << " " << state_right(1) << " " << state_right(2) << " " << state_right(3) << " " << state_right(4) << " " << state_right(5) << std::endl << std::endl;

    std::cout << "Propagation time [days]: " << (state_right(7) - state_left(7)) / 86400.0 << std::endl << std::endl;

#ifndef AD_INSTRUMENTATION
    /*
    EMTG::math::Matrix <doubleType> state_right_forward_1(state_vector_size, 1, 0.0);
    EMTG::math::Matrix <doubleType> state_right_backward_1(state_vector_size, 1, 0.0);
    EMTG::math::Matrix <doubleType> state_right_forward_2(state_vector_size, 1, 0.0);
    EMTG::math::Matrix <doubleType> state_right_backward_2(state_vector_size, 1, 0.0);
    EMTG::math::Matrix <doubleType> state_right_forward_3(state_vector_size, 1, 0.0);
    EMTG::math::Matrix <doubleType> state_right_backward_3(state_vector_size, 1, 0.0);
    std::vector<double> pert = { 10.0, 10.0, 10.0, 1.0, 1.0e-1, 1.0e-1, 1.0e-1, 10.0 };
    for (size_t pert_index = 0; pert_index < 7; ++pert_index)
    {
        state_left.assign_zeros();
        state_right.assign_zeros();
        for (size_t k = 0; k < 8; ++k)
        {
            if (k == pert_index)
            {
                state_left(k) = x[k] + pert[pert_index];
            }
            else
            {
                state_left(k) = x[k];
            }
        }
        state_left(8) = 0.0;
        state_left(9) = 0.0;
        
        test_acceleration_model.setLaunchEpoch(test_launch_epoch);
        test_acceleration_model.setEpoch(state_left(7));
        integratedPropagator->setCurrentEpoch(state_left(7));
        if (options.mission_type == EMTG::PhaseType::SundmanCoastPhase)
        {
            integratedPropagator->setCurrentIndependentVariable(0.0);
        }
        else
        {
            integratedPropagator->setCurrentIndependentVariable(test_state(7));
        }
        integratedPropagator->propagate(TOF, false);
        state_right_forward_1 = state_right;

        state_left.assign_zeros();
        state_right.assign_zeros();
        for (size_t k = 0; k < 8; ++k)
        {
            if (k == pert_index)
            {
                state_left(k) = x[k] - pert[pert_index];
            }
            else
            {
                state_left(k) = x[k];
            }
        }
        state_left(8) = 0.0;
        state_left(9) = 0.0;

        test_acceleration_model.setLaunchEpoch(test_launch_epoch);
        test_acceleration_model.setEpoch(state_left(7));
        integratedPropagator->setCurrentEpoch(state_left(7));
        if (options.mission_type == EMTG::PhaseType::SundmanCoastPhase)
        {
            integratedPropagator->setCurrentIndependentVariable(0.0);
        }
        else
        {
            integratedPropagator->setCurrentIndependentVariable(test_state(7));
        }
        integratedPropagator->propagate(TOF, false);
        state_right_backward_1 = state_right;

        state_left.assign_zeros();
        state_right.assign_zeros();
        for (size_t k = 0; k < 8; ++k)
        {
            if (k == pert_index)
            {
                state_left(k) = x[k] + 2.0*pert[pert_index];
            }
            else
            {
                state_left(k) = x[k];
            }
        }
        state_left(8) = 0.0;
        state_left(9) = 0.0;

        test_acceleration_model.setLaunchEpoch(test_launch_epoch);
        test_acceleration_model.setEpoch(state_left(7));
        integratedPropagator->setCurrentEpoch(state_left(7));
        if (options.mission_type == EMTG::PhaseType::SundmanCoastPhase)
        {
            integratedPropagator->setCurrentIndependentVariable(0.0);
        }
        else
        {
            integratedPropagator->setCurrentIndependentVariable(test_state(7));
        }
        integratedPropagator->propagate(TOF, false);
        state_right_forward_2 = state_right;

        state_left.assign_zeros();
        state_right.assign_zeros();
        for (size_t k = 0; k < 8; ++k)
        {
            if (k == pert_index)
            {
                state_left(k) = x[k] - 2.0*pert[pert_index];
            }
            else
            {
                state_left(k) = x[k];
            }
        }
        state_left(8) = 0.0;
        state_left(9) = 0.0;

        test_acceleration_model.setLaunchEpoch(test_launch_epoch);
        test_acceleration_model.setEpoch(state_left(7));
        integratedPropagator->setCurrentEpoch(state_left(7));
        if (options.mission_type == EMTG::PhaseType::SundmanCoastPhase)
        {
            integratedPropagator->setCurrentIndependentVariable(0.0);
        }
        else
        {
            integratedPropagator->setCurrentIndependentVariable(test_state(7));
        }
        integratedPropagator->propagate(TOF, false);
        state_right_backward_2 = state_right;

        state_left.assign_zeros();
        state_right.assign_zeros();
        for (size_t k = 0; k < 8; ++k)
        {
            if (k == pert_index)
            {
                state_left(k) = x[k] + 3.0*pert[pert_index];
            }
            else
            {
                state_left(k) = x[k];
            }
        }
        state_left(8) = 0.0;
        state_left(9) = 0.0;

        test_acceleration_model.setLaunchEpoch(test_launch_epoch);
        test_acceleration_model.setEpoch(state_left(7));
        integratedPropagator->setCurrentEpoch(state_left(7));
        if (options.mission_type == EMTG::PhaseType::SundmanCoastPhase)
        {
            integratedPropagator->setCurrentIndependentVariable(0.0);
        }
        else
        {
            integratedPropagator->setCurrentIndependentVariable(test_state(7));
        }
        integratedPropagator->propagate(TOF, false);
        state_right_forward_3 = state_right;

        state_left.assign_zeros();
        state_right.assign_zeros();
        for (size_t k = 0; k < 8; ++k)
        {
            if (k == pert_index)
            {
                state_left(k) = x[k] - 3.0*pert[pert_index];
            }
            else
            {
                state_left(k) = x[k];
            }
        }
        state_left(8) = 0.0;
        state_left(9) = 0.0;

        test_acceleration_model.setLaunchEpoch(test_launch_epoch);
        test_acceleration_model.setEpoch(state_left(7));
        integratedPropagator->setCurrentEpoch(state_left(7));
        if (options.mission_type == EMTG::PhaseType::SundmanCoastPhase)
        {
            integratedPropagator->setCurrentIndependentVariable(0.0);
        }
        else
        {
            integratedPropagator->setCurrentIndependentVariable(test_state(7));
        }
        integratedPropagator->propagate(TOF, false);
        state_right_backward_3 = state_right;

        for (size_t k = 0; k < 8; ++k)
        {
            //std::cout << ((2.0/3.0)*state_right_forward_1(k) - (1.0/12.0)*state_right_forward_2(k) - (2.0/3.0)*state_right_backward_1(k) + (1.0 / 12.0)*state_right_backward_2(k)) / (pert[pert_index]) << " ";
            std::cout << ((3.0 / 4.0)*state_right_forward_1(k) - (3.0 / 20.0)*state_right_forward_2(k) + (1.0 / 60.0)*state_right_forward_3(k) - (3.0 / 4.0)*state_right_backward_1(k) + (3.0 / 20.0)*state_right_backward_2(k) - (1.0 / 60.0)*state_right_backward_3(k)) / (pert[pert_index]) << " ";
        }

        std::cout << std::endl << std::endl;
    }
    */
#endif

    getchar();
}