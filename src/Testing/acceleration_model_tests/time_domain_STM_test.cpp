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
#include "file_utilities.h"
#include "IntegratedAdaptiveStepPropagator.h"
#include "IntegratedFixedStepPropagator.h"
#include "mission.h"
#include "missionoptions.h"
#include "PropagatorFactory.h"
#include "RungeKutta8.h"
#include "SpacecraftAccelerationModel.h"
#include "SpiceUsr.h" 
#include "TimeDomainSpacecraftEOM.h"
#include "universe.h"


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
    size_t STM_start_index = 10;
    size_t num_STM_rows = 14;
    size_t num_STM_columns = 14;
    std::vector<std::string> Xdescriptions;
    EMTG::Astrodynamics::SpacecraftAccelerationModel test_acceleration_model(&myMission.options,
        &myMission.options.Journeys[0],
        &TheUniverse[0],
        &Xdescriptions,
        &mySpacecraft,
        num_STM_rows - 1,
        num_STM_columns - 1,
        STM_start_index);

    size_t GSADindex = 0;

    doubleType test_launch_epoch = myMission.options.launch_window_open_date; // seconds past MJD J2000
    doubleType TOFprevious = 0.0 * 86400.0;
    doubleType TOF = 9.8* 86400.0;
    doubleType test_epoch;

    test_launch_epoch.setDerivative(GSADindex++, 1.0);
    TOFprevious.setDerivative(GSADindex++, 1.0);
    TOF.setDerivative(GSADindex++, 1.0);

    test_epoch = test_launch_epoch + TOFprevious;

    //initialize state vector
    double LT_dump;
    std::vector<double> x(8);
    double shit = test_epoch _GETVALUE - (51544.5 * 86400.0);
    //spkez_c(399, shit, "J2000", "NONE", 10, x.data(), &LT_dump);
    spkez_c(504, test_epoch _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 599, x.data(), &LT_dump);
    //spkez_c(499, test_epoch _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 4, x.data(), &LT_dump);
    //spkez_c(602, test_epoch _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 699, x.data(), &LT_dump);
    x[0] = -1169636.86875674;
    x[1] = -781158.68641361;
    x[2] = -337861.83290009;
    x[3] = 4.67341199;
    x[4] = 3.04126260;
    x[5] = 1.31500098;
    x[6] = myMission.options.maximum_mass;
    //x[7] = test_epoch _GETVALUE;

    //create perturbed initial state for TOF derivative use
    std::vector <double> x_pert(7, 0.0);
    std::vector <double> dx_dTOFprevious(7, 0.0);
    std::vector <double> dx_dTOFcurrent(7, 0.0);
    //spkez_c(399, shit + 10.0, "J2000", "NONE", 10, x_pert.data(), &LT_dump);
    spkez_c(504, test_epoch _GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", 599, x_pert.data(), &LT_dump);
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
        test_state(k).setDerivative(GSADindex++, 1.0);
        //test_state(k).setDerivative(0, dx_dTOFprevious[k]);
        //test_state(k).setDerivative(1, dx_dTOFprevious[k]);
        //test_state(k).setDerivative(2, dx_dTOFcurrent[k]);
    }

    // current epoch
    test_state(7) = test_epoch;
    test_state(7).setDerivative(GSADindex++, 1.0);

    // virtual chemical thruster propellant tank
    test_state(8) = 0.0;
    test_state(8).setDerivative(GSADindex++, 1.0);

    // virtual electric thruster propellant tank 
    test_state(9) = 0.0;
    test_state(9).setDerivative(GSADindex++, 1.0);

    // control decision variables
    EMTG::math::Matrix<doubleType> control(4, 1, 0.0);
    for (size_t i = 0; i < 3; ++i)
    {
#ifdef AD_INSTRUMENTATION
        control(i).setValue(DoubleDistribution(RNG));
        //control(i).setValue(0.0);
        //segmentControls[3 * segment + i].setValue(0.0);
        //control(i).setDerivative(GSADindex++, 1.0);
#else
        //segmentControls[3 * segment + i] = DoubleDistribution(RNG);
        segmentControls[3 * segment + i] = 0.0;
#endif
    }

    test_acceleration_model.setDutyCycle(myMission.options.engine_duty_cycle);
    test_acceleration_model.setLaunchEpoch(test_launch_epoch);
    test_acceleration_model.setEpoch(test_state(7));
    //test_acceleration_model.computeAcceleration(test_state, control, true);
    //EMTG::math::Matrix<doubleType> acceleration = test_acceleration_model.getAcceleration();

    // create the spacecraft equations of motion
    EMTG::Astrodynamics::TimeDomainSpacecraftEOM FiniteBurnEOM = EMTG::Astrodynamics::TimeDomainSpacecraftEOM();
    FiniteBurnEOM.setSpacecraftAccelerationModel(&test_acceleration_model);

    // create the integrand pointer
    EMTG::Integration::Integrand * thing = &FiniteBurnEOM;

    // create the integration scheme
    // {x, y, z, vx, vy, vz, m, ux, uy, uz}
    size_t total_states_to_integrate = STM_start_index + (num_STM_rows - 1) * (num_STM_columns - 1);
    size_t number_of_prop_var_derivative_states = 10;
    EMTG::Integration::RungeKutta8 rk8(thing, total_states_to_integrate, number_of_prop_var_derivative_states);

    // specify the current number of states that you want to integrate 
    // basically do you want states+STM or just states?
    rk8.setNumStatesToIntegratePtr(total_states_to_integrate);

    // create the propagator
    EMTG::math::Matrix <doubleType> state_left(state_vector_size, 1, 0.0);
    EMTG::math::Matrix <double> dstate_leftdTOF(10, 2, 0.0); // initial state partials w.r.t. previous and current TOF
    EMTG::math::Matrix <doubleType> state_right(state_vector_size, 1, 0.0);
    EMTG::math::Matrix <double> dstate_rightdTOF(10, 2, 0.0); // final state partials w.r.t. previous and current TOF
    EMTG::math::Matrix <double> STM(num_STM_rows, num_STM_columns, 0.0);
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
        num_STM_rows - 1,
        num_STM_columns - 1,
        STM_start_index,
        state_left,
        state_right,
        STM,
        dstate_leftdTOF,
        (EMTG::Integration::Integrand*) &FiniteBurnEOM,
        &rk8,
        &BoundaryTarget_dStepSizedPropVar,
        options.Journeys[0].override_integration_step_size
        ? options.Journeys[0].integration_step_size
        : options.integration_time_step_size);

    integratedPropagator->setCurrentEpoch(test_state(7));
    integratedPropagator->setIndexOfEpochInStateVec(7); // where is epoch located in the state vector?
    integratedPropagator->setCurrentIndependentVariable(test_state(7));


    //options.X_scale_factors[0] = 1.0;
    //options.X_scale_factors[1] = 1.0;
    //options.X_scale_factors[2] = 1.0;
    //integratedFSProp.propagate(TOF / options.num_timesteps, control, true);
    integratedPropagator->propagate(TOF / options.num_timesteps, true);

    dstate_rightdTOF = dstate_leftdTOF;

    std::cout << std::setprecision(16);

    std::cout << "Numerically integrated STM:" << std::endl;
    for (size_t i = 0; i < num_STM_rows; ++i)
    {
        for (size_t j = 0; j < num_STM_columns; ++j)
        {
            std::cout << STM(i, j) << "     ";
        }
        std::cout << std::endl << std::endl;
    }

    std::cout << std::endl << std::endl;

    std::cout << "AD computed STM:" << std::endl;
    EMTG::math::Matrix <double> AD_STM(num_STM_rows, num_STM_columns, 0.0);
    size_t derivIndex = 3;
    for (size_t i = 0; i < state_vector_size; ++i)
    {
        for (size_t derivIndex = 3; derivIndex < num_STM_columns + 2; ++derivIndex)
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
    for (size_t i = 0; i < num_STM_rows; ++i)
    {
        for (size_t j = 0; j < num_STM_columns; ++j)
        {
            std::cout << AD_STM(i, j) - STM(i, j) << "     ";
        }
        std::cout << std::endl << std::endl;
    }

    std::cout << std::endl << std::endl;

    std::cout << "Relative STM error:" << std::endl;
    for (size_t i = 0; i < num_STM_rows; ++i)
    {
        for (size_t j = 0; j < num_STM_columns; ++j)
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

    std::cout << "Time of Flight derivatives" << std::endl;
    std::cout << "Format: analytical, GSAD, absolute error, relative error" << std::endl;
    std::cout << "dxdTOFcurrent: " << dstate_rightdTOF(0, 1) << " " << state_right(0).getDerivative(2) << " " << dstate_rightdTOF(0, 1) - state_right(0).getDerivative(2) << " " << (dstate_rightdTOF(0, 1) - state_right(0).getDerivative(2)) / state_right(0).getDerivative(2) << std::endl;
    std::cout << "dxdTOFprevious: " << dstate_rightdTOF(0, 0) << " " << state_right(0).getDerivative(1) << " " << dstate_rightdTOF(0, 0) - state_right(0).getDerivative(1) << " " << (dstate_rightdTOF(0, 0) - state_right(0).getDerivative(1)) / state_right(0).getDerivative(1) << std::endl;
    std::cout << "dydTOFcurrent: " << dstate_rightdTOF(1, 1) << " " << state_right(1).getDerivative(2) << " " << dstate_rightdTOF(1, 1) - state_right(1).getDerivative(2) << " " << (dstate_rightdTOF(1, 1) - state_right(1).getDerivative(2)) / state_right(1).getDerivative(2) << std::endl;
    std::cout << "dydTOFprevious: " << dstate_rightdTOF(1, 0) << " " << state_right(1).getDerivative(1) << " " << dstate_rightdTOF(1, 0) - state_right(1).getDerivative(1) << " " << (dstate_rightdTOF(1, 0) - state_right(1).getDerivative(1)) / state_right(1).getDerivative(1) << std::endl;
    std::cout << "dzdTOFcurrent: " << dstate_rightdTOF(2, 1) << " " << state_right(2).getDerivative(2) << " " << dstate_rightdTOF(2, 1) - state_right(2).getDerivative(2) << " " << (dstate_rightdTOF(2, 1) - state_right(2).getDerivative(2)) / state_right(2).getDerivative(2) << std::endl;
    std::cout << "dzdTOFprevious: " << dstate_rightdTOF(2, 0) << " " << state_right(2).getDerivative(1) << " " << dstate_rightdTOF(2, 0) - state_right(2).getDerivative(1) << " " << (dstate_rightdTOF(2, 0) - state_right(2).getDerivative(1)) / state_right(2).getDerivative(1) << std::endl;
    std::cout << "dvxdTOFcurrent: " << dstate_rightdTOF(3, 1) << " " << state_right(3).getDerivative(2) << " " << dstate_rightdTOF(3, 1) - state_right(3).getDerivative(2) << " " << (dstate_rightdTOF(3, 1) - state_right(3).getDerivative(2)) / state_right(3).getDerivative(2) << std::endl;
    std::cout << "dvxdTOFprevious: " << dstate_rightdTOF(3, 0) << " " << state_right(3).getDerivative(1) << " " << dstate_rightdTOF(3, 0) - state_right(3).getDerivative(1) << " " << (dstate_rightdTOF(3, 0) - state_right(3).getDerivative(1)) / state_right(3).getDerivative(1) << std::endl;
    std::cout << "dvydTOFcurrent: " << dstate_rightdTOF(4, 1) << " " << state_right(4).getDerivative(2) << " " << dstate_rightdTOF(4, 1) - state_right(4).getDerivative(2) << " " << (dstate_rightdTOF(4, 1) - state_right(4).getDerivative(2)) / state_right(4).getDerivative(2) << std::endl;
    std::cout << "dvydTOFprevious: " << dstate_rightdTOF(4, 0) << " " << state_right(4).getDerivative(1) << " " << dstate_rightdTOF(4, 0) - state_right(4).getDerivative(1) << " " << (dstate_rightdTOF(4, 0) - state_right(4).getDerivative(1)) / state_right(4).getDerivative(1) << std::endl;
    std::cout << "dvzdTOFcurrent: " << dstate_rightdTOF(5, 1) << " " << state_right(5).getDerivative(2) << " " << dstate_rightdTOF(5, 1) - state_right(5).getDerivative(2) << " " << (dstate_rightdTOF(5, 1) - state_right(5).getDerivative(2)) / state_right(5).getDerivative(2) << std::endl;
    std::cout << "dvzdTOFprevious: " << dstate_rightdTOF(5, 0) << " " << state_right(5).getDerivative(1) << " " << dstate_rightdTOF(5, 0) - state_right(5).getDerivative(1) << " " << (dstate_rightdTOF(5, 0) - state_right(5).getDerivative(1)) / state_right(5).getDerivative(1) << std::endl;
    std::cout << "dmdTOFcurrent: " << dstate_rightdTOF(6, 1) << " " << state_right(6).getDerivative(2) << " " << dstate_rightdTOF(6, 1) - state_right(6).getDerivative(2) << " " << (dstate_rightdTOF(6, 1) - state_right(6).getDerivative(2)) / state_right(6).getDerivative(2) << std::endl;
    std::cout << "dmdTOFprevious: " << dstate_rightdTOF(6, 0) << " " << state_right(6).getDerivative(1) << " " << dstate_rightdTOF(6, 0) - state_right(6).getDerivative(1) << " " << (dstate_rightdTOF(6, 0) - state_right(6).getDerivative(1)) / state_right(6).getDerivative(1) << std::endl;
    std::cout << "depochdIndVarcurrent: " << dstate_rightdTOF(7, 1) << " " << state_right(7).getDerivative(2) << " " << dstate_rightdTOF(7, 1) - state_right(7).getDerivative(2) << " " << (dstate_rightdTOF(7, 1) - state_right(7).getDerivative(2)) / state_right(7).getDerivative(2) << std::endl;
    std::cout << "depochdIndVarprevious: " << dstate_rightdTOF(7, 0) << " " << state_right(7).getDerivative(1) << " " << dstate_rightdTOF(7, 0) - state_right(7).getDerivative(1) << " " << (dstate_rightdTOF(7, 0) - state_right(7).getDerivative(1)) / state_right(7).getDerivative(1) << std::endl;
    std::cout << "dchemdTOFcurrent: " << dstate_rightdTOF(8, 1) << " " << state_right(8).getDerivative(2) << " " << dstate_rightdTOF(8, 1) - state_right(8).getDerivative(2) << " " << (dstate_rightdTOF(8, 1) - state_right(8).getDerivative(2)) / state_right(8).getDerivative(2) << std::endl;
    std::cout << "dchemdTOFprevious: " << dstate_rightdTOF(8, 0) << " " << state_right(8).getDerivative(1) << " " << dstate_rightdTOF(8, 0) - state_right(8).getDerivative(1) << " " << (dstate_rightdTOF(8, 0) - state_right(8).getDerivative(1)) / state_right(8).getDerivative(1) << std::endl;
    std::cout << "delecdTOFcurrent: " << dstate_rightdTOF(9, 1) << " " << state_right(9).getDerivative(2) << " " << dstate_rightdTOF(9, 1) - state_right(9).getDerivative(2) << " " << (dstate_rightdTOF(9, 1) - state_right(9).getDerivative(2)) / state_right(9).getDerivative(2) << std::endl;
    std::cout << "delecdTOFprevious: " << dstate_rightdTOF(9, 0) << " " << state_right(9).getDerivative(1) << " " << dstate_rightdTOF(9, 0) - state_right(9).getDerivative(1) << " " << (dstate_rightdTOF(9, 0) - state_right(9).getDerivative(1)) / state_right(9).getDerivative(1) << std::endl;

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
    std::cout << "Final Epoch [JD TDB]: " << state_right(7) / 86400.0 << std::endl;


    getchar();


}