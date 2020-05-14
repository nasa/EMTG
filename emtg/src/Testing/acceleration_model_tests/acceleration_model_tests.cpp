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

    EMTG::missionoptions options(options_file_name);
    std::cout << options.error_message << std::endl;
    if (!(options.file_status == 0))
    {
        std::cout << "Aborting program run." << std::endl;
        std::cin.ignore();
        return 0;
    }

    //configure the LaunchVehicleOptions and SpacecraftOptions objects
    EMTG::HardwareModels::LaunchVehicleOptions myLaunchVehicleOptions = EMTG::HardwareModels::CreateLaunchVehicleOptions(options);
    EMTG::HardwareModels::SpacecraftOptions mySpacecraftOptions = EMTG::HardwareModels::CreateSpacecraftOptions(options);

    EMTG::HardwareModels::LaunchVehicle myLaunchVehicle(myLaunchVehicleOptions);
    EMTG::HardwareModels::Spacecraft mySpacecraft(mySpacecraftOptions);

    //Load appropriate SPICE kernels
    //load DE430, leap seconds kernel, frame kernel

    //load all SPICE ephemeris data
    std::vector<fs::path> SPICE_files;
    EMTG::filesystem::get_all_files_with_extension(fs::path(options.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files);
    EMTG::filesystem::get_all_files_with_extension(fs::path(options.universe_folder + "/ephemeris_files/"), ".cmt", SPICE_files);

    std::string filestring;
    for (size_t k = 0; k < SPICE_files.size(); ++k)
    {
        filestring = options.universe_folder + "/ephemeris_files/" + SPICE_files[k].string();
        furnsh_c(filestring.c_str());
    }
    std::string leapsecondstring = options.universe_folder + "/ephemeris_files/" + options.SPICE_leap_seconds_kernel;
    std::string referenceframestring = options.universe_folder + "/ephemeris_files/" + options.SPICE_reference_frame_kernel;
    furnsh_c(leapsecondstring.c_str());
    furnsh_c(referenceframestring.c_str());

    //disable quit-on-SPICE-error so that we can see what happens if the leap second and/or frame kernels don't load properly
    erract_c((SpiceChar*)"SET", 100, (SpiceChar*)"RETURN");

    //disable SPICE error printing. This is because we can, and will often, go off the edge of an ephemeris file.
    //errprt_c((SpiceChar*)"SET", 100, (SpiceChar*)"NONE");

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
    //if (options.ephemeris_source == 2)
    {
        for (size_t j = 0; j < options.number_of_journeys; ++j)
        {
            std::vector<int> body_index_array;

            //first boundary point
            if (options.Journeys[j].destination_list[0] > 0)
                body_index_array.push_back(options.Journeys[j].destination_list[0] - 1);

            //last boundary point
            if (options.Journeys[j].destination_list[1] > 0)
                body_index_array.push_back(options.Journeys[j].destination_list[1] - 1);

            //sequence
            for (size_t b = 0; b < options.Journeys[j].sequence_input.size(); ++b)
                if (options.Journeys[j].sequence_input[b] > 0)
                    body_index_array.push_back(options.Journeys[j].sequence_input[b] - 1);

            //perturbation list
            for (size_t b = 0; b < TheUniverse[j].perturbation_menu.size(); ++b)
                body_index_array.push_back(TheUniverse[j].perturbation_menu[b]);

            //distance constraint list
            for (size_t constraintIndex = 0; constraintIndex < options.PhaseDistanceConstraintDefinitions.size(); ++constraintIndex)
            {
                if (options.PhaseDistanceConstraintDefinitions[constraintIndex].find("j" + std::to_string(j) + "p") < 1024)
                {
                    std::vector<std::string> ConstraintDefinitionCell;
                    boost::split(ConstraintDefinitionCell,
                        options.PhaseDistanceConstraintDefinitions[constraintIndex],
                        boost::is_any_of("_"),
                        boost::token_compress_on);

                    if (boost::to_lower_copy(ConstraintDefinitionCell[1]) != "cb")
                    {
                        int bodyIndex = std::stoi(ConstraintDefinitionCell[1]) - 1;

                        body_index_array.push_back(bodyIndex);
                    }
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
                        360.0,
                        TheUniverse[j].mu));
                }
            }//end loop over bodies in the universe

             //is this universe's central body the sun? If not, let's add this body with respect to the sun. Let's add extra ephemeris points, too.
            if (!(TheUniverse[j].central_body_SPICE_ID == 10))
            {
                SplineUniverse_keyList.push_back(std::make_tuple(
                    TheUniverse[j].central_body_SPICE_ID,
                    10,
                    1000,
                    1.32712440018e+11));
            }
        }
        SplineUniverse.reinitialize(SplineUniverse_keyList,
            options.launch_window_open_date + options.Journeys.front().journey_wait_time_bounds[0] - 10.0 * 86400.0,
            options.launch_window_open_date + options.Journeys.front().journey_wait_time_bounds[0] + 18262.5 * 86400.0);
    }
#endif

    {
        //do this the simple way
        options.description.clear();

        for (size_t j = 0; j < options.number_of_journeys; ++j)
        {
            options.Journeys[j].sequence.clear();
            options.Journeys[j].impulses_per_phase.clear();
            options.Journeys[j].phase_type.clear();

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
            for (size_t p = 0; p < options.max_phases_per_journey; ++p)
            {
                size_t bodyIndex = options.Journeys[j].sequence_input[p];
                if (bodyIndex > 0 && bodyIndex < (TheUniverse[j].size_of_flyby_menu / 2) + 1) //this is a legitimate flyby
                {
                    if (bodyIndex - 1 > TheUniverse[j].flyby_menu.size())
                    {
                        std::cout << "ERROR: Journey " << j << " phase " << p << " body index " << bodyIndex << " exceeds size of flyby menu. Aborting." << std::endl;
                        throw 888;
                    }

                    options.Journeys[j].sequence.push_back(bodyIndex);

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

            //then phase types and number of impulses
            for (size_t p = 0; p < options.Journeys[j].number_of_phases + 1; ++p)
            {
                if (options.mission_type == EMTG::PhaseType::VARIABLE_PHASE_TYPE)
                    options.Journeys[j].phase_type.push_back(options.Journeys[j].phase_type_input[p]);
                else
                    options.Journeys[j].phase_type.push_back(options.mission_type);

                options.Journeys[j].impulses_per_phase.push_back(options.Journeys[j].impulses_per_phase_input[p]);
            }
        }
    }

    EMTG::Mission myMission(options,
                            TheUniverse,
                            myLaunchVehicle,
                            mySpacecraft);

     
    size_t state_vector_size = 10;
    size_t STM_start_index = 10;
    size_t num_STM_rows = 13;
    size_t num_STM_columns = 13;
    EMTG::Astrodynamics::SpacecraftAccelerationModel test_acceleration_model(&myMission.options, 
                                                                             &myMission.options.Journeys[0], 
                                                                             &TheUniverse[0], 
                                                                             &mySpacecraft, 
                                                                             num_STM_rows, 
                                                                             num_STM_columns,
                                                                             STM_start_index);
    
    size_t GSADindex = 0;

    doubleType test_launch_epoch = myMission.options.launch_window_open_date; // seconds past MJD J2000
    doubleType TOFprevious = 0.0 * 86400.0;
    doubleType TOF = 1.0 * 86400.0;
    doubleType test_epoch;

    test_launch_epoch.setDerivative(GSADindex++, 1.0);
    TOFprevious.setDerivative(GSADindex++, 1.0);
    TOF.setDerivative(GSADindex++, 1.0);

    test_epoch = test_launch_epoch + TOFprevious;

    //initialize state vector
    double LT_dump;
    std::vector<double> x(8);
    double shit = test_epoch _GETVALUE - (51544.5 * 86400.0);
    spkez_c(399, shit, "J2000", "NONE", 10, x.data(), &LT_dump);
    //spkez_c(501, test_epoch _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 599, x.data(), &LT_dump);
    //spkez_c(499, test_epoch _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 4, x.data(), &LT_dump);
    //spkez_c(602, test_epoch _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 699, x.data(), &LT_dump);
    x[6] = myMission.options.maximum_mass;
    //x[7] = test_epoch _GETVALUE;

    //create perturbed initial state for TOF derivative use
    std::vector <double> x_pert(7, 0.0);
    std::vector <double> dx_dTOFprevious(7, 0.0);
    std::vector <double> dx_dTOFcurrent(7, 0.0);
    spkez_c(399, shit + 10.0, "J2000", "NONE", 10, x_pert.data(), &LT_dump);
    //spkez_c(501, test_epoch _GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", 599, x_pert.data(), &LT_dump);
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
        //control(i).setValue(DoubleDistribution(RNG));
        control(i).setValue(0.0);
        //segmentControls[3 * segment + i].setValue(0.0);
        control(i).setDerivative(GSADindex++, 1.0);
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
    size_t total_states_to_integrate = STM_start_index + num_STM_rows * num_STM_columns;
    size_t number_of_prop_var_derivative_states = 10;
    EMTG::Integration::RungeKutta8 rk8(thing, total_states_to_integrate, number_of_prop_var_derivative_states);

    // specify the current number of states that you want to integrate 
    // basically do you want states+STM or just states?
    rk8.setNumStatesToIntegratePtr(total_states_to_integrate);

    // create the propagator
    EMTG::Astrodynamics::IntegratedFixedStepPropagator integratedFSProp(myMission.options, TheUniverse[0], num_STM_rows, num_STM_columns, STM_start_index);
    integratedFSProp.setIntegrand(thing);
    integratedFSProp.setCurrentEpoch(test_state(7));
    integratedFSProp.setIndexOfEpochInStateVec(7); // where is epoch located in the state vector?
    integratedFSProp.setCurrentIndependentVariable(test_state(7));
    integratedFSProp.setIntegrationScheme(&rk8);

    
    EMTG::math::Matrix <doubleType> state_left (state_vector_size, 1, 0.0);
    EMTG::math::Matrix <double> dstate_leftdTOF(10, 2, 0.0); // initial state partials w.r.t. previous and current TOF
    EMTG::math::Matrix <doubleType> state_right (state_vector_size, 1, 0.0);
    EMTG::math::Matrix <double> dstate_rightdTOF(10, 2, 0.0); // final state partials w.r.t. previous and current TOF
    EMTG::math::Matrix <double> STM (num_STM_rows, num_STM_columns, 0.0);

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

    integratedFSProp.setStateLeft(state_left);
    integratedFSProp.setStateRight(state_right);
    integratedFSProp.setSTM(STM);
    integratedFSProp.setdStatedIndependentVariable(dstate_leftdTOF);
    integratedFSProp.setdCurrentIndVardPropVarPrevious(1.0); // previous phase flight times and launch epoch
    integratedFSProp.setdCurrentIndVardPropVar(0.0); // current phase flight time
    double BoundaryTarget_dStepSizedPropVar = 1.0 / options.num_timesteps;
    integratedFSProp.setBoundaryTarget_dStepSizedPropVar(&BoundaryTarget_dStepSizedPropVar);

    //options.X_scale_factors[0] = 1.0;
    //options.X_scale_factors[1] = 1.0;
    //options.X_scale_factors[2] = 1.0;
    integratedFSProp.propagate(TOF / options.num_timesteps, control, true);
    //integratedFSProp.propagate(TOF, true);

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
        for (size_t derivIndex = 3; derivIndex < num_STM_columns + 3; ++derivIndex)
        {
            AD_STM(i, derivIndex - 3) = state_right(i).getDerivative(derivIndex);
            std::cout << AD_STM(i, derivIndex - 3) << "     ";
        }
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
    std::cout << "dxdTOFcurrent: "     << dstate_rightdTOF(0, 1) << " " << state_right(0).getDerivative(2) << " " << dstate_rightdTOF(0, 1) - state_right(0).getDerivative(2) << " " << (dstate_rightdTOF(0, 1) - state_right(0).getDerivative(2)) / state_right(0).getDerivative(2) << std::endl;
    std::cout << "dxdTOFprevious: "    << dstate_rightdTOF(0, 0) << " " << state_right(0).getDerivative(1) << " " << dstate_rightdTOF(0, 0) - state_right(0).getDerivative(1) << " " << (dstate_rightdTOF(0, 0) - state_right(0).getDerivative(1)) / state_right(0).getDerivative(1) << std::endl;
    std::cout << "dydTOFcurrent: "     << dstate_rightdTOF(1, 1) << " " << state_right(1).getDerivative(2) << " " << dstate_rightdTOF(1, 1) - state_right(1).getDerivative(2) << " " << (dstate_rightdTOF(1, 1) - state_right(1).getDerivative(2)) / state_right(1).getDerivative(2) << std::endl;
    std::cout << "dydTOFprevious: "    << dstate_rightdTOF(1, 0) << " " << state_right(1).getDerivative(1) << " " << dstate_rightdTOF(1, 0) - state_right(1).getDerivative(1) << " " << (dstate_rightdTOF(1, 0) - state_right(1).getDerivative(1)) / state_right(1).getDerivative(1) << std::endl;
    std::cout << "dzdTOFcurrent: "     << dstate_rightdTOF(2, 1) << " " << state_right(2).getDerivative(2) << " " << dstate_rightdTOF(2, 1) - state_right(2).getDerivative(2) << " " << (dstate_rightdTOF(2, 1) - state_right(2).getDerivative(2)) / state_right(2).getDerivative(2) << std::endl;
    std::cout << "dzdTOFprevious: "    << dstate_rightdTOF(2, 0) << " " << state_right(2).getDerivative(1) << " " << dstate_rightdTOF(2, 0) - state_right(2).getDerivative(1) << " " << (dstate_rightdTOF(2, 0) - state_right(2).getDerivative(1)) / state_right(2).getDerivative(1) << std::endl;
    std::cout << "dvxdTOFcurrent: "    << dstate_rightdTOF(3, 1) << " " << state_right(3).getDerivative(2) << " " << dstate_rightdTOF(3, 1) - state_right(3).getDerivative(2) << " " << (dstate_rightdTOF(3, 1) - state_right(3).getDerivative(2)) / state_right(3).getDerivative(2) << std::endl;
    std::cout << "dvxdTOFprevious: "   << dstate_rightdTOF(3, 0) << " " << state_right(3).getDerivative(1) << " " << dstate_rightdTOF(3, 0) - state_right(3).getDerivative(1) << " " << (dstate_rightdTOF(3, 0) - state_right(3).getDerivative(1)) / state_right(3).getDerivative(1) << std::endl;
    std::cout << "dvydTOFcurrent: "    << dstate_rightdTOF(4, 1) << " " << state_right(4).getDerivative(2) << " " << dstate_rightdTOF(4, 1) - state_right(4).getDerivative(2) << " " << (dstate_rightdTOF(4, 1) - state_right(4).getDerivative(2)) / state_right(4).getDerivative(2) << std::endl;
    std::cout << "dvydTOFprevious: "   << dstate_rightdTOF(4, 0) << " " << state_right(4).getDerivative(1) << " " << dstate_rightdTOF(4, 0) - state_right(4).getDerivative(1) << " " << (dstate_rightdTOF(4, 0) - state_right(4).getDerivative(1)) / state_right(4).getDerivative(1) << std::endl;
    std::cout << "dvzdTOFcurrent: "    << dstate_rightdTOF(5, 1) << " " << state_right(5).getDerivative(2) << " " << dstate_rightdTOF(5, 1) - state_right(5).getDerivative(2) << " " << (dstate_rightdTOF(5, 1) - state_right(5).getDerivative(2)) / state_right(5).getDerivative(2) << std::endl;
    std::cout << "dvzdTOFprevious: "   << dstate_rightdTOF(5, 0) << " " << state_right(5).getDerivative(1) << " " << dstate_rightdTOF(5, 0) - state_right(5).getDerivative(1) << " " << (dstate_rightdTOF(5, 0) - state_right(5).getDerivative(1)) / state_right(5).getDerivative(1) << std::endl;
    std::cout << "dmdTOFcurrent: "     << dstate_rightdTOF(6, 1) << " " << state_right(6).getDerivative(2) << " " << dstate_rightdTOF(6, 1) - state_right(6).getDerivative(2) << " " << (dstate_rightdTOF(6, 1) - state_right(6).getDerivative(2)) / state_right(6).getDerivative(2) << std::endl;
    std::cout << "dmdTOFprevious: "    << dstate_rightdTOF(6, 0) << " " << state_right(6).getDerivative(1) << " " << dstate_rightdTOF(6, 0) - state_right(6).getDerivative(1) << " " << (dstate_rightdTOF(6, 0) - state_right(6).getDerivative(1)) / state_right(6).getDerivative(1) << std::endl;
    std::cout << "dchemdTOFcurrent: "  << dstate_rightdTOF(8, 1) << " " << state_right(8).getDerivative(2) << " " << dstate_rightdTOF(8, 1) - state_right(8).getDerivative(2) << " " << (dstate_rightdTOF(8, 1) - state_right(8).getDerivative(2)) / state_right(8).getDerivative(2) << std::endl;
    std::cout << "dchemdTOFprevious: " << dstate_rightdTOF(8, 0) << " " << state_right(8).getDerivative(1) << " " << dstate_rightdTOF(8, 0) - state_right(8).getDerivative(1) << " " << (dstate_rightdTOF(8, 0) - state_right(8).getDerivative(1)) / state_right(8).getDerivative(1) << std::endl;
    std::cout << "delecdTOFcurrent: "  << dstate_rightdTOF(9, 1) << " " << state_right(9).getDerivative(2) << " " << dstate_rightdTOF(9, 1) - state_right(9).getDerivative(2) << " " << (dstate_rightdTOF(9, 1) - state_right(9).getDerivative(2)) / state_right(9).getDerivative(2) << std::endl;
    std::cout << "delecdTOFprevious: " << dstate_rightdTOF(9, 0) << " " << state_right(9).getDerivative(1) << " " << dstate_rightdTOF(9, 0) - state_right(9).getDerivative(1) << " " << (dstate_rightdTOF(9, 0) - state_right(9).getDerivative(1)) / state_right(9).getDerivative(1) << std::endl;
   
    getchar();
    

}