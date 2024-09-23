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

#include "SpiceUsr.h" 
#include "missionoptions.h"
#include "EMTG_Matrix.h"
#include "Kepler_Lagrange_Laguerre_Conway_Der.h"
#include "STM.h"

#include "file_utilities.h"
#include "universe.h"

#include "rk8_fixed.h"
#include "FBLT_EOM.h"

#include "ConstructHardware.h"

#include "boost/random/uniform_real.hpp"
#include "boost/random/mersenne_twister.hpp"

#include "doubleType.h"

//#define FiniteDifference
//#define AD_INSTRUMENTATION

// needed for GSAD 4B
//size_t GSAD::adouble::reserve_size = 10;

// needed for GSAD A4
GSAD::adouble GSAD::adouble::temp;
std::vector<size_t>::size_type GSAD::adouble::point = 0;
//GSAD::adouble::derivative_t GSAD::adouble::pair_temp = std::make_pair(0, 0);
//size_t GSAD::adouble::reserve_size = 10;


int main(int argc, char* argv[])
{

    int GSADindex = 0;
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

    // configure the LaunchVehicleOptions and SpacecraftOptions objects
    EMTG::HardwareModels::LaunchVehicleOptions myLaunchVehicleOptions;
    EMTG::HardwareModels::SpacecraftOptions mySpacecraftOptions;
    ConstructHardware(options, myLaunchVehicleOptions, mySpacecraftOptions);

    // create spacecraft and launch vehicle objects
#ifdef AD_INSTRUMENTATION
    EMTG::HardwareModels::LaunchVehicle myLaunchVehicle(myLaunchVehicleOptions);
    EMTG::HardwareModels::Spacecraft mySpacecraft(mySpacecraftOptions);
#else
    EMTG::HardwareModels::LaunchVehicle<double> myLaunchVehicle(myLaunchVehicleOptions);
    EMTG::HardwareModels::Spacecraft<double> mySpacecraft(mySpacecraftOptions);
#endif

    //add an Earth-Mars sequence
    options.Journeys.front().sequence.push_back(3);
    options.Journeys.front().sequence.push_back(4);

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
            if (options.Journeys[j].destination_list[0] > 0)
                body_index_array.push_back(options.Journeys[j].destination_list[0] - 1);
            if (options.Journeys[j].destination_list[1] > 0)
                body_index_array.push_back(options.Journeys[j].destination_list[1] - 1);
            for (size_t b = 0; b < options.Journeys[j].sequence_input.size(); ++b)
                if (options.Journeys[j].sequence_input[b] > 0)
                    body_index_array.push_back(options.Journeys[j].sequence_input[b] - 1);
            for (size_t b = 0; b < TheUniverse[j].perturbation_menu.size(); ++b)
                body_index_array.push_back(TheUniverse[j].perturbation_menu[b]);

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
                        100,
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
            100000.0 * 86400.0);
    }
#endif

    // double mu_normalized = 1.0;
    // double LU = TheUniverse.front().LU;
    // double TU = TheUniverse.front().TU;
    double mu_normalized = TheUniverse[0].mu;
    double LU = 1.0;
    double TU = 1.0;

    // this determines the size of the STMs as well as the length of the state vector
    int statecount = 7;
#ifdef AD_INSTRUMENTATION
    int ns = statecount + 10 * 10;
    //int ns = statecount;
    bool computeDerivatives = true;
#else    
    int ns = statecount + 10 * 10;
    bool computeDerivatives = true;
#endif
    int STMrows = 10;
    int STMcolumns = 10;


    
    size_t numSegments = 50;
    double specified_launch_epoch = 56774.0; // MJD past J2000
    double previous_phase_flight_time = 0.5; // days
    double current_phase_flight_time = 10.0; // days
    double max_rk_step_size = 3600.0; // seconds
    
    size_t derivIndex = 0;

#ifdef AD_INSTRUMENTATION
    GSAD::adouble launch_epoch;
    GSAD::adouble TOFprevious;
    GSAD::adouble TOF;
    GSAD::adouble t_left;

    // seed the derivative directions for the phase TOF variables
    launch_epoch.setDerivative(derivIndex, 1.0);
    TOFprevious.setDerivative(derivIndex + 1, 1.0);
    TOF.setDerivative(derivIndex + 2, 1.0);
    derivIndex = derivIndex + 3;

    // set values for the time of flight variables
    launch_epoch.setValue((specified_launch_epoch * 86400.0));
    TOFprevious.setValue((previous_phase_flight_time * 86400.0));
    TOF.setValue((current_phase_flight_time * 86400.0));
#else
    double launch_epoch;
    double TOFprevious;
    double TOF;
    double t_left;

    launch_epoch = (specified_launch_epoch * 86400.0);
    TOFprevious = (previous_phase_flight_time * 86400.0);
    TOF = (current_phase_flight_time * 86400.0);
#endif
    t_left = launch_epoch + TOFprevious;

#ifdef AD_INSTRUMENTATION
    GSAD::adouble float_rk_steps_per_segment = TOF / numSegments / max_rk_step_size;
    double wholeStepsPerSegment = floor(float_rk_steps_per_segment _GETVALUE);
    GSAD::adouble segment_final_step_fraction = float_rk_steps_per_segment - wholeStepsPerSegment;
#else
    double float_rk_steps_per_segment = TOF / numSegments / max_rk_step_size;
    double wholeStepsPerSegment = floor(float_rk_steps_per_segment);
    double segment_final_step_fraction = float_rk_steps_per_segment - wholeStepsPerSegment;
#endif
    size_t rk_steps_per_segment = (ceil(float_rk_steps_per_segment)) _GETVALUE;
    
#ifdef AD_INSTRUMENTATION
    std::vector< GSAD::adouble > rk_step_lengths(rk_steps_per_segment, max_rk_step_size); // seconds
#else
    std::vector<double> rk_step_lengths(rk_steps_per_segment, max_rk_step_size); // seconds
#endif
    std::vector<double> drk_step_lengthsdTOF(rk_steps_per_segment, 0.0);

    for (size_t k = 0; k < rk_step_lengths.size(); ++k) {
        if (k % (rk_steps_per_segment - 1) == 0 && k != 0) {
            rk_step_lengths[k] = segment_final_step_fraction*max_rk_step_size;
            drk_step_lengthsdTOF[k] = 1.0 / (numSegments);
        }
    }


#ifdef AD_INSTRUMENTATION
    std::vector <GSAD::adouble> state0(ns, 0.0); //initial state
    std::vector <GSAD::adouble> state1(ns, 0.0); //state at the end of the half phase
    std::vector <GSAD::adouble> state_current(ns, 0.0);
    std::vector <GSAD::adouble> state_next(ns, 0.0);

    //set derivative of each state variable with respect to itself equal to 1 (to seed the auto-diff routine)    
    for (size_t i = 0; i < statecount; ++i) {
        state0[i].setDerivative(derivIndex++, 1.0);
    }
#else
    std::vector <double> state0(ns, 0.0); //initial state
    std::vector <double> state1(ns, 0.0); //state at the end of the half phase
    std::vector <double> state_current(ns, 0.0);
    std::vector <double> state_next(ns, 0.0);
#endif
    EMTG::math::Matrix <double> dstate0dTOF(7, 2, 0.0); //initial state partials w.r.t. TOF
    EMTG::math::Matrix <double> dstate1dTOF(7, 2, 0.0); //final state partials w.r.t. TOF
    EMTG::math::Matrix <double> dstate_currentdTOF(7, 2, 0.0); //partials of the left hand state vector w.r.t. TOF
    EMTG::math::Matrix <double> dstate_nextdTOF(7, 2, 0.0); //partials of the right hand state vector w.r.t. TOF




#ifdef AD_INSTRUMENTATION
    std::vector <GSAD::adouble> segmentControls(3 * numSegments, 0.0);
    std::vector< std::vector <GSAD::adouble> > spacecraft_state(numSegments, std::vector<GSAD::adouble>(7, 0.0));
    std::vector <GSAD::adouble> control(3, 0.0);
#else
    std::vector <double> segmentControls(3 * numSegments, 0.0);
    std::vector< std::vector <double> > spacecraft_state(numSegments, std::vector<double>(7, 0.0));
    std::vector <double> control(3, 0.0);
#endif

    boost::mt19937 RNG;
    boost::uniform_real<> DoubleDistribution;
    DoubleDistribution = boost::uniform_real<>(0.0, 0.55);
    RNG.seed(time(NULL));

    //create a vector of STM's.....one STM per FBLT segment
    std::vector < EMTG::math::Matrix<double> > Phi_archive(numSegments, EMTG::math::Matrix<double>(STMrows, STMcolumns, 0.0));
    std::vector < EMTG::math::Matrix<double> > Phi_archive_stripped(numSegments, EMTG::math::Matrix<double>(STMrows, STMcolumns, 0.0));

    //initialize state vector
    double LT_dump;
    std::vector<double> x(7);
    spkez_c(399, t_left _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 10, x.data(), &LT_dump);
    //spkez_c(504, t_left _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 599, x.data(), &LT_dump);
    //spkez_c(499, t_left _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 4, x.data(), &LT_dump);
    //spkez_c(602, t_left _GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", 699, x.data(), &LT_dump);

    //create perturbed initial state for TOF derivative use
    std::vector <double> x_pert(7, 0.0);
    std::vector <double> dx_dTOFprevious(7, 0.0);
    std::vector <double> dx_dTOFcurrent(7, 0.0);
    spkez_c(399, t_left _GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", 10, x_pert.data(), &LT_dump);
    //spkez_c(504, t_left _GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", 599, x_pert.data(), &LT_dump);
    //spkez_c(602, t_left _GETVALUE - (51544.5 * 86400.0) + 10.0, "J2000", "NONE", 699, x_pert.data(), &LT_dump);
    dx_dTOFprevious[0] = (x_pert[0] - x[0]) / 10.0;
    dx_dTOFprevious[1] = (x_pert[1] - x[1]) / 10.0;
    dx_dTOFprevious[2] = (x_pert[2] - x[2]) / 10.0;
    dx_dTOFprevious[3] = (x_pert[3] - x[3]) / 10.0;
    dx_dTOFprevious[4] = (x_pert[4] - x[4]) / 10.0;
    dx_dTOFprevious[5] = (x_pert[5] - x[5]) / 10.0;
    dx_dTOFprevious[6] = 0.0;
    x[6] = options.maximum_mass;

    for (size_t i = 0; i < 7; ++i)
    {
#ifdef AD_INSTRUMENTATION
        state0[i].setValue(x[i]);
        state0[i].setDerivative(2, dx_dTOFcurrent[i]);
        state0[i].setDerivative(1, dx_dTOFprevious[i]);
        state0[i].setDerivative(0, dx_dTOFprevious[i]);
#else
        state0[i] = x[i];
#endif
    }

    // STM is initialized to the identity matrix
    for (size_t i = statecount; i < ns; ++i) {
        state0[i] = 0.0;
    }
    for (size_t i = statecount; i < ns; i = i + STMrows + 1) {
        state0[i] = 1.0;
    }

    // seed derivatives for spacecraft state with respect to current and previous phase flight times
    for (size_t i = 0; i < 3; i++) {
        dstate0dTOF(i, 0) = 0.0;
        dstate0dTOF(i + 3, 0) = 0.0;
        dstate0dTOF(i, 1) = dx_dTOFprevious[i];
        dstate0dTOF(i + 3, 1) = dx_dTOFprevious[i + 3];
    }

    // boundary point mass is unaffected by phase TOFs because it is a decision variable
   


    // Place random throttle values into the master control vector
    std::vector<std::vector<int>> control_active_indices(numSegments, std::vector<int>(3, -1));
    for (size_t segment = 0; segment < numSegments; ++segment)
    {
        for (size_t i = 0; i < 3; ++i)
        {
#ifdef AD_INSTRUMENTATION
            segmentControls[3 * segment + i].setValue(DoubleDistribution(RNG));
            //segmentControls[3 * segment + i].setValue(0.0);
            control_active_indices[segment][i] = derivIndex;
            segmentControls[3 * segment + i].setDerivative(derivIndex++, 1.0);
#else
            //segmentControls[3 * segment + i] = DoubleDistribution(RNG);
            segmentControls[3 * segment + i] = 0.0;
#endif
        }
    }

    //seed for derivatives of current epoch with respect to current and previous phase flight times
    std::vector <double> dt_leftdTOF(2, 0.0);
    dt_leftdTOF[0] = 0.0;
    dt_leftdTOF[1] = 1.0;

    //seed for derivative of an FBLT segment length with respect to the current phase flight time
    //previous phase flight times have no influence on this

#ifdef AD_INSTRUMENTATION
    GSAD::adouble segmentTime = TOF / numSegments;
#else
    double segmentTime = TOF / numSegments;
#endif
    double dsegmentTimedTOF = 1.0 / numSegments;

    //dummy variables
#ifdef AD_INSTRUMENTATION
    GSAD::adouble available_thrust;
    GSAD::adouble mdot;
    GSAD::adouble Isp;
    GSAD::adouble power;
    GSAD::adouble active_power;
    GSAD::adouble u_command = 0.0;
#else
    double available_thrust;
    double mdot;
    double Isp;
    double power;
    double active_power;
    double u_command = 0.0;
#endif

    int number_of_active_engines;
    void * Universepointer;
    void * Controllerpointer;

    //instantiate an integrator
#ifdef AD_INSTRUMENTATION
    EMTG::integration::rk8 * integrator;
    integrator = new EMTG::integration::rk8 (ns);
#else
    EMTG::integration::rk8<double> * integrator;
    integrator = new EMTG::integration::rk8<double>(ns);
#endif

    state_current = state0;
    dstate_currentdTOF = dstate0dTOF;
    double elapsedTime = 0.0;

    size_t num_states_to_integrate = ns;

    auto start = std::chrono::steady_clock::now();
    //FBLT propagation loop
    for (size_t segment = 0; segment < numSegments; ++segment)
    {
        control[0] = segmentControls[3 * segment];
        control[1] = segmentControls[3 * segment + 1];
        control[2] = segmentControls[3 * segment + 2];

        //state_current[7] = segmentControls[3 * segment];
        //state_current[8] = segmentControls[3 * segment + 1];
        //state_current[9] = segmentControls[3 * segment + 2];
        

       
        // propagate for the number of RK steps required to get to the next segment
        // we need to take one extra partial RK step to hit the segment boundary exactly
        // (hence the +1)
        // for (size_t k = 0; k < wholeStepsPerSegment + 1; ++k) {
        for (size_t k = 0; k < rk_step_lengths.size(); ++k) {

            integrator->rk8_step(
                state_current,
                dstate_currentdTOF,
                state_next,
                dstate_nextdTOF,
                control,
                u_command,
                options.engine_duty_cycle,
                t_left,
                dt_leftdTOF,
                launch_epoch,
                rk_step_lengths[k],
                drk_step_lengthsdTOF[k],
                EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                available_thrust,
                mdot,
                Isp,
                power,
                active_power,
                number_of_active_engines,
                STMrows,
                STMcolumns,
                options,
                TheUniverse[0],
                &mySpacecraft,
                computeDerivatives,
                0,//j
                0,//p
                0,//step
                mu_normalized,
                LU,
                TU,
                1.0,
                num_states_to_integrate);

            t_left += rk_step_lengths[k]; //advance the mission clock
            
            //state at the end of the segment becomes the new current state
            state_current = state_next;
            //state TOF partial derivatives at the right of the current segment become the partials for the left of the next 
            dstate_currentdTOF = dstate_nextdTOF;
            elapsedTime += rk_step_lengths[k] _GETVALUE;
            std::cout << "Time: " << elapsedTime / 86400.0 << " days" << std::endl;
        }

        // for fixed step, only the epochs of the segment boundaries
        // shift with a change in TOF, not the boundaries of the individual
        // RK steps
        //dt_leftdTOF[0] += dsegmentTimedTOF;

        spacecraft_state[segment] = state_next;
        
        if (computeDerivatives) {
            int index = statecount;
            for (size_t i = 0; i < STMrows; ++i)
            {
                for (size_t j = 0; j < STMcolumns; ++j)
                {
                    Phi_archive[segment](i, j) = state_next[index] _GETVALUE;
                    ++index;
                }
            }

            //STM for the next step is initialized to the identity matrix
            for (size_t i = statecount; i < ns; ++i)
                state_current[i] = 0.0;

            for (size_t i = statecount; i < ns; i = i + STMrows + 1)
                state_current[i] = 1.0;
        }

        std::cout << "Step " << segment << std::endl;

    } //end FBLT propagation loop
    auto end = std::chrono::steady_clock::now();
    auto execution_time = end - start;

    state1 = state_next;
    dstate1dTOF = dstate_nextdTOF;

    EMTG::math::Matrix <double> derivatives(STMrows, STMcolumns, 0.0);
    EMTG::math::Matrix <double> derivatives_stripped(STMrows, STMcolumns, 0.0);

    //initialize the derivatives matrix to the identity
    for (size_t i = 0; i < STMrows; ++i)
    {
        derivatives(i, i) = 1.0;
        derivatives_stripped(i, i) = 1.0;
    }

    Phi_archive_stripped = Phi_archive;
    for (int segment = numSegments - 1; segment >= 0; --segment)
    {
        for (size_t row = 0; row < 7; ++row)
        {
            for (size_t column = 7; column < 10; ++column)
            {
                Phi_archive_stripped[segment](row, column) = 0.0;
            }
        }
    }


    //build the derivatives matrix through successive STM multiplication
    //compute the cumulative STM
    for (int segment = numSegments - 1; segment > 0; --segment)
    {
        derivatives_stripped *= Phi_archive_stripped[segment];
    }

    //multiply in the first step's STM
    derivatives = derivatives_stripped * Phi_archive[0];


    std::cout << "Numerically integrated resultant STM:" << std::endl;
    std::cout << std::setprecision(16);
    for (size_t i = 0; i < STMrows; ++i)
    {
        for (size_t j = 0; j < STMcolumns; ++j)
            std::cout << derivatives(i, j) << "     ";

        std::cout << std::endl << std::endl;
    }

    std::cout << std::endl << std::endl;

#ifdef AD_INSTRUMENTATION
    //AD calculated derivatives
    EMTG::math::Matrix<double> PhiAD(STMrows, STMcolumns, 0.0);
    EMTG::math::Matrix<double> PhiERROR(STMrows, STMcolumns, 0.0);
    EMTG::math::Matrix<double> PhiERRORrel(STMrows, STMcolumns, 0.0);

    for (size_t i = 0; i < statecount; ++i)
    {
        for (size_t j = 0; j < statecount; ++j)
        {
            PhiAD(i, j) = state1[i].getDerivative(j+3);
        }
        PhiAD(i, 7) = state1[i].getDerivative(10);
        PhiAD(i, 8) = state1[i].getDerivative(11);
        PhiAD(i, 9) = state1[i].getDerivative(12);
    }

    for (size_t i = 0; i < STMrows; ++i)
        for (size_t j = 0; j < STMcolumns; ++j)
            PhiERROR(i, j) = PhiAD(i, j) - derivatives(i, j);

    for (size_t i = 0; i < STMrows; ++i)
        for (size_t j = 0; j < STMcolumns; ++j)
            PhiERRORrel(i, j) = (PhiAD(i, j) - derivatives(i, j)) / PhiAD(i, j);

    std::cout << "AD calculated resultant STM:" << std::endl;
    for (size_t i = 0; i < STMrows; ++i)
    {
        for (size_t j = 0; j < STMcolumns; ++j)
        {
            std::cout << PhiAD(i, j) << "     ";
        }

        std::cout << std::endl << std::endl;
    }

    std::cout << std::endl;

    std::cout << "Error between AD STM and integrated STM:" << std::endl;
    for (size_t i = 0; i < STMrows; ++i)
    {
        for (size_t j = 0; j < STMcolumns; ++j)
        {
            std::cout << PhiERROR(i, j) << "     ";
        }

        std::cout << std::endl << std::endl;
    }

    std::cout << std::endl;

    std::cout << "Relative error between AD STM and integrated STM:" << std::endl;
    for (size_t i = 0; i < STMrows; ++i)
    {
        for (size_t j = 0; j < STMcolumns; ++j)
        {
            std::cout << PhiERRORrel(i, j) << "     ";
        }

        std::cout << std::endl << std::endl;
    }

    std::cout << std::endl;

    std::cout << "Check final states" << std::endl;
    std::cout << "x: " << state1[0] << std::endl;
    std::cout << "y: " << state1[1] << std::endl;
    std::cout << "z: " << state1[2] << std::endl;
    std::cout << "vx: " << state1[3] << std::endl;
    std::cout << "vy: " << state1[4] << std::endl;
    std::cout << "vz: " << state1[5] << std::endl;
    std::cout << "m: " << state1[6] << std::endl;

    std::cout << std::endl;

    std::cout << "Time of Flight derivatives" << std::endl;
    std::cout << "Format: analytical, GSAD, absolute error, relative error" << std::endl;
    std::cout << "dxdTOFcurrent: " << dstate1dTOF(0, 0) << " " << state1[0].getDerivative(2) << " " << dstate1dTOF(0, 0) - state1[0].getDerivative(2)   << " " << (dstate1dTOF(0, 0) - state1[0].getDerivative(2)) / state1[0].getDerivative(2) << std::endl;
    std::cout << "dxdTOFprevious: " << dstate1dTOF(0, 1) << " " << state1[0].getDerivative(1) << " " << dstate1dTOF(0, 1) - state1[0].getDerivative(1)  << " " << (dstate1dTOF(0, 1) - state1[0].getDerivative(1)) / state1[0].getDerivative(1) << std::endl;
    std::cout << "dydTOFcurrent: " << dstate1dTOF(1, 0) << " " << state1[1].getDerivative(2) << " " << dstate1dTOF(1, 0) - state1[1].getDerivative(2)   << " " << (dstate1dTOF(1, 0) - state1[1].getDerivative(2)) / state1[1].getDerivative(2) << std::endl;
    std::cout << "dydTOFprevious: " << dstate1dTOF(1, 1) << " " << state1[1].getDerivative(1) << " " << dstate1dTOF(1, 1) - state1[1].getDerivative(1)  << " " << (dstate1dTOF(1, 1) - state1[1].getDerivative(1)) / state1[1].getDerivative(1) << std::endl;
    std::cout << "dzdTOFcurrent: " << dstate1dTOF(2, 0) << " " << state1[2].getDerivative(2) << " " << dstate1dTOF(2, 0) - state1[2].getDerivative(2)   << " " << (dstate1dTOF(2, 0) - state1[2].getDerivative(2)) / state1[2].getDerivative(2) << std::endl;
    std::cout << "dzdTOFprevious: " << dstate1dTOF(2, 1) << " " << state1[2].getDerivative(1) << " " << dstate1dTOF(2, 1) - state1[2].getDerivative(1)  << " " << (dstate1dTOF(2, 1) - state1[2].getDerivative(1)) / state1[2].getDerivative(1) << std::endl;
    std::cout << "dvxdTOFcurrent: " << dstate1dTOF(3, 0) << " " << state1[3].getDerivative(2) << " " << dstate1dTOF(3, 0) - state1[3].getDerivative(2)  << " " << (dstate1dTOF(3, 0) - state1[3].getDerivative(2)) / state1[3].getDerivative(2) << std::endl;
    std::cout << "dvxdTOFprevious: " << dstate1dTOF(3, 1) << " " << state1[3].getDerivative(1) << " " << dstate1dTOF(3, 1) - state1[3].getDerivative(1) << " " << (dstate1dTOF(3, 1) - state1[3].getDerivative(1)) / state1[3].getDerivative(1) << std::endl;
    std::cout << "dvydTOFcurrent: " << dstate1dTOF(4, 0) << " " << state1[4].getDerivative(2) << " " << dstate1dTOF(4, 0) - state1[4].getDerivative(2)  << " " << (dstate1dTOF(4, 0) - state1[4].getDerivative(2)) / state1[4].getDerivative(2) << std::endl;
    std::cout << "dvydTOFprevious: " << dstate1dTOF(4, 1) << " " << state1[4].getDerivative(1) << " " << dstate1dTOF(4, 1) - state1[4].getDerivative(1) << " " << (dstate1dTOF(4, 1) - state1[4].getDerivative(1)) / state1[4].getDerivative(1) << std::endl;
    std::cout << "dvzdTOFcurrent: " << dstate1dTOF(5, 0) << " " << state1[5].getDerivative(2) << " " << dstate1dTOF(5, 0) - state1[5].getDerivative(2)  << " " << (dstate1dTOF(5, 0) - state1[5].getDerivative(2)) / state1[5].getDerivative(2) << std::endl;
    std::cout << "dvzdTOFprevious: " << dstate1dTOF(5, 1) << " " << state1[5].getDerivative(1) << " " << dstate1dTOF(5, 1) - state1[5].getDerivative(1) << " " << (dstate1dTOF(5, 1) - state1[5].getDerivative(1)) / state1[5].getDerivative(1) << std::endl;
    std::cout << "dmdTOFcurrent: " << dstate1dTOF(6, 0) << " " << state1[6].getDerivative(2) << " " << dstate1dTOF(6, 0) - state1[6].getDerivative(2)   << " " << (dstate1dTOF(6, 0) - state1[6].getDerivative(2)) / state1[6].getDerivative(2) << std::endl;
    std::cout << "dmdTOFprevious: " << dstate1dTOF(6, 1) << " " << state1[6].getDerivative(1) << " " << dstate1dTOF(6, 1) - state1[6].getDerivative(1)  << " " << (dstate1dTOF(6, 1) - state1[6].getDerivative(1)) / state1[6].getDerivative(1) << std::endl;

    (PhiAD.element_divide(derivatives)).print_to_file("phi_ratio.txt");
    PhiAD.print_to_file("PhiAD.txt");
    PhiERROR.print_to_file("PhiERROR.txt");
    PhiERRORrel.print_to_file("PhiERRORrel.txt");
#endif
    std::cout << "Propagation wallclock time: " << std::chrono::duration <double, std::milli>(execution_time).count() << " ms" << std::endl;
    getchar();
    delete integrator;
}