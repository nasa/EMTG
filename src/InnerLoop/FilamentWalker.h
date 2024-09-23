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

//header file for filament walker
//math by Arnold Englander
//code by Jacob Englander

#pragma once

#include "problem.h"
#include "EMTG_enums.h"
#include "SNOPT_interface.h"
#include "EMTG_innerloop_solution.h"

#include <vector>
#include <string>
#include <tuple>
#include <random>

namespace EMTG
{
    namespace Solvers
    {
        class FilamentWalker
        {
        public:
            //constructors
            FilamentWalker();
            FilamentWalker(problem* myProblem,
                           SNOPT_interface* mySNOPT);

            //get/set

            //execute
            void walk();


        private:
            //internal methods
            void output(std::tuple< EMTG_innerloop_solution, double, double, double, time_t >& solution);
            bool remedialStart(); //return true if succeeded, false if failed
            void remedialHop(bool& firstHop);
            bool routineStart(); //return true if succeeded, false if failed
            void routineHop();
            bool slide(); //return true if succeeded, false if failed


            //fields
            problem* myProblem;
            SNOPT_interface* mySNOPT;
            std::ofstream outputfile;
            size_t nX;
            bool WalkSoftly;

            //random number generator
            std::mt19937 RNG;
            size_t myRandomSeed;
            std::uniform_real_distribution<> DoubleDistribution;

            //timers and counters
            time_t FilamentWalker_start_time;
            size_t current_remedial_start_index;
            size_t current_routine_start_index;
            size_t routine_start_failures_since_reversal;
            bool reversal_occurred;
            
            //tuning parameters
            double pareto_alpha;
            size_t max_remedial_starts;
            size_t max_routine_start_attempts;
            size_t max_routine_start_failures_before_reversal;
            size_t max_slide_steps;
            double max_routine_start_hop_size;
            double equality_constraint_tolerance;

            //container for where we currently live
            std::vector<double> X;
            std::vector<double> X_before_hop;
            std::vector<double> X_after_hop;
            std::vector<double> X_after_slide;
            std::vector<double> dX_hop;
            std::vector<double> dX_slide;
            std::vector<double> dX_iteration;

            //lists
            //in each case we store a tuple:
            //<0> the solution itself
            //<1> the distance since the previous solution (negative for remedial start)
            //<2> the distance hopped (negative for remedial start)
            //<3> the distance slid (negative for remedial start)
            //<4> the time in seconds since we started
            //2D list of solutions found, indexed by remedial start and then by routine start
            std::vector< std::vector< std::tuple< EMTG_innerloop_solution, double, double, double, time_t > > > list_of_filaments;
            //1D list of solutions found in the current remedial start, which we hope means the current filament
            std::vector< std::tuple< EMTG_innerloop_solution, double, double, double, time_t > > list_of_points_on_current_filament;
        };
    }//close namespace solvers
}//close namespace EMTG