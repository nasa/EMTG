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

//header file for Monotonic Basin Hopping
//for EMTGv9
//Jacob Englander 6-22-2017

#include "problem.h"

#include "randutils.hpp"

#include "SNOPT_interface.h"
#include "EMTG_innerloop_solution.h"

#pragma once

namespace EMTG 
{
    namespace Solvers
    {

        class MBH
        {
        public:
            //constructor
            MBH();
            MBH(EMTG::problem* myProblem,
                SNOPT_interface* mySNOPT);

            //destructor
            virtual ~MBH() {};

            //methods
            void initialize(EMTG::problem* Problem_input,
                SNOPT_interface* mySNOPT);
            void reset_point();
            void seed(std::vector<double>& seed_vector);
            void slide();
            int run();

            double getJGlobalIncumbent() const { return this->JGlobalIncumbent; }
            void setJGlobalIncumbent(const double& JGlobalIncumbent) { this->JGlobalIncumbent = JGlobalIncumbent; }

            std::vector<double> getX_most_feasible() const { return this->X_most_feasible; }
    
        private:
            void hop();
            void time_hop();
            void print_archive_header(const std::string& filename);
            void print_archive_line(const std::string& filename, std::tuple< EMTG_innerloop_solution, time_t > solution);

            //fields

            //pointer to problem object
            EMTG::problem* myProblem;

            //archive of feasible solutions
            std::vector< std::tuple< EMTG_innerloop_solution, time_t > > archive;

            //counters
            time_t MBH_start_time;

            //container for where we currently live
            std::vector<double> X_incumbent;
            std::vector<double> X_global_incumbent;
            std::vector<double> X_after_hop;
            std::vector<double> X_after_hop_unscaled;
            std::vector<double> X_after_slide;
            std::vector<double> X_after_slide_unscaled;
            std::vector<double> X_most_feasible;
            std::vector<double> F_after_slide;
            

            double Jincumbent;
            double JGlobalIncumbent;
            double FeasibilityIncumbent;
            double FeasibilityGlobalIncumbent;

            //helper arrays
            std::vector<int> time_variable_indices;
            std::vector<int> significant_variable_indices;

            //counters
            size_t number_of_solutions; //how many feasible solutions found so far
            size_t number_of_improvements; //how many times have we improved within a basin
            size_t number_of_resets;
            size_t number_of_failures_since_last_improvement;
            size_t number_of_attempts;

            //track the worst constraint violation
            size_t worst_constraint;
            double max_constraint_violation;

            //track the step size
            double step_size;

            //track whether or not the sparsity file and XFfile have been printed
            bool printed_sparsity;

            //track whether or not the feasible point finder is active
            bool feasible_point_finder_active;

            //random number generator
            randutils::mt19937_rng RNG;

            //SNOPT object
            SNOPT_interface* mySNOPT;
        };


    }//close namespace Solvers
}//close namespace EMTG