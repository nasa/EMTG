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

//header file for NLP options class
#pragma once

#include "problem.h"
#include "EMTG_enums.h"

#include <string.h>

#ifndef PORTABLE_SOLVER
#include "missionoptions.h"
#endif

namespace EMTG
{
    namespace Solvers
    {
        class NLPoptions
        {

        public:
            //constructor
            NLPoptions();

#ifndef PORTABLE_SOLVER
            NLPoptions(const missionoptions& options);
#endif

            //get/set
            NLPMode get_SolverMode() { return this->SolverMode; }
            bool get_enable_auto_scale() { return this->enable_auto_scale; }
            bool get_enable_NLP_chaperone() { return this->enable_NLP_chaperone; }
            bool get_check_derivatives() { return this->check_derivatives; }
            bool get_quiet_NLP() { return this->quiet_NLP; }
            bool get_print_NLP_movie_frames() { return this->print_NLP_movie_frames; }
            bool get_stop_on_goal_attain() { return this->stop_on_goal_attain; }

            size_t get_major_iterations_limit() { return this->major_iterations_limit; }
            size_t get_minor_iterations_limit() { return this->minor_iterations_limit; }
            double get_max_step() { return this->max_step; }
            time_t get_max_run_time_seconds() { return this->max_run_time_seconds; }

            double get_feasibility_tolerance() { return this->feasibility_tolerance; }
            double get_optimality_tolerance() { return this->optimality_tolerance; }

            double get_objective_goal() { return this->objective_goal; }

            std::string get_specs_file_path() { return this->specs_file_path; }
            std::string get_output_file_path() { return this->output_file_path; }

            void set_SolverMode(const NLPMode& SolverMode) { this->SolverMode = SolverMode; }
            void set_enable_auto_scale(const bool& enable_auto_scale) { this->enable_auto_scale = enable_auto_scale; }
            void set_enable_NLP_chaperone(const bool& enable_NLP_chaperone) { this->enable_NLP_chaperone = enable_NLP_chaperone; }
            void set_check_derivatives(const bool& check_derivatives) { this->check_derivatives = check_derivatives; }
            void set_print_NLP_movie_frames(const bool& print_NLP_movie_frames) { this->print_NLP_movie_frames = print_NLP_movie_frames; }
            void set_major_iterations_limit(const size_t& major_iterations_limit) { this->major_iterations_limit = major_iterations_limit; }
            void set_minor_iterations_limit(const size_t& minor_iterations_limit) { this->minor_iterations_limit = minor_iterations_limit; }
            void set_max_run_time_seconds(const time_t& max_run_time_seconds) { this->max_run_time_seconds = max_run_time_seconds; }
            void set_max_step(const double& max_step) { this->max_step = max_step; }
            void set_feasibility_tolerance(const double& feasibility_tolerance) { this->feasibility_tolerance = feasibility_tolerance; }
            void set_optimality_tolerance(const double& optimality_tolerance) { this->optimality_tolerance = optimality_tolerance; }
            void set_specs_file_path(const std::string& specs_file_path) { this->specs_file_path = specs_file_path; }
            void set_output_file_path(const std::string& output_file_path) { this->output_file_path = output_file_path; }

        protected:
            //fields
            NLPMode SolverMode;
            bool enable_auto_scale;
            bool enable_NLP_chaperone;
            bool check_derivatives;
            bool quiet_NLP;
            bool print_NLP_movie_frames;
            bool stop_on_goal_attain;

            size_t major_iterations_limit;
            size_t minor_iterations_limit;
            time_t max_run_time_seconds;

            double max_step;
            double feasibility_tolerance;
            double optimality_tolerance;
            double objective_goal;
            std::string specs_file_path;
            std::string output_file_path;            
        };//end class NLPOptions
    }//end namespace Solvers
}//end namespace EMTG