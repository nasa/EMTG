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

//Monotonic Basin Hopping
//for EMTG version 8
//Jacob Englander 7-27-2012
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <math.h>

#include "problem.h"
#include "monotonic_basin_hopping.h"
#include "EMTG_math.h"

#include "snoptProblemExtension.h"

namespace EMTG { namespace Solvers {
    //constructors
    MBH::MBH() {}

    MBH::MBH(EMTG::problem* myProblem,
        SNOPT_interface* mySNOPT)
    {
        //initialize the MBH variables
        initialize(myProblem, mySNOPT);

        //search through the problem object and identify which decision variables are flight time variables
        if (this->myProblem->options.MBH_time_hop_probability > 0.0)
        {
            for (int entry = 0; entry < this->myProblem->total_number_of_NLP_parameters; ++entry)
                if (this->myProblem->Xdescriptions[entry].find("flight time") < 1024)
                {
                    this->time_variable_indices.push_back(entry);
                }
        }

        //search through the problem object and identify which decision variables are significant (i.e. non-control)
        if (this->myProblem->options.MBH_time_hop_probability > 0.0)
        {
            for (int entry = 0; entry < this->myProblem->total_number_of_NLP_parameters; ++entry)
                if ( !(myProblem->Xdescriptions[entry].find("u_x") < 1024 
                    || myProblem->Xdescriptions[entry].find("u_y") < 1024 
                    || myProblem->Xdescriptions[entry].find("u_z") < 1024 
                    || myProblem->Xdescriptions[entry].find("u_command") < 1024) )
                {
                    this->significant_variable_indices.push_back(entry);
                }
        }
    }

    
    //method to initialize the MBH solver
    //resets all storage fields
    void MBH::initialize(EMTG::problem* myProblem,
        SNOPT_interface* mySNOPT)
    {
        this->myProblem = myProblem;
        this->mySNOPT = mySNOPT;
        
        //size the storage vectors
        this->archive.clear();
        this->X_incumbent.resize(this->myProblem->total_number_of_NLP_parameters, 0.0);
        this->X_after_hop.resize(this->myProblem->total_number_of_NLP_parameters, 0.0);
        this->X_after_hop_unscaled.resize(this->myProblem->total_number_of_NLP_parameters, 0.0);
        this->X_after_slide.resize(this->myProblem->total_number_of_NLP_parameters, 0.0);
        this->X_after_slide_unscaled.resize(this->myProblem->total_number_of_NLP_parameters, 0.0);

        //clear the scores
        this->Jincumbent = EMTG::math::LARGE;
        this->JGlobalIncumbent = EMTG::math::LARGE;
        this->FeasibilityIncumbent = EMTG::math::LARGE;

        //reset the Jacobian printed flag
        this->printed_sparsity = false;

        //seed the random number generator based on the node's clock
        if (this->myProblem->options.MBH_RNG_seed < 0)
            this->RNG.seed();
        else
            this->RNG.seed(this->myProblem->options.MBH_RNG_seed);
    }

    //function to reset to a new point, either at the beginning of the myProblem or when a certain number of failures is reached
    void MBH::reset_point()
    {
        //reset the number of failures
        this->number_of_failures_since_last_improvement = 0;

        //increment the number of global resets
        ++this->number_of_resets;

        //turn on the feasible point finder if it is enabled
        if (this->myProblem->options.ACE_feasible_point_finder)
            this->feasible_point_finder_active = true;
        else
            this->feasible_point_finder_active = false;

        //generate a new random trial point
        for (int k = 0; k < myProblem->total_number_of_NLP_parameters; ++k)
            this->X_after_hop[k] =
                (this->RNG.uniform(0.0, 1.0) * (this->myProblem->Xupperbounds[k] - this->myProblem->Xlowerbounds[k])
                    + this->myProblem->Xlowerbounds[k])
                    / this->myProblem->X_scale_factors[k];

        this->Jincumbent = EMTG::math::LARGE;

        if (this->number_of_resets == 1)
            this->X_most_feasible = this->X_after_hop;

		for (size_t Xindex = 0; Xindex < this->myProblem->total_number_of_NLP_parameters; ++Xindex)
			this->X_after_hop_unscaled[Xindex] = this->X_after_hop[Xindex] * this->myProblem->X_scale_factors[Xindex];
    }

    //function to seed a new point based on an input trial point
    void MBH::seed(std::vector<double>& seed_vector)
    {
        //reset the number of failures
        number_of_failures_since_last_improvement = 0;

        //read a trial point
        if (seed_vector.size() > this->myProblem->total_number_of_NLP_parameters)
        {
            std::cout << "Seed vector is longer than myProblem decision vector. Truncating." << std::endl;
        }
        for (size_t k = 0; k < std::min(myProblem->Xlowerbounds.size(), seed_vector.size()); ++k)
            this->X_after_hop[k] = seed_vector[k] / this->myProblem->X_scale_factors[k];
        if (seed_vector.size() < myProblem->total_number_of_NLP_parameters)
        {
            std::cout << "Seed vector is shorter than myProblem decision vector. Extending with random values." << std::endl;
            for (size_t k = seed_vector.size(); k < myProblem->total_number_of_NLP_parameters; ++k)
                this->X_after_hop[k] = this->RNG.uniform(0.0, 1.0) * ((this->myProblem->Xupperbounds[k] - this->myProblem->Xlowerbounds[k])
                    + this->myProblem->Xlowerbounds[k])
                        / this->myProblem->X_scale_factors[k];
        }
        
        for (size_t Xindex = 0; Xindex < this->myProblem->total_number_of_NLP_parameters; ++Xindex)
            this->X_after_hop_unscaled[Xindex] = this->X_after_hop[Xindex] * this->myProblem->X_scale_factors[Xindex];

        try
        {
            this->myProblem->evaluate(this->X_after_hop_unscaled, this->myProblem->F, this->myProblem->G, true);
        }
        catch (std::exception &error)
        {
            if (!myProblem->options.quiet_basinhopping)
            {
                std::cout << "EMTG::Invalid initial point or failure in objective function while evaluating initial point." << std::endl;
                std::cout << error.what() << std::endl;
            }
        }
        this->Jincumbent = EMTG::math::LARGE;
    }

    //function to perform a "hop" operation
    void MBH::hop()
    {

        if (myProblem->options.MBH_hop_distribution == 0)
        {
            //perform a uniform "hop"
            for (size_t k = 0; k < myProblem->total_number_of_NLP_parameters; ++k)
            {
                double r = this->RNG.uniform(0.0, 1.0);

                this->step_size = 2 * (r - 0.5) * myProblem->options.MBH_max_step_size;

                this->X_after_hop[k] = this->X_incumbent[k] + this->step_size
                    * (this->myProblem->Xupperbounds[k] - this->myProblem->Xlowerbounds[k]) / this->myProblem->X_scale_factors[k];
            }
        }
        else if (myProblem->options.MBH_hop_distribution == 1)
        {
            //perform a Cauchy "hop"            
            for (size_t k = 0; k < myProblem->total_number_of_NLP_parameters; ++k)
            {
                //Cauchy generator v2 as per instructions from Arnold Englander 1-12-2014 using Cauchy CDF

                double r = 0.0;

                while (r < math::SMALL)
                    r = this->RNG.uniform(0.0, 1.0);

                this->step_size = myProblem->options.MBH_max_step_size * tan(math::PI*(r - 0.5));

                this->X_after_hop[k] = this->X_incumbent[k] + this->step_size
                    * (this->myProblem->Xupperbounds[k] - this->myProblem->Xlowerbounds[k]) / this->myProblem->X_scale_factors[k];
            }
        }
        else if (myProblem->options.MBH_hop_distribution == 2)
        {
            double alpha = myProblem->options.MBH_Pareto_alpha;
            //perform a Pareto hop
            
            for (size_t k = 0; k < myProblem->total_number_of_NLP_parameters; ++k)
            {
                //Pareto distribution from x = 1.0 with alpha, modified to start from x = 0.0

                //s*(((alpha-1)/xU_min)/power((xU_min/u),-alpha))
                double r = ( (alpha - 1.0) / math::SMALL ) / pow((math::SMALL / (math::SMALL + this->RNG.uniform(0.0, 1.0))), -alpha);
                int s = this->RNG.uniform(0.0, 1.0) > 0.5 ? 1 : -1;

                this->step_size = s * myProblem->options.MBH_max_step_size * r;

                if (isfinite(step_size))
                    this->X_after_hop[k] = this->X_incumbent[k] + this->step_size
                    * (this->myProblem->Xupperbounds[k] - this->myProblem->Xlowerbounds[k]) / this->myProblem->X_scale_factors[k];
            }
        }
        else if (myProblem->options.MBH_hop_distribution == 3)
        {
            //perform a Gaussian hop
            double sigma = myProblem->options.MBH_max_step_size;
            double sigma2 = sigma * sigma;
            
            for (size_t k = 0; k < myProblem->total_number_of_NLP_parameters; ++k)
            {
                double r = this->RNG.uniform(0.0, 1.0);
                int s = this->RNG.uniform(0.0, 1.0) > 0.5 ? 1 : -1;

                this->step_size = s / (sigma * sqrt(math::TwoPI)) * exp(-r*r / (2*sigma2));

                this->X_after_hop[k] = this->X_incumbent[k] + this->step_size
                    * (this->myProblem->Xupperbounds[k] - this->myProblem->Xlowerbounds[k]) / this->myProblem->X_scale_factors[k];
            }
        }

        //MBH clipping
        for (size_t k = 0; k < this->myProblem->total_number_of_NLP_parameters; ++k)
        {
            double Xmax = this->myProblem->Xupperbounds[k] / this->myProblem->X_scale_factors[k];
            double Xmin = this->myProblem->Xlowerbounds[k] / this->myProblem->X_scale_factors[k];
            if (this->X_after_hop[k] > Xmax)
                this->X_after_hop[k] = Xmax;
            else if (this->X_after_hop[k] < Xmin)
                this->X_after_hop[k] = Xmin;
        }
    }

    //function to perform a "time hop" operation
    void MBH::time_hop()
    {
        //loop through any time variables and if (uniform random < threshold) then add/subtract a synodic period
        for (size_t timeindex = 0; timeindex < time_variable_indices.size(); ++timeindex)
        {
            if (this->RNG.uniform(0.0, 1.0) < myProblem->options.MBH_time_hop_probability)
            {
                int k = time_variable_indices[timeindex];
                int s = this->RNG.uniform(0.0, 1.0) > 0.5 ? 1 : -1;
                this->X_after_hop[k] = this->X_incumbent[k] 
                    + s * myProblem->synodic_periods[timeindex] / this->myProblem->X_scale_factors[k];

                //MBH clipping
                double Xmax = this->myProblem->Xupperbounds[k] / this->myProblem->X_scale_factors[k];
                double Xmin = this->myProblem->Xlowerbounds[k] / this->myProblem->X_scale_factors[k];
                if (this->X_after_hop[k] > Xmax)
                    this->X_after_hop[k] = Xmax;
                else if (this->X_after_hop[k] < Xmin)
                    this->X_after_hop[k] = Xmin;
            }
        }
    }

    //function to perform a "slide" operation, i.e. run SNOPT
    void MBH::slide()
    {
        //Step 1: set the current state equal to the initial guess
        for (size_t Xindex = 0; Xindex < this->myProblem->total_number_of_NLP_parameters; ++Xindex)
            this->X_after_hop_unscaled[Xindex] = this->X_after_hop[Xindex] * this->myProblem->X_scale_factors[Xindex];
        
        this->mySNOPT->setX0_unscaled(this->X_after_hop_unscaled);

        //print the sparsity file and XF files if this is the first pass, otherwise don't to save time and hard drive cycles
        if (!this->printed_sparsity)
        {
            this->printed_sparsity = true;
            if (!myProblem->options.quiet_basinhopping || myProblem->options.run_inner_loop == InnerLoopSolverType::NLP)
            {
                if (myProblem->options.short_output_file_names)
                {

                    //myProblem->output_Jacobian_sparsity_information(myProblem->options.working_directory + "//" + "SparsityDescriptions.csv");
                    myProblem->output_problem_bounds_and_descriptions(myProblem->options.working_directory + "//" + "XFfile.csv");
                }
                else
                {
                    //myProblem->output_Jacobian_sparsity_information(myProblem->options.working_directory + "//" + myProblem->options.mission_name + "_" + myProblem->options.description + "_SparsityDescriptions.csv");
                    myProblem->output_problem_bounds_and_descriptions(myProblem->options.working_directory + "//" + myProblem->options.mission_name + "_" + myProblem->options.description + "XFfile.csv");
                }
            }
        }

        if (this->myProblem->options.NLP_solver_type == 1)
        {
            std::cout << "WORHP interface is deprecated" << std::endl;
        }
        else
        {
            //run SNOPT
            
            try
            {
				if (this->myProblem->options.enable_Scalatron)
				{
					this->myProblem->myScalatron->initialize(this->X_after_hop_unscaled,
						this->myProblem->Xlowerbounds,
						this->myProblem->Xupperbounds,
						this->myProblem->F,
						this->myProblem->Flowerbounds,
						this->myProblem->Fupperbounds,
						this->myProblem->iGfun,
						this->myProblem->jGvar,
						this->myProblem->G,
						0);
				}

                this->mySNOPT->run_NLP(false);

                this->X_after_slide_unscaled = this->mySNOPT->getX_unscaled();
                for (size_t Xindex = 0; Xindex < this->myProblem->total_number_of_NLP_parameters; ++Xindex)
                {
                    this->X_after_slide[Xindex] = this->X_after_slide_unscaled[Xindex] / this->myProblem->X_scale_factors[Xindex];
                }

            }
            catch (std::exception &error)
            {
                std::cout << error.what() << std::endl;
                //prevent a crash, yay
                this->X_after_slide_unscaled = this->mySNOPT->getX_unscaled();
                for (size_t Xindex = 0; Xindex < this->myProblem->total_number_of_NLP_parameters; ++Xindex)
                {
                    this->X_after_slide[Xindex] = this->X_after_slide_unscaled[Xindex] / this->myProblem->X_scale_factors[Xindex];
                }

                std::cout << "SNOPT has crashed on mission " << myProblem->options.description << ". Creating dumpfile." << std::endl;
                std::stringstream dumpstream;
                dumpstream << myProblem->options.working_directory << "//" << myProblem->options.mission_name << "_" << myProblem->options.description << ".SNOPTcrash";
                std::ofstream dumpfile(dumpstream.str().c_str(), std::ios::out | std::ios::trunc);
                dumpfile << "SNOPT crashed, caught by MBH try-catch block, on mission " << myProblem->options.mission_name << "_" << myProblem->options.description << std::endl;
                dumpfile << "After " << this->number_of_solutions << std::endl;
                dumpfile.precision(20);
                dumpfile << "decision vector follows" << std::endl;
                dumpfile << std::endl;
                dumpfile << this->X_after_slide_unscaled[0];
                for (size_t k = 1; k < myProblem->total_number_of_NLP_parameters; ++k)
                    dumpfile << " " << this->X_after_slide_unscaled[k];
                dumpfile << std::endl;
                dumpfile << "Initial guess was:" << std::endl;
                dumpfile << this->X_after_hop_unscaled[0];
                for (size_t k = 1; k < myProblem->total_number_of_NLP_parameters; ++k)
                    dumpfile << " " << this->X_after_hop_unscaled[k];
                dumpfile << std::endl;
                /*dumpfile << "F was:" << std::endl;
                dumpfile << F[0];
                for (size_t k = 1; k < myProblem->total_number_of_constraints; ++k)
                    dumpfile << " " << F[k];
                dumpfile << std::endl;

                this->myProblem->unscale(x);
                myProblem->evaluate(this->myProblem->X, myProblem->F, myProblem->G, 0, myProblem->iGfun, myProblem->jGvar);
                dumpfile << "G was:" << std::endl;
                dumpfile << this->myProblem->G[0];
                for (size_t k = 1; k < this->myProblem->G.size(); ++k)
                    dumpfile << " " << this->myProblem->G[k];
                dumpfile.close();*/
            }
        }

        //Step 4: set the myProblem state to where the solver left off
        try
        {
            myProblem->evaluate(this->X_after_slide_unscaled, myProblem->F, myProblem->G, false);
        }
        catch (std::exception &error)
        {
            if (!myProblem->options.quiet_basinhopping) 
            {
                std::cout << error.what() << std::endl;
                std::cout << "EMTG::Invalid initial point or failure in objective function while evaluating final point. MBH" << std::endl;
            }
        }
    }

    
    //function to run MBH
    int MBH::run()
    {
        bool new_point = true;
        bool seeded_step = false;
        double best_feasibility = math::LARGE;

        //If we have seeded MBH, start from that seed. Otherwise we will need to generate a new random point.
        if (myProblem->options.seed_MBH)
        {
            seeded_step = true;
            new_point = false;
            this->feasible_point_finder_active = this->myProblem->options.ACE_feasible_point_finder;
            if (myProblem->options.skip_first_nlp_run)
            {
                seeded_step = false;

                this->X_incumbent = this->X_after_hop;

                double feasibility, normalized_feasibility, distance_from_equality_filament, decision_variable_infeasibility;
                size_t worst_decision_variable;

                this->myProblem->check_feasibility(this->X_incumbent,
                    this->myProblem->F,
                    worst_decision_variable,
                    worst_constraint,
                    feasibility,
                    normalized_feasibility,
                    distance_from_equality_filament,
                    decision_variable_infeasibility);

                this->Jincumbent = this->myProblem->F[0];
                best_feasibility = fmax(normalized_feasibility, decision_variable_infeasibility);
            }
        }

        this->number_of_solutions = 0;
        this->number_of_resets = 0;
        this->number_of_attempts = 1;
        bool continue_flag = true;

        //print the archive header
        std::string archive_file;
        if (myProblem->options.short_output_file_names)
            archive_file = myProblem->options.working_directory + "//" + myProblem->options.mission_name + "archive.emtg_archive";
        else
            archive_file = myProblem->options.working_directory + "//" + myProblem->options.mission_name + "_" + myProblem->options.description + "archive.emtg_archive";

        if (!myProblem->options.quiet_basinhopping || myProblem->options.MBH_always_write_archive)
            print_archive_header(archive_file);

        this->MBH_start_time = time(NULL);

        do
        {
            ++this->number_of_attempts;
            if (new_point)
            {
                //Step 1: generate a random new point
                this->reset_point();

                this->number_of_failures_since_last_improvement = 0;
            }
            else if (!seeded_step)
            {                
                //Step 1 (alternate): perturb the existing point
                this->hop();

                if (myProblem->options.MBH_time_hop_probability > 0.0  && best_feasibility >= this->myProblem->options.snopt_feasibility_tolerance)
                    this->time_hop();
            }
			//evaluate the case once
			try
			{
				this->myProblem->evaluate(this->X_after_hop_unscaled, this->myProblem->F, this->myProblem->G, true);
			}
            catch (std::exception &error)
			{
                if (!myProblem->options.quiet_basinhopping)
                {
                    std::cout << "EMTG::Invalid initial point or failure in objective function while evaluating initial point." << std::endl;
                    std::cout << error.what() << std::endl;
                }
			}
            
            //if seeding MBH, only the first step runs from the seed. After that hopping occurs.
            seeded_step = false;

            //Step 2: apply the slide operator      
            this->slide();
            double ObjectiveFunctionValue = this->myProblem->F[0];
            std::vector<double> Xunscaled = this->mySNOPT->getX_unscaled();
            std::vector<double> F = this->myProblem->F;

            //Step 3: determine if the new trial point is feasible and if so, operate on it
            double feasibility, normalized_feasibility, distance_from_equality_filament, decision_variable_infeasibility;
            size_t worst_decision_variable;

            this->myProblem->check_feasibility(Xunscaled,
                this->myProblem->F,
                worst_decision_variable,
                worst_constraint,
                feasibility,
                normalized_feasibility,
                distance_from_equality_filament,
                decision_variable_infeasibility);

            //print worst constraint violation
            if (!myProblem->options.quiet_basinhopping)
            {
                std::cout << "Worst constraint is F[" << this->worst_constraint << "]: " << this->myProblem->Fdescriptions[this->worst_constraint] << std::endl;
                std::cout << "with violation " << feasibility << std::endl;
            }

            //if we ALWAYS want our archive file, write it out
            if (myProblem->options.MBH_always_write_archive)
            {
                std::string solution_name = "R" + std::to_string(this->number_of_resets) + "r" + std::to_string(this->number_of_attempts);

                this->archive.push_back(std::make_tuple(
                    EMTG_innerloop_solution(this->myProblem,
                        this->X_after_slide_unscaled,
                        this->myProblem->F,
                        solution_name),
                        time(NULL) - this->MBH_start_time));

                this->print_archive_line(archive_file, this->archive.back());
            }

            //note: I do not trust SNOPT's "requested accuracy could not be achieved" return state - I prefer my own feasibility check
            bool isFeasible = false;
            if (this->myProblem->options.enable_NLP_chaperone)
            {
                if (normalized_feasibility < myProblem->options.snopt_feasibility_tolerance && decision_variable_infeasibility < myProblem->options.snopt_feasibility_tolerance)
                    isFeasible = true;
            }
            else
            {
                if ((normalized_feasibility < myProblem->options.snopt_feasibility_tolerance && decision_variable_infeasibility < myProblem->options.snopt_feasibility_tolerance)
                    || this->mySNOPT->getInform() < 10)
                    isFeasible = true;
            }

            if (isFeasible)
            {
                //Step 3.1: if the trial point is feasible, add it to the archive. Also disable the feasible point finder
                this->feasible_point_finder_active = false;

                if (!myProblem->options.quiet_basinhopping)
                {
                    if (!myProblem->options.MBH_always_write_archive) //because if MBH_always_write_archive was on, then we would have already written to the archive a few lines ago
                    {
                        std::string solution_name = "R" + std::to_string(this->number_of_resets) + "r" + std::to_string(this->number_of_attempts);

                        this->archive.push_back(std::make_tuple(
                            EMTG_innerloop_solution(this->myProblem,
                                this->X_after_slide_unscaled,
                                this->myProblem->F,
                                solution_name),
                                time(NULL) - this->MBH_start_time));

                        this->print_archive_line(archive_file, this->archive.back());
                    }
                    std::cout << "Hop evaluated mission " << myProblem->options.description << " with fitness " << ObjectiveFunctionValue << std::endl;
                }

                ++number_of_solutions;

                //Step 3.2: if the point came from a hop operation and is superior to the current point, adopt it as the new current point
                if (!(new_point))
                {
                    if (ObjectiveFunctionValue < this->Jincumbent)
                    {
                        this->Jincumbent = ObjectiveFunctionValue;
                        this->X_incumbent = this->X_after_slide;

                        if (!myProblem->options.quiet_basinhopping)
                            std::cout << "New local best" << std::endl;

                        ++this->number_of_improvements;

                        this->number_of_failures_since_last_improvement = 0;

                        //if we have chosen to always write every improvement for later animation, then we need to print out a "hop" file
                        if (this->myProblem->options.MBH_write_every_improvement)
                        {
                            //I think that the following line is not needed. It also introduces a risk - the final solution file is written out using myProblem->Xopt and so
                            //we really don't want it to update except when the global incumbent updates
                            //this->myProblem->Xopt = this->X_after_slide_unscaled;
                            this->myProblem->what_the_heck_am_I_called(SolutionOutputType::HOP, this->number_of_attempts);
                            this->myProblem->output(this->myProblem->options.outputfile);
                        }
                    }
                    else
                        ++this->number_of_failures_since_last_improvement;
                }
                else //if this is from a reset, adopt it as the current point and turn off the "new point" flag
                {
                    this->Jincumbent = ObjectiveFunctionValue;
                    this->X_incumbent = this->X_after_slide;
                    new_point = false;
                    
                    //if we have chosen to always write every improvement for later animation, then we need to print out a "hop" file
                    if (this->myProblem->options.MBH_write_every_improvement)
                    {

                        //I think that the following line is not needed. It also introduces a risk - the final solution file is written out using myProblem->Xopt and so
                        //we really don't want it to update except when the global incumbent updates
                        //this->myProblem->Xopt = this->X_after_slide_unscaled;
                        this->myProblem->what_the_heck_am_I_called(SolutionOutputType::HOP, this->number_of_attempts);
                        this->myProblem->output(this->myProblem->options.outputfile);
                    }
                }

                //Step 3.3 update the global best if applicable
                if (ObjectiveFunctionValue < this->JGlobalIncumbent)
                {
                    this->JGlobalIncumbent = ObjectiveFunctionValue;
                    this->X_global_incumbent = this->X_after_slide;
                    myProblem->Xopt = this->X_after_slide_unscaled;
                    myProblem->best_cost = this->JGlobalIncumbent;
                    

                    if (!myProblem->options.quiet_basinhopping)
                        std::cout << "New global best" << std::endl;

                    //Write out a results file for the current global best
                    bool successfully_evaluated = false;
                    try
                    {
                        this->myProblem->evaluate(myProblem->Xopt, myProblem->F, myProblem->G, false);
                        successfully_evaluated = true;
                    }
                    catch (std::exception &error)
                    {
                        std::cout << "Failure to evaluate " << myProblem->options.description << std::endl;
                        std::cout << error.what() << std::endl;
                        ObjectiveFunctionValue = EMTG::math::LARGE;
                        break;
                    }

                    if (successfully_evaluated)
                    {
                        this->myProblem->what_the_heck_am_I_called(SolutionOutputType::SUCCESS);
                        this->myProblem->output(this->myProblem->options.outputfile);
                    }
                }

            }
            else if (this->feasible_point_finder_active && best_feasibility >= this->myProblem->options.snopt_feasibility_tolerance)
            {
                //if we have not yet found our first feasible point and the ACE feasible point finder is enabled
                //then we should see if this point is "more feasible" than best one we have so far
                if (fmax(normalized_feasibility, decision_variable_infeasibility) < best_feasibility)
                {
                    if (!myProblem->options.quiet_basinhopping)
                    {
                        std::cout << "Acquired slightly less infeasible point with constraint feasibility " << normalized_feasibility << std::endl;
                        std::cout << "and decision variable feasibility " << decision_variable_infeasibility << std::endl;
                        //std::std::cout << "Worst constraint is F[" << this->worst_constraint << "]: " << this->myProblem->Fdescriptions[this->worst_constraint] << std::endl;
                        //std::std::cout << "with abs(violation) " << feasibility << std::endl;
                    }
                    this->Jincumbent = ObjectiveFunctionValue;
                    this->X_incumbent = this->X_after_slide;
                    myProblem->Xopt = this->X_after_slide_unscaled; //we store the unscaled Xcurrent
                    best_feasibility = fmax(normalized_feasibility, decision_variable_infeasibility);
                    this->X_most_feasible = this->X_after_slide_unscaled;
                    new_point = false;

                    if (this->myProblem->options.MBH_write_every_improvement)
                    {
                        this->myProblem->what_the_heck_am_I_called(SolutionOutputType::HOP, this->number_of_attempts);
                        this->myProblem->output(this->myProblem->options.outputfile);
                    }
                }
                ++this->number_of_failures_since_last_improvement;
            }
            else
            {
                ++this->number_of_failures_since_last_improvement;

                if (fmax(normalized_feasibility, decision_variable_infeasibility) < this->FeasibilityIncumbent)
                {
                    this->FeasibilityIncumbent = fmax(normalized_feasibility, decision_variable_infeasibility);
                    this->X_most_feasible = this->X_after_slide_unscaled;
                }
            }

            if (this->number_of_failures_since_last_improvement >= myProblem->options.MBH_max_not_improve)
                new_point = true;

        } while (this->number_of_attempts < myProblem->options.MBH_max_trials && (time(NULL) - this->MBH_start_time) < myProblem->options.MBH_max_run_time);

        if (!myProblem->options.quiet_basinhopping)
        {
            std::cout << std::endl;
            if (this->number_of_solutions > 0)
                std::cout << "Best value found was " << this->JGlobalIncumbent << std::endl;
            else
                std::cout << "No feasible solutions found." << std::endl;
        }

        return this->number_of_solutions;
    }

    //function to print the archive header to a text file
    void MBH::print_archive_header(const std::string& filename)
    {
        //print the archive in CSV format
        //the final entry in each line is the fitness value

        std::ofstream outputfile (filename.c_str(), std::ios::trunc);

        //archive column headers
        if (this->myProblem->options.MBH_archive_state_vector)
            for (int entry = 0; entry < myProblem->total_number_of_NLP_parameters; ++entry)
                outputfile << myProblem->Xdescriptions[entry] << ",";
        outputfile << "reset count,";
        outputfile << "step count,";
        outputfile << "solution timestamp,";
        outputfile << "Objective function" << std::endl;

        outputfile.close();
    }

    //function to print the archive line to a text file
    void MBH::print_archive_line(const std::string& filename, std::tuple< EMTG_innerloop_solution, time_t > solution)
    {
        //print the archive in CSV format
        //the final entry in each line is the fitness value

        std::ofstream outputfile (filename.c_str(), std::ios::app);

        //archive lines
        if (this->myProblem->options.MBH_archive_state_vector)
        {
            for (size_t entry = 0; entry < myProblem->total_number_of_NLP_parameters; ++entry)
            {
                if (this->myProblem->Xdescriptions[entry].find("epoch") < 1024 || this->myProblem->Xdescriptions[entry].find("time") < 1024)
                {
                    outputfile << std::get<0>(solution).getX()[entry] / 86400.0 << ",";
                }
                else
                    outputfile << std::get<0>(solution).getX()[entry] << ",";
            }
        }

        outputfile << this->number_of_resets << ",";
        outputfile << this->number_of_attempts << ",";
        outputfile << std::get<1>(solution) << ",";
        outputfile << std::get<0>(solution).getobjective_function_value() << ",";
        outputfile << std::endl;
        outputfile.close();
    }
    

}} //close namespace