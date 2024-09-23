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


//EMTGv9 Problem class

#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

#include "problem.h"
#include "monotonic_basin_hopping.h"
#include "SNOPT_interface.h"
#include "NLPoptions.h"
#include "FilamentWalker.h"
#include "EMTG_math.h"

#include "boost/date_time.hpp"
#include "boost/date_time/local_time/local_date_time.hpp"

#include "ScalatronFactory.h"


namespace EMTG
{

    problem::problem() :
        total_number_of_constraints(0),
        total_number_of_NLP_parameters(0),
        number_of_solutions(0),
        best_cost(1.0e+100)
    {
		if (this->options.enable_Scalatron)
			this->myScalatron = Scalatron::Create_Scalatron(Scalatron::ScalatronType::PJRN);
    }


	problem::~problem()
	{
		if (this->options.enable_Scalatron)
			delete this->myScalatron;
	}

    //evaluate stub
    void problem::evaluate(const std::vector<doubleType>& X)
    {
        this->evaluate(X, this->F, this->G, false);
    }

    bool problem::optimize() 
    {
        bool successfully_optimized = false;

        switch (options.run_inner_loop)
        {
            case InnerLoopSolverType::RUN_TRIALX: //run trialX
            {
                this->options.current_trialX = this->construct_initial_guess();

                this->Xopt.resize(this->total_number_of_NLP_parameters);
                for (size_t Xindex = 0; Xindex < this->total_number_of_NLP_parameters; ++Xindex)
                    Xopt[Xindex] = options.current_trialX[Xindex];

                if (!(options.current_trialX.size() == Xdescriptions.size()))
                {
                    std::cout << "Invalid initial guess. trialX vector is the wrong length " << options.current_trialX.size() << ", expected " << Xdescriptions.size() << "." << std::endl;
                    if (options.current_trialX.size() < Xdescriptions.size())
                    {
                        throw std::invalid_argument("trialX (your initial guess) is too small. Aborting. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }
                    else
                    {
                        std::cout << "trialX (your initial guess)  is too big. This probably means you have the initial guess to a different problem but this error is not fatal so we shall proceed." << std::endl;
                    }
                }

				double normalized_feasibility = 1e100;
				double decision_variable_feasibility = 1e100;

                try
                {
                    this->evaluate(this->Xopt, this->F, this->G, true);

                    successfully_optimized = true;
                }
                catch (std::exception &error)
                {
                    std::cout << error.what() << std::endl;
                    std::cout << "Failure to evaluate " << this->options.mission_name << std::endl;
                }
				
                this->X = this->Xopt;

                if (successfully_optimized)
                {
                    double feasibility, normalized_feasibility, distance_from_equality_filament, decision_variable_infeasibility;
                    size_t worst_constraint, worst_decision_variable;

                    this->check_feasibility(this->X,
                        this->F,
                        worst_decision_variable,
                        worst_constraint,
                        feasibility,
                        normalized_feasibility,
                        distance_from_equality_filament,
                        decision_variable_infeasibility);

                    std::cout << "J = " << this->F[0] << std::endl;
                    if (normalized_feasibility > options.snopt_feasibility_tolerance || decision_variable_infeasibility > options.snopt_feasibility_tolerance)
                    {
                        this->what_the_heck_am_I_called(SolutionOutputType::FAILURE);
                        std::cout << "Acquired infeasible point ";
                    }
                    else
                    {
                        this->what_the_heck_am_I_called(SolutionOutputType::SUCCESS);
                        std::cout << "Acquired feasible point ";
                    }

                    std::cout << "with feasibility " << normalized_feasibility << std::endl;
                    std::cout << "Worst constraint is F[" << worst_constraint << "]: " << this->Fdescriptions[worst_constraint] << std::endl;
                    std::cout << "with violation " << feasibility << std::endl;

                    if (decision_variable_infeasibility < options.snopt_feasibility_tolerance)
                    {
                        std::cout << "Decision vector is feasible." << std::endl;
                    }
                    else
                    {
                        std::cout << "Decision vector is infeasible." << std::endl;
                        std::cout << "Worst decision variable is X[" << worst_decision_variable << "]: " << this->Xdescriptions[worst_decision_variable] << std::endl;
                        std::cout << "with violation " << decision_variable_infeasibility << std::endl;
                    }

                    this->output(this->options.outputfile);
                }

                this->output_problem_bounds_and_descriptions(this->options.working_directory + "//" + "XFfile.csv");

                break;
            }

#ifndef AD_INSTRUMENTATION
            case InnerLoopSolverType::MBH: //run MBH
            {
                Solvers::NLPoptions myNLPoptions(this->options);

                Solvers::SNOPT_interface mySNOPT(this, myNLPoptions);
                EMTG::Solvers::MBH solver(this, &mySNOPT);

                if (options.seed_MBH)
                {
                    this->options.current_trialX = this->construct_initial_guess();

                    solver.seed(options.current_trialX);
                }

                this->number_of_solutions = solver.run();

                if (this->number_of_solutions > 0)
                    this->what_the_heck_am_I_called(SolutionOutputType::SUCCESS);
                else
                    this->what_the_heck_am_I_called(SolutionOutputType::FAILURE);

                if (this->number_of_solutions == 0 || solver.getJGlobalIncumbent() > 1.0e+10)
                {
                    this->Xopt = solver.getX_most_feasible(); //we store the unscaled Xcurrent

                    options.outputfile = options.working_directory + "//FAILURE_" + options.mission_name + ".emtg";
                }                

                try
                {
                    this->evaluate(this->Xopt, this->F, this->G, false);
                    successfully_optimized = true;

                    this->output(this->options.outputfile);

                    //print an XFfile
                    this->X = this->Xopt;
                    this->output_problem_bounds_and_descriptions(this->options.working_directory + "//" + "XFfile.csv");
                }
                catch (std::exception &error)
                {
                    std::cout << error.what() << std::endl;
                    std::cout << "Failure to evaluate " << options.mission_name + "_" + options.description << std::endl;
                    F[0] = EMTG::math::LARGE;
                    this->number_of_solutions = 0;

                }

                break;
            }
            case InnerLoopSolverType::ACDE:
            {
                std::cout << "ACDE is deprecated" << std::endl;

                break;
            }
            case InnerLoopSolverType::NLP: //run NLP on an existing initial guess
            {
                { //boring stuff, not specific

                    //set up outputs
                    if (this->options.short_output_file_names)
                    {
                        options.outputfile = options.working_directory + "//" + options.mission_name + ".emtg";
                        //this->output_Jacobian_sparsity_information(this->options.working_directory + "//" + "SparsityDescriptions.csv");
                        this->output_problem_bounds_and_descriptions(this->options.working_directory + "//" + "XFfile.csv");
                    }
                    else
                    {
                        options.outputfile = options.working_directory + "//" + options.mission_name;
						if (std::strcmp("",this->options.description.c_str()) != 0)
							 options.outputfile += "_" + options.description;
						options.outputfile += ".emtg";
                        //this->output_Jacobian_sparsity_information(this->options.working_directory + "//" + this->options.mission_name + "_" + this->options.description + "_SparsityDescriptions.csv");
                        this->output_problem_bounds_and_descriptions(this->options.working_directory + "//" + this->options.mission_name + "_" + this->options.description + "XFfile.csv");
                    }

                    //make an initial guess, 'cuz that's important and stuff
                    this->options.current_trialX = this->construct_initial_guess();
                }

                Solvers::NLPoptions myNLPoptions(this->options);

                Solvers::SNOPT_interface mySNOPT(this, myNLPoptions);

                mySNOPT.setX0_unscaled(this->options.current_trialX);
                
				try
				{
					this->evaluate(this->options.current_trialX, this->F, this->G, true);
                    successfully_optimized = true;
				}
                catch (std::exception &error)
                {
                    std::cout << error.what() << std::endl;
					std::cout << "Failure to evaluate initial guess for " << this->options.mission_name << std::endl;
				}

				if (this->options.enable_Scalatron)
				{
					this->myScalatron->initialize(this->options.current_trialX,
						this->Xlowerbounds,
						this->Xupperbounds,
						this->F,
						this->Flowerbounds,
						this->Fupperbounds,
						this->iGfun,
						this->jGvar,
						this->G,
						0);
				}

                mySNOPT.run_NLP(false);

                this->Xopt = mySNOPT.getX_unscaled();
                try
                {
                    this->evaluate(this->Xopt, this->F, this->G, false);
                    successfully_optimized = true;
                }
                catch (std::exception &error)
                {
                    std::cout << error.what() << std::endl;
                    std::cout << "Could not evaluate NLP solution." << std::endl;
                }

                this->X = this->Xopt;
                this->output_problem_bounds_and_descriptions(this->options.working_directory + "//" + this->options.mission_name + "_" + this->options.description + "XFfile.csv");



                double feasibility, normalized_feasibility, distance_from_equality_filament, decision_variable_infeasibility;
                size_t worst_constraint, worst_decision_variable;

                if (successfully_optimized)
                {
                    this->check_feasibility(this->X,
                        this->F,
                        worst_decision_variable,
                        worst_constraint,
                        feasibility,
                        normalized_feasibility,
                        distance_from_equality_filament,
                        decision_variable_infeasibility);

                    std::cout << "J = " << this->F[0] << std::endl;
                    if (normalized_feasibility > options.snopt_feasibility_tolerance || decision_variable_infeasibility > options.snopt_feasibility_tolerance)
                    {
                        this->what_the_heck_am_I_called(SolutionOutputType::FAILURE);
                        std::cout << "Acquired infeasible point ";
                    }
                    else
                    {
                        this->what_the_heck_am_I_called(SolutionOutputType::SUCCESS);
                        std::cout << "Acquired feasible point ";
                    }

                    std::cout << "with feasibility " << normalized_feasibility << std::endl;
                    std::cout << "Worst constraint is F[" << worst_constraint << "]: " << this->Fdescriptions[worst_constraint] << std::endl;
                    std::cout << "with violation " << feasibility << std::endl;

                    if (decision_variable_infeasibility < options.snopt_feasibility_tolerance)
                    {
                        std::cout << "Decision vector is feasible." << std::endl;
                    }
                    else
                    {
                        std::cout << "Decision vector is infeasible." << std::endl;
                        std::cout << "Worst decision variable is X[" << worst_decision_variable << "]: " << this->Xdescriptions[worst_decision_variable] << std::endl;
                        std::cout << "with violation " << decision_variable_infeasibility << std::endl;
                    }

                    this->output(this->options.outputfile);

                    //print an XFfile
                    this->X = this->Xopt;
                    this->output_problem_bounds_and_descriptions(this->options.working_directory + "//" + "XFfile.csv");
                }

                break;
            }
            case InnerLoopSolverType::FilamentWalker: //filament walker
            {
                { //boring stuff, not specific
                    //set up outputs
                    if (this->options.short_output_file_names)
                    {
                        options.outputfile = options.working_directory + "//" + options.mission_name + "_" + ".emtg";
                        //this->output_Jacobian_sparsity_information(this->options.working_directory + "//" + "SparsityDescriptions.csv");
                        this->output_problem_bounds_and_descriptions(this->options.working_directory + "//" + "XFfile.csv");
                    }
                    else
                    {
                        options.outputfile = options.working_directory + "//" + options.mission_name + "_" + options.description + ".emtg";
                        //this->output_Jacobian_sparsity_information(this->options.working_directory + "//" + this->options.mission_name + "_" + this->options.description + "_SparsityDescriptions.csv");
                        this->output_problem_bounds_and_descriptions(this->options.working_directory + "//" + this->options.mission_name + "_" + this->options.description + "XFfile.csv");
                    }
                }

                //do filament walker things
                Solvers::NLPoptions myNLPoptions(this->options);

                Solvers::SNOPT_interface mySNOPT(this, myNLPoptions);
                Solvers::FilamentWalker myFilamentWalker(this, &mySNOPT);
                myFilamentWalker.walk();

                //write the output - but right now filament walkers don't really have output
                //we have to fake Xopt so we don't get a crash
                this->Xopt = this->Xupperbounds;


                try
                {
                    this->evaluate(this->Xopt, this->F, this->G, false);
                    successfully_optimized = true;
                }
                catch (std::exception &error)
                {
                    std::cout << error.what() << std::endl;
                    std::cout << "Could not evaluate Filament Walker solution." << std::endl;
                }

                break;
            }
#endif
        }

        return successfully_optimized;
    }//end optimize()


    //function to output X and F bounds, descriptions
    void problem::output_problem_bounds_and_descriptions()
    {
        this->output_problem_bounds_and_descriptions("XFfile.csv");
    }//end output_problem_bounds_and_descriptions()

    void problem::output_problem_bounds_and_descriptions(std::string filestring)
    {
        std::ofstream outputfile(filestring.c_str(), std::ios::trunc);
        outputfile.precision(13);
        outputfile << "Xindex, Description, Lowerbound, Upperbound, ScaleFactor, Value" << std::endl;
        for (size_t k = 0; k < Xdescriptions.size(); ++k)
            outputfile << k << "," << Xdescriptions[k] << "," << Xlowerbounds[k] << "," << Xupperbounds[k] << "," << X_scale_factors[k] << "," << X[k] _GETVALUE<< std::endl;

        outputfile << "Findex, Description, Lowerbound, Upperbound, Value" << std::endl;
        for (size_t k = 0; k < Fdescriptions.size(); ++k)
            outputfile << k << "," << Fdescriptions[k] << "," << Flowerbounds[k] << "," << Fupperbounds[k] << "," << F[k] _GETVALUE << std::endl;

        outputfile.close();
    }//end output_problem_bounds_and_descriptions()

    //function to output the Jacobian sparsity information
    void problem::output_Jacobian_sparsity_information(std::string filestring)
    {
        std::ofstream outputfile(filestring.c_str(), std::ios::trunc);

        outputfile <<  "Linear Constraints (A):" << std::endl;
        for (size_t k = 0; k < Adescriptions.size(); ++k)
            outputfile << iAfun[k] << "," << jAvar[k] << "," << k << "," << Adescriptions[k] << std::endl;

        outputfile << std::endl;

        outputfile << "Nonlinear Constraints (G):" << std::endl;
        for (size_t k = 0; k < Gdescriptions.size(); ++k)
            outputfile << iGfun[k] << "," << jGvar[k] << "," << k << "," << Gdescriptions[k] << std::endl;

        outputfile.close();
    }//end output_Jacobian_sparsity_information()

    //function to check feasibility of a solution
    void problem::check_feasibility(const std::vector<doubleType>& X,
        const std::vector<doubleType>& F,
        size_t& worst_decision_variable,
        size_t& worst_constraint,
        double& max_constraint_violation,
        double& normalized_max_constraint_violation,
        double& distance_from_equality_filament,
        double& decision_variable_infeasibility,
        bool calledFromNLP)
    {
        worst_decision_variable = 0;
        max_constraint_violation = 0.0;
        distance_from_equality_filament = 0.0;
        decision_variable_infeasibility = 0.0;


        //check decision vector feasibility
        for (size_t Xindex = 0; Xindex < this->total_number_of_NLP_parameters; ++Xindex)
        {
            if (X[Xindex] != X[Xindex]) //check for solutions that have NaN in their constraint vectors and make sure they get thrown out!
            {

                if (calledFromNLP)
                {
                    throw std::runtime_error("nan alert: X[" + std::to_string(Xindex) + "]: " + this->Xdescriptions[Xindex]);
                }
                else
                {
                    std::cout << "decision vector nan alert --------------------------------------------------------------" << std::endl;
                    std::cout << "X[" << Xindex << "]: " << this->Xdescriptions[Xindex] << std::endl;
                    decision_variable_infeasibility = math::LARGE;
                    worst_decision_variable = Xindex;
                }
            }
            else if (X[Xindex] > Xupperbounds[Xindex]) //exceed upper bound
            {
                if ((X[Xindex] - Xupperbounds[Xindex]) / this->X_scale_factors[Xindex] > fabs(decision_variable_infeasibility))
                {
                    decision_variable_infeasibility = fabs((X[Xindex]) _GETVALUE - Xupperbounds[Xindex]) / this->X_scale_factors[Xindex];
                    worst_decision_variable = Xindex;
                }
            }
            else if (X[Xindex] < Xlowerbounds[Xindex]) //smaller than lower bound
            {
                if ((Xlowerbounds[Xindex] - X[Xindex]) / this->X_scale_factors[Xindex] > fabs(decision_variable_infeasibility))
                {
                    decision_variable_infeasibility = fabs((X[Xindex]) _GETVALUE - Xlowerbounds[Xindex]) / this->X_scale_factors[Xindex];
                    worst_decision_variable = Xindex;
                }
            }
        }

        //check constraint feasibility
        for (size_t Findex = 1; Findex < this->total_number_of_constraints; ++Findex)
        {
            if (F[Findex] != F[Findex]) //check for solutions that have NaN in their constraint vectors and make sure they get thrown out!
            {
                if (calledFromNLP)
                {
                    throw std::runtime_error("nan alert: F[" + std::to_string(Findex) + "]: " + this->Fdescriptions[Findex]);
                }
                else
                {
                    std::cout << "constraint nan alert --------------------------------------------------------------" << std::endl;
                    std::cout << "F[" << Findex << "]: " << this->Fdescriptions[Findex] << std::endl;
                    max_constraint_violation = math::LARGE;
                    worst_constraint = Findex;
                }
            }
            else if (F[Findex] > Fupperbounds[Findex]) //exceed upper bound
            {
                if (F[Findex] - Fupperbounds[Findex] > fabs(max_constraint_violation))
                {
                    max_constraint_violation = (F[Findex]) _GETVALUE - Fupperbounds[Findex];
                    worst_constraint = Findex;
                }

                if (this->F_equality_or_inequality[Findex - 1])
                {
                    distance_from_equality_filament += ((F[Findex]) _GETVALUE - Fupperbounds[Findex]) * ((F[Findex]) _GETVALUE - Fupperbounds[Findex]);
                }
            }
            else if (F[Findex] < Flowerbounds[Findex]) //smaller than lower bound
            {
                if (Flowerbounds[Findex] - F[Findex] > fabs(max_constraint_violation))
                {
                    max_constraint_violation = (F[Findex]) _GETVALUE - Flowerbounds[Findex];
                    worst_constraint = Findex;
                }

                if (this->F_equality_or_inequality[Findex - 1])
                {
                    distance_from_equality_filament += ((Flowerbounds[Findex] - F[Findex]) * (Flowerbounds[Findex] - F[Findex]))  _GETVALUE;
                }
            }
        }

        //we're not actually going to normalize
        normalized_max_constraint_violation = fabs(max_constraint_violation);
        
        distance_from_equality_filament = sqrt(distance_from_equality_filament);
    }

    //functions to check inputs (X) and outputs (F) for NaN, inf
    bool problem::check_inputs_for_invalid_entries(const std::vector<doubleType>& X, const bool& PrintAnyway)
    {
        int failure_entry = -1;
        std::string failure_type = "";


        //look for failures
        for (size_t Xentry = 0; Xentry < this->Xdescriptions.size(); ++Xentry)
        {
            if (std::isnan(X[Xentry] _GETVALUE))
            {
                failure_entry = Xentry;
                failure_type = "NaN";
                break;
            }
            else if (!std::isfinite(X[Xentry] _GETVALUE))
            {
                failure_entry = Xentry;
                failure_type = "inf";
                break;
            }
        }

        //if a failure is found, write out a dump and throw an error
        if (failure_entry >= 0 || PrintAnyway)
        {
            std::stringstream outfilestream;
            boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
            std::stringstream timestream;
            timestream << static_cast<int>(now.date().month()) << now.date().day() << now.date().year() << "_" << now.time_of_day().hours() << now.time_of_day().minutes() << now.time_of_day().seconds();
            outfilestream << this->options.working_directory << "/" << "bad_decision_vector_" << timestream.str() << ".Xcrash";
            std::ofstream dumpfile(outfilestream.str().c_str(), std::ios::out | std::ios::trunc);

            dumpfile << "Bad entry in X" << std::endl;
            dumpfile << "Bad entry is #" << failure_entry << std::endl;
            dumpfile << this->options.mission_name << std::endl;
            dumpfile << std::endl;
            for (size_t Xentry = 0; Xentry < this->Xdescriptions.size(); ++Xentry)
            {
                dumpfile << "X[" << Xentry << "] (" << this->Xdescriptions[Xentry] << ") = " << X[Xentry] << std::endl;
            }
            dumpfile.close();

            return true;
        }

        //if you got this far there was no problem
        return false;
    }//end check_inputs_for_invalid_entries()

    bool problem::check_outputs_for_invalid_entries(const std::vector<doubleType>& X, const std::vector<doubleType>& F, const bool& PrintAnyway)
    {
        int failure_entry = -1;
        std::string failure_type = "";


        //look for failures
        for (size_t Fentry = 0; Fentry < this->Fdescriptions.size(); ++Fentry)
        {
            if (std::isnan(F[Fentry] _GETVALUE))
            {
                failure_entry = Fentry;
                failure_type = "NaN";
                break;
            }
            else if (!std::isfinite(F[Fentry] _GETVALUE))
            {
                failure_entry = Fentry;
                failure_type = "inf";
                break;
            }
        }

        //if a failure is found, write out a dump and throw an error
        if (failure_entry >= 0 || PrintAnyway)
        {
            std::stringstream outfilestream;
            boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
            std::stringstream timestream;
            timestream << static_cast<int>(now.date().month()) << now.date().day() << now.date().year() << "_" << now.time_of_day().hours() << now.time_of_day().minutes() << now.time_of_day().seconds();
            outfilestream << this->options.working_directory << "/" << "bad_objective_constraint_vector_" << timestream.str() << ".Fcrash";
            std::ofstream dumpfile(outfilestream.str().c_str(), std::ios::out | std::ios::trunc);
            std::vector<double> X_which_caused_problem;

            dumpfile << "Bad entry in F in process rank" << std::endl;
            dumpfile << "Bad entry is #" << failure_entry << std::endl;
            dumpfile << this->options.mission_name << std::endl;
            dumpfile << std::endl;
            for (size_t Xentry = 0; Xentry < this->Xdescriptions.size(); ++Xentry)
            {
                
                if (this->Xdescriptions[Xentry].find("epoch") < 1024 || this->Xdescriptions[Xentry].find("time") < 1024)
                {
                    dumpfile << "X[" << Xentry << "] (" << this->Xdescriptions[Xentry] << ") = " << X[Xentry] _GETVALUE / 86400.0 << std::endl;
                    X_which_caused_problem.push_back(X[Xentry] _GETVALUE / 86400.0);
                }
                else
                {
                    dumpfile << "X[" << Xentry << "] (" << this->Xdescriptions[Xentry] << ") = " << X[Xentry] _GETVALUE << std::endl;
                    X_which_caused_problem.push_back(X[Xentry] _GETVALUE);
                }
            }
            dumpfile << std::endl;
            for (size_t Fentry = 0; Fentry < this->Fdescriptions.size(); ++Fentry)
            {
                dumpfile << "F[" << Fentry << "] (" << this->Fdescriptions[Fentry] << ") = " << F[Fentry] _GETVALUE << std::endl;
            }
            dumpfile.close();

            missionoptions optionscopy = this->options;
            for (size_t Xindex = 0; Xindex < this->total_number_of_NLP_parameters; ++Xindex)
                optionscopy.trialX.push_back({ this->Xdescriptions[Xindex], X_which_caused_problem[Xindex] });

            optionscopy.run_inner_loop = InnerLoopSolverType::NLP;
            std::stringstream optionsfilenamestream;
            optionsfilenamestream << this->options.working_directory << "/" << this->options.mission_name << "_bad_objective_constraint_vector.emtgopt";
            optionscopy.write(optionsfilenamestream.str(), optionscopy.print_only_non_default_options);

            return true;
        }

        //if you got this far there was no problem
        return false;
    }//end check_outputs_for_invalid_entries()

    bool problem::check_derivatives_for_invalid_entries(const std::vector<doubleType>& X, const std::vector<double>& G, const bool& PrintAnyway)
    {
        int failure_entry = -1;
        std::string failure_type = "";


        //look for failures
        for (size_t Gentry = 0; Gentry < this->Gdescriptions.size(); ++Gentry)
        {
            if (std::isnan(G[Gentry]))
            {
                failure_entry = Gentry;
                failure_type = "NaN";
                break;
            }
            else if (!std::isfinite(G[Gentry]))
            {
                failure_entry = Gentry;
                failure_type = "inf";
                break;
            }
        }

        //if a failure is found, write out a dump and throw an error
        if (failure_entry >= 0 || PrintAnyway)
        {
            std::stringstream outfilestream;
            boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
            std::stringstream timestream;
            timestream << static_cast<int>(now.date().month()) << now.date().day() << now.date().year() << "_" << now.time_of_day().hours() << now.time_of_day().minutes() << now.time_of_day().seconds();
            outfilestream << this->options.working_directory << "/" << "bad_objective_constraint_vector_" << timestream.str() << ".Gcrash";
            std::ofstream dumpfile(outfilestream.str().c_str(), std::ios::out | std::ios::trunc);
            std::vector<double> X_which_caused_problem;

            dumpfile << "Bad entry in G" << std::endl;
            dumpfile << "Bad entry is #" << failure_entry << std::endl;
            dumpfile << this->options.mission_name << std::endl;
            dumpfile << std::endl;
            for (size_t Xentry = 0; Xentry < this->Xdescriptions.size(); ++Xentry)
            {

                if (this->Xdescriptions[Xentry].find("epoch") < 1024 || this->Xdescriptions[Xentry].find("time") < 1024)
                {
                    dumpfile << "X[" << Xentry << "] (" << this->Xdescriptions[Xentry] << ") = " << X[Xentry] _GETVALUE / 86400.0 << std::endl;
                    X_which_caused_problem.push_back(X[Xentry] _GETVALUE / 86400.0);
                }
                else
                {
                    dumpfile << "X[" << Xentry << "] (" << this->Xdescriptions[Xentry] << ") = " << X[Xentry] _GETVALUE << std::endl;
                    X_which_caused_problem.push_back(X[Xentry] _GETVALUE);
                }
            }
            dumpfile << std::endl;
            for (size_t Gentry = 0; Gentry < this->Gdescriptions.size(); ++Gentry)
            {
                dumpfile << "G[" << Gentry << "] (" << this->Gdescriptions[Gentry] << ") = " << G[Gentry] << std::endl;
            }
            dumpfile.close();


            missionoptions optionscopy = this->options;
            for (size_t Xindex = 0; Xindex < this->total_number_of_NLP_parameters; ++Xindex)
                optionscopy.trialX.push_back({ this->Xdescriptions[Xindex], X_which_caused_problem[Xindex] });

            optionscopy.run_inner_loop = InnerLoopSolverType::NLP;
            std::stringstream optionsfilenamestream;
            optionsfilenamestream << this->options.working_directory << "/" << this->options.mission_name << "_bad_Jacobian_vector.emtgopt";
            optionscopy.write(optionsfilenamestream.str(), !optionscopy.print_only_non_default_options);

            return true;
        }

        //if you got this far there was no problem
        return false;
    }//end check_derivatives_for_invalid_entries()

    //function to determine which constraints are equalities vs inequalities
    void problem::locate_equality_constraints()
    {
        this->F_equality_or_inequality.clear();
        for (size_t Findex = 1; Findex < this->Fdescriptions.size(); ++Findex)
            this->F_equality_or_inequality.push_back(this->Fupperbounds[Findex] - this->Flowerbounds[Findex] <= this->options.snopt_feasibility_tolerance ? true : false);
    }//end locate_equality_constraints()

    //function to scale to [0, 1] hypercube
    void problem::scale_to_unit_hypercube()
    {
        //first rescale A
        for (size_t Aindex = 0; Aindex < this->A.size(); ++Aindex)
        {
            size_t Xindex = this->jAvar[Aindex];
            A[Aindex] *= (this->Xupperbounds[Xindex] - this->Xlowerbounds[Xindex]) / this->X_scale_factors[Xindex];
        }

        //then rescale the decision variables, which also rescales G
        for (size_t Xindex = 0; Xindex < this->total_number_of_NLP_parameters; ++Xindex)
        {
            this->X_scale_factors[Xindex] = this->Xupperbounds[Xindex] - this->Xlowerbounds[Xindex];
        }

    }//end scale_to_unit_hypercube

    void problem::what_the_heck_am_I_called(const SolutionOutputType& mySolutionType, const size_t number_of_hops)
    {
        std::string path = options.working_directory + "//";

        std::string base_name; 
        if (this->options.override_default_output_file_name)
        {
            base_name = this->options.forced_output_file_name;
        }
        else if (this->options.short_output_file_names)
        {
            base_name = this->options.mission_name;
        }
        else
        {
            base_name = this->options.mission_name + "_" + this->options.description;
        }

        switch (mySolutionType)
        {
            case SolutionOutputType::FAILURE:
            {
                this->options.outputfile = path + "FAILURE_" + base_name + ".emtg";
                break;
            }
            case SolutionOutputType::HOP:
            {
                this->options.outputfile = path + "Hop" + std::to_string(number_of_hops) + ".emtg";
                break;
            }
            case SolutionOutputType::SUCCESS:
            {
                this->options.outputfile = path + base_name + ".emtg";
                break;
            }
        }
    }//end what_the_heck_am_I_called()
} /* namespace EMTG */

