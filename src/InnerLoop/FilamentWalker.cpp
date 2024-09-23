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

//filament walker
//math by Arnold Englander
//code by Jacob Englander

#include "FilamentWalker.h"
#include "EMTG_math.h"

#include <iostream>
#include <fstream>

namespace EMTG
{
    namespace Solvers
    {
        //constructors
        FilamentWalker::FilamentWalker() :
            max_remedial_starts(10),
            max_routine_start_attempts(10000),
            max_routine_start_failures_before_reversal(10),
            max_slide_steps(1000),
            max_routine_start_hop_size(1.0e-4),
            equality_constraint_tolerance(1.0e-5),
            WalkSoftly(false),
            myRandomSeed(0)
        {
            this->RNG.seed(myRandomSeed);

            /*std::string seedfilename = this->myProblem->options.working_directory + "/" + this->myProblem->options.mission_name + "random_seed";
            std::fstream seedstream = std::fstream(seedfilename, std::ios::trunc);
            seedstream << this->myRandomSeed << std::endl;
            seedstream.close();*/
        }

        FilamentWalker::FilamentWalker(problem* myProblem,
                                       SNOPT_interface* mySNOPT) :
            FilamentWalker()
        {
            this->myProblem = myProblem;
            this->pareto_alpha = myProblem->options.MBH_Pareto_alpha;
            this->nX = this->myProblem->total_number_of_NLP_parameters;
            this->X_before_hop.resize(this->nX, 0.0);
            this->X_after_hop.resize(this->nX, 0.0);
            this->X_after_slide.resize(this->nX, 0.0);
            this->dX_hop.resize(this->nX, 0.0);
            this->dX_slide.resize(this->nX, 0.0);
            this->dX_iteration.resize(this->nX, 0.0);

            this->mySNOPT = mySNOPT;
        }

        //execute
        void FilamentWalker::walk()
        {
            //initialize
            this->current_remedial_start_index = 0;
            this->FilamentWalker_start_time = time(NULL);

            //create the output file
            std::string outputfilename = this->myProblem->options.working_directory + "/" + this->myProblem->options.mission_name + "_filament.emtg_filamentarchive";
            this->outputfile = std::ofstream(outputfilename, std::ios::trunc);
            outputfile.width(19); outputfile << "start_type";
            outputfile.width(19); outputfile << "remedial_index";
            outputfile.width(19); outputfile << "routine_index";
            outputfile.width(19); outputfile << "distance_stepped";
            outputfile.width(19); outputfile << "distance_hopped";
            outputfile.width(19); outputfile << "distance_slid";
            outputfile.width(19); outputfile << "eq_violation";
            outputfile.width(19); outputfile << "ineq_violation";
            outputfile.width(19); outputfile << "objective";
            outputfile.width(19); outputfile << "timestamp";

            /*for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
            {
                outputfile.width(19); outputfile << this->myProblem->Xdescriptions[Xindex];
            }*/
            outputfile << std::endl;

            //do walking things
            while (this->current_remedial_start_index < this->max_remedial_starts)
            {
                //conduct a remedial start, and if successful, do the rest
                if (this->remedialStart())
                {
                    while (this->current_routine_start_index < this->max_routine_start_attempts)
                    {
                        this->routineStart();
                    }

                    this->list_of_filaments.push_back(this->list_of_points_on_current_filament);
                    this->list_of_points_on_current_filament.clear();
                }

                //increment the remedial start counter
                ++this->current_remedial_start_index;
            }

            //close the output file
            this->outputfile.close();
        }//end walk()

        //output
        void FilamentWalker::output(std::tuple< EMTG_innerloop_solution, double, double, double, time_t >& solution)
        {
            EMTG_innerloop_solution& mySolution = std::get<0>(solution);
            if (this->current_routine_start_index)
            {
                if (this->reversal_occurred)
                {
                    outputfile.width(19); outputfile << "reverse";
                }
                else
                {
                    outputfile.width(19); outputfile << "routine";
                }
            }
            else
            {
                outputfile.width(19); outputfile << "remedial";
            }
            outputfile.width(19); outputfile << this->current_remedial_start_index;
            outputfile.width(19); outputfile << this->current_routine_start_index;


            double distance_since_previous_solution = std::get<1>(solution);
            double distance_hopped = std::get<2>(solution);
            double distance_slid = std::get<3>(solution);

            outputfile.width(19); outputfile << distance_since_previous_solution;
            outputfile.width(19); outputfile << distance_hopped;
            outputfile.width(19); outputfile << distance_slid;

            outputfile.width(19); outputfile << mySolution.getequality_constraint_violation();
            outputfile.width(19); outputfile << mySolution.getinequality_constraint_violation();
            outputfile.width(19); outputfile << mySolution.getobjective_function_value();
            outputfile.width(19); outputfile << std::get<4>(solution); //timestamp


            for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
            {
                outputfile.width(19); outputfile << mySolution.getX(Xindex);
            }
            outputfile << std::endl;
        }//end output()

        //internal methods
        bool FilamentWalker::remedialStart()
        {
            bool firstHop = true;
            double best_feasibility = 1.0e+100;
            size_t start_attempts = 0;
            this->current_routine_start_index = 0;
            this->routine_start_failures_since_reversal = 0;
            ++this->current_remedial_start_index;
            this->dX_iteration.assign(this->nX, 0.0);

            while (best_feasibility > this->equality_constraint_tolerance
                && start_attempts < this->max_routine_start_attempts)
            {
                //increment counter
                ++start_attempts;

                //hop
                this->remedialHop(firstHop);

                //slide
                if (this->slide())
                {
                    //get feasibility
                    double feasibility = this->mySNOPT->getF().front();

                    if (feasibility < best_feasibility)
                    {
                        best_feasibility = feasibility;
                        this->X = this->X_after_slide;
                    }
                }
            }

            if (best_feasibility < this->equality_constraint_tolerance)
            {
                std::string solution_name = "R" + std::to_string(this->current_remedial_start_index) + "r0";

                this->list_of_points_on_current_filament.push_back(std::make_tuple(
                                                                   EMTG_innerloop_solution(myProblem,
                                                                                           this->X_after_slide,
                                                                                           this->myProblem->F,
                                                                                           solution_name),
                                                                   -1,
                                                                   -1,
                                                                   -1,
                                                                   time(NULL) - this->FilamentWalker_start_time));

                this->output(this->list_of_points_on_current_filament.back());

                return true;
            }
            else
                return false;
        }//end remedialStart

        void FilamentWalker::remedialHop(bool& firstHop)
        {
            //if this is the first hop of the current remedial start, then use a uniform distribution in the whole space
            //otherwise use a Pareto distribution (Englander and Englander 2014)
            if (firstHop)
            {
                for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
                    this->X_after_hop[Xindex] = this->DoubleDistribution(this->RNG);

                this->X = this->X_after_hop;

                firstHop = false;
            }
            else
            {
                for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
                {
                    double r = ((this->pareto_alpha - 1.0) / math::SMALL) / pow((math::SMALL / (math::SMALL + DoubleDistribution(RNG))), -this->pareto_alpha);
                    int s = DoubleDistribution(RNG) > 0.5 ? 1 : -1;

                    double step_size = s * r;

                    if (isfinite(step_size))
                        this->X_after_hop[Xindex] = this->X[Xindex] + step_size;
                }
            }
        }

        bool FilamentWalker::routineStart()
        {
            //increment index
            ++this->current_routine_start_index;

            //hop
            this->routineHop();

            //slide - if it works, add a point to the current filament
            if (this->slide())
            {
                this->reversal_occurred = false;
                this->routine_start_failures_since_reversal = 0;

                for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
                    this->dX_slide[Xindex] = this->X_after_slide[Xindex] - this->X_after_hop[Xindex];

                double problem_size = sqrt((double)this->nX);

                double distance_slid = math::norm(this->dX_slide);

                double distance_hopped = math::norm(this->dX_hop);

                for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
                    this->dX_iteration[Xindex] = this->X_after_slide[Xindex] - this->X_before_hop[Xindex];
                double distance_since_previous_solution = math::norm(this->dX_iteration);

                std::string solution_name = "R" + std::to_string(this->current_remedial_start_index) + "r" + std::to_string(this->current_routine_start_index);

                this->list_of_points_on_current_filament.push_back(std::make_tuple(
                                                                                   EMTG_innerloop_solution(myProblem,
                                                                                   this->X_after_slide, 
                                                                                   this->myProblem->F,
                                                                                   solution_name),
                                                                   distance_since_previous_solution,
                                                                   distance_hopped,
                                                                   distance_slid,
                                                                   time(NULL) - this->FilamentWalker_start_time));

                this->output(this->list_of_points_on_current_filament.back());

                this->X = this->X_after_slide;
                return true;
            }//end if succcessful slide
            else//unsuccessful slide
            {
                ++this->routine_start_failures_since_reversal;
            }

            return false;
        }//end routineStart

        void FilamentWalker::routineHop()
        {
            //start from the current point
            this->X_before_hop = this->X;

            int direction_bias;
            if (this->routine_start_failures_since_reversal > this->max_routine_start_failures_before_reversal)
            {
                this->reversal_occurred = true;
                this->routine_start_failures_since_reversal = 0;
                direction_bias = -1;
            }
            else
            {
                direction_bias = 1;
            }

            //generate the hop
            for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
            {
                double dx_bias = (this->X_before_hop[Xindex] > (1.0 - this->max_routine_start_hop_size) ? 
                    0.0 : (this->X_before_hop[Xindex] < this->max_routine_start_hop_size ? 
                        0.0 : this->dX_iteration[Xindex]));

                //clear the bias if an inequality constraint is active that uses this variable
                for (size_t Findex = 1; Findex < this->mySNOPT->getnF(); ++Findex)
                {
                    if (this->mySNOPT->getF()[Findex] < this->mySNOPT->getFlowerbounds()[Findex]
                        || this->mySNOPT->getF()[Findex] > this->mySNOPT->getFupperbounds()[Findex])
                    {
                        for (size_t Gindex = 0; Gindex < this->mySNOPT->getnG(); ++Gindex)
                        {
                            if (this->mySNOPT->getiGfun()[Gindex] == Findex
                                && this->mySNOPT->getjGvar()[Gindex] == Xindex)
                            {
                                dx_bias = 0.0;
                                break;
                            }
                        }
                    }
                }

                this->dX_hop[Xindex] = dx_bias * direction_bias + (this->DoubleDistribution(this->RNG) * 2.0 - 1.0) * this->max_routine_start_hop_size;
                this->X_after_hop[Xindex] = this->X_before_hop[Xindex] + this->dX_hop[Xindex];
            }
        }//end hop()

        bool FilamentWalker::slide()
        {
            this->mySNOPT->setX0_scaled(this->X_after_hop);
            this->mySNOPT->run_NLP();
            this->X_after_slide = this->mySNOPT->getX_scaled();

            double current_equality_constraint_violation = this->mySNOPT->getF().front();

            if (!this->WalkSoftly)
                std::cout << current_equality_constraint_violation << std::endl;

            //store this point if it satisfies the equality constraints
            if (current_equality_constraint_violation < this->equality_constraint_tolerance)
            {
                try
                {
                    myProblem->evaluate(this->mySNOPT->getX_unscaled());
                    return true;
                }
                catch (std::exception &error)
                {
                    std::cout << error.what() << std::endl;
                    return false;
                }
            }
            else
                return false;
        }//end slide()
    }//close namespace solvers
}//close namespace EMTG