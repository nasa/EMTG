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

//SNOPT interface

#include "SNOPT_interface.h"
#include <time.h>
#include <cmath>

namespace EMTG
{
    namespace Solvers
    {
        SNOPT_interface::SNOPT_interface(problem* myProblem,
            const NLPoptions& myOptions) :
            NLP_interface::NLP_interface(myProblem, myOptions)
        {
            //set up SNOPT-specific stuff
            //take advantage of what the base constructor has already done for us, but recognize that we have to change some things
            //obnoxiously, we have to change from size_t to SNOPT_INT_TYPE for iGfun/jGvar/iAfun/iAvar
            this->iAfun.resize(this->nA);
            this->jAvar.resize(this->nA);
            this->A = this->myProblem->A;
            this->iGfun.resize(this->nG);
            this->jGvar.resize(this->nG);
            this->lenA = this->nA;
            this->lenG = this->nG;
            if (this->myOptions.get_SolverMode() == NLPMode::FilamentFinder)
            {
                for (size_t Aindex = 0; Aindex < this->nA; ++Aindex)
                {
                    this->iAfun[Aindex] = this->myProblem->iAfun[Aindex];
                    this->jAvar[Aindex] = this->myProblem->jAvar[Aindex];
                }
                for (size_t Gindex = 0; Gindex < this->nG; ++Gindex)
                {
                    this->iGfun[Gindex] = NLP_interface::iGfun[Gindex];
                    this->jGvar[Gindex] = NLP_interface::jGvar[Gindex];
                }
            }
            else
            {
                for (size_t Aindex = 0; Aindex < this->nA; ++Aindex)
                {
                    this->iAfun[Aindex] = this->myProblem->iAfun[Aindex];
                    this->jAvar[Aindex] = this->myProblem->jAvar[Aindex];
                }
                for (size_t Gindex = 0; Gindex < this->nG; ++Gindex)
                {
                    this->iGfun[Gindex] = this->myProblem->iGfun[Gindex];
                    this->jGvar[Gindex] = this->myProblem->jGvar[Gindex];
                }
            }
            
            this->xmul.resize(this->nX, 0.0);
            this->xstate.resize(this->nX, 0);
            this->Fstate.resize(this->nF, 0);

#ifdef SNOPT72
            this->nxnames = 1;
            this->nFnames = 1;
            this->xnames = new char[nxnames * 8];
            this->Fnames = new char[nFnames * 8];
#endif

            this->ObjRow = 0;
            this->ObjAdd = 0;

            //we don't have any real intuition for what the intial "states" of the F functions will be
            //the same applies for the constraint Lagrange multipliers, so set all of the Fstates = 0 and Fmuls to 0.0
            this->Fstate = std::vector<SNOPT_INT_TYPE>(this->nF, 0);
            this->Fmul = std::vector<SNOPT_DOUBLE_TYPE>(this->nF, 0.0);
#ifndef SNOPT76
            this->mySNOPT.setProblemSize(this->nX, this->nF);
            this->mySNOPT.setObjective(ObjRow, ObjAdd);
#else
			this->mySNOPT.initialize(this->myOptions.get_output_file_path().c_str(),1);
#endif 
			this->mySNOPT.setUserspace((SNOPT_INT_TYPE*)&NLP_start_time, 500, (SNOPT_DOUBLE_TYPE*)this, 500);

#if defined SNOPT75 //in SNOPT 7.5 and later, NeA and NeG are set along with A and G
            this->mySNOPT.setA(this->nA, this->nA, this->iAfun.data(), this->jAvar.data(), this->A.data());
            this->mySNOPT.setG(this->nG, this->nG, this->iGfun.data(), this->jGvar.data());
#elif defined SNOPT72
            this->mySNOPT.setNeA(this->nA);
            this->mySNOPT.setNeG(this->nG);
            this->mySNOPT.setA(lenA, iAfun.data(), jAvar.data(), A.data());
            this->mySNOPT.setG(lenG, iGfun.data(), jGvar.data());
#endif

            this->mySNOPT.setIntParameter("Total real workspace", 500 * (this->nF * this->nF));
            this->mySNOPT.setIntParameter("Total integer workspace", 500 * (this->nF * this->nF));
            //this->mySNOPT.setIntParameter("Total character workspace", lenA + lenG);
#ifdef SNOPT72
            this->mySNOPT.setXNames(xnames, nxnames);
            this->mySNOPT.setFNames(Fnames, nFnames);
#endif
            this->mySNOPT.setProbName((char*)"EMTG");
#ifndef SNOPT76
            this->mySNOPT.setUserFun(SNOPT_user_function);
#endif
            this->mySNOPT.setIntParameter((char*)"Iterations limit", 10000 * this->myOptions.get_major_iterations_limit());
            this->mySNOPT.setIntParameter((char*)"Major iterations limit", this->myOptions.get_major_iterations_limit());
            this->mySNOPT.setIntParameter((char*)"Minor iterations limit", this->myOptions.get_minor_iterations_limit());
            this->mySNOPT.setRealParameter((char*)"Major step limit", this->myOptions.get_max_step());
            //this->mySNOPT.setRealParameter((char*)"Function precision", 1.0e-8);

            this->mySNOPT.setIntParameter((char*)"Derivative option", 1);
            this->mySNOPT.setIntParameter((char*)"Minor print level", 0);
            this->mySNOPT.setRealParameter((char*)"Major feasibility tolerance", this->myOptions.get_feasibility_tolerance());

            //recommended settings from Michael Saunders to avoid crash in SNOPT 7.4
            //this->mySNOPT.setRealParameter((char*)"LU factor tolerance", 2.0);
            //this->mySNOPT.setRealParameter((char*)"LU update tolerance", 2.0);
            //this->mySNOPT.setIntParameter((char*)"Factorization frequency", 25);
            //this->mySNOPT.setParameter((char*)"LU complete pivoting"); //this seems to cause problems

            this->mySNOPT.setIntParameter((char*)"Major Print Level", 1);
            this->mySNOPT.setRealParameter((char*)"Major optimality tolerance", this->myOptions.get_optimality_tolerance());
            if (this->myOptions.get_check_derivatives())
            {
                this->mySNOPT.setIntParameter((char*)"Print file", 1);
                this->mySNOPT.setIntParameter((char*)"Summary file", 1);
                this->mySNOPT.setIntParameter((char*)"System information", 1);
                this->mySNOPT.setIntParameter((char*)"Verify level", 3); //0 = cheap test 1 = individual gradients checked (OK or BAD) 2 = Individual columns of the Jacobian are checked 3 = 1 and 2 happen -1 = Derivative checking is disabled
            }
            if (this->myOptions.get_quiet_NLP())
            {
                this->mySNOPT.setIntParameter((char*)"Major Print Level", 0);
                this->mySNOPT.setIntParameter((char*)"Minor print level", 0);
                this->mySNOPT.setIntParameter((char*)"Print No", 0);
                this->mySNOPT.setIntParameter((char*)"Summary file", 0);
                this->mySNOPT.setParameter((char*)"Suppress parameters");
            }

            if (this->myOptions.get_SolverMode())
                this->mySNOPT.setParameter((char*)"Minimize");
            else
                this->mySNOPT.setParameter((char*)"Feasible point");

            //read a specs file if there is one
            std::string specs_file_path = this->myOptions.get_specs_file_path();
            if (specs_file_path != "")
                this->mySNOPT.setSpecsFile(specs_file_path.c_str());

            //set output file path if there is one
            std::string output_file_path = this->myOptions.get_output_file_path();
            if (output_file_path != "")
                this->mySNOPT.setPrintFile(output_file_path.c_str());

            //set first feasibility flag
            this->first_feasibility = false;
        }//end constructor

        void SNOPT_interface::run_NLP(const bool& X0_is_scaled)
        {
            //compute the scaled decision vector
            if (!X0_is_scaled)
                this->scaleX0();
            else
                this->unscaleX0();

            this->X_scaled = this->X0_scaled;

            //check that the user function works, and return values of F
            this->myProblem->evaluate(this->X0_unscaled, this->F, this->myProblem->G, false);

            //set the SNOPT bounds equal to the problem bounds
            for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
            {
                this->Xupperbounds[Xindex] = (this->myProblem->Xupperbounds[Xindex] - this->myProblem->Xlowerbounds[Xindex]) / this->myProblem->X_scale_factors[Xindex];
                this->Xlowerbounds[Xindex] = 0.0;
            }

            //reset multipliers
            for (size_t Findex = 0; Findex < this->nF; ++Findex)
            {
                this->Fstate[Findex] = 0.0;
                this->Fmul[Findex] = 0.0;
                //this->F[Findex] = 0.0;

                ////rules on setting Fstate from page 17 of sndoc7
                //if (this->F[Findex] <= this->myProblem->Flowerbounds[Findex])
                //    this->Fstate[Findex] = 4;
                //else if (this->F[Findex] >= this->myProblem->Fupperbounds[Findex])
                //    this->Fstate[Findex] = 5;
                //else
                //    this->Fstate[Findex] = 3;
            }
            for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
            {
                this->xstate[Xindex] = 0.0;
                this->xmul[Xindex] = 0.0;

                ////rules on setting xstate from page 17 of sndoc7
                //if (this->X_unscaled[Xindex] <= this->myProblem->Xlowerbounds[Xindex])
                //    this->xstate[Xindex] = 4;
                //else if (this->X_unscaled[Xindex] >= this->myProblem->Xupperbounds[Xindex])
                //    this->xstate[Xindex] = 5;
                //else
                //    this->xstate[Xindex] = 3;
            }

#ifdef AD_INSTRUMENTATION
            std::vector<double> X_scaled_double(this->nX);
            std::vector<double> Fdouble(this->nF);
            for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
                X_scaled_double[Xindex] = this->X_scaled[Xindex] _GETVALUE;

            this->mySNOPT.setX(X_scaled_double.data(), this->Xlowerbounds.data(), this->Xupperbounds.data(), this->xmul.data(), this->xstate.data());
            this->mySNOPT.setF(Fdouble.data(), this->Flowerbounds.data(), this->Fupperbounds.data(), this->Fmul.data(), this->Fstate.data());
#else
#ifndef SNOPT76
            this->mySNOPT.setX(this->X_scaled.data(), this->Xlowerbounds.data(), this->Xupperbounds.data(), this->xmul.data(), this->xstate.data());
            this->mySNOPT.setF(this->F.data(), this->Flowerbounds.data(), this->Fupperbounds.data(), this->Fmul.data(), this->Fstate.data());
#endif
#endif

            //run SNOPT
            this->NLP_start_time = time(NULL);
            this->feasibility_metric_NLP_incumbent = 1.0e+101;
            this->J_NLP_incumbent = math::LARGE;
            this->movie_frame_count = 0;

#ifdef SNOPT76
			int nS = 0;
			int nInf = 0;
			double sInf = 0.0;
			this->inform = mySNOPT.solve(0, this->nF, this->X_scaled.size(), this->ObjAdd,
	                                     this->ObjRow, SNOPT_user_function,
	     								 this->iAfun.data(), this->jAvar.data(), this->A.data(), this->nA,
	    								 this->iGfun.data(), this->jGvar.data(), this->nG,
	     								 this->Xlowerbounds.data(), this->Xupperbounds.data(), this->Flowerbounds.data(), this->Fupperbounds.data(),
	     							     this->X_scaled.data(), this->xstate.data(), this->xmul.data(),
	     							     this->F.data(), this->Fstate.data(), this->Fmul.data(),
	     							     nS, nInf, sInf);
#else
            this->inform = mySNOPT.solve(0);			
#endif
            //unscale the various things that might need to be uncaled
            this->unscaleX();

            //perform a feasibility check
            this->myProblem->check_feasibility(this->X_unscaled,
                this->F,
                this->worst_decision_variable,
                this->worst_constraint,
                this->feasibility_metric,
                this->normalized_feasibility_metric,
                this->distance_from_equality_filament,
                this->decision_vector_feasibility_metric);

            double worst_feasibility = fmax(normalized_feasibility_metric, decision_vector_feasibility_metric);

            //adopt the incumbent if appropriate
            if (this->myOptions.get_enable_NLP_chaperone())
            {
                this->unscaleX_NLP_incumbent();
                if (worst_feasibility < this->myOptions.get_feasibility_tolerance()
                    && this->feasibility_metric_NLP_incumbent < this->myOptions.get_feasibility_tolerance())
                {
                    // Both incumbent point and NLP exit point are feasible
                    if (this->J_NLP_incumbent < this->F.front())
                    {
                        this->X_unscaled = this->X_NLP_incumbent_unscaled;
                        this->X_scaled = this->X_NLP_incumbent_scaled;
                        this->F = this->F_NLP_incumbent;

                        std::cout << "NLP incumbent point and exit point are feasible and incumbent point is superior to exit point." << std::endl;
                    }
                    else
                    {
                        std::cout << "NLP incumbent point and exit point are feasible and exit point is superior to incumbent point." << std::endl;
                    }
                }
                //if the incumbent point is feasible but the exit point is not
                else if (this->feasibility_metric_NLP_incumbent < worst_feasibility
                    && this->feasibility_metric_NLP_incumbent < this->myOptions.get_feasibility_tolerance())
                {
                    this->X_unscaled = this->X_NLP_incumbent_unscaled;
                    this->X_scaled = this->X_NLP_incumbent_scaled;
                    this->F = this->F_NLP_incumbent;

                    std::cout << "NLP incumbent point was feasible and the exit point was not." << std::endl;
                }
                else if (this->feasibility_metric_NLP_incumbent < worst_feasibility)
                {
                    this->X_unscaled = this->X_NLP_incumbent_unscaled;
                    this->X_scaled = this->X_NLP_incumbent_scaled;
                    this->F = this->F_NLP_incumbent;

                    std::cout << "NLP incumbent point was infeasible but less infeasible than the exit point." << std::endl;
                }
                else
                {
                    std::cout << "NLP exit point was superior to incumbent point." << std::endl;
                }
            }//end NLP chaperone
        }

        //SNOPT user function
#ifdef SNOPT72
        int
#else
        void
#endif
        SNOPT_interface::SNOPT_user_function(SNOPT_INT_TYPE    *Status, SNOPT_INT_TYPE *n,   SNOPT_DOUBLE_TYPE SNOPT_X[],
                                                         SNOPT_INT_TYPE    *needF,  SNOPT_INT_TYPE *neF, SNOPT_DOUBLE_TYPE SNOPT_F[],
                                                         SNOPT_INT_TYPE    *needG,  SNOPT_INT_TYPE *neG, SNOPT_DOUBLE_TYPE SNOPT_G[],
                                                         char               *cu,    SNOPT_INT_TYPE *lencu,
                                                         SNOPT_INT_TYPE    iu[],    SNOPT_INT_TYPE *leniu,
                                                         SNOPT_DOUBLE_TYPE ru[],    SNOPT_INT_TYPE *lenru)
        {
            //Step 1: get pointer to self class
            SNOPT_interface* self = (SNOPT_interface*) ru;

            //Step 2: unwrap the inputs
            for (size_t Xindex = 0; Xindex < self->nX; ++Xindex)
                self->X_scaled[Xindex] = SNOPT_X[Xindex];

            //Step 3: process the SNOPT decision vector
            self->unscaleX();

            //Step 4: Call the objective and constraint function
#ifdef SAFE_SNOPT
            try
            {
#endif
                if (self->myOptions.get_SolverMode() == NLPMode::FilamentFinder) //single-objective unconstrained, special Jacobian
                {

                    self->myProblem->evaluate(self->X_unscaled, self->myProblem->F, self->myProblem->G, *needG);
                    
                    //then sum up the squares of the equality constraint values
                    self->F.front() = 0.0;
                    for (size_t Findex = 1; Findex < self->myProblem->total_number_of_constraints; ++Findex)
                    {
                        if (self->myProblem->F_equality_or_inequality[Findex - 1])
                            self->F.front() += self->myProblem->F[Findex] * self->myProblem->F[Findex];
                    }


                    //compute dfdX
                    for (size_t Gindex = 0; Gindex < self->nG; ++Gindex)
                    {
                        self->G[Gindex] = 0.0;
                    }

                    size_t Problem_nG = self->myProblem->Gdescriptions.size();
                    for (size_t Gindex = 0; Gindex < Problem_nG; ++Gindex)
                    {
                        size_t Findex = self->myProblem->iGfun[Gindex];
                        size_t Xindex = self->myProblem->jGvar[Gindex];

                        if (Findex > 0)
                        {
                            if (self->myProblem->F_equality_or_inequality[Findex - 1])
                                self->G[Xindex] += 2.0 * self->myProblem->F[Findex] * self->myProblem->G[Gindex];
                        }
                    }

                    //inequality constraints
                    for (size_t Findex = 1; Findex <= self->myProblem->F_indices_of_filament_critical_inequality_constraints.size(); ++Findex)
                    {
                        self->F[Findex] = self->myProblem->F[self->myProblem->F_indices_of_filament_critical_inequality_constraints[Findex - 1]];
                    }
                    for (size_t Gindex = self->nX; Gindex < self->nG; ++Gindex)
                    {
                        self->G[Gindex] = self->myProblem->G[self->original_G_indices_of_filament_critical_inequality_constraints[Gindex - self->nX]];
                    }

                }
                else //regular problem
                {
                    self->myProblem->evaluate(self->X_unscaled, self->F, self->G, *needG);
                }
#ifdef SAFE_SNOPT
            }
            catch (int errorcode)
            {
                if (errorcode == 13)  //integration step error
                    *Status = -1;
            }
            catch (std::runtime_error runtime_error)
            {
                
                if (!self->myOptions.get_quiet_NLP())
                {
                    std::cout << runtime_error.what() << std::endl;
                    std::cout << std::endl;
                }
                *Status = -1;
            }
#endif

            //Step 5: wrap the outputs
            for (size_t Findex = 0; Findex < self->nF; ++Findex)
                SNOPT_F[Findex] = self->F[Findex];
            for (size_t Gindex = 0; Gindex < self->nG; ++Gindex)
                SNOPT_G[Gindex] = self->G[Gindex];

            //Step 6: NLP chaperone
            if (*needG && self->myOptions.get_enable_NLP_chaperone())
            {
                //first compute the current feasibility
                double f_current;
                double f_abs_current;
                double distance_from_equality_filament;
                size_t worst_constraint;
                size_t worst_decision_variable;
                double decision_variable_feasibility_metric;

                try
                {
                    self->myProblem->check_feasibility(self->X_unscaled,
                        self->F,
                        worst_decision_variable,
                        worst_constraint,
                        f_abs_current,
                        f_current,
                        distance_from_equality_filament,
                        decision_variable_feasibility_metric,
                        true);
                }
                catch (std::runtime_error runtime_error)
                {

                    if (!self->myOptions.get_quiet_NLP())
                    {
                        std::cout << runtime_error.what() << std::endl;
                        std::cout << std::endl;
                    }
                    *Status = -1;
                }


                //stop immediately if goal met
                if (self->myOptions.get_stop_on_goal_attain())
                {
					std::cout << "snoptf0:"<<self->myProblem->getUnscaledObjective() << std::endl;
					std::cout << "goal:"<<self->myOptions.get_objective_goal() << std::endl;
                    if (f_current < self->myOptions.get_feasibility_tolerance() && self->myProblem->getUnscaledObjective() < self->myOptions.get_objective_goal())
                    {						
                        if (!self->myOptions.get_quiet_NLP())
                            std::cout << "NLP goal satisfied, exiting NLP" << std::endl;
                        *Status = -2;
                    }
                }

                //if the current point is feasible and better than the incumbent
                //raw pointer math for speed
                if (fmax(f_current, decision_variable_feasibility_metric) < self->myOptions.get_feasibility_tolerance()
                    && self->feasibility_metric_NLP_incumbent < self->myOptions.get_feasibility_tolerance()
                    && decision_variable_feasibility_metric == 0.0)
                {
                    // Both the incumbent and the current point are feasible
                    if (self->F.front() < self->J_NLP_incumbent)
                    {
                        // The current point is more optimal than the incumbent, so update the incumbent
                        self->J_NLP_incumbent = self->F.front();
                        self->feasibility_metric_NLP_incumbent = fmax(f_current, decision_variable_feasibility_metric);

                        self->X_NLP_incumbent_unscaled = self->X_unscaled;
                        self->X_NLP_incumbent_scaled = self->X_scaled;
                        self->F_NLP_incumbent = self->F;
                    }
					
					// Check if this is the first time things have been feasible. It shouldnt be, because the incumbent was already feasible
					// But do it anyway, even though this code should never be called.
					if (self->first_feasibility == false)
					{
						self->first_feasibility = true;
						self->myProblem->Xopt = self->X_NLP_incumbent_unscaled;
						self->myProblem->F = self->F_NLP_incumbent;
	                    self->myProblem->evaluate(self->X_NLP_incumbent_unscaled, self->F_NLP_incumbent, self->G, false);

                        self->myProblem->what_the_heck_am_I_called(SolutionOutputType::SUCCESS);
						self->myProblem->output(self->myProblem->options.outputfile);
					}
                }
                //if the current point is more feasible than the incumbent
                else if (fmax(f_current, decision_variable_feasibility_metric) < self->feasibility_metric_NLP_incumbent)
                {
                    self->J_NLP_incumbent = self->F.front();
                    self->feasibility_metric_NLP_incumbent = fmax(f_current, decision_variable_feasibility_metric);

                    self->X_NLP_incumbent_unscaled = self->X_unscaled;
                    self->X_NLP_incumbent_scaled = self->X_scaled;
                    self->F_NLP_incumbent = self->F;
					
					if (fmax(f_current, decision_variable_feasibility_metric) < self->myOptions.get_feasibility_tolerance() && self->first_feasibility == false)
					{
						self->first_feasibility = true;
						self->myProblem->Xopt = self->X_NLP_incumbent_unscaled;
						self->myProblem->F = self->F_NLP_incumbent;
	                    self->myProblem->evaluate(self->X_NLP_incumbent_unscaled, self->F_NLP_incumbent, self->G, false);

                        self->myProblem->what_the_heck_am_I_called(SolutionOutputType::SUCCESS);
                        self->myProblem->output(self->myProblem->options.outputfile);
					}
                }
            }//end NLP chaperone

            //Step 7: print movie frames
            if (*needG && self->myOptions.get_print_NLP_movie_frames())
            {
                self->myProblem->X = self->X_unscaled;
                self->myProblem->F = self->F;
                self->myProblem->output_problem_bounds_and_descriptions(self->myProblem->options.working_directory + "//" + "NLP_frame_" + std::to_string(self->movie_frame_count++) + ".csv");
            }

            //Step 8: check the timer and kill SNOPT if necessary
            time_t now = time(NULL);
            if (now - (time_t)*iu > self->myOptions.get_max_run_time_seconds())
            {
                if (!self->myOptions.get_quiet_NLP())
                    std::cout << "Exceeded NLP time limit of " << self->myOptions.get_max_run_time_seconds() << " seconds. Aborting NLP run." << std::endl;
                *Status = -2;
            }

            //Step 9: do debug things if appropriate
#ifdef DEBUG_CHECK_VARIABLES_AND_CONSTRAINTS_IN_USER_FUNCTION
            if (self->myProblem->check_inputs_for_invalid_entries(x, false))
            {
                std::cout << "Invalid entry in decision vector" << endl;
                *Status = -1;
            }
            else if (self->myProblem->check_outputs_for_invalid_entries(Problem->X.data(), F, false))
            {
                std::cout << "Invalid entry in objective/constraint vector" << endl;
                *Status = -1;
            }
            else if (self->myProblem->check_derivatives_for_invalid_entries(Problem->X.data(), G, false))
            {
                std::cout << "Invalid entry in Jacobian" << endl;
                *Status = -1;
            }
#endif

#ifdef SNOPT72
            return 0;
#else
            return;
#endif
        }
    }//end namespace Solvers
}//end namespace EMTG