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

//Class interface to WORHP, to tidy-up global solvers which need to use it
//note this won't compile without WORHP so the whole thing lives in an #ifdef guard
//Jacob Englander 11/8/2013

#ifdef _use_WORHP
#include <cmath>
#include <float.h>

#include "missionoptions.h"
#include "problem.h"
#include "EMTG_WORHP_interface.h"



namespace EMTG { namespace Solvers {
    //default constructor is never used
    EMTG_WORHP_interface::EMTG_WORHP_interface()
    {
    }

    //constructor when a problem object is supplied
    EMTG_WORHP_interface::EMTG_WORHP_interface(problem* Problem_input)
    {
        Problem = Problem_input;

        //this code adapted from the WORHP C++ example
        this->WORHP_Workspace_pointer = &this->WORHP_Workspace;
        this->WORHP_OptVar.initialised = false;
        this->WORHP_Workspace.initialised = false;
        this->WORHP_Params.initialised = false;
        this->WORHP_Control.initialised = false;

        /* Check Version of library and header files */
        CHECK_WORHP_VERSION

        /* Uncomment this to get more info on data structures */
        // WorhpDiag(&o, &w, &p, &c);

        /* Read Worhp XML parameter file */
        ReadParams(&this->status, (char*)"param.xml", &this->WORHP_Params);
        this->status = OK; // ignore errors: WorhpReadParams falls back to standard values.
        this->WORHP_Params.MatrixCC = false;

        /* Specify number of variables and constraints. */
        this->WORHP_OptVar.n = this->Problem->total_number_of_NLP_parameters;
        this->WORHP_OptVar.m = this->Problem->total_number_of_constraints - 1; //first entry is the objective function

        /* Specify nonzeros of derivative matrices. */
        //how many derivative entries are there for the objective function?
        this->nDF = 0;
        this->linear_objective_function = false;
        for (int k = 0; k < this->Problem->iAfun.size(); ++k)
        {
            if (this->Problem->iAfun[k] == 0)
            {
                ++nDF;
                this->iDFfun.push_back(0);
                this->jDFvar.push_back(Problem->jAvar[k] + 1);
                this->linear_objective_function = true;
                this->DFindices.push_back(k);
            }
            else
            {
                this->iDGfun.push_back(Problem->iAfun[k]);
                this->jDGvar.push_back(Problem->jAvar[k] + 1);
                this->DGindices.push_back(k);
            }
        }
        for (int k = 0; k < Problem->iGfun.size(); ++k)
        {
            if (Problem->iGfun[k] == 0)
            {
                ++nDF;
                this->iDFfun.push_back(0);
                this->jDFvar.push_back(Problem->jGvar[k] + 1);
                this->DFindices.push_back(k);
            }
            else
            {
                this->iDGfun.push_back(Problem->iGfun[k]);
                this->jDGvar.push_back(Problem->jGvar[k] + 1);
                this->DGindices.push_back(k);
            }
        }
        this->nDG = Problem->iAfun.size() + Problem->iGfun.size() - nDF;

        for (int k = 0; k < this->nDG; ++k)
        {
            EMTG_sparse_matrix_entry NewEntry(this->iDGfun[k], this->jDGvar[k], k);
            this->DGpattern.push_back(NewEntry);
        }
        std::sort(this->DGpattern.begin(), this->DGpattern.end(), [](const EMTG_sparse_matrix_entry &a, const EMTG_sparse_matrix_entry &b)
        {
            if (a.ColIndex == b.ColIndex) return a.RowIndex < b.RowIndex;
            else return a.ColIndex < b.ColIndex;
        });

        //which entries of the gradient and Jacobian are not specified analytically?
        //discover this by running the evaluate function and seeing which entries are unchanged
        this->Problem->G = vector<double> (this->Problem->G.size(), -1.01e+100);
        vector<double> centerX(this->Problem->total_number_of_NLP_parameters, 0.5);
        Problem->unscale(centerX.data());
        Problem->evaluate(this->Problem->X.data(), this->Problem->F.data(), this->Problem->G.data(), 1, this->Problem->iGfun, this->Problem->jGvar);

        for (int k = 0; k < this->nDF; ++k)
        {
            if (this->Problem->G[this->DFindices[k]] < -1.0e+100)
                this->DFfidif_indices.push_back(this->DFindices[k]);
        }
        for (int k = 0; k < this->nDG; ++k)
        {
            if (this->Problem->G[this->DGindices[k]] < -1.0e+100)
                this->DGfidif_indices.push_back(this->DGindices[k]);
        }

        this->WORHP_Workspace.DF.nnz = this->nDF;        /* dense, structure init by solver */
        this->WORHP_Workspace.DG.nnz = this->nDG;    /* dense, structure init by solver */
        this->WORHP_Workspace.HM.nnz = this->WORHP_OptVar.n;        /* diagonal, structure init by user */

        /* Data structure initialisation. */
        WorhpInit(&this->WORHP_OptVar, &this->WORHP_Workspace, &this->WORHP_Params, &this->WORHP_Control);
        if (this->WORHP_Control.status != FirstCall)
            cout << "Main: Initialisation failed." << endl;
        //return EXIT_FAILURE;

        /* Define bounds for X */
        for (int k = 0; k < Problem->total_number_of_NLP_parameters; ++k)
        {
            this->WORHP_OptVar.XL[k] = 0.0;
            this->WORHP_OptVar.XU[k] = 1.0;

        }

        /* Define bounds for G */
        for (int k = 0; k < Problem->total_number_of_constraints - 1; ++k)
        {
            this->WORHP_OptVar.GL[k] = Problem->Flowerbounds[k+1];
            this->WORHP_OptVar.GU[k] = Problem->Fupperbounds.[k+1];
        }

        /*
        * Specify matrix structures in CS format, using Fortran indexing, 
        * i.e. 1...N instead of 0...N-1, to describe the matrix structure.
        */

        //define the gradient DF structure in CS format
        if (this->WORHP_Workspace.DF.NeedStructure)
        {
            for (int i = 0; i < this->nDF; ++i)
            {
                this->WORHP_Workspace.DF.row[i] = this->jDFvar[i];
            }
        }

        //define the gradient DG structure in CS format
        if (this->WORHP_Workspace.DG.NeedStructure)
        {
            for (int i = 0; i < this->nDG; ++i)
            {
                this->WORHP_Workspace.DG.row[i] = this->DGpattern[i].RowIndex;
                this->WORHP_Workspace.DG.col[i] = this->DGpattern[i].ColIndex;
            }
        }

        /*
        * Define HM as a diagonal LT-CS-matrix, but only if needed by WORHP
        */
        if (this->WORHP_Workspace.HM.NeedStructure) 
        {
            for(int i = 0; i < this->WORHP_Workspace.HM.nnz; ++i) 
            {
                this->WORHP_Workspace.HM.row[i] = i + 1;
                this->WORHP_Workspace.HM.col[i] = i + 1;
            }
        }
    }

    //destructor
    EMTG_WORHP_interface::~EMTG_WORHP_interface()
    {
        WorhpFree(&this->WORHP_OptVar, &this->WORHP_Workspace, &this->WORHP_Params, &this->WORHP_Control);
    }

    //function to intialize the solver
    int EMTG_WORHP_interface::SetInitialGuess(double* X)
    {
        for (int i = 0; i < Problem->total_number_of_NLP_parameters; ++i)
        {
            this->WORHP_OptVar.X[i] = X[i];
            this->WORHP_OptVar.Lambda[i] = 0.0;
        }

        for (int j = 0; j < Problem->total_number_of_constraints; ++j)
            this->WORHP_OptVar.Mu[j] = 0.0;

        return 0;
    }

    //function to run the solver
    int EMTG_WORHP_interface::Solve()
    {
        int count = 0;
          /* Reverse Communication loop */
          while(this->WORHP_Control.status < TerminateSuccess &&  this->WORHP_Control.status > TerminateError) {

            if (GetUserAction(&this->WORHP_Control, callWorhp)) {
              Worhp(&this->WORHP_OptVar, &this->WORHP_Workspace, &this->WORHP_Params, &this->WORHP_Control);
              /* DO NOT call DoneUserAction here! */
            }

            if (GetUserAction(&this->WORHP_Control, iterOutput)) {
              IterationOutput(&this->WORHP_OptVar, &this->WORHP_Workspace, &this->WORHP_Params, &this->WORHP_Control);
              DoneUserAction(&this->WORHP_Control, iterOutput);
            }

            //new call structure
            //if WORHP wants DF, or DG
            //run the objective function with derivatives enabled
            //depending on what information is needed, copy information to appropriate vectors
            //fill in linear derivatives
            //if WORHP does not want DF or DG then call objective function without derivatives enabled
            if (GetUserAction(&this->WORHP_Control, evalDF) || GetUserAction(&this->WORHP_Control, evalDG))
            {
                Problem->unscale(this->WORHP_OptVar.X);
                Problem->evaluate(Problem->X.data(), Problem->F.data(), Problem->G.data(), 1, Problem->iGfun, Problem->jGvar);

                if (GetUserAction(&this->WORHP_Control, evalF)) 
                {
                    this->WORHP_OptVar.F = Problem->F[0];
                    DoneUserAction(&this->WORHP_Control, evalF);
                }

                if (GetUserAction(&this->WORHP_Control, evalG))
                {
                    for (int i = 0; i < this->WORHP_OptVar.m; ++i)
                        this->WORHP_OptVar.G[i] = Problem->F[i+1];
                    DoneUserAction(&this->WORHP_Control, evalG);
                }

                if (GetUserAction(&this->WORHP_Control, evalDF)) 
                {
                    if (this->linear_objective_function)
                    {
                        //process any linear derivatives of the objective function
                        for (int i = 0; i < this->nDF; ++i)
                            this->WORHP_Workspace.DF.val[i] = this->Problem->A[i] * this->WORHP_OptVar.X[this->WORHP_Workspace.DF.row[i] - 1];
                    }
                    else
                    {
                        //process any nonlinear derivatives of the objective function

                        //compute via finite differencing any derivatives which were not computed analytically

                        //copy the nonlinear gradient entries to WORHP
                    }
                    DoneUserAction(&this->WORHP_Control, evalDF);
                }

                if (GetUserAction(&this->WORHP_Control, evalDG)) 
                {
                    //compute via finite differencing any derivatives which were not computed analytically

                    //copy the nonlinear gradient entries to WORHP

                    DoneUserAction(&this->WORHP_Control, evalDG);
                }
            }
            else if (GetUserAction(&this->WORHP_Control, evalF) || GetUserAction(&this->WORHP_Control, evalG))
            {
                Problem->unscale(this->WORHP_OptVar.X);
                Problem->evaluate(Problem->X.data(), Problem->F.data(), Problem->G.data(), 0, Problem->iGfun, Problem->jGvar);

                if (GetUserAction(&this->WORHP_Control, evalF)) 
                {
                    this->WORHP_OptVar.F = Problem->F[0];
                    //std::cout << "F" << std::endl;
                    DoneUserAction(&this->WORHP_Control, evalF);
                }

                if (GetUserAction(&this->WORHP_Control, evalG))
                {
                    for (int i = 0; i < this->WORHP_OptVar.m; ++i)
                        this->WORHP_OptVar.G[i] = Problem->F[i+1];
                    //std::cout << "G" << std::endl;
                    DoneUserAction(&this->WORHP_Control, evalG);
                }
            }

            if (GetUserAction(&this->WORHP_Control, evalHM)) {
              /* Only scale the F part of HM */
              DoneUserAction(&this->WORHP_Control, evalHM);
            }
    
            if (GetUserAction(&this->WORHP_Control, fidif)) {
                //std::cout << "Finite difference" << std::endl;
                //std::cout << count++ << std::endl;
              WorhpFidif(&this->WORHP_OptVar, &this->WORHP_Workspace, &this->WORHP_Params, &this->WORHP_Control);
            }
          }
  
          StatusMsg(&this->WORHP_OptVar, &this->WORHP_Workspace, &this->WORHP_Params, &this->WORHP_Control);

        return 0;
    }


    //methods for the EMTG_sparse_matrix_entry class
    EMTG_sparse_matrix_entry::EMTG_sparse_matrix_entry(const int& R, const int& C, const int& O)
    {
        this->RowIndex = R;
        this->ColIndex = C;
        this->OriginalIndex = O;
    }

    EMTG_sparse_matrix_entry::~EMTG_sparse_matrix_entry(){}

}} //close namespace


#endif //_use_WORHP