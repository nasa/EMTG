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

#ifndef _EMTG_WORHP_INTERFACE
#define _EMTG_WORHP_INTERFACE

#ifdef _use_WORHP
#include "worhp\worhp.h"
#include "missionoptions.h"
#include "problem.h"

namespace EMTG { namespace Solvers {

    class EMTG_sparse_matrix_entry
    {
    public:
        //constructor
        EMTG_sparse_matrix_entry(const int& R, const int& C, const int& O);

        //destructor
        ~EMTG_sparse_matrix_entry();

        //members
        int RowIndex;
        int ColIndex;
        int OriginalIndex;
    };//close class EMTG_sparse_matrix_entry

    class EMTG_WORHP_interface
    {
    public:
        //constructor
        EMTG_WORHP_interface();
        EMTG_WORHP_interface(problem* Problem_Input);

        //destructor
        ~EMTG_WORHP_interface();

        //methods
        int SetInitialGuess(double* X);
        int Solve();

        //fields
        OptVar    WORHP_OptVar;
        Workspace WORHP_Workspace;
        Workspace* WORHP_Workspace_pointer;
        Params    WORHP_Params;
        Control   WORHP_Control;
        problem* Problem;
        int status;
        int nDF; //number of derivative entries for the gradient
        int nDG; //number of derivative entries for the Jacobian
        int nHM; //number of derivative entries for the Hessian
        vector<EMTG_sparse_matrix_entry> DGpattern;
        vector<int> iDFfun;
        vector<int> jDFvar;
        vector<int> DFindices;
        vector<int> DFfidif_indices;
        vector<int> iDGfun;
        vector<int> jDGvar;
        vector<int> DGindices;
        vector<int> DGfidif_indices;
        bool linear_objective_function;
    };//close class EMTG_WORHP_interface

}}//close namespace


#endif
#endif
