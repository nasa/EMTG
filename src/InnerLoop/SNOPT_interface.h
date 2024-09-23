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

//header file for SNOPT interface
#pragma once

#include "NLP_interface.h"
#include "snoptProblemExtension.h"

namespace EMTG
{
    namespace Solvers
    {

        class SNOPT_interface : public NLP_interface
        {
        public:
            //constructor
            SNOPT_interface() : NLP_interface::NLP_interface() {};
            SNOPT_interface(problem* myProblem,
                const NLPoptions& myOptions);

            //run NLP
            virtual void run_NLP(const bool& X0_is_scaled = true);
#ifdef SNOPT72
            static int
#else
            static void
#endif
                SNOPT_user_function(SNOPT_INT_TYPE    *Status, SNOPT_INT_TYPE *n, SNOPT_DOUBLE_TYPE x[],
                                    SNOPT_INT_TYPE    *needF, SNOPT_INT_TYPE *neF, SNOPT_DOUBLE_TYPE F[],
                                    SNOPT_INT_TYPE    *needG, SNOPT_INT_TYPE *neG, SNOPT_DOUBLE_TYPE G[],
                                    char       *cu, SNOPT_INT_TYPE *lencu,
                                    SNOPT_INT_TYPE    iu[], SNOPT_INT_TYPE *leniu,
                                    SNOPT_DOUBLE_TYPE ru[], SNOPT_INT_TYPE *lenru);

            //get
            SNOPT_INT_TYPE getInform() const { return this->inform; }

        private:
            //SNOPT user function
            //internal feasibility check
            //fields
            snoptProblemExtension mySNOPT;

            SNOPT_INT_TYPE inform;
            SNOPT_INT_TYPE neF;
            SNOPT_INT_TYPE lenA;

            std::vector<SNOPT_INT_TYPE> iAfun;
            std::vector<SNOPT_INT_TYPE> jAvar;
            std::vector<SNOPT_DOUBLE_TYPE> A;

            SNOPT_INT_TYPE lenG;
            std::vector<SNOPT_INT_TYPE> iGfun;
            std::vector<SNOPT_INT_TYPE> jGvar;
            std::vector<SNOPT_DOUBLE_TYPE> xmul;
            std::vector<SNOPT_INT_TYPE>    xstate;

            std::vector<SNOPT_DOUBLE_TYPE> Fmul;
            std::vector<SNOPT_INT_TYPE>    Fstate;
			
#ifdef SNOPT72
            SNOPT_INT_TYPE nxnames;
            SNOPT_INT_TYPE nFnames;
            char *xnames;
            char *Fnames;
#endif

            SNOPT_INT_TYPE    ObjRow;
            SNOPT_DOUBLE_TYPE ObjAdd;
        };//end class SNOPT_interface

        
    }//end namespace Solvers
}//end namespace EMTG