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

//* ------- SNOPT PROBLEM EXTENSION
//
//        This header is built to extend a standard snoptProblem (as per snopt v7.x)
//
//        WARNING: This header OVERRIDES some snopt functions.  This means polymorphism for those functions will break
//            To use this properly DO NOT assume polymorphism is functional on snjac or solve().
//            i.e. DO NOT assume a pointer to an original snoptProblem which has been assigned a snoptProblemExtension will
//            call the new child solve().  IT WILL CALL THE PARENT.
//        TO AVOID THIS PROBLEM: only use pointers to snoptProblemExtension or extension objects directly.
//
//
// ---------------------------------------------*/

#pragma once

#include <cstring>
#include <iostream>
#include <cassert>

#ifdef SNOPT72
    #include "snoptProblem.hh"
    #include "snopt.hh"
    #include "snfilewrapper.hh"
    #define SNOPT_INT_TYPE integer
    #define SNOPT_DOUBLE_TYPE doublereal
    #define SNOPT_BASE_TYPE public::snoptProblem
#else
    #include "snoptProblem.hpp"
    #include "snopt.h"
    #define SNOPT_INT_TYPE int
    #define SNOPT_DOUBLE_TYPE double
    #define SNOPT_BASE_TYPE public::snoptProblemA
#endif

namespace EMTG
{
    namespace Solvers
    {

#ifndef SNPT_EXT_H
#define SNPT_EXT_H
        class snoptProblemExtension : SNOPT_BASE_TYPE {

        public:
#ifdef SNOPT72
            snoptProblemExtension(bool const & suppress_input = false) : snoptProblem() { //extended contructor for enabling Alex's snopt silencer.  Parent default constructor is fine in most cases
#else
            snoptProblemExtension(bool const & suppress_input = false) : snoptProblemA() { //extended contructor for enabling Alex's snopt silencer.  Parent default constructor is fine in most cases
#endif
#ifdef QUIET_SNOPT
                if (suppress_input) {
                    initCalled = 1;
                    this->setIntParameter((char*) "Print No", 1);
                };
#endif        
#ifdef SNOPT72    
                leniu = 0;
                lenru = 0;
#endif
            };
#ifndef SNOPT76
            SNOPT_INT_TYPE getNeA() { return neA; };
            SNOPT_INT_TYPE getNeG() { return neG; };
#endif 
			
#ifdef SNOPT72    
            void setSummaryFile(char asummaryname[]) {
                assert(initCalled = 1);
                if (iSumm != 0) {
                    snclose_(&iSumm);
                }
                iSumm = 6;
#ifdef _STONEAGECplusplus
                strcpy(summaryname, asummaryname);
#else
                strcpy_s(summaryname, asummaryname);
#endif
                prnt_len = strlen(summaryname);
                snopenappend_(&iSumm, summaryname, &inform, prnt_len);
                this->setIntParameter((char*)"Summary file", iSumm);
            }
#endif


#ifdef SNOPT72    
            void setUserspace(integer *iu_in, integer leniu_in, double *ru_in, integer lenru_in) {

                leniu = leniu_in;
                lenru = lenru_in;

                iu = iu_in;
                ru = ru_in;

            }
#endif

#ifdef SNOPT72
            //this is a very similar, but overridden version of computeJac from the parent.  This is NOT polymorphism safe.
            //this method actually passes in the userspace correctly.
            integer computeJac() {
                //Ensures all user data has been initialized.
                userDataSet();
                this->snmema(mincw, miniw, minrw);
                if (mincw > lencw | miniw > leniw | minrw > lenrw) {

                    this->realloc(mincw, miniw, minrw);

                    this->setIntParameter((char*)"Total real workspace   ", lenrw);
                    this->setIntParameter((char*)"Total integer workspace", leniw);
                }
                snjac_(&inform, &neF, &n, usrfun, iAfun, jAvar, &lenA, &neA, A, iGfun, jGvar, &lenG, &neG, x, xlow, xupp, &mincw, &miniw, &minrw, cu, &lencu, iu, &leniu, ru, &lenru, cw, &lencw, iw, &leniw, rw, &lenrw, 8 * 500, 8 * 500);

                fortranStyleAG = 1;
                return inform;

            }
#else
            int computeJac()
            {
#ifdef SNOPT75
                snoptProblemA::computeJac(this->neA, this->neG);
                return inform;
#elif SNOPT76
				return 0;
#else
                snoptProblemA::computeJac();
                return inform;
#endif //SNOPT75
            }
#endif


#ifdef SNOPT72
            //this is a very similar, but overridden version of solve from the parent.  This is NOT polymorphism safe.
            //this method actually passes in the userspace correctly.
            integer solve(integer starttype) {
                assert(initCalled == 1);

                userDataSet();

                if (neA == -1 | neG == -1) {
                    std::cout << "Warning: neA and neG must be set before calling" << "snoptProblem::solve()\n";
                    exit(1);
                }
                integer npname = strlen(Prob);
                integer nS, nInf;
                doublereal sInf;
                this->increment(); //Convert array entries to Fortran style
                this->setMemory();
                snopta_(&starttype, &neF, &n, &nxnames, &nFnames, &ObjAdd, &ObjRow, Prob, usrfun, iAfun, jAvar, &lenA, &neA, A, iGfun, jGvar, &lenG, &neG, xlow, xupp, xnames, Flow, Fupp, Fnames, x, xstate, xmul, F, Fstate, Fmul, &inform, &mincw, &miniw, &minrw, &nS, &nInf, &sInf, cw, &lencw, iu, &leniu, ru, &lenru, cw, &lencw, iw, &leniw, rw, &lenrw, npname, 8 * (nxnames), 8 * (nFnames), 8 * 500, 8 * 500);
                this->decrement();  //Convert array entries to C style

                return inform;
            }
#endif
        protected:
#ifdef SNOPT72
            integer    lenru, leniu, lencu;
            doublereal *ru;
            integer    *iu;
            char       *cu;
#endif    

            char summaryname[200];



        };//end class snoptProblemExtension
    }//end namespace Solvers
}//end namespace EMTG


#endif