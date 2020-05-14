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

//header file for EMTG innerloop solution
//for use with Filament Walker, maybe also MBH

#include "EMTG_innerloop_solution.h"

namespace EMTG
{
    EMTG_innerloop_solution::EMTG_innerloop_solution(const problem* myProblem, 
                                                     const std::vector<double>& X,
                                                     const std::vector<double>& F,
                                                     const std::string& name) :
        X(X),
        F(F),
        name(name)
    {
        this->compute_constraint_violations(myProblem);
    }

    void EMTG_innerloop_solution::compute_constraint_violations(const problem* myProblem)
    {
        this->objective_function_value = this->F[0];
        this->equality_constraint_violation = 0.0;
        this->inequality_constraint_violation = 0.0;
        this->overall_constraint_violation = 0.0;

        for (size_t Findex = 1; Findex < myProblem->total_number_of_constraints; ++Findex)
        {
            double violation = this->F[Findex] < myProblem->Flowerbounds[Findex] ?
                myProblem->Flowerbounds[Findex] - this->F[Findex] : (this->F[Findex] > myProblem->Fupperbounds[Findex] ? 
                    this->F[Findex] - myProblem->Fupperbounds[Findex] : 0.0);
            double F2 = violation * violation;
            this->overall_constraint_violation += F2;

            if (myProblem->F_equality_or_inequality[Findex - 1])
                this->equality_constraint_violation += F2;
            else
                this->inequality_constraint_violation += F2;
        }

        this->equality_constraint_violation = sqrt(this->equality_constraint_violation);
        this->inequality_constraint_violation = sqrt(this->inequality_constraint_violation);
        this->overall_constraint_violation = sqrt(this->overall_constraint_violation);
    }
}