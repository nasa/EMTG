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

#pragma once
#include "problem.h"

#include <vector>
#include <string>

namespace EMTG
{
    class EMTG_innerloop_solution
    {
    public:
        //constructors
        EMTG_innerloop_solution(const problem* myProblem,
                                const std::vector<double>& X = std::vector<double>(1, 0.0),
                                const std::vector<double>& F = std::vector<double>(1, 0.0),
                                const std::string& name = std::string("default_name"));

        //get/set
        inline void setX(const std::vector<double>& X) { this->X = X; }
        inline void setF(const std::vector<double>& F) { this->F = F; }
        inline void setname(const std::vector<double>& X) { this->X = X; }

        inline std::vector<double> getX() const { return this->X; }
        inline double getX(const size_t& Xindex) const { return this->X[Xindex]; }
        inline std::vector<double> getF() const { return this->F; }
        inline double getF(const size_t& Findex) const { return this->F[Findex]; }
        inline double getobjective_function_value() const { return this->objective_function_value; }
        inline double getequality_constraint_violation() const { return this->equality_constraint_violation; }
        inline double getinequality_constraint_violation() const { return this->inequality_constraint_violation; }
        inline double getoverall_constraint_violation() const { return this->overall_constraint_violation; }
        inline std::string getname() const { return this->name; }

        //other
        void compute_constraint_violations(const problem* myProblem);

    private:
        //fields
        std::vector<double> X;
        std::vector<double> F;
        time_t timestamp;
        size_t iteration;
        double objective_function_value;
        double equality_constraint_violation;
        double inequality_constraint_violation;
        double overall_constraint_violation;
        std::string name;

    };
}