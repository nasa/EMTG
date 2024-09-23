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

#pragma once

#include "doubleType.h"

#include <string>
#include <vector>

#include "missionoptions.h"
#include "EMTG_solver_utilities.h"

namespace EMTG
{
    class sparsey_thing
    {
    public:
        //constructor
        sparsey_thing() {};

        virtual void setup_calcbounds(
            std::vector<double>* Xupperbounds,
            std::vector<double>* Xlowerbounds,
            std::vector<double>* X_scale_factors,
            std::vector<double>* Fupperbounds,
            std::vector<double>* Flowerbounds,
			std::vector<double>* F_scale_factors,
            std::vector<std::string>* Xdescriptions,
            std::vector<std::string>* Fdescriptions,
            std::vector<size_t>* iGfun,
            std::vector<size_t>* jGvar,
            std::vector<std::string>* Gdescriptions,
            std::vector<size_t>* iAfun,
            std::vector<size_t>* jAvar,
            std::vector<std::string>* Adescriptions,
            std::vector<double>* A);
    protected:

        //utilities
        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const std::string& variable_name,
            size_t& sparsity_index_container);

        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const std::string& variable_name,
            std::vector<size_t>& sparsity_index_container);

        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xindex,
            size_t& sparsity_index_container);

        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xindex,
            std::vector<size_t>&  sparsity_index_container);

        void create_sparsity_entry_vector(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const int& number_of_entries,
            const std::string& variable_name,
            std::vector<size_t>& sparsity_index_container);
        
        //calcbounds fields
        std::string prefix;
        std::vector<double>* Xupperbounds;
        std::vector<double>* Xlowerbounds;
        std::vector<double>* X_scale_factors;
        std::vector<double>* Fupperbounds;
        std::vector<double>* Flowerbounds;
		std::vector<double>* F_scale_factors;
        std::vector<std::string>* Xdescriptions;
        std::vector<std::string>* Fdescriptions;
        std::vector<size_t>* iGfun;
        std::vector<size_t>* jGvar;
        std::vector<std::string>* Gdescriptions;
        std::vector<size_t>* iAfun;
        std::vector<size_t>* jAvar;
        std::vector<std::string>* Adescriptions;
        std::vector<double>* A;
    };//end class sparsey_thing
}