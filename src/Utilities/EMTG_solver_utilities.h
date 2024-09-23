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

//EMTG solver utilities header
//Jacob Englander 8-31-2016

#pragma once

#include <vector>
#include <string>

namespace EMTG
{
    namespace solver_utilities
    {
        int detect_duplicate_Jacobian_entries(const std::vector<size_t>& iGfun,
            const std::vector<size_t>& jGvar,
            const size_t& Findex,
            const size_t& Xindex);

        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const std::string& variable_name,
            const std::vector<std::string>& Fdescriptions,
            const std::vector<std::string>& Xdescriptions,
            std::vector<std::string>& Gdescriptions,
            std::vector<size_t>& iGfun,
            std::vector<size_t>& jGvar,
            size_t& sparsity_index_container);

        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const std::string& variable_name,
            const std::vector<std::string>& Fdescriptions,
            const std::vector<std::string>& Xdescriptions,
            std::vector<std::string>& Gdescriptions,
            std::vector<size_t>& iGfun,
            std::vector<size_t>& jGvar,
            std::vector<size_t>& sparsity_index_container);

        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xindex,
            const std::vector<std::string>& Fdescriptions,
            const std::vector<std::string>& Xdescriptions,
            std::vector<std::string>& Gdescriptions,
            std::vector<size_t>& iGfun,
            std::vector<size_t>& jGvar,
            size_t& sparsity_index_container);

        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xindex,
            const std::vector<std::string>& Fdescriptions,
            const std::vector<std::string>& Xdescriptions,
            std::vector<std::string>& Gdescriptions,
            std::vector<size_t>& iGfun,
            std::vector<size_t>& jGvar,
            std::vector<size_t>& sparsity_index_container);

        void create_sparsity_entry_vector(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const int& number_of_entries,
            const std::string& variable_name,
            const std::vector<std::string>& Fdescriptions,
            const std::vector<std::string>& Xdescriptions,
            std::vector<std::string>& Gdescriptions,
            std::vector<size_t>& iGfun,
            std::vector<size_t>& jGvar,
            std::vector<size_t>& sparsity_index_container);
    }//end namespace solver_utilities
}//end namespace EMTG