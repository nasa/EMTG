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

//EMTG solver utilities
//Jacob Englander 8-31-2016

#include "EMTG_solver_utilities.h"
#include <sstream>
#include <iostream>

namespace EMTG
{
    namespace solver_utilities
    {
        int detect_duplicate_Jacobian_entries(const std::vector<size_t>& iGfun,
            const std::vector<size_t>& jGvar,
            const size_t& Findex,
            const size_t& Xindex)
        {
            for (size_t Gindex = 0; Gindex < iGfun.size(); ++Gindex)
            {
                if (iGfun[Gindex] == Findex && jGvar[Gindex] == Xindex) 
                {
                    return Gindex;
                }
            }

            return -1;
        }//end detect_duplicate_Jacobian_entries()


        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const std::string& variable_name,
            const std::vector<std::string>& Fdescriptions,
            const std::vector<std::string>& Xdescriptions,
            std::vector<std::string>& Gdescriptions,
            std::vector<size_t>& iGfun,
            std::vector<size_t>& jGvar,
            size_t& sparsity_index_container)
        {
            size_t Xsize = Xdescriptions.size();
            // std::cout<<variable_name<<std::endl;
            if (ForwardPass)
            {
                for (size_t Xindex = Xstart; Xindex < Xsize; ++Xindex)
                {
                    if (Xdescriptions[Xindex].find(variable_name) < 1024)
                    {
                        std::stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << Fdescriptions[Findex] << " F[" << Findex << "] with respect to X[" << Xindex << "]: " << Xdescriptions[Xindex];
                        // std::cout<<EntryNameStream<<std::endl;
                        int tempIdx = detect_duplicate_Jacobian_entries(iGfun, jGvar, Findex, Xindex);
                        if (tempIdx < 0)
                        {
                            iGfun.push_back(Findex);
                            jGvar.push_back(Xindex);
                            Gdescriptions.push_back(EntryNameStream.str());
                            sparsity_index_container = Gdescriptions.size() - 1;
                        }
                        else
                            sparsity_index_container = tempIdx;
                        return Xindex;
                    }
                }
            }
            else
            {
                for (int Xindex = Xstart; Xindex >= 0; --Xindex)
                {
                    if (Xdescriptions[Xindex].find(variable_name) < 1024)
                    {
                        std::stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << Fdescriptions[Findex] << " F[" << Findex << "] with respect to X[" << Xindex << "]: " << Xdescriptions[Xindex];

                        int tempIdx = detect_duplicate_Jacobian_entries(iGfun, jGvar, Findex, Xindex);
                        if (tempIdx < 0)
                        {
                            iGfun.push_back(Findex);
                            jGvar.push_back(Xindex);
                            Gdescriptions.push_back(EntryNameStream.str());
                            sparsity_index_container = Gdescriptions.size() - 1;
                        }
                        else
                            sparsity_index_container = tempIdx;
                        return Xindex;
                    }
                }
            }
            throw std::invalid_argument("Decision variable '" + variable_name + "' does not exist. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            return 32767;
        }//end create_sparsity_entry

        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const std::string& variable_name,
            const std::vector<std::string>& Fdescriptions,
            const std::vector<std::string>& Xdescriptions,
            std::vector<std::string>& Gdescriptions,
            std::vector<size_t>& iGfun,
            std::vector<size_t>& jGvar,
            std::vector<size_t>& sparsity_index_container)
        {
            sparsity_index_container.push_back(SIZE_MAX);
            return create_sparsity_entry(Findex,
                                         Xstart,
                                         ForwardPass,
                                         variable_name,
                                         Fdescriptions,
                                         Xdescriptions,
                                         Gdescriptions,
                                         iGfun,
                                         jGvar,
                                         sparsity_index_container.back());

            if (sparsity_index_container.back() == SIZE_MAX)
                sparsity_index_container.pop_back();
        }



        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xindex,
            const std::vector<std::string>& Fdescriptions,
            const std::vector<std::string>& Xdescriptions,
            std::vector<std::string>& Gdescriptions,
            std::vector<size_t>& iGfun,
            std::vector<size_t>& jGvar,
            size_t& sparsity_index_container)
        {
            std::stringstream EntryNameStream;
            EntryNameStream << "Derivative of " << Fdescriptions[Findex] << " F[" << Findex << "] with respect to X[" << Xindex << "]: " << Xdescriptions[Xindex];

            int tempIdx = detect_duplicate_Jacobian_entries(iGfun, jGvar, Findex, Xindex);
            if (tempIdx < 0)
            {
                iGfun.push_back(Findex);
                jGvar.push_back(Xindex);
                Gdescriptions.push_back(EntryNameStream.str());
                sparsity_index_container = Gdescriptions.size() - 1;
            }
            else
                sparsity_index_container = tempIdx;
            return Xindex;
        }

        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xindex,
            const std::vector<std::string>& Fdescriptions,
            const std::vector<std::string>& Xdescriptions,
            std::vector<std::string>& Gdescriptions,
            std::vector<size_t>& iGfun,
            std::vector<size_t>& jGvar,
            std::vector<size_t>&  sparsity_index_container)
        {
            sparsity_index_container.push_back(SIZE_MAX);
            return create_sparsity_entry(Findex,
                Xindex,
                Fdescriptions,
                Xdescriptions,
                Gdescriptions,
                iGfun,
                jGvar,
                sparsity_index_container.back());

            if (sparsity_index_container.back() == SIZE_MAX)
                sparsity_index_container.pop_back();
        }

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
            std::vector<size_t>& sparsity_index_container)
        {
            size_t Xsize = Xdescriptions.size();
            int counter = 0;
            if (ForwardPass)
            {
                for (size_t Xindex = Xstart; Xindex < Xsize; ++Xindex)
                {
                    if (Xdescriptions[Xindex].find(variable_name) < 1024)
                    {
                        std::stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << Fdescriptions[Findex] << " F[" << Findex << "] with respect to X[" << Xindex << "]: " << Xdescriptions[Xindex];

                        int tempIdx = detect_duplicate_Jacobian_entries(iGfun, jGvar, Findex, Xindex);
                        if (tempIdx < 0)
                        {
                            iGfun.push_back(Findex);
                            jGvar.push_back(Xindex);
                            Gdescriptions.push_back(EntryNameStream.str());
                            sparsity_index_container.push_back(Gdescriptions.size() - 1);
                        }
                        else
                            sparsity_index_container.push_back(tempIdx);

                        ++counter;
                        if (counter == number_of_entries)
                            break;
                    }
                }
            }
            else
            {
                for (int Xindex = Xstart; Xindex >= 0; --Xindex)
                {
                    if (Xdescriptions[Xindex].find(variable_name) < 1024)
                    {
                        std::stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << Fdescriptions[Findex] << " F[" << Findex << "] with respect to X[" << Xindex << "]: " << Xdescriptions[Xindex];

                        int tempIdx = detect_duplicate_Jacobian_entries(iGfun, jGvar, Findex, Xindex);
                        if (tempIdx < 0)
                        {
                            iGfun.push_back(Findex);
                            jGvar.push_back(Xindex);
                            Gdescriptions.push_back(EntryNameStream.str());
                            sparsity_index_container.push_back(Gdescriptions.size() - 1);
                        }
                        else
                            sparsity_index_container.push_back(tempIdx);

                        ++counter;
                        if (counter == number_of_entries)
                            break;
                    }
                }
            }
        }//end create_sparsity_entry_vector
    }//end namespace solver_utilities
}//end namespace EMTG