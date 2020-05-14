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

#include <fstream>
#include <iostream>

#include <boost/filesystem.hpp>

#include "file_utilities.h"

namespace fs = ::boost::filesystem;
namespace EMTG 
{
    namespace file_utilities
    {
        // return the filenames of all files that have the specified extension
        // in the specified directory and all subdirectories

        void get_all_files_with_extension(const fs::path& root, const std::string& ext, std::vector<fs::path>& ret)
        {  
          if (!fs::exists(root)) return;

          if (fs::is_directory(root))
          {
            fs::recursive_directory_iterator it(root);
            fs::recursive_directory_iterator endit;
            while(it != endit)
            {
                fs::path file(*it);
        
                if (fs::is_regular_file(file) && file.extension() == ext)
                ret.push_back(file.filename());
              ++it;
            }
          }
        }

        std::istream& safeGetline(std::istream& is, std::string& t)
        {
            t.clear();

            // The characters in the stream are read one-by-one using a std::streambuf.
            // That is faster than reading them one-by-one using the std::istream.
            // Code that uses streambuf this way must be guarded by a sentry object.
            // The sentry object performs various tasks,
            // such as thread synchronization and updating the stream state.

            std::istream::sentry se(is, true);
            std::streambuf* sb = is.rdbuf();

            for (;;) {
                int c = sb->sbumpc();
                switch (c) {
                case '\n':
                    return is;
                case '\r':
                    if (sb->sgetc() == '\n')
                        sb->sbumpc();
                    return is;
                case std::streambuf::traits_type::eof():
                    // Also handle the case when the last line has no line ending
                    if (t.empty())
                        is.setstate(std::ios::eofbit);
                    return is;
                default:
                    t += (char)c;
                }
            }
        }
    }//close namespace file_utilities
}//close namespace EMTG