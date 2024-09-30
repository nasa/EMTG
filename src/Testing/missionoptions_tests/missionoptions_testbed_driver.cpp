// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2024 United States Government as represented by the
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

//test driver for new missionoptions
//Jacob Englander 1/9/2019

#include "missionoptions.h"
#include "journeyoptions.h"

#include <iostream>

#include <exception>

int main(int argc, char* argv[])
{
    try
    {
        //startup stuff
        std::cout << "program starting" << std::endl;

        //parse the options file
        std::string options_file_name;
        if (argc == 1)
            options_file_name = "default.emtgopt";
        else if (argc == 2)
            options_file_name.assign(argv[1]);

        std::cout << options_file_name << std::endl;

        EMTG::missionoptions options(options_file_name);

        options.write("tests/newdefault.emtgopt", !options.print_only_non_default_options);
    }
    catch (std::exception myException)
    {
        std::cout << "Oops!" << std::endl;
        std::cout << myException.what() << std::endl;
    }
}

