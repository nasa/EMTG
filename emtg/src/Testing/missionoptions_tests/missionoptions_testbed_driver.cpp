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

