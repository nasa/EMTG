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

#include "EMTG_string_utilities.h"

namespace EMTG {
    namespace string_utilities {

        //adapted from from http://stackoverflow.com/questions/7132957/c-scientific-notation-format-number
        std::string convert_number_to_formatted_string(const double& number, int expSize)
        {
            std::ostringstream oss;
            oss.precision(14);
            std::string output;
            oss.precision(14);
            oss << std::scientific << (fabs(number) > 1.0e-99 ? number : 0.0);
            std::string numstring = oss.str();
            size_t ePos = numstring.find("e");
            numstring.replace(ePos, 1, "E");
            size_t dPos = numstring.find(".");
            if (ePos == 0)
            {
                //no exponent
                return numstring;
            }
            else if (dPos == 0)
            {
                //not decimal
                return numstring;
            }
            else
            {
                if (number < 0.0)
                    output = "";
                else
                    output = " ";
                output += numstring.substr(0, ePos);

                while (output.size() < 17)
                    output += "0";

                output += numstring.substr(ePos, 2);
                if (numstring.size() - expSize > ePos + 1)
                    output += numstring.substr(numstring.size() - expSize, numstring.size());
                else
                {
                    //expSize too big (or bug -> e used but no exponent?)
                }
                return output;
            }
        }//end convert_number_to_formatted_string()

        bool string_contains_substring(const std::string& myString, const std::string& mySubString)
        {
            return (myString.find(mySubString) < 1024);
        }//end string_contains_substring()
    }//end namespace string_utilities
}//end namespace EMTG
