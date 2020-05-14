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

#include "covariance_reader.h"
#include "file_utilities.h"

#include "boost/algorithm/string.hpp"

#include "SpiceUsr.h"

#include <fstream>

namespace EMTG
{
    covariance_reader::covariance_reader() :
        fitSplines(false)
    {
    }

    covariance_reader::covariance_reader(const std::string& inputfilename, const std::string& leap_seconds_path)
        : covariance_reader()
    {
        this->initialize(inputfilename, leap_seconds_path);
    }

    covariance_reader::~covariance_reader()
    {
        if (this->fitSplines)
        {
            this->fitSplines = false;

            for (size_t headerIndex = 1; headerIndex < 50; ++headerIndex)
            {
                gsl_spline_free(std::get<1>(this->data[this->MatrixHeaders[headerIndex]]));
            }

            gsl_interp_accel_free(this->SplineAccelerator);
        }
    }

    void covariance_reader::initialize(const std::string& inputfilename, const std::string& leap_seconds_path)
    {
        //Step 0: load a leap-seconds kernel. We need it to do time conversions
        furnsh_c(leap_seconds_path.c_str());

        //Step 1: clear the data object
        this->data.clear();
        
        //Step 2: open the input file
        std::ifstream inputfile(inputfilename);
        std::string line;

        //Step 3: ingest the column headers
        //Step 3.1: read the line
        file_utilities::safeGetline(inputfile, line);

        //Step 3.2: remove the comment character
        line.erase(line.begin());
        
        //Step 3.3: split the line
        boost::split(this->MatrixHeaders,
            line,
            boost::is_any_of(","),
            boost::token_compress_on);

        for (std::string& column : this->MatrixHeaders)
            boost::trim(column);

        //Step 3.4: create the (empty) columns
        for (std::string column : this->MatrixHeaders)
            std::get<0>(this->data[column]) = std::vector<double>();

        //Step 4: ingest the data

        this->dataRows = 0;
        while (file_utilities::safeGetline(inputfile, line))
        {
            if (line.size() > 0)
            {
                if (!(line.front() == *"#"))
                {
                    ++this->dataRows;

                    //split the line
                    std::vector<std::string> linecell;
                    boost::split(linecell,
                        line,
                        boost::is_any_of(","),
                        boost::token_compress_on);

                    //loop over columns
                    for (size_t columnIndex = 0; columnIndex < this->MatrixHeaders.size(); ++columnIndex)
                    {
                        std::string column = this->MatrixHeaders[columnIndex];

                        if (column == "epoch")
                        {
                            //transform to epoch
                            double SecondsPastJ2000;
                            str2et_c(linecell[0].c_str(), &SecondsPastJ2000);
                            std::get<0>(this->data[column]).push_back(SecondsPastJ2000 + 51544.5 * 86400.0);
                        }
                        else
                        {
                            //transform to double
                            std::get<0>(this->data[column]).push_back(std::stod(linecell[columnIndex]));
                        }
                    }//end loop over columns
                }
            }
        }//end loop over lines

        //Step 5: close the input file
        inputfile.close();

        //Step 6: unload the SPICE file
        unload_c(leap_seconds_path.c_str());

        //Step 7: fit splines!
        this->fit_splines();
    }//end initialize()

    void covariance_reader::fit_splines()
    {
        this->fitSplines = true;

        //allocate accelerator
        this->SplineAccelerator = gsl_interp_accel_alloc();

        for (size_t columnIndex = 1; columnIndex < 50; ++columnIndex)
        {
            std::string header = this->MatrixHeaders[columnIndex];
            //allocate
            std::get<1>(this->data[header]) = gsl_spline_alloc(gsl_interp_cspline, this->dataRows);

            //initialize
            gsl_spline_init(std::get<1>(this->data[this->MatrixHeaders[columnIndex]]), std::get<0>(this->data["epoch"]).data(), std::get<0>(this->data[header]).data(), this->dataRows);
        }
    }//end fit_splines()

}//end namespace EMTG