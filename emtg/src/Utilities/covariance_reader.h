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

//class that reads covariance-inverse matrices
//assume row ordering

#pragma once

#include <string>
#include <vector>
#include <map>

#include "gsl/gsl_spline.h"
#include "gsl/gsl_errno.h" //I'm not really sure what this is but the GSL examples all have it, probably error handling

#include "EMTG_Matrix.h"

namespace EMTG
{
    class covariance_reader
    {
    public:
        //constructor
        covariance_reader();

        covariance_reader(const std::string& inputfilename, const std::string& leap_seconds_path);

        virtual ~covariance_reader();

        void initialize(const std::string& inputfilename, const std::string& leap_seconds_path);

        void fit_splines();

        inline void getPinv(const double& epoch, math::Matrix<double>& Pinv)
        {
            if (epoch > std::get<0>(this->data["epoch"]).back() || epoch < std::get<0>(this->data["epoch"]).front())
            {
                throw std::runtime_error("Covariance reader does not have data for epoch " + std::to_string(epoch / 86400.0) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else
            {
                for (size_t i = 0; i < 7; ++i)
                {
                    for (size_t j = 0; j < 7; ++j)
                    {
                        size_t headerIndex = i * 7 + j + 1;
                        Pinv(i, j) = gsl_spline_eval(std::get<1>(this->data[this->MatrixHeaders[headerIndex]]), epoch, this->SplineAccelerator);
                    }
                }
            }
        }

        inline void getPinvDerivative(const double& epoch, math::Matrix<double>& PinvDerivative)
        {
            if (epoch > std::get<0>(this->data["epoch"]).back() || epoch < std::get<0>(this->data["epoch"]).front())
            {
                throw std::runtime_error("Covariance reader does not have data for epoch " + std::to_string(epoch / 86400.0) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else
            {
                for (size_t i = 0; i < 7; ++i)
                {
                    for (size_t j = 0; j < 7; ++j)
                    {
                        size_t headerIndex = i * 7 + j + 1;
                        PinvDerivative(i, j) = gsl_spline_eval_deriv(std::get<1>(this->data[this->MatrixHeaders[headerIndex]]), epoch, this->SplineAccelerator);
                    }
                }
            }
        }

    protected:

        //fields
        size_t dataRows;
        std::vector<std::string> MatrixHeaders;
        std::map< std::string, std::tuple<std::vector<double>, gsl_spline* > > data;

        //GSL spline structs
        bool fitSplines;
        gsl_interp_accel *SplineAccelerator;
    };//end clss writey_thing
}