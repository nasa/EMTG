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

//class that reads EMTG-MIRAGE ".ephemeris" files
//flexible, adjusts to what headers are or are not present

#pragma once

#include <string>
#include <vector>
#include <map>
#include <stdexcept>

#include "gsl/gsl_spline.h"
#include "gsl/gsl_errno.h" //I'm not really sure what this is but the GSL examples all have it, probably error handling

namespace EMTG
{
    class ephemeris_reader
    {
    public:
        //constructor
        ephemeris_reader();

        ephemeris_reader(const std::string& inputfilename, const std::string& leap_seconds_path);

        virtual ~ephemeris_reader();

        void initialize(const std::string& inputfilename, const std::string& leap_seconds_path);

        void fit_splines();

        inline void getPosition(const double& epoch, double* PositionArray)
        {
            if (epoch > this->data["epoch"].back() || epoch < this->data["epoch"].front())
            {
                throw std::runtime_error("Ephemeris reader does not have data for epoch " + std::to_string(epoch / 86400.0) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else
            {
                PositionArray[0] = gsl_spline_eval(this->Spline_x, epoch, this->SplineAccelerator);
                PositionArray[1] = gsl_spline_eval(this->Spline_y, epoch, this->SplineAccelerator);
                PositionArray[2] = gsl_spline_eval(this->Spline_z, epoch, this->SplineAccelerator);
            }
        }

        inline void getVelocity(const double& epoch, double* VelocityArray)
        {
            if (epoch > this->data["epoch"].back() || epoch < this->data["epoch"].front())
            {
                throw std::runtime_error("Ephemeris reader does not have data for epoch " + std::to_string(epoch / 86400.0) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else
            {
                VelocityArray[0] = gsl_spline_eval(this->Spline_vx, epoch, this->SplineAccelerator);
                VelocityArray[1] = gsl_spline_eval(this->Spline_vy, epoch, this->SplineAccelerator);
                VelocityArray[2] = gsl_spline_eval(this->Spline_vz, epoch, this->SplineAccelerator);
            }
        }

        inline void getMass(const double& epoch, double* MassPointer)
        {
            if (epoch > this->data["epoch"].back() || epoch < this->data["epoch"].front())
            {
                throw std::runtime_error("Ephemeris reader does not have data for epoch " + std::to_string(epoch / 86400.0) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else
            {
                *MassPointer = gsl_spline_eval(this->Spline_mass, epoch, this->SplineAccelerator);
            }
        }

        void get7State(const double& epoch, double* StateArray)
        {
            this->getPosition(epoch, StateArray);
            this->getVelocity(epoch, StateArray + 3);
            this->getMass(epoch, StateArray + 7);
        }

        inline void getPositionDerivative(const double& epoch, double* PositionDerivativeArray)
        {
            if (epoch > this->data["epoch"].back() || epoch < this->data["epoch"].front())
            {
                throw std::runtime_error("Ephemeris reader does not have data for epoch " + std::to_string(epoch / 86400.0) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else
            {
                PositionDerivativeArray[0] = gsl_spline_eval_deriv(this->Spline_x, epoch, this->SplineAccelerator);
                PositionDerivativeArray[1] = gsl_spline_eval_deriv(this->Spline_y, epoch, this->SplineAccelerator);
                PositionDerivativeArray[2] = gsl_spline_eval_deriv(this->Spline_z, epoch, this->SplineAccelerator);
            }
        }

        inline void getVelocityDerivative(const double& epoch, double* VelocityDerivativeArray)
        {
            if (epoch > this->data["epoch"].back() || epoch < this->data["epoch"].front())
            {
                throw std::runtime_error("Ephemeris reader does not have data for epoch " + std::to_string(epoch / 86400.0) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else
            {
                VelocityDerivativeArray[0] = gsl_spline_eval_deriv(this->Spline_vx, epoch, this->SplineAccelerator);
                VelocityDerivativeArray[1] = gsl_spline_eval_deriv(this->Spline_vy, epoch, this->SplineAccelerator);
                VelocityDerivativeArray[2] = gsl_spline_eval_deriv(this->Spline_vz, epoch, this->SplineAccelerator);
            }
        }

        inline void getMassDerivative(const double& epoch, double* MassDerivativePointer)
        {
            if (epoch > this->data["epoch"].back() || epoch < this->data["epoch"].front())
            {
                throw std::runtime_error("Ephemeris reader does not have data for epoch " + std::to_string(epoch / 86400.0) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else
            {
                *MassDerivativePointer = gsl_spline_eval_deriv(this->Spline_mass, epoch, this->SplineAccelerator);
            }
        }

        void get7StateDerivative(const double& epoch, double* StateDerivativeArray)
        {
            this->getPositionDerivative(epoch, StateDerivativeArray);
            this->getVelocityDerivative(epoch, StateDerivativeArray + 3);
            this->getMassDerivative(epoch, StateDerivativeArray + 6);
        }

        void get7StateAndDerivative(const double& epoch, double* StateAndDerivativeArray)
        {
            this->get7State(epoch, StateAndDerivativeArray);
            this->get7StateDerivative(epoch, StateAndDerivativeArray + 7);
        }

    protected:

        //fields
        size_t dataRows;
        std::map< std::string, std::vector<double> > data;

        //GSL spline structs
        bool fitSplines;
        gsl_interp_accel *SplineAccelerator;
        gsl_spline *Spline_x;
        gsl_spline *Spline_y;
        gsl_spline *Spline_z;
        gsl_spline *Spline_vx;
        gsl_spline *Spline_vy;
        gsl_spline *Spline_vz;
        gsl_spline *Spline_mass;
    };//end clss writey_thing
}