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

//header for SplineEphem Universe container
//all this does is hold all of your bodies
//Jacob Englander 11-8-2016

#ifdef SPLINE_EPHEM
#ifndef SPLINEPEHEM_UNIVERSE
#define SPLINEPEHEM_UNIVERSE

#include "SplineEphem_body.h"

#include <vector>
#include <tuple>
#include <memory>

namespace SplineEphem
{
    class universe
    {
    public:
        //constructors
        universe();
        universe(const std::vector< std::tuple<int, int, int, double> >& keyList, 
            const double& ephemeris_window_open = 30000.0*86400.0,
            const double& ephemeris_window_close = 100000.0*86400.0);

        //destructors
        ~universe();

        //methods
        void reinitialize(const std::vector< std::tuple<int, int, int, double> >& keyList,
            const double& ephemeris_window_open = 30000.0*86400.0,
            const double& ephemeris_window_close = 100000.0*86400.0);
        void getBodyPosition(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* PositionArray);
        void getBodyVelocity(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* VelocityArray);
        void getBody6State(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* StateArray);
        void getBodyPositionDerivative(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* PositionDerivativeArray);
        void getBodyVelocityDerivative(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* VelocityDerivativeArray);
        void getBody6StateDerivative(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* StateDerivativeArray);
        void getBody6StateAndDerivative(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* StateAndDerivativeArray);

        double getEphemerisWindowOpen(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID);
        double getEphemerisWindowClose(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID);

    private:
        //private methods
        void initialize(const std::vector< std::tuple<int, int, int, double> >& keyList,
            const double& ephemeris_window_open = 30000.0*86400.0,
            const double& ephemeris_window_close = 100000.0*86400.0);
        void clear();
        int getBodyIndex(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID);

        //fields
        std::vector< std::unique_ptr<SplineEphem::body> > bodies;
        size_t number_of_bodies;

    };//end class universe
}//close namespace SplineEphem

#endif//SPLINEPEHEM_UNIVERSE
#endif//SPLINE_EPHEM