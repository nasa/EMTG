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

//EMTGv8 thruster library
#include "missionoptions.h"

#pragma once

namespace EMTG 
{ 
    namespace Astrodynamics 
    {
                                        
        inline void get_thruster_coefficients_from_library(const EMTG::missionoptions& options,
        double& minP,
        double& maxP,
        double& at,
        double& bt,
        double& ct,
        double& dt,
        double& et,
        double& gt,
        double& ht,
        double& af,
        double& bf,
        double& cf,
        double& df,
        double& ef,
        double& gf,
        double& hf)
        {
            static const double min_power[] = { 2.7, 0.055, 0.0963, 0.0963, 0.64, 0.64, 2.5, 2.5 };
            static const double max_power[] = { 13.1, 0.075, 0.544, 0.544, 7.33, 7.33, 5, 4.5 };
            
            // need to adjust the index to account for non-engine options (0-5 & 14-17)
            // options.engine_type = 6 corresponds to the first item in each of these arrays
            // options.engine_type > 13 means you are using options beyond these default thrusters
            int k = 0;
            if (options.engine_type < 13 && options.engine_type > 5)
            {
                k = options.engine_type - 6;
            }

            /* For these arrays :
                0: AEPS_High_Thrust_and_Isp
                1: BIT3_High_Thrust_and_Isp
                2: Halo12_High_Thrust
                3: Halo12_High_Isp
                4: NEXTC_High_Thrust
                5: NEXTC_High_Isp
                6: PPS5000_High_Thrust
                7: PPS5000_High_Isp
            */

            //first, the coefficients for thrust (thrust in mN, power in W)
            static const double Ht[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            static const double Gt[] = { -1.147443E-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            static const double Et[] = { 4.845195E-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            static const double Dt[] = { -7.553864E+00, 0.000000E+00, -4.330000E-03, 2.841500E-01, 3.457600E-01, 4.940000E-03, 0.0, 0.0 };
            static const double Ct[] = { 5.088113E+01, -1.142860E-04, 2.965300E-01, -2.644810E+00, -5.863490E+00, 6.023900E-01, -6.525970E-06, -5.000000E-06 };
            static const double Bt[] = { -8.491459E+01, 3.685710E-02, -2.374800E-01, 9.136610E+00, 5.829300E+01, 2.616480E+01, 1.027110E-01, 8.750000E-02 };
            static const double At[] = { 1.504258E+02, -1.021140E+00, 1.458520E+00, -6.576110E+00, -1.339446E+01, 9.190150E+00, -6.545450E+01, -5.750000E+01 };

            //next, the coefficients for mass flow rate (mdot in mg/s, power in W)
            static const double Hf[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            static const double Gf[] = { 1.375908E-04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            static const double Ef[] = { -1.125028E-02, 8.943030E-08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            static const double Df[] = { 3.064583E-01, -2.331060E-05, -2.531211E+01, -3.472538E+01, 1.796000E-02, -1.390000E-02, 0.0, 0.0 };
            static const double Cf[] = { -3.758711E+00, 2.269710E-03, 7.530300E+00, 3.477985E+01, -2.894700E-01, 2.355000E-01, -4.480490E-07, -5.801260E-08 };
            static const double Bf[] = { 2.162108E+01, -9.783410E-02, 1.065476E+01, -8.527540E+00, 1.876640E+00, -4.756800E-01, 6.133670E-03, 3.122140E-03 };
            static const double Af[] = { -2.672510E+01, 1.631260E+00, -6.436300E-01, 1.095510E+00, 4.869900E-01, 2.088050E+00, -2.942510E+00, 2.198290E-01 };

            minP = min_power[k];
            maxP = max_power[k];
            at = At[k];
            bt = Bt[k];
            ct = Ct[k];
            dt = Dt[k];
            et = Et[k];
            gt = Gt[k];
            ht = Ht[k];

            af = Af[k];
            bf = Bf[k];
            cf = Cf[k];
            df = Df[k];
            ef = Ef[k];
            gf = Gf[k];
            hf = Hf[k];
        }//end get_thruster_coefficients_from_library()


    }//end namespace Astrodynamics
}//end namespace EMTG