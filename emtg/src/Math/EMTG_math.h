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

//math functions for EMTG

#pragma once

#include <math.h>
#include <stdlib.h>
#include <vector>
#include "doubleType.h"

namespace EMTG 
{ 
    namespace math
    {
        //constants
        const double PI = 3.14159265358979323;
        const double PIover2 = 3.14159265358979323 / 2.0;
        const double TwoPI = 3.14159265358979323 * 2.0;
        const double SMALL  = 1.0e-13;
        const double LARGE  = 1.0e+30;
		const double deg2rad = PI / 180.0;

        template <typename T> int sgn(T val) 
        {
            return (T(0) < val) - (val < T(0));
        }

        template<class T> inline T acosh(const T& x) {return 2*log(sqrt((x+1.0)/2.0) + sqrt((x-1.0)/2.0));}

        template<class T> inline T asinh(const T& x) {    return log (x+sqrt(1+x*x));}

        inline doubleType safe_acos(const doubleType& x)
        {
#ifdef AD_INSTRUMENTATION
            //preserve the derivative information in acos even when we clip the argument
            if (x > 1.0)
            {
                doubleType x_prime = x;
                x_prime.setValue(1.0 - SMALL);

                return acos(x_prime);
            }
            else if (x < -1.0)
            {
                doubleType x_prime = x;
                x_prime.setValue(-1.0 + SMALL);

                return acos(x_prime);
        }
            else
            {
                return acos(x);
            }
#else
            return x > 1.0 ? 0.0 : (x < -1.0 ? PI : acos(x));
#endif
        }

        template<class T> inline T fmod(const T& x, const double& modArgument)
        {
            if (x > 0.0)
            {
                T temp = x;

                while (temp > modArgument)
                {
                    temp -= modArgument;
                }

                return temp;
            }
            else if (x < 0.0)
            {
                T temp = x;

                while (temp < 0.0)
                {
                    temp += modArgument;
                }

                return temp;
            }
            else
            {
                return x;
            }
        }

        inline doubleType safe_asin(const doubleType& x)
        {
#ifdef AD_INSTRUMENTATION
            //preserve the derivative information in acos even when we clip the argument
            if (x > 1.0)
            {
                doubleType x_prime = x;
                x_prime.setValue(1.0 - SMALL);

                return asin(x_prime);
            }
            else if (x < -1.0)
            {
                doubleType x_prime = x;
                x_prime.setValue(-1.0 + SMALL);

                return asin(x_prime);
            }
            else
            {
                return asin(x);
            }
#else
            return x > 1.0 ? -PIover2 : (x < -1.0 ? PIover2 : asin(x));
#endif
        }

        template<class T> inline T norm(const std::vector<T>& A)
        {
            size_t dim = A.size();
            T norm2 = 0.0;
            for (size_t k = 0; k < dim; ++k)
                norm2 += A[k] * A[k];

            return sqrt(norm2);
        }

        inline doubleType absclip(const doubleType& X, const double& maxabs)
        {
            doubleType Xstar = X;
            if (X < -maxabs)
            {
#ifdef AD_INSTRUMENTATION
                Xstar.setValue(-maxabs);
#else
                Xstar = -maxabs;
#endif
            }
            else if (X > maxabs)
            {
#ifdef AD_INSTRUMENTATION
                Xstar.setValue(maxabs);
#else
                Xstar = maxabs;
#endif
            }

            return Xstar;
        }

    }//end namespace math
}//end namespace EMTG
