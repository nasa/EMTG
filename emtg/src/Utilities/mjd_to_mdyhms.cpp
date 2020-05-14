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

//code to convert from MJD to month, day, year, hours, minutes, seconds
//Brent Barbee
//added to EMTG 7-24-2012


#include <cmath>

namespace EMTG
{

/* Begin function. */
void mjd_to_mdyhms(double MJD, int *month, int *day, int *year, int *hrs, int *mins, double *secs)
{
    /* Local variable declarations. */
    double JD;
    int    LMonth[12];
    double year_doub;
    double T1900;
    double leapYrs;
    int    dayofYr;
    double days;
    int    k;
    int    num_days;
    double temp;

    JD = MJD + 2400000.5;

    LMonth[0]  = 31;
    LMonth[1]  = 28;
    LMonth[2]  = 31;
    LMonth[3]  = 30;
    LMonth[4]  = 31;
    LMonth[5]  = 30;
    LMonth[6]  = 31;
    LMonth[7]  = 31;
    LMonth[8]  = 30;
    LMonth[9]  = 31;
    LMonth[10] = 30;
    LMonth[11] = 31;

    T1900 = (JD - 2415019.5)/365.25;
   
    year_doub = 1900.0 + floor(T1900);
   
    leapYrs = floor((year_doub - 1900.0 - 1.0)*0.25);
   
    days = (JD - 2415019.5) - (((year_doub-1900.0)*365.0) + leapYrs);
   
    if(days < 1.0)
    {
      
        year_doub -= 1.0;
      
        leapYrs = floor((year_doub - 1900.0 - 1.0)*0.25);
      
        days = (JD - 2415019.5) - (((year_doub-1900.0)*365.0) + leapYrs);
      
    }

    *year = (int)(year_doub);
   
    if(fmod(year_doub, 4.0) == 0.0)
    {
      
        LMonth[1] = 29;
      
    }
   
    dayofYr = (int)(floor(days));
   
    num_days = 0;
    
    k = 0;
   
    while(1)
    {
        num_days += LMonth[k];
   
        if((num_days+1) > dayofYr)
        {
            break;
        }
        else
        {
            k++;
        }
    }
   
    *month = k + 1;
   
    *day = dayofYr - (num_days - LMonth[k]);
   
    temp = (days - (double)(dayofYr))*24.0;
    *hrs  = (int)(floor(temp));
    *mins = (int)(floor((temp-(double)(*hrs))*60.0));
    *secs = (temp - (double)(*hrs) - ((double)(*mins)/60.0))*3600.0;

    /* Return from this function. */
    return;

} /* End of function. */

}//end namespace EMTG