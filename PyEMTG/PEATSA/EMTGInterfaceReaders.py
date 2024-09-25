#EMTG: Evolutionary Mission Trajectory Generator
#An open-source global optimization tool for preliminary mission design
#Provided by NASA Goddard Space Flight Center
#
#Copyright (c) 2014 - 2024 United States Government as represented by the
#Administrator of the National Aeronautics and Space Administration.
#All Other Rights Reserved.
#
#Licensed under the NASA Open Source License (the "License"); 
#You may not use this file except in compliance with the License. 
#You may obtain a copy of the License at:
#https://opensource.org/license/nasa1-3-php
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
#express or implied.   See the License for the specific language
#governing permissions and limitations under the License.

#Code to read EMTG interface files
#Jeremy Knittel 2-27-2018

# Input is [year,month,date,time of day in fractional hours]
def date_list_to_JD(date_list):
    from numpy import sign, floor
    
    # Extract the year
    K = date_list[0]
    # Extract the month
    M = date_list[1]
    # Extract the day
    I = date_list[2]
    # Extract the time of day, in hours
    UT = date_list[3]
    
    # Calculate the Julian Date. 
    # From: http://aa.usno.navy.mil/faq/docs/JD_Formula.php
    return 367*K - floor( ( 7 * ( K + floor((M+9.0)/12.0) ) ) / 4.0 ) + floor( (275*M) / 9.0 ) + I + 1721013.5 + UT/24 - 0.5*sign(100 * K + M - 190002.5) + 0.5

def convert_month_string_to_int(input):
    
    # I dont think its necessary to comment this function. I think a fifth grader can understand this
    if input == "jan" or input == "Jan" or input == "JAN" or input == "JANUARY" or input == "january" or input == "January" or input == "1" or input == "01":
        return 1
    if input == "feb" or input == "Feb" or input == "FEB" or input == "FEBRUARY" or input == "february" or input == "February" or input == "2" or input == "02":
        return 2
    if input == "mar" or input == "Mar" or input == "MAR" or input == "MARCH" or input == "march" or input == "March" or input == "3" or input == "03":
        return 3
    if input == "apr" or input == "Apr" or input == "APR" or input == "APRIL" or input == "april" or input == "April" or input == "4" or input == "04":
        return 4
    if input == "may" or input == "May" or input == "MAY" or input == "5" or input == "05":
        return 5
    if input == "jun" or input == "Jun" or input == "JUN" or input == "JUNE" or input == "june" or input == "June" or input == "6" or input == "06":
        return 6
    if input == "jul" or input == "Jul" or input == "JUL" or input == "JULY" or input == "july" or input == "July" or input == "7" or input == "07":
        return 7
    if input == "aug" or input == "Aug" or input == "AUG" or input == "AUGUST" or input == "august" or input == "August" or input == "8" or input == "08":
        return 8
    if input == "sep" or input == "Sep" or input == "SEP" or input == "SEPTEMBER" or input == "september" or input == "September" or input == "9" or input == "09":
        return 9
    if input == "oct" or input == "Oct" or input == "OCT" or input == "OCTOBER" or input == "october" or input == "October" or input == "10" or input == "10":
        return 10
    if input == "nov" or input == "Nov" or input == "NOV" or input == "NOVEMBER" or input == "november" or input == "November" or input == "11" or input == "11":
        return 11
    if input == "dec" or input == "Dec" or input == "DEC" or input == "DECEMBER" or input == "december" or input == "December" or input == "12" or input == "12":
        return 12

def convert_time_string_to_hours_float(input):
    # Split the line by colon, assuming it is hours:minutes:seconds.seconds_fraction
    input_split = input.split(":")
    
    # Calculate the fractional hours past midnight
    return float(input_split[0]) + float(input_split[1]) / 60.0 + float(input_split[2]) / 3600.0

def read_spk_data_to_list_of_lists(spkin_file):
    
    # Get a handle to the file
    file_handle = open(spkin_file)
    
    # Initialize the output data array
    data_out = [[] for idx in range(8)]
    
    # Loop through all lines
    for line in file_handle:
        
        if line[0] == "#":
            continue
        
        # Split the line based on spaces, removing the newlines
        linesplit = line.rstrip(" \r\n").split(",")
        
        # Remove any blank entries in the list
        try:
            linesplit.remove("")
        except:
            pass
            
        date_split = linesplit[0].replace("  "," ").split(" ")
                        
        # Get the epoch
        JD = date_list_to_JD([int(date_split[0]), # year
                              convert_month_string_to_int(date_split[1]), # month
                              int(date_split[2]), # days
                              convert_time_string_to_hours_float(date_split[3])]) # hours
        
        # Add the julian date
        data_out[0].append(JD)
        
        # Loop through the rest of the data
        for idx in range(1,len(linesplit)):
            data_out[idx].append(float(linesplit[idx]))
                    
    # Return the ephemeris data
    return data_out
        
         