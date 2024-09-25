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

"""
combine_PEATSA_output_csv_files.py
=========================================

Contains function that combines the contents of multiple
PEATSA output csv files (the Iteration#.csv files).

It only makes sense to use this function if all the input .csv files
have the same columns.

"""

def combine_PEATSA_output_csv_files(inputFiles, outputFile):
    """
    Function that combines the contents of multiple
    PEATSA output csv files (the Iteration#.csv files).
    
    It only makes sense to use this function if all the input .csv files
    have the same columns.

    Parameters
    ----------
    inputFiles : List of strings
        List of full file paths to all input files to be combined.
        Column headers for the output file are taken from the first element in the list
    outputFile : String
        Full file path to the csv file to which the combined contents are to be written.

    Returns
    -------
    None.

    """
    
    # open the output file
    fop = open(outputFile, "w+")
    
    # loop through the input files
    for i in range(len(inputFiles)):
        inputFile = inputFiles[i]
        fip = open(inputFile, "r")
        count = 0
        while True:
            count += 1
            line = fip.readline()
            
            if not line:
                fip.close() # close file and break loop if we have reached the end of the file
                break
            if count == 1:
                if i == 0: # only write the header line if we are reading from the first file
                    fop.write(line)
            else:
                fop.write(line) # write the data line
        
    
    # close the output file
    fop.close()
    return

def find_latest_PEATSA_output_csv_file_in_directory(directory):
    """
    Function that finds the latest "Iteration#.csv" file in a given directory
    
    Parameters
    ----------
    directory : String
        The directory in which to look for files
        
    Returns
    -------
    fileName : String
        The latest PEATSA output csv file in the directory. None if there are none
        
    """
    import os
    
    # find all relevant files in directory
    files = []
    for file in os.listdir(directory):
        if file.endswith(".csv") and file.startswith("Iteration"):
            files.append(file)
            
    # now that we have all the files, find the one with the largest number in it
    fileInt = -1
    if files:
        for file in files:
            tempInt = int(file.lstrip("Iteration").rstrip(".csv"))
            if tempInt > fileInt:
                fileInt = tempInt
                fileName = file
        fileName = os.path.join(directory, fileName)
    else:
        fileName = None
            
    
    return fileName

if __name__ == "__main__":
    basePath = 'C:/Eagle/peatsa/ResultsToCheck/maxdrymass/'
    inputPaths = [basePath + 'PEATSA_Eagle_MaxDryMass_FHE_VVVE_EVVE_2031_2036_1DaySteps',
                  basePath + 'PEATSA_Eagle_MaxDryMass_FHE_VVVEJ_VVEJ_VEEJ_2033_2034_1DaySteps_2',
                  basePath + 'PEATSA_Eagle_MaxDryMass_FHE_VVVEJ_VVEJ_VEEJ_2036_1DaySteps',
                  basePath + 'PEATSA_Eagle_VGA_MaxDryMass_FHE_VVE_VEE_2031_2036_1DaySteps_2']
    inputFiles = []
    for inputPath in inputPaths:
        inputFiles.append(find_latest_PEATSA_output_csv_file_in_directory(inputPath))
    print(inputFiles)
    outputFile = basePath + 'combined_files/attempt3.csv'

    combine_PEATSA_output_csv_files(inputFiles, outputFile)
    print('done')