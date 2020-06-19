#EMTG: Evolutionary Mission Trajectory Generator
#An open-source global optimization tool for preliminary mission design
#Provided by NASA Goddard Space Flight Center
#
#Copyright (c) 2014 - 2018 United States Government as represented by the
#Administrator of the National Aeronautics and Space Administration.
#All Other Rights Reserved.
#
#Licensed under the NASA Open Source License (the "License"); 
#You may not use this file except in compliance with the License. 
#You may obtain a copy of the License at:
#https://opensource.org/licenses/NASA-1.3
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
#express or implied.   See the License for the specific language
#governing permissions and limitations under the License.

#base class for high-fidelity trajectory driver
#Jacob Englander 1-23-2018

import os
import HighFidelityTrajectory
import Universe
import sys

if len(sys.argv) < 4:
        raise Exception("Unknown number of command line options!\n\
                         Syntax: python HighFidelityDriver.py source_emtgopt source_emtg outputFilePath\n")
                         
originalOptionsFile = sys.argv[1]
originalMissionFile = sys.argv[2]
outputFilePath = sys.argv[3]

myHighFidelityTrajectory = HighFidelityTrajectory.HighFidelityTrajectory(originalMissionPath=originalMissionFile,
                                                                         originalOptionsPath=originalOptionsFile,
                                                                         useSundman=False)

myHighFidelityTrajectory.CreateMission()
myHighFidelityTrajectory.CreateInitialGuess()

newOptions = myHighFidelityTrajectory.getMissionOptions()
newOptions.short_output_file_names = 1
newOptions.mission_type = 9 #variable phase type
newOptions.MBH_time_hop_probability = 0.0
newOptions.mission_name = originalMissionFile.replace('.emtg', '') + '_HighFidelity'
newOptions.DisassembleMasterConstraintVectors()
newOptions.DisassembleMasterDecisionVector()

# create a universe object for each central body in the new options
central_body_names = []
universe_list = []
for JourneyIndex in range(0, len(newOptions.Journeys)):
    central_body_name = newOptions.Journeys[JourneyIndex].journey_central_body
    if central_body_name not in central_body_names:
        central_body_names.append(central_body_name)
        universe_file = central_body_name + ".emtg_universe"
        universe = Universe.Universe(os.path.join(newOptions.universe_folder, universe_file))
        universe_list.append(universe)
        

# add all of the bodies in each Journey's central body universe file as perturbers
for JourneyIndex in range(0, len(newOptions.Journeys)):
    central_body_name = newOptions.Journeys[JourneyIndex].journey_central_body
    for universe in universe_list:
        if universe.central_body_name == central_body_name:
            perturbation_body_IDs = [index + 1 for index in universe.perturbation_indices]
            newOptions.Journeys[JourneyIndex].perturbation_bodies = perturbation_body_IDs

newOptions.write_options_file(outputFilePath + '/' + newOptions.mission_name + '.emtgopt')