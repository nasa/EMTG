# EMTG: Evolutionary Mission Trajectory Generator
# An open-source global optimization tool for preliminary mission design
# Provided by NASA Goddard Space Flight Center
#
# Copyright (c) 2013 - 2020 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Other Rights Reserved.
#
# Licensed under the NASA Open Source License (the "License"); 
# You may not use this file except in compliance with the License. 
# You may obtain a copy of the License at:
# https://opensource.org/licenses/NASA-1.3
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
# express or implied.   See the License for the specific language
# governing permissions and limitations under the License.

#doohicky that makes missionoptions and journeyoptions C++ files based on JSON inputs
#Jacob Englander 1/8/2019

import csv
import time
from optionValidator import validate
 
#get the current epoch
now = time.strftime("%c")

EMTG_path = 'C:/EMTG/'

EMTG_MissionOptionsStructure_file = EMTG_path + 'OptionsOverhaul/list_of_missionoptions.csv'
EMTG_JourneyOptionsStructure_file = EMTG_path + 'OptionsOverhaul/list_of_journeyoptions.csv'

#we need lists of all of the options of interest
missionOptionsDefinitions = []
journeyOptionsDefinitions = []

with open(EMTG_JourneyOptionsStructure_file) as csvFile:
    reader = csv.reader(csvFile)
    for row in reader:
        if reader.line_num == 1:
            header = row
        else:
            optionDictionary = {}
            for key, cell in zip(header, row):
                if cell != '':
                    optionDictionary[key] = cell
            validate(optionDictionary)
            journeyOptionsDefinitions.append(optionDictionary)

with open(EMTG_MissionOptionsStructure_file) as csvFile:
    reader = csv.reader(csvFile)
    for row in reader:
        if reader.line_num == 1:
            header = row
        else:
            optionDictionary = {}
            for key, cell in zip(header, row):
                if cell != '':
                    optionDictionary[key] = cell
            validate(optionDictionary)
            missionOptionsDefinitions.append(optionDictionary)

#create journeyoptions.h
from make_journeyoptions_header import *
make_journeyoptions_header(journeyOptionsDefinitions, now, path=EMTG_path + 'src/Core')
#make_journeyoptions_header(journeyOptionsDefinitions, now)

#create journeyoptions.cpp
from make_journeyoptions_source import *
make_journeyoptions_source(journeyOptionsDefinitions, now, path=EMTG_path + 'src/Core')
#make_journeyoptions_source(journeyOptionsDefinitions, now)

#create missionoptions.h
from make_missionoptions_header import *
make_missionoptions_header(missionOptionsDefinitions, now, path=EMTG_path + 'src/Core')
#make_missionoptions_header(missionOptionsDefinitions, now)

#create missionoptions.cpp
from make_missionoptions_source import *
make_missionoptions_source(missionOptionsDefinitions, now, path=EMTG_path + 'src/Core')
#make_missionoptions_source(missionOptionsDefinitions, now)

#build PyMissionOptions.h
from make_PyMissionOptions import *
#make_PyMissionOptions(journeyOptionsDefinitions, missionOptionsDefinitions, now, path=EMTG_path + 'PyEMTG')

#build PyEMTG JourneyOptions
from make_journeyoptions_python import *
make_PyEMTG_JourneyOptions(journeyOptionsDefinitions, now, path=EMTG_path + 'PyEMTG')
#make_PyEMTG_JourneyOptions(journeyOptionsDefinitions, now)

#build PyEMTG MissionOptions
from make_missionoptions_python import *
make_PyEMTG_MissionOptions(missionOptionsDefinitions, now, path=EMTG_path + 'PyEMTG')
#make_PyEMTG_MissionOptions(missionOptionsDefinitions, now)
