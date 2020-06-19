#helper script to convert a mission from regular CoastPhase to SundmanCoastPhase

import sys
import os
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(currentdir + "/..")
from copy import deepcopy
from MissionOptions import MissionOptions
from Mission import Mission
from math import pi
from kepler import cart2kep

# Ensure the correct number of command line options were provided
if len(sys.argv) < 3:
    raise Exception("Unknown number of command line options!\n\
                     Syntax: python convert_CoastPhase_to_SundmanCoastPhase.py source_emtgopt destination_directory")

# expand username if a relative path was used
source_emtgopt = os.path.expanduser(sys.argv[1]).replace('\\','/')
destination_dir = os.path.expanduser(sys.argv[2]).replace('\\','/')

if ".emtgopt" not in source_emtgopt:
    raise Exception("You did not pass in a valid .emtgopt file")

print('parsing ', source_emtgopt)
if not os.path.exists(source_emtgopt.replace('.emtgopt','.emtg')):
    raise Exception("cannot find " + source_emtgopt.replace('.emtgopt','.emtg'))

#create the output directory if it doesn't exist
if not os.path.exists(destination_dir):
    os.makedirs(destination_dir)

#ingest the old mission and script
originalOptions = MissionOptions(source_emtgopt)
originalMission = Mission(source_emtgopt.replace('.emtgopt','.emtg'))

#create a new options object
newOptions = deepcopy(originalOptions)

#seed it from the old mission
newOptions.trialX = []
for Xindex in range(0, len(originalMission.DecisionVector)):
    newOptions.trialX.append([originalMission.Xdescriptions[Xindex].replace('\n','').replace('\r',''), originalMission.DecisionVector[Xindex]])

newOptions.DisassembleMasterDecisionVector()
newOptions.ConvertDecisionVector()
newOptions.AssembleMasterDecisionVector()

#step through journeys, converting as-needed
for journeyIndex in range(0, newOptions.number_of_journeys):
    thisJourneyOptions = newOptions.Journeys[journeyIndex]
    thisJourney = originalMission.Journeys[journeyIndex]

    if thisJourneyOptions.phase_type == 7: #regular CoastPhase

        #Step 1: turn this journey into a SundmanCoastPhase
        thisJourneyOptions.phase_type = 8

        #Step 2: compute the true anomaly change for this journey
        initialState = thisJourney.missionevents[0].SpacecraftState
        finalState = thisJourney.missionevents[-1].SpacecraftState
        
        mu = thisJourney.mu

        initial_OE = cart2kep(initialState[0:3], initialState[3:6], mu)
        final_OE = cart2kep(finalState[0:3], finalState[3:6], mu)

        #compute the change in true anomaly
        deltaTA = final_OE[5] - initial_OE[5]

        print(initial_OE[5], final_OE[5], deltaTA)
        if deltaTA < 0.0:
            deltaTA = 2 * pi + deltaTA

        #Step 4: edit the journey's trialX
        newtrialX = []
        for entry in thisJourneyOptions.trialX:
            newtrialX.append([entry[0].replace('CoastPhase','SundmanCoastPhase'), entry[1]])
        newtrialX.append(['p0SundmanCoastPhase: phase Sundman independent variable',deltaTA])

        thisJourneyOptions.trialX = newtrialX

#assemble the master decision vector
newOptions.AssembleMasterDecisionVector()

#write the new file
newOptions.mission_name = newOptions.mission_name + '_SundmanCoast'
newOptions.write_options_file(destination_dir + '/' + newOptions.mission_name + '.emtgopt')