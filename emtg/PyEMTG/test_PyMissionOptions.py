import sys
sys.path.append('c:/EMTG/bin')

import PyMissionOptions

myPyMissionOptions = PyMissionOptions.MissionOptions('c:/emtg/scripts/DefaultNEXTHT.emtgopt')

print(myPyMissionOptions.mission_name)
print(myPyMissionOptions.Journeys)