#test MissionOptions.py
from sys import path
path.append('c:/EMTG/OptionsOverhaul')
import MissionOptions

myOptions = MissionOptions.MissionOptions('c:/EMTG/build/tests/newdefault.emtgopt')
myOptions.write_options_file('c:/EMTG/build/tests/newdefaultPython.emtgopt')

mySecondOptions = MissionOptions.MissionOptions('c:/EMTG/build/tests/newdefaultPython.emtgopt')
mySecondOptions.write_options_file('c:/EMTG/build/tests/newerdefaultPython.emtgopt')

myOptions = MissionOptions.MissionOptions()
myOptions.write_options_file('c:/EMTG/build/tests/default.emtgopt')