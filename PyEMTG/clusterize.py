#script to make a .emtgopt cluster-compatible

import MissionOptions
import sys


if len(sys.argv) != 2:
    print("call syntax is clusterize SCRIPTNAME.emtgopt")
else:
    myOptions = MissionOptions.MissionOptions(sys.argv[1])
    if myOptions.success == 1:
        myOptions.universe_folder = '/archive/Utilities/Universe/'
        myOptions.HardwarePath = '/archive/Utilities/HardwareModels/'
        
        myOptions.write_options_file(sys.argv[1])
    else:
        print(sys.argv[1] + ' not found')