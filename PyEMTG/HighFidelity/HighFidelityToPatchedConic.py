import sys
sys.path.append("/Users/jknittel/EMTG/PyEMTG")
sys.path.append("/Users/jknittel/EMTG/PyEMTG/PEATSA")
import MissionOptions

MO = MissionOptions.MissionOptions("CAESAR_3952kg_launch_open_307kg_SRC_CustomReturnSC.emtgopt")

import PEATSAchef

chef = PEATSAchef.PEATSAchef()

newMO = chef.ThreeD_Flybys_to_PatchedConic(MO)

newMO.mission_name = "PatchedConicLaunchcOpen"

newMO.write_options_file("PatchedConicLaunchOpen.emtgopt")