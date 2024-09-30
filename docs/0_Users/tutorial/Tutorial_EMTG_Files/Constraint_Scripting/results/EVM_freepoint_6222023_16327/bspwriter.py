import os
mkspk_path = 'C:/cspice/exe//mkspk'
brief_path = 'C:/cspice/exe//brief'
import sys 
sys.path.append('C:/Users/ts34254/windows_dev/emtg4/PyEMTG/')
sys.path.append('C:/Users/ts34254/windows_dev/emtg4/PyEMTG//SimpleMonteCarlo')
sys.path.append('C:/Users/ts34254/windows_dev/emtg4/PyEMTG//SpiceyPy_Utilities') 
import clean_spiceicles 
clean_spiceicles.do_the_stuff(['C:/EMTG/Tutorials/Constraint_Scripting/results/EVM_freepoint_6222023_16327','EVM_freepoint.ephemeris','C:/EMTG/Tutorials/EVM_universe/ephemeris_files/', 600])
if os.path.exists('EVM_freepoint.bsp') :
    os.remove('EVM_freepoint.bsp')
os.system(mkspk_path + ' -setup EVM_freepoint.cmd -input EVM_freepoint_clean.ephemeris -output EVM_freepoint.bsp')
os.system(brief_path + ' EVM_freepoint.bsp')
