import sys
import os
os.environ["OMP_NUM_THREADS"] = "4" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "4" # export OPENBLAS_NUM_THREADS=4 

if len(sys.argv) != 5:
    raise Exception("Unknown number of command line options!\n\
                     Syntax: python make_many_bsps.py options_file working_directory destination_dir path_to_EMTG_executable")

thisdir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(thisdir)
sys.path.append(thisdir + '/../')
sys.path.append(thisdir + '/../SpiceyPy_Utilities')

options_file = sys.argv[1]
working_directory = sys.argv[2]
destination_dir = sys.argv[3]
EMTGpath = sys.argv[4]

import MissionOptions
import Mission

#does a seed exist?
seedfilename = os.path.join(working_directory, options_file.replace('.emtgopt','.emtg'))
if os.path.exists(seedfilename):    
    #load
    myOptions = MissionOptions.MissionOptions(os.path.join(working_directory, options_file))
        
    #seed
    seedMission = Mission.Mission(seedfilename)

    myOptions.trialX = []
    for Xindex in range(0, len(seedMission.DecisionVector)):
        myOptions.trialX.append([seedMission.Xdescriptions[Xindex].replace('\n','').replace('\r',''), seedMission.DecisionVector[Xindex]])

    myOptions.DisassembleMasterDecisionVector()
    myOptions.ConvertDecisionVector()    
        
    #redirect output
    myOptions.forced_working_directory = destination_dir + '/'
    myOptions.override_mission_subfolder = 0
        
    #enable ephemeris printing
    myOptions.generate_forward_integrated_ephemeris = 1

    #set up SPICE writer
    myOptions.spice_utilities_path = '/Utilities/cspice/exe/'
    myOptions.spice_utility_extension ='' #Linux
    myOptions.pyemtg_path = '/home/jaengla2/emtg/PyEMTG/'

    #trialX mode
    myOptions.run_inner_loop = 0
        
    #save
    myOptions.write_options_file(destination_dir + '/' + options_file)
        
    #run!
    print('running case ' + destination_dir + '/' + options_file)
    os.system(EMTGpath + ' ' + destination_dir + '/' + options_file + ' > ' + destination_dir + '/' + options_file + '.out &')