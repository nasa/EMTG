#bulletproof EMTG
#automatically warm-starts NSGAII runs if they die
#needs to be put in the directory where you are running your EMTG script from

import subprocess
import sys
import copy
import os
import glob
import datetime
import time
import MissionOptions

def running_process(process):
    "check if process is running. < process > is the name of the process."

    proc = subprocess.Popen(["if pgrep " + process + " >/dev/null 2>&1; then echo 'True'; else echo 'False'; fi"], stdout=subprocess.PIPE, shell=True)

    (Process_Existance, err) = proc.communicate()
    return Process_Existance

#for later:
if len(sys.argv) != 2:
    print("call syntax is bulletproof SCRIPTNAME.emtgopt")
else:
    reference_script = sys.argv[1]
    current_script = reference_script
    runcount = 0
    #run the initial script
    commandstring = 'mpirun -np 64 ~/EMTG/emtg ' './' + current_script + ' > ' + current_script + '.out 2>&1 &'
    print(commandstring)
    os.system(commandstring)
    time.sleep(15)
    while True:
        #check if there is an EMTG running
        if running_process("emtg") == 'False\n':
            print('EMTG has crashed', datetime.datetime.now())
            #figure out what the most recent results directory was
            most_recent_results_directory = max([os.path.join('../EMTG_v8_results',d) for d in os.listdir('../EMTG_v8_results')], key=os.path.getmtime)
            
            #grab the archive and latest generation file from the most recently created results folder
            print('copying archive file')
            os.system('cp ' + most_recent_results_directory + '/NSGAII_archive.NSGAII .')
            populationList = []
            for populationfile in glob.glob(most_recent_results_directory + '/NSGAII_population_gen*.NSGAII'):
                populationList.append(int(populationfile.rstrip('.NSGAII').split('_')[-1]))
            sortedPopulationList = sorted(populationList)
            
            crashpopulation = most_recent_results_directory + '/NSGAII_population_gen_' + str(sortedPopulationList[-1]) + '.NSGAII'
            print('copying population file ', crashpopulation)
            os.system('cp ' + crashpopulation + ' .')
            
            #figure out what generation we quit on
            crashed_generation = sortedPopulationList[-1]
            print('EMTG crashed on generation ', crashed_generation)
            
            #load the reference script, set the appropriate warm start value and archive and population files, then save it with a new name. If the generation we quit on is the last generation, STOP
            OptionsStructure = MissionOptions.MissionOptions(reference_script)
            if crashed_generation == OptionsStructure.outerloop_genmax:
                #we're done, stop
                print('final generation complete, stopping')
                break
            #otherwise, make the new options script
            runcount += 1
            OptionsStructure.outerloop_warmstart = crashed_generation - 1
            OptionsStructure.outerloop_warm_population = crashpopulation.split('/')[-1]
            OptionsStructure.outerloop_warm_archive = 'NSGAII_archive.NSGAII'
            current_script = reference_script.rstrip('.emtgopt') + '_Restart' + str(runcount) + '.emtgopt'
            print('saving options file ', current_script)
            OptionsStructure.write_options_file(current_script)
            
            #run the new script
            commandstring = 'mpirun -np 64 ~/EMTG/emtg ' './' + current_script + ' > ' + current_script + '.out 2>&1 &'
            print(commandstring)
            os.system(commandstring)
        else:
            print('EMTG is alive', datetime.datetime.now())
        #sleep for a while
        #set this to 600 seconds for now
        sys.stdout.flush()
        time.sleep(600)