#script to convert an entire directory of .emtg and .emtgopt files all the way from MGALT to PSFB high fidelity


import shutil
import sys
import os
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(currentdir + "/..")
sys.path.append(currentdir + "/../HighFidelity/")

import Mission
import MissionOptions
import Universe
import Convert_TwoPointShootingLowThrust_to_PSFB
import HighFidelityTrajectory


#a handy function
def run_options_file(optionsFile, destination_dir, path_to_EMTG):
    file = os.path.abspath(destination_dir + optionsFile)
    command = path_to_EMTG + ' ' + file + ' > ' + file + '.out'
    print(command)
    os.system(command)

# Ensure the correct number of command line options were provided
if len(sys.argv) < 4:
    raise Exception("Unknown number of command line options!\n\
                     Syntax: python batch_convert_MGALT_to_PSFB_HiFi.py source_directory destination_directory path_to_EMTG [useSundman]\n\
                     where useSundman is True if you want Sundman gravity assist, False if you don't, and defaults to false.")

# expand username if a relative path was used
source_dir = os.path.expanduser(sys.argv[1])
destination_dir = os.path.expanduser(sys.argv[2])
path_to_EMTG = os.path.expanduser(sys.argv[3])
useSundman = False
if len(sys.argv) > 4:
    useSundman = os.path.expanduser(sys.argv[4])

#create the output directory if it doesn't exist
if not os.path.exists(destination_dir):
    os.makedirs(destination_dir)

#create a temporary directory for high-fidelity MGALT outputs
if not os.path.exists(destination_dir + '/temp_high_fidelity/'):
    os.makedirs(destination_dir + '/temp_high_fidelity/')

#convert to high-fidelity
print('Converting from patched-conic to high-fidelity flybys.')
for root,dirs,files in os.walk(source_dir):
    for file in files:
        if file.endswith('.emtg') and 'FAILURE' not in file:
            missionFile = file
            optionsFile = file.replace('.emtg','.emtgopt')
            missionPath = os.path.join(root, missionFile)
            optionsPath = os.path.join(root, optionsFile)

            myHighFidelityTrajectory = HighFidelityTrajectory.HighFidelityTrajectory(originalMissionPath=missionPath,
                                                                                     originalOptionsPath=optionsPath,
                                                                                     useSundman=useSundman)

            myHighFidelityTrajectory.CreateMission()
            myHighFidelityTrajectory.CreateInitialGuess()

            newOptions = myHighFidelityTrajectory.getMissionOptions()
            newOptions.forced_post_launch_coast = 0.0 #to avoid conflict with the constraint for journeys
            newOptions.forced_pre_flyby_coast = 0.0 #to avoid conflict with the constraint for journeys
            newOptions.forced_post_flyby_coast = 0.0 #to avoid conflict with the constraint for journeys
            newOptions.short_output_file_names = 1
            newOptions.mission_type = 9 #variable phase type
            newOptions.MBH_time_hop_probability = 0.0
            newOptions.mission_name = missionFile.replace('.emtg', '') + '_HighFidelity'
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
                    universe.filename = central_body_name
                    universe_list.append(universe)
        

            # add all of the bodies in each Journey's central body universe file as perturbers
            for JourneyIndex in range(0, len(newOptions.Journeys)):
                central_body_name = newOptions.Journeys[JourneyIndex].journey_central_body
                for universe in universe_list:
                    if universe.central_body_name == central_body_name or universe.filename == central_body_name:
                        newOptions.Journeys[JourneyIndex].perturbation_bodies = []
                        for bodyIndex in universe.perturbation_indices:
                            if universe.bodies[bodyIndex].mu > 100.0:
                                newOptions.Journeys[JourneyIndex].perturbation_bodies.append(bodyIndex + 1)
                        break
                print(newOptions.Journeys[JourneyIndex].journey_name, newOptions.Journeys[JourneyIndex].perturbation_bodies)

            #write the script in trialX mode
            newOptions.background_mode = 1
            newOptions.run_inner_loop = 0
            newOptions.override_working_directory = 1
            newOptions.forced_working_directory = destination_dir + '/temp_high_fidelity/'
            newOptions.override_mission_subfolder = 0
            newOptions.write_options_file(destination_dir + '/temp_high_fidelity/' + newOptions.mission_name + '.emtgopt')
            
            print('wrote ' + destination_dir + '/temp_high_fidelity/' + newOptions.mission_name)

#run all of the new high-fidelity scripts in trialX mode to produce new output files
print('Running high-fidelity flyby cases to produce initial guesses')
optionsFiles = [file for file in os.listdir(destination_dir + '/temp_high_fidelity/') if (file.endswith(".emtgopt"))]
#do this in parallel if you can
try:
    from joblib import Parallel, delayed
    Parallel(n_jobs=-1)(delayed(run_options_file)(optionsFile, destination_dir + '/temp_high_fidelity/', path_to_EMTG) for optionsFile in optionsFiles)
except ImportError:
    for optionsFile in optionsFiles:
        run_options_file(optionsFile, destination_dir + '/temp_high_fidelity/', path_to_EMTG)

#now we have a bunch of .emtg and .emtgopt files, so let's make them into PSFB
print('Creating PSFB cases')
for root,dirs,files in os.walk(destination_dir + '/temp_high_fidelity/'):
    for file in files:
        if file.endswith('.emtg'):
            missionFile = file
            optionsFile = file.replace('.emtg','.emtgopt').replace('FAILURE_','')
            missionPath = os.path.join(root, missionFile)
            optionsPath = os.path.join(root, optionsFile)

            myPSFB_Converter = Convert_TwoPointShootingLowThrust_to_PSFB.MissionConverter_TPSLT_to_PSFB(originalMissionPath=missionPath,
                                                                                            originalOptionsPath=optionsPath,
                                                                                            OriginalTranscription=['FBLT', 3])

            myPSFB_Converter.CreateMission()

            myPSFB_Converter.CreateInitialGuess()

            newOptions = myPSFB_Converter.getMissionOptions()

            newOptions.trialX = myPSFB_Converter.CreateInitialGuess()
            
            newOptions.DisassembleMasterDecisionVector()

            newOptions.ConvertDecisionVector()

            newOptions.short_output_file_names = 1

            newOptions.mission_type = 9 #variable phase type
            
            newOptions.propagatorType = 1#integrator

            newOptions.mission_name = missionFile.replace('.emtg','').replace('FBLT','PSFB').replace('FAILURE_','')
            
            #some special purpose code! might want to remove this at some later date
            newOptions.Journeys[-1].CoastPhaseMatchPointFraction = 0.9

            #write the script in MBH mode for regional exploitation
            #newOptions.run_inner_loop = 1
            #newOptions.MBH_Pareto_alpha = 1.5
            #newOptions.MBH_max_step_size = 0.1
            newOptions.background_mode = 1
            newOptions.forced_working_directory = destination_dir
            newOptions.override_mission_subfolder = 1
            newOptions.forced_mission_subfolder = '/'
            newOptions.write_options_file(destination_dir + '/' + newOptions.mission_name + '.emtgopt')
            
            
            print('wrote ' + destination_dir + '/' + newOptions.mission_name)

#delete the temporary directory?
#TODO add that after testing
print('Running high-fidelity flyby, PSFB cases to produce initial guesses')
optionsFiles = [file for file in os.listdir(destination_dir) if (file.endswith(".emtgopt"))]
#do this in parallel if you can
try:
    from joblib import Parallel, delayed
    Parallel(n_jobs=-1)(delayed(run_options_file)(optionsFile, destination_dir + '/', path_to_EMTG) for optionsFile in optionsFiles)
except ImportError:
    for optionsFile in optionsFiles:
        run_options_file(optionsFile, destination_dir + '/', path_to_EMTG)
    
#rename output files if they have "FAILURE" in the name
missionFiles = [file for file in os.listdir(destination_dir) if (file.endswith(".emtg"))]

print('Renaming "failure" files')
for file in missionFiles:
    if 'FAILURE_' in file:
        os.rename(destination_dir + file, destination_dir + file.replace('FAILURE_',''))