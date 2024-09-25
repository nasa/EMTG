#EMTG: Evolutionary Mission Trajectory Generator
#An open-source global optimization tool for preliminary mission design
#Provided by NASA Goddard Space Flight Center
#
#Copyright (c) 2014 - 2024 United States Government as represented by the
#Administrator of the National Aeronautics and Space Administration.
#All Other Rights Reserved.
#
#Licensed under the NASA Open Source License (the "License"); 
#You may not use this file except in compliance with the License. 
#You may obtain a copy of the License at:
#https://opensource.org/license/nasa1-3-php
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
#express or implied.   See the License for the specific language
#governing permissions and limitations under the License.

"""
PEATSA.py
============
This file contains the main PEATSA executable script.

Calling sequence: python PEATSA.py <path/to/options/file.py>

PEATSA nomenclature:
    
   PEATSAbox: sortable container for a mission and its options script
   PEATSAchef: creates .emtgopt files
   PEATSAcrust: output (.emtg files)
   PEATSAdelivery: plots
   PEATSAdough: input (.emtgopt files)
   PEATSAmeal: a PEATSA iteration
   PEATSAmenu: holds options for PEATSA run (i.e., the PEATSAtoppings)
   PEATSAoven: executes EMTG cases
   PEATSAtopping: the input options
   PEATSAwaiter: pre/post-processor

Jeremy Knittel 1-23-2018
"""

def main():
    """
    The PEATSA main() function requires a single command-line input argument: A PEATSA options script.
    
    The function performs the following steps:
    
    #1) Reads a PEATSA options script (command-line argument)
    
    #2) Creates a set of .emtgopt files
    
    #3) Executes the .emtgopt files to produce .emtg files
    
    #4) Analyzes the resulting .emtg files and oganizes them 
    
    #4) Creates plots
    
    #5) Determines which .emtgopts made poor .emtgs and need to be rerun
    
    #6) Picks better initial guesses for the cases that need to be rerun
    
    #7) Returns to #3)
    """
    
    print("--------------------- Starting PEATSA ---------------------")
    
    # Import system modules
    import sys
    import logging
    import os
    import time
    import shutil
    
    dir_path = os.path.dirname(os.path.realpath(__file__)) # get path to PEATSA.py
    sys.path.append(dir_path + "/..") # add PyETMG directory to path
    
    # Import PEATSA modules
    try:
        import PyHardware
        pyhardware_warning = False
    except:
        pyhardware_warning = True
        print ("WARNING: PyHardware not available. Might cause errors. You have been warned")
    import PEATSAmenu
    import PEATSAchef
    import PEATSAwaiter
    import PEATSAbusboy
    import PEATSAoven
    import PEATSAdelivery
    import PEATSAbox
    
    # Ensure that an options script was specified
    if len(sys.argv) == 1:
        raise Exception('Specify the options script! Calling sequence: python PEATSA.py <path/to/options/file.py>')
    elif len(sys.argv) > 2:
        raise Exception('Too many options specified. Just give one options script! Calling sequence: python PEATSA.py <path/to/options/file.py>')
    
    print("Attempting to load options script: " + sys.argv[1])
    
    # Create the options class
    PEATSAorder = PEATSAmenu.PEATSAmenu(sys.argv[1])
    
    print("Attempting to open logfile: " + PEATSAorder.logfile)
    
    # Set up a log file script
    try:
        logging.basicConfig(filename=PEATSAorder.logfile,filemode='w',level=logging.INFO,format='%(asctime)s %(message)s')
        logging.info("Successfully opened logfile")
        print("Successfully opened logfile, all status updates will now be logged instead of printed to console.")
    except:
        raise Exception("Error opening logfile")
        
    def myError(excType, excValue, traceback):
        """
        Define a custom error handler so that the errors go to the logfile, not the console
        """
        logging.error("*************** PEATSA encountered an error. ****************",
                     exc_info=(excType, excValue, traceback))
    sys.excepthook = myError
    
    # start peatsa
    logging.info("Welcome to PEATSA.")
    logging.info("Running version " + sys.version)
    
    # check if PyHardware class was imported
    if pyhardware_warning:
        logging.info("WARNING: PyHardware not available. Might cause errors. You have been warned")
        
    logging.info("Making working directories")
    
    PEATSAorder.working_directory = os.path.expanduser(PEATSAorder.working_directory)
    
    root_working_directory = PEATSAorder.working_directory
    
    # keep_only_current_and_previous requires a fresh start and trade study study type
    if PEATSAorder.keep_only_current_and_previous == True:
        if PEATSAorder.start_type != "Fresh" or PEATSAorder.PEATSA_type != 2:
            raise Exception("keep_only_current_and_previous == True is currently only compatible with a Fresh start and PEATSA_type == 2.")
            
    
    # If warm-starting, then a new folder is not needed
    if PEATSAorder.start_type == "Warm":
        PEATSAorder.working_directory = PEATSAorder.working_directory + "/" + PEATSAorder.run_name + "/"
        
        # Get a path of the results and case folders
        base_results_dir = PEATSAorder.working_directory + 'results/'
        base_cases_dir = PEATSAorder.working_directory + 'cases/'
        base_docs_dir = PEATSAorder.working_directory + 'docs/'
        base_images_dir = PEATSAorder.working_directory + 'images/'
    else:
        # A new working folder is needed. But, first we need to make sure we arent
        # writing over previous results
        
        # Get a base name for the working directory without anything appended
        base_name = PEATSAorder.working_directory + "/" + PEATSAorder.run_name
        new_name = base_name
        counter = 1
        
        # While the name exists, try increasing the counter. This way the working
        # folder will be the first unused working_directory_#
        while os.path.isdir(new_name):
            counter += 1
            new_name = base_name + "_" + str(counter)
            
        # Found it. Update the peatsa order
        PEATSAorder.working_directory = new_name + "/"
        os.mkdir(PEATSAorder.working_directory)
        
        # Get a path of the results and case folders
        base_results_dir = PEATSAorder.working_directory + 'results/'
        base_cases_dir = PEATSAorder.working_directory + 'cases/'
        base_docs_dir = PEATSAorder.working_directory + 'docs/'
        base_images_dir = PEATSAorder.working_directory + 'images/'
        if PEATSAorder.start_type == "Fresh":
            base_hardware_dir = PEATSAorder.working_directory + 'HardwareModels/'
        
        # Make the main four working subdirectories
        os.mkdir(base_results_dir)
        os.mkdir(base_cases_dir)
        os.mkdir(base_docs_dir)
        os.mkdir(base_images_dir)
        if PEATSAorder.start_type == "Fresh":
            os.mkdir(base_hardware_dir)
        
        # Since this is a new peatsa run, set the iteration
        PEATSAorder.iteration = 0
    
    if os.path.isdir(root_working_directory + "/running"):
        os.remove(root_working_directory + "/running")
    os.symlink(PEATSAorder.working_directory,root_working_directory + "/running")
    
    # Copy the options to the working directory
    logging.info("Copying options script to working directory")
    PEATSAorder.write_to_file(PEATSAorder.working_directory + 'PEATSA_ExecutedOptions.py')
    
    # Loop through built in plotter files
    plot_file_ct = 0
    new_plot_files = []
    for plotfile in PEATSAorder.built_in_plotter_files:
        # Load the plot options
        PO = PEATSAmenu.PEATSAinstagram(plotfile)
        
        # Write out the plot options
        PO.write_options(PEATSAorder.working_directory + "Plot" + str(plot_file_ct) + "Options.py")
        
        # Add the new plot file name to the list
        new_plot_files.append(PEATSAorder.working_directory + "Plot" + str(plot_file_ct) + "Options.py")
        
        # Update the counter
        plot_file_ct += 1
    
    # Point the peatsa options at the updated plot files
    PEATSAorder.built_in_plotter_files = new_plot_files
    
    # Write out the midbake options
    PEATSAorder.write_to_file(PEATSAorder.working_directory + 'PEATSA_MidBake_Options.py',True,True)
    
    # Put the docs directory on the PEATSAorder object
    PEATSAorder.docs_dir = base_docs_dir
    PEATSAorder.images_dir = base_images_dir
    if PEATSAorder.start_type == "Fresh":
        PEATSAorder.hardware_dir = base_hardware_dir
    
    # Fix home folder relative paths if they exist
    PEATSAorder.restart_run_root_directory = [(os.path.expanduser(restart_dir[0]),restart_dir[1]) for restart_dir in PEATSAorder.restart_run_root_directory]
    PEATSAorder.emtg_root_directory = os.path.expanduser(PEATSAorder.emtg_root_directory)
    
    # Create the oven (running emtg) class
    oven = PEATSAoven.PEATSAoven()
    # Create the waiter (re-seeder) class
    waiter = PEATSAwaiter.PEATSAwaiter()
    # Create the busboy (harvester) class
    busboy = PEATSAbusboy.PEATSAbusboy()
    # Create the delivery (plotting) class
    delivery = PEATSAdelivery.PEATSAdelivery()
    
    #create an empty "boxes2run" variable
    boxes2run = []
    boxes = []
    
    # create a dictionary to store all of the used seeds
    history = {'seed':[],'run':[],'objective':[], 'firstSolveFeasible':[], 'bestSolutionIndex':[], 'solutionAttempts':[]}
    
    # Check if we are running things in parallel
    # if PEATSAorder.parse_in_parallel:
        # # I dont really know how this stuff works. I got it from here:
        # # https://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods
        # # But, it is necessary to run things in parallel
        
        # try:
            # def _pickle_method(method):
                # func_name = method.im_func.__name__
                # obj = method.im_self
                # cls = method.im_class
                # return _unpickle_method, (func_name, obj, cls)
            
            # def _unpickle_method(func_name, obj, cls):
                # for cls in cls.mro():
                    # try:
                        # func = cls.__dict__[func_name]
                    # except KeyError:
                        # pass
                    # else:
                        # break
                # return func.__get__(obj, cls)
            
            
            # from joblib import Parallel, delayed
            # import copy_reg
            # import types
        
            # copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
        
        # except:
            # logging.info("Cant setup parallel, so running in serial")
            # PEATSAorder.parse_in_parallel = 0
          
    # Load the extenral seed boxes. Do so now, so that if there is a problem, it doesnt take down
    # the peatsa run  
    external_seed_boxes = busboy.load_external_seed_boxes(PEATSAorder)
    
    # Set flags regarding if midbake options have changed
    update_external_seed_cases_flag = False
    update_seed_criteria_or_fingerprint_flag = False
    
    # Set a flag to indicate this is our first warm iteration if it is
    if PEATSAorder.start_type == "Warm":
        if_first_warm_iteration = 1
         
    # Iterate until we reach a break
    notdone = 1
    while notdone:
        
        # Add a new dictionary to the history lists
        history['seed'].append({})
        history['run'].append({})
        history['objective'].append({})
        history['firstSolveFeasible'].append({})
        history['bestSolutionIndex'].append({})
        history['solutionAttempts'].append({})
        
        if PEATSAorder.iteration > 0:
            for key in history['objective'][-2].keys():
                history["objective"][-1].update({key:(history["objective"][-2][key][0],history["objective"][-2][key][1],0,0)})
        
        # Get the paths to the current results and case folders
        if PEATSAorder.keep_only_current_and_previous == True:
            PEATSAorder.results_directory = base_results_dir + "IterationCurrent/"
            PEATSAorder.cases_directory   = base_cases_dir   + "IterationCurrent/"
        else:
            PEATSAorder.results_directory = base_results_dir + "Iteration" + str(PEATSAorder.iteration) + "/"
            PEATSAorder.cases_directory   = base_cases_dir   + "Iteration" + str(PEATSAorder.iteration) + "/"
        PEATSAorder.csv_file          = base_docs_dir    + "Iteration" + str(PEATSAorder.iteration) + ".csv"
        
        # If this isnt the first entry of a warm start, then we need to create the results and case directories
        # if PEATSAorder.keep_only_current_and_previous, this will be IterationCurrent
        if not PEATSAorder.start_type == "Warm":
            # if this isn't the first iteration of a fresh start
            # that is using keep_only_current_and_previous, then the
            # IterationCurrent directory already exists, so we don't
            # want to recreate it.
            if (not PEATSAorder.keep_only_current_and_previous) or  PEATSAorder.iteration == 0:
                os.mkdir(PEATSAorder.cases_directory)
                os.mkdir(PEATSAorder.results_directory)
        
        # Handle the starting conditions
        if PEATSAorder.iteration == 0 or (PEATSAorder.start_type == "Warm" and if_first_warm_iteration == 1):
            if PEATSAorder.start_type == "Fresh":
                # We are starting a brand new peatsa run    
                
                # Create a chef to create all of the .emtgopt files
                chef = PEATSAchef.PEATSAchef()
                # Create the emtgopt files
                boxes2run = chef.create_PEATSAs(PEATSAorder)
                # Initialize the boxes list
                boxes = []
                # Check which cases actually need to be run
                boxes2run = waiter.check_only_run_if(PEATSAorder,boxes2run,True)
                
                # If user wants to go ahead and seed these cases from the external cases
                if PEATSAorder.seed_from_seed_folders_on_fresh_start and len(external_seed_boxes):
                    
                    # If the objective type is not those, then the seed cases need to be found now
                    boxes2run = waiter.find_all_seed_cases(PEATSAorder,boxes2run,boxes,external_seed_boxes)
            
                    # Create new .emtgopt files 
                    boxes2run = waiter.reseed(PEATSAorder,boxes2run,history)
                
                
            elif PEATSAorder.start_type == "Hot":
                # At least one iteration of results already exists.
                    
                boxes = []
                for result_dir in PEATSAorder.restart_run_root_directory:
                    if result_dir[1] == 'Best':
                        if not os.path.isdir(PEATSAorder.working_directory + "/BestFromPrevious/"):
                            os.mkdir(PEATSAorder.working_directory + "/BestFromPrevious/")
                            
                        delivery.CopyBestPeatsas(result_dir[0],PEATSAorder.working_directory + "/BestFromPrevious/",-1)
                        
                        # Call the busboy to parse the results
                        boxes, history = busboy.parse_results(PEATSAorder,PEATSAorder.working_directory + "/BestFromPrevious/",boxes,0,history)
                    
                    elif result_dir[1] == 'Full':
                        # Call the waiter to parse the results_directory
                        boxes_out, history = busboy.fresh_parse_results(PEATSAorder,result_dir[0],history)
                        boxes += boxes_out
                    else:
                        raise Exception("Unknown parse type")
                    
                # Sort the results and eliminate lesser cases
                boxes = busboy.sort_and_filter(PEATSAorder,boxes)
                                    
                # Write the 0th iteration results to file
                busboy.write_to_csv(PEATSAorder,boxes)
        
                # Write the history file
                busboy.write_history_file(PEATSAorder,history)
                
                # If the options say so, copy the results files into the iteration 0 folder
                if PEATSAorder.copy_previous_results == 1:
                    boxes = busboy.copy_results(PEATSAorder,boxes)
        
                # Call the delivery
                delivery.post_process(PEATSAorder,boxes)
                
                # No cases to be run yet. Need to loop back around now
                PEATSAorder.iteration += 1
                continue
                
            elif PEATSAorder.start_type == "Warm":
                # Flip the switch
                if_first_warm_iteration = 0
        else:
            # We are not on the first iteration, so we need to determine which cases need to be re-run
            
            # Check which cases need to be re-run
            if PEATSAorder.PEATSA_type == 5:
                boxes2run,boxes = waiter.propulate(PEATSAorder,boxes)
            else:
                boxes2run = waiter.check_optimality(PEATSAorder,boxes)
            
            # If the objective type is 1, the seed cases were already found. Either way, if we 
            # have new external seed cases, then we need to refind seed cases
            if PEATSAorder.objective_type != 1 or update_external_seed_cases_flag == True:
                
                # Check if we need to update the external boxes
                if update_external_seed_cases_flag == True:
                    external_seed_boxes = busboy.load_external_seed_boxes(PEATSAorder)
                
                # If the objective type is not those, then the seed cases need to be found now
                boxes2run = waiter.find_all_seed_cases(PEATSAorder,boxes2run,boxes,external_seed_boxes)
            
            # Create new .emtgopt files 
            boxes2run = waiter.reseed(PEATSAorder,boxes2run,history)
         
        # Check to make sure there are still cases to run 
        if len(boxes2run) == 0:
            notdone = False
            continue
          
        # Check if we can fill all cores with the current PEATSA run
        # We must first check how many cores we have available:
        net_cores = 0
        # Loop through all cores
        for host in PEATSAorder.nCores:
            net_cores += host[1]
        # Check if have more cores than cases
        if net_cores > len(boxes2run) and len(boxes2run) > 0:
            import MissionOptions
            # We have more cores than cases. lets create some copies of cases so that all cores are filled
            loop_number = 1
            new_boxes = []
            while net_cores > len(boxes2run) + len(new_boxes):
                loop_number += 1
                # Loop through all boxes
                for box in boxes2run:
                    # Load the options into a mission options objet
                    MO = MissionOptions.MissionOptions(box.PEATSApath + "/" + box.PEATSAdough_path)
                    # Give the object the new mission name
                    MO.mission_name += "_v" + str(loop_number)
                    # Save the new mission options object to file with the new mission name
                    MO.write_options_file(box.PEATSApath + "/" + MO.mission_name + ".emtgopt")
                    # Create a new box for this case
                    new_box = PEATSAbox.PEATSAbox(PEATSAorder,box.PEATSApath,MO.mission_name + ".emtgopt")
                    # Append the new box
                    new_boxes.append(new_box)
                    
                    # Check if we have created enough cases yet
                    if net_cores <= len(boxes2run) + len(new_boxes):
                        break
            # Add the new boxes into the boxes to run
            for box in new_boxes:
                boxes2run.append(box)
        
        if PEATSAorder.if_run_cases == 0:
            logging.info("Debugging has been completed, go check what happened")
            raise Exception("Options say not to run cases. If you did want to run cases set if_run_cases = 1")
        # Before running the cases, update the seed history
        for box in boxes2run:
            # Get the base mission name, as this will be used as the key in the dictionary
            key = box.mission_name.split("_seeded_by")[0].split("_v")[0]
            # Check if this key is already in the dictionary
            if key in history['seed'][-1].keys():
                # It is. so multiple versions of this case are being run. Add the case to list
                history['seed'][-1][key].append(box.used_seed)
            else:
                # This case is not yet in the history. Add it, along with the seed details as the first entry in a list
                history['seed'][-1].update({key:[box.used_seed]}) 
                   
        # Run the cases
        logging.info('made it to the oven with ' + str(len(boxes2run)))
        if PEATSAorder.execution_type == 1:
            oven.executeV2(PEATSAorder, boxes2run)
        elif PEATSAorder.execution_type == 2 or PEATSAorder.execution_type == 3:
            oven.executeWithPebble(PEATSAorder, boxes2run)
        else:
            oven.cook(PEATSAorder,boxes2run)
        
        logging.info('made it through the oven')
        
        # Store the relevent options so that we know if they have changed when we load the midbake options
        local_seed_folders = PEATSAorder.seed_folders
        local_fingerprint = PEATSAorder.fingerprint
        local_seed_criteria = PEATSAorder.seed_criteria
        local_extra_csv_columns = PEATSAorder.extra_csv_column_definitions
        
        # Update the mid-bake options in case the user changed something
        try:
            PEATSAorder.load_options_from_file(PEATSAorder.working_directory + 'PEATSA_MidBake_Options.py',True)
        except Exception as e:
            logging.info("WARNING: Unable to update options mid-bake:" + str(e))
            if "text" in dir(e):
                logging.info("*******************************************")
                logging.info(e.text)
                logging.info(' ' * e.offset + "^")
                logging.info("*******************************************")
        
        # Check if local_seed_folders have changed
        if local_seed_folders != PEATSAorder.seed_folders:
            # They did change, so update the flag
            update_external_seed_cases_flag = True
        else:
            # They did not change, but make sure the flag is off
            update_external_seed_cases_flag = False
        
        # Check if the seed criteria or fingerprint formulas changed    
        if local_fingerprint != PEATSAorder.fingerprint or local_seed_criteria != PEATSAorder.seed_criteria or local_extra_csv_columns != PEATSAorder.extra_csv_column_definitions:
            # They did change. update the flag
            update_seed_criteria_or_fingerprint_flag = True
            
            # need to update the fingerprint of the external seeds, too
            update_external_seed_cases_flag = True
        else:
            # They did not change, but make sure the flag is off
            update_seed_criteria_or_fingerprint_flag = False
    
        # Call the busboy to parse the results
        logging.info('made it to parse_results')
        boxes, history = busboy.parse_results(PEATSAorder,PEATSAorder.results_directory,boxes,update_seed_criteria_or_fingerprint_flag,history)
        #logging.info(boxes[0].PEATSApath)
        #logging.info("had " + str(len(boxes)) + " cases")
        
        # Sort the results and eliminate lesser cases
        boxes = busboy.sort_and_filter(PEATSAorder,boxes)
        #logging.info("kept " + str(len(boxes)) + " cases")
        #logging.info("After filtering:")
        #for box in boxes:
            #string = box.PEATSApath + box.PEATSAcrust_path
            #logging.info(string)
        
        # Combine stuff from "IterationCurrent" to "IterationPrevious"
        # Also need to update the information in the spreadsheet about where
        # the results are.
        if PEATSAorder.keep_only_current_and_previous:
            if PEATSAorder.iteration == 0:
                # if iteration 0, then we need to create the "IterationPrevious" directories
                os.mkdir(base_cases_dir   + "IterationPrevious/")
                os.mkdir(base_results_dir   + "IterationPrevious/")
                
                # update the results directories for each "box"
                # they were IterationCurrent, but they need to be changed to IterationPrevious
                for box in boxes:
                    box.PEATSApath = base_results_dir   + "IterationPrevious/"
                
                # also, there will not actually be anything in IterationPrevious
                # to start with, so "combining" will actually just be moving
                # Current into Previous
                
                # move cases
                sourceDir = PEATSAorder.cases_directory
                destinationDir = base_cases_dir   + "IterationPrevious/"
                filesToMove = os.listdir(sourceDir)
                for f in filesToMove:
                    os.rename(sourceDir + f, destinationDir + f)
                    
                # move results
                sourceDir = PEATSAorder.results_directory
                destinationDir = base_results_dir   + "IterationPrevious/"
                filesToMove = os.listdir(sourceDir)
                for f in filesToMove:
                    os.rename(sourceDir + f, destinationDir + f)
                    
            else:
                # IterationPrevious directories already exist and contain stuff
                # we want to keep everything that is in boxes and get rid of everything that is
                # not in boxes
                sourceDirCases = PEATSAorder.cases_directory
                destinationDirCases = base_cases_dir   + "IterationPrevious/"
                sourceDirResults = PEATSAorder.results_directory
                destinationDirResults = base_results_dir   + "IterationPrevious/"
                for box in boxes:
                    if "Previous" not in box.PEATSApath:
                        # if Previous IS in PEATSApath then we don't have to do anything
                        # because the one in previous is the one we want to keep and it is
                        # already there.
                        # on the other hand, if Previous is NOT in PEATSApath, then Current
                        # is the one we want to keep. We need to get rid of any of the same
                        # fingerprint in cases/IterationPrevious and results/IterationPrevious
                        # before moving over the one from Current.
                        # update the directories for each "box"
                        # they were IterationCurrent, but they need to be changed to IterationPrevious
                        string = box.PEATSApath + box.PEATSAcrust_path
                        box.PEATSApath = base_results_dir   + "IterationPrevious/"
                        
                        # if an emtgopt file with the same fingerprint (case number) exists in IterationPrevious already, delete it
                        # remove "seeded_by" and "v" because these still reference the same case,
                        # just different versions of it
                        key = box.mission_name.split("_seeded_by")[0].split("_v")[0]
                        #caseNumberString = box.PEATSAdough_path.lstrip("TradeStudy_").rstrip(".emtgopt")
                        #underscoreIndex = caseNumberString.find("_")
                        #if underscoreIndex > -1:
                        #    caseNumberString[0:underscoreIndex]
                            
                        # at this point, caseNumberString is just "Case###"
                        import glob
                        # filesToDelete will find everything that starts with our TradeStudy case index, which is what we want
                        # will remove all the files, regardless of file name
                        #filesToDelete = glob.glob(destinationDirCases + "TradeStudy_" + caseNumberString + "*")
                        
                        # key gives us TradeStudy_Case0 (e.g.)
                        # delete everything that starts with that, followed immediately by and underscore
                        filesToDelete = glob.glob(destinationDirCases + key + "_*")
                        for f in filesToDelete:
                            os.remove(f)
                            
                        # delete from old results as well as from old cases
                        filesToDelete = glob.glob(destinationDirResults + key + "_*")
                        for f in filesToDelete:
                            os.remove(f)
                            
                        # also delete any FAILURE emtg files
                        filesToDelete = glob.glob(destinationDirResults + "FAILURE_*")
                        for f in filesToDelete:
                            os.remove(f)
                            
                        # also delete any snopt crash files
                        filesToDelete = glob.glob(destinationDirResults + "*.SNOPTcrash")
                        for f in filesToDelete:
                            os.remove(f)
                        
                        # move the case to Previous
                        if os.path.exists(sourceDirCases + box.PEATSAdough_path):
                            os.rename(sourceDirCases + box.PEATSAdough_path, destinationDirCases + box.PEATSAdough_path)
                        
                        #logging.info("moved " + sourceDirCases + box.PEATSAdough_path + " to " + destinationDirCases + box.PEATSAdough_path)
                        
                        
                        
                        # move the result to Previous
                        # also need to move the emtg_spacecraftopt file
                        # also need to move the _probe.emtg file if it exists
                        # the .emtgopt needs to be in both places
                        spacecraftoptFile = box.PEATSAcrust_path.replace('.emtg', '.emtg_spacecraftopt')
                        probeFile = box.PEATSAcrust_path.replace('.emtg', '_probe.emtg')
                        
                        #logging.info("moving " + sourceDirResults + box.PEATSAdough_path + " to " + destinationDirResults + box.PEATSAdough_path)
                        #logging.info("moving " + sourceDirResults + box.PEATSAcrust_path + " to " + destinationDirResults + box.PEATSAcrust_path)
                        #logging.info("moving " + sourceDirResults + spacecraftoptFile + " to " + destinationDirResults + spacecraftoptFile)
                                                
                        
                        if os.path.exists(sourceDirResults + box.PEATSAdough_path):
                            os.rename(sourceDirResults + box.PEATSAdough_path, destinationDirResults + box.PEATSAdough_path)
                        if os.path.exists(sourceDirResults + box.PEATSAcrust_path):
                            os.rename(sourceDirResults + box.PEATSAcrust_path, destinationDirResults + box.PEATSAcrust_path)
                        if os.path.exists(sourceDirResults + spacecraftoptFile):
                            os.rename(sourceDirResults + spacecraftoptFile, destinationDirResults + spacecraftoptFile)
                        if os.path.exists(sourceDirResults + probeFile):
                            os.rename(sourceDirResults + probeFile, destinationDirResults + probeFile)
                        
                        
                
            # finish by empyting out the IterationCurrent directories
            # in preparation for the next iteration
            import glob
            filesToDelete = glob.glob(PEATSAorder.cases_directory + "*")
            for f in filesToDelete:
                os.remove(f)
            filesToDelete = glob.glob(PEATSAorder.results_directory + "*")
            for f in filesToDelete:
                os.remove(f)
                
        # Write this iteration's results to file
        busboy.write_to_csv(PEATSAorder,boxes,history)
        
        # Write the history file
        PEATSAorder = busboy.write_history_file(PEATSAorder,history)
        
        # Call the delivery
        delivery.post_process(PEATSAorder,boxes)
            
        
        # We are all done 
        PEATSAorder.iteration += 1
        if PEATSAorder.iteration >= PEATSAorder.max_iterations:
            break
        
    logging.info("All done! The meal is over. Closing the restaurant. See ya when you want your next slice")
        

    return

if __name__ == "__main__":
    main()
    
    