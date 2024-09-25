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
PEATSAbusboy.py
==================

Class to analyze all of the .emtg files from a PEATSA iteration and put them into a list of PEATSAbox objects

"""


def parallelParser(instance, name, i, nCases, args=(), kwargs=None):
    """
    Non-class method required for parallelization.
    We use this to call parse_fun instead of calling parse_fun directly.
    See torek's response from April 2018 at https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map/7309686#7309686
    """
    import logging
    if kwargs is None:
        kwargs = {}
    logging.info(str(i) + '/' + str(nCases))
    return getattr(instance, name)(*args, **kwargs)
    
def parallelParserRegions(instance, name, region, PEATSAorder, emtgopt_files, results_folder, kwargs=None):
    """
    Non-class method required for parallelization.
    We use this to call parse_fun instead of calling parse_fun directly.
    See torek's response from April 2018 at https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map/7309686#7309686
    """
    import logging
    if kwargs is None:
        kwargs = {}
        
    #i = region[0] # lower bound index
    results = []
    #logging.info(str(i))
    #while i <= region[1]: # upper bound index
    #    args = (PEATSAorder, emtgopt_files[i], results_folder,)
    #    results.append(getattr(instance, name)(*args, **kwargs))
        #logging.info(str(i))
    logging.info('started parsing a region')
    for emtgopt_file in emtgopt_files:
        args = (PEATSAorder, emtgopt_file, results_folder,)
        results.append(getattr(instance, name)(*args, **kwargs))
    logging.info('finished parsing a region')
    return results

def parallelParserNotReadingRegions(instance, name, PEATSAorder, PEATSAboxes, kwargs=None):
    """
    Non-class method required for parallelization.
    We use this to call parse_fun instead of calling parse_fun directly.
    See torek's response from April 2018 at https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map/7309686#7309686
    """
    import logging
    if kwargs is None:
        kwargs = {}
        
    results = []

    logging.info('started post-read parsing a region')
    for box in PEATSAboxes:
        args = (PEATSAorder, box,)
        results.append(getattr(instance, name)(*args, **kwargs))
    logging.info('finished post-read parsing a region')
    return results
    
def parallelCallback(x):
    """
    Callback is required for parallelization. This doesn't actually do anything.

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    return

class PEATSAbusboy(object):
    """
    Used to analyze all of the .emtg files from a PEATSA iteration and put them into a list of PEATSAbox objects
    
    Parameters
    ----------
    None.

    Returns
    -------
    None.
    """
        
    def load_external_seed_boxes(self,PEATSAorder):
        """
        Description.

        Parameters
        ----------
        PEATSAorder : TYPE
            DESCRIPTION.

        Returns
        -------
        external_seed_boxes : TYPE
            DESCRIPTION.

        """
        import os
        
        # Initialize an output list
        external_seed_boxes = []
        
        # loop through all seed folders and parse them
        for seed_folder in PEATSAorder.seed_folders:
    
            # Check what type of folder this is
            if seed_folder[0] == 0:
        
                # Get a list of all the subfolders in the main results directory. Each of these should be an iteration directory
                subfolders = [seed_folder[1] + "/" + folder + "/" for folder in os.listdir(seed_folder[1]) if os.path.isdir(seed_folder[1] + "/" + folder)]
    
                # Loop through all of the subfolders (iteration directories)
                for folder in subfolders:
            
                    # Parse the results for this subfolder and add it to the rest of the results
                    external_seed_boxes = self.parse_results(PEATSAorder,folder,external_seed_boxes)
            
            elif seed_folder[0] == 1:
        
                # Parse the results for this subfolder and add it to the rest of the results
                external_seed_boxes = self.parse_results(PEATSAorder,seed_folder[1],external_seed_boxes)
    
        if len(external_seed_boxes):
            # Sort the seed results and eliminate lesser cases
            external_seed_boxes = self.sort_and_filter(PEATSAorder,external_seed_boxes)    
        
        external_ctr = 0
        for box in external_seed_boxes:
            external_ctr += 1
            box.mission_name = "External_Case" + str(external_ctr) + "_" + box.mission_name
            box.in_study = False
        
        # Return the list of boxes
        return external_seed_boxes
        
    def fresh_parse_results(self,PEATSAorder,restart_dir,history=None):
        """
        Description

        Parameters
        ----------
        PEATSAorder : TYPE
            DESCRIPTION.
        restart_dir : TYPE
            DESCRIPTION.
        history : TYPE, optional
            DESCRIPTION. The default is None.

        Raises
        ------
        Exception
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        import sys
        import os
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import Mission
        import logging
        logging.info("Parsing a full results root folder")
        
        # Initialize the output array
        boxes_out = []
                
        # Get a list of all the subfolders in the main results directory. Each of these should be an iteration directory
        subfolders = [restart_dir + "/" + folder + "/" for folder in os.listdir(restart_dir) if os.path.isdir(restart_dir + "/" + folder)]
    
        # Loop through all of the subfolders (iteration directories)
        for folder in subfolders:
                            
            if history == None:
                # Parse the results for this subfolder and add it to the rest of the results
                boxes_out = self.parse_results(PEATSAorder,folder,boxes_out,False,history)
            else:
                # Parse the results for this subfolder and add it to the rest of the results
                boxes_out, history = self.parse_results(PEATSAorder,folder,boxes_out,False,history)
        
        # Check to mkae sure we got some results
        if len(boxes_out) == 0:
            # Maybe the user forgot the results folder on top of the peatsa path   
            
            # Get a list of all the subfolders in the main results directory. Each of these should be an iteration directory
            subfolders = [restart_dir + "/results/" + folder + "/" for folder in os.listdir(restart_dir + "/results") if os.path.isdir(restart_dir + "/results/" + folder)]
    
            # Loop through all of the subfolders (iteration directories)
            for folder in subfolders:
                
                if history == None:
                    # Parse the results for this subfolder and add it to the rest of the results
                    boxes_out = self.parse_results(PEATSAorder,folder,boxes_out,False,history)
                else:
                    # Parse the results for this subfolder and add it to the rest of the results
                    boxes_out, history = self.parse_results(PEATSAorder,folder,boxes_out,False,history)
                    
        # Check to see if we have results now
        if len(boxes_out) == 0:
            raise Exception("Not able to find any emtg results. Check to make sure that the 'restart_run_root_directory' path is correct")
        
        # Return the list of all results
        if history == None:
            return boxes_out
        else:
            return boxes_out, history
        
    def parse_results(self,PEATSAorder,results_folder,PEATSAboxes,update_seed_criteria_or_fingerprint_flag = False,history = None):
        """
        Description.

        Parameters
        ----------
        PEATSAorder : TYPE
            DESCRIPTION.
        results_folder : TYPE
            DESCRIPTION.
        PEATSAboxes : TYPE
            DESCRIPTION.
        update_seed_criteria_or_fingerprint_flag : TYPE, optional
            DESCRIPTION. The default is False.
        history : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        import sys
        import os
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        import PEATSAbox
        logging.info("Parsing folder and updating results: " + results_folder)
        
        # Check if we need to reload the boxes
        if update_seed_criteria_or_fingerprint_flag:
            # we do
            tempBoxes = []
            
            # Loop through all the boxes
            for box in PEATSAboxes:
                # Append a new box after reloading the data
                tempBoxes.append(PEATSAbox.PEATSAbox(PEATSAorder,box.PEATSApath,box.PEATSAdough_path,box.PEATSAcrust_path))
                
            # Update the list of boxes with the new data
            PEATSAboxes = tempBoxes
            
        #logging.info("Before parsing:")
        #for box in PEATSAboxes:
        #    string = box.PEATSApath + "/" + box.PEATSAcrust_path
        #    logging.info(string)
        # Get a listing of all the emtgopt files in this directory.
        emtgopt_files = [file for file in os.listdir(results_folder) if file.endswith('.emtgopt')]
        
        # way faster to do this now and save it than do it in parse_fun
        self.emtg_files = [file for file in os.listdir(results_folder) if file.endswith('.emtg')]

        # Check if parsing will be done in parallel
        if PEATSAorder.parse_in_parallel:
            import multiprocessing as mp

            ######################################
            # new way of parallel post-processing:
            nCases = len(emtgopt_files) # number of things we need to post-process
            if nCases == 0:
                if history == None:
                    return PEATSAboxes
                else:
                    return PEATSAboxes, history
                # no cases: parallelism will break if we allow it to happen
                # there are no new_boxes to add to peatsa boxes
                # just break
            
            # read the cases in serial
            #logging.info('Reading cases prior to parallel processing.')
            #new_boxes = []
            #for emtgoptfile in emtgopt_files:
            #    # Call the parse function to get the new box
            #    # Add the new case to the list
            #    new_boxes.append(self.parse_fun_read_only(PEATSAorder,emtgoptfile,results_folder))

            #logging.info('Finished reading cases prior to parallel processing.')
            # now, do the parallelizable stuff in parallel
            nProcessors = mp.cpu_count() # get the actual number of processors
            
            processorsDesired = min(nProcessors, PEATSAorder.parse_in_parallel_nCores) # number of processors: the minimum of the number of processors we have and the number the user requested
            logging.info('Starting parallel post-processing')
            nWorkers = min(processorsDesired, nCases)
            logging.info('Using ' + str(nWorkers) + ' processes for ' + str(nCases) + ' cases')
            pool = mp.Pool(nWorkers) # create the parallel pool
                
            regions = []
            nCasesPerWorker = int(nCases / nWorkers)
            for i in range(nWorkers):
                lb = i * nCasesPerWorker
                if i == nWorkers - 1:
                    ub = nCases - 1
                else:
                    ub = (i + 1) * nCasesPerWorker - 1
                regions.append((lb, ub))
                
            for i in range(len(regions)):
                logging.info('Region ' + str(i) + ': ' + str(regions[i][0]) + ' to ' + str(regions[i][1]))
                
            #results = [pool.apply_async(parallelParserNotReadingRegions, args = (self, 'parse_fun_after_reading', PEATSAorder, new_boxes[regions[i][0]:regions[i][1]+1],), callback = parallelCallback) for i in range(len(regions))]
            results = [pool.apply_async(parallelParserRegions, args = (self, 'parse_fun', regions[i], PEATSAorder, emtgopt_files[regions[i][0]:regions[i][1]+1], results_folder,), callback = parallelCallback) for i in range(len(regions))]
            pool.close()
            map(mp.pool.ApplyResult.wait, results)
            pool.join()
            
            logging.info('Finished parallel post-processing')
            new_boxes = [r.get() for r in results] # grab the results of the parsing
            list_of_list_of_boxes = [r.get() for r in results] # results is a list of lists. we want just a single list of all those items
            new_boxes = [item for sublist in list_of_list_of_boxes for item in sublist]
            logging.info('Put processed results in new_boxes')
            
            # do the parsing:
            #results = [pool.apply_async(parallelParser, args = (self, 'parse_fun', i, nCases, (PEATSAorder, emtgopt_files[i], results_folder,)), callback = parallelCallback) for i in range(nCases)]
            #pool.close()
            #pool.join()
            #logging.info('Finished parallel post-processing')
            #new_boxes = [r.get() for r in results] # grab the results of the parsing
            #logging.info('Put processed results in new_boxes')
            
            ##################################################
            
            if history != None:
                for new_box in new_boxes:
                    history = self.update_history(PEATSAorder,new_box,history)
            # Add in the new boxes
            PEATSAboxes += new_boxes        

        else:
            # Loop through all of the emtgopt files in serial
            for emtgoptfile in emtgopt_files:
                
                # Call the parse function to get the new box
                new_box = self.parse_fun(PEATSAorder,emtgoptfile,results_folder)
                
                # Add the new case to the list
                PEATSAboxes.append(new_box)
                
                if history != None:
                    history = self.update_history(PEATSAorder,new_box,history)
        
        #logging.info("After parsing:")
        #for box in PEATSAboxes:
        #    string = box.PEATSApath + "/" + box.PEATSAcrust_path
        #    logging.info(string)
            
        # Return the full list
        if history == None:
            return PEATSAboxes
        else:
            return PEATSAboxes, history
        


    def update_history(self,PEATSAorder,new_box,history):
        """
        Description
        """
        
        # Get the base mission name, as this will be used as the key in the dictionary
        key = new_box.mission_name.split("_seeded_by")[0].split("_v")[0]

        # Check if the case finished:
        if new_box.PEATSAcrust_path == "FAILURE_fake.emtg":
            ifFinished = 0
        else:
            ifFinished = 1
        
        seed = None
        if PEATSAorder.iteration > 0:    
            # Find the seed from the seed history
            seed = -1
            if "seeded_by" not in new_box.mission_name:
                # There isnt actually a seed
                seed = None
            else:
                # Excract hte seeded by string
                seed_target = new_box.mission_name.split("seeded_by")[1]
                # Loop through the seeds in the history for this case
                for possible_seed in history['seed'][-1][key]:
                    # Check if this seed name is in the seeded by string
                    if possible_seed != None and possible_seed[0] in seed_target:
                        # It is, we found it
                        seed = possible_seed
                        break
                    
            # Make sure we found the seed, and throw an error if not
            if seed == -1:
                raise Exception("Cant find this case in the seed history: " + new_box.mission_name)
        
        # Did this case converge?
        if new_box.converged:
            # Check if this case was seedeed
            if seed == None:
                # It was not, so yes, it improved from its "seed"
                ifImprovedFromSeed = 1
            else:
                # It was seeded. Get the seeds' objective value
                seed_obj_val = seed[1]
                
                # Check if the seed has a better objective value or the case does
                if PEATSAorder.max_or_min == "max" and new_box.objective_value > seed_obj_val:
                    ifImprovedFromSeed = 1
                elif PEATSAorder.max_or_min == "min" and new_box.objective_value < seed_obj_val:
                    ifImprovedFromSeed = 1
                else:
                    ifImprovedFromSeed = 0                        
        else:
            # The case did not converge, so clearly it did not do as well as the seed
            ifImprovedFromSeed = 0    
        
        
        # Check if this key is already in the dictionary
        if key in history['run'][-1].keys():
            # It is. so multiple versions of this case are being run. Add the case to list
            history['run'][-1][key].append((ifFinished,new_box.converged,ifImprovedFromSeed,seed))
        else:
            # This case is not yet in the history. Add it, along with the seed details as the first entry in a list
            history['run'][-1].update({key:[(ifFinished,new_box.converged,ifImprovedFromSeed,seed)]}) 
        
        # Check if this key is the objetive dictionary
        if key not in history['objective'][-1].keys():
            # Add this key to the history
            history['objective'][-1].update({key:(0,0,0,0)}) 
            
        # Check if this key is in the firstSolveFeasible dictionary
        if key not in history['firstSolveFeasible'][-1].keys():
            # Add this key to the history
            history['firstSolveFeasible'][-1].update({key:(0)}) # just one element in the tuple: 1 if first attempt converged, 0 if it did not
            
        # Check if this key is in the bestSolutionIndex dictionary
        if key not in history['bestSolutionIndex'][-1].keys():
            # Add this key to the history
            history['bestSolutionIndex'][-1].update({key:(0)}) # just one element in the tuple: 0 if no converged solutions, int for best attempt otherwise
            
         # Check if this key is in the solutionAttempts dictionary
        if key not in history['solutionAttempts'][-1].keys():
            # Add this key to the history
            history['solutionAttempts'][-1].update({key:(0)}) # just one element in the tuple: number of solution attempts
            
        # how many solution attempts? relevant whether feasible or not
        history['solutionAttempts'][-1][key] = (new_box.number_of_solution_attempts)
        
        # Check if this case converged
        if new_box.converged:
            # did the first solution attempt converge?
            history['firstSolveFeasible'][-1][key] = (new_box.first_nlp_solve_feasible)
            
            # what was the solution attempt index of the best solution found?
            history['bestSolutionIndex'][-1][key] = (new_box.solution_attempt_index_that_produced_best_feasible_solution)
    
            # Check if this is the first iteration
            if PEATSAorder.iteration == 0:
                # It is. 
            
                # Check if there is already a converged case
                if history['objective'][-1][key][0]:
                    # There is. 
                    
                    #Check if this one did better
                    if PEATSAorder.max_or_min == "min" and history['objective'][-1][key][1] > new_box.objective_value:
                        # It did. Update the tuple
                        history['objective'][-1][key] = (1,new_box.objective_value,1,1)
                    elif PEATSAorder.max_or_min == "max" and history['objective'][-1][key][1] < new_box.objective_value:
                        # It did. Update the tuple
                        history['objective'][-1][key] = (1,new_box.objective_value,1,1)
                        
                # Check if this case converged
                else:
                    # It did. Update the history
                    history['objective'][-1][key] = (1,new_box.objective_value,1,1)
            # This is not the first iteration, so we have to compare to previous iterations
            else:
                # Check if this case has been converged before
                if history['objective'][-2][key][0] == 0:
                    # It has not.                 
                    
                    # Make sure something better form this iteration isnt already in there
                    if history['objective'][-1][key][0]:     
                        #Check if this one did better
                        if PEATSAorder.max_or_min == "min" and history['objective'][-1][key][1] > new_box.objective_value:
                            # It did. Update the tuple
                            history['objective'][-1][key] = (1,new_box.objective_value,1,1)
                        elif PEATSAorder.max_or_min == "max" and history['objective'][-1][key][1] < new_box.objective_value:
                            # It did. Update the tuple
                            history['objective'][-1][key] = (1,new_box.objective_value,1,1)           
                    else:
                        # There is not anything in there, so this is certainly better    
                        history['objective'][-1][key] = (1,new_box.objective_value,1,1)
                else:
                    # It has been converged before. Check if it has improved
                    if PEATSAorder.max_or_min == "min" and history['objective'][-2][key][1] > new_box.objective_value:
                        # Make sure something better form this iteration isnt already in there
                        if history['objective'][-1][key][0]:     
                            if history['objective'][-1][key][1] > new_box.objective_value:
                                # It did. Update the tuple
                                history['objective'][-1][key] = (1,new_box.objective_value,0,1)           
                        else:
                            # There is not anything in there, so this is certainly better    
                            history['objective'][-1][key] = (1,new_box.objective_value,0,1)
                    elif PEATSAorder.max_or_min == "max" and history['objective'][-2][key][1] < new_box.objective_value:
                        # Make sure something better form this iteration isnt already in there
                        if history['objective'][-1][key][0]:     
                            if history['objective'][-1][key][1] < new_box.objective_value:
                                # It did. Update the tuple
                                history['objective'][-1][key] = (1,new_box.objective_value,0,1)           
                        else:
                            # There is not anything in there, so this is certainly better    
                            history['objective'][-1][key] = (1,new_box.objective_value,0,1)
        
        return history

    def parse_fun(self,PEATSAorder,emtgoptfile,results_folder):
        """
        Description
        
        """
        import PEATSAbox
        import os
        
        # Get a listing of all the emtgopt files in this directory.
        #if PEATSAorder.parse_in_parallel:            
        # Check if this file is in the directory
        if emtgoptfile.rstrip("opt") in self.emtg_files:
            # It is! Hurray
            emtgfile = emtgoptfile.rstrip("opt")
        elif "FAILURE_" + emtgoptfile.rstrip("opt") in self.emtg_files:
            # This case didnt converge, but its in here
            emtgfile = "FAILURE_" + emtgoptfile.rstrip("opt")
        else:
            # Huh. Weird. Couldnt find the file. Lets just create fake one
            emtgfile = "FAILURE_fake.emtg"

        # else:
            # emtg_files = [file for file in os.listdir(results_folder) if file.endswith('.emtg')]
            
            # # Check if this file is in the directory
            # if emtgoptfile.rstrip("opt") in emtg_files:
                # # It is! Hurray
                # emtgfile = emtgoptfile.rstrip("opt")
            # elif "FAILURE_" + emtgoptfile.rstrip("opt") in emtg_files:
                # # This case didnt converge, but its in here
                # emtgfile = "FAILURE_" + emtgoptfile.rstrip("opt")
            # else:
                # # Huh. Weird. Couldnt find the file. Lets just create fake one
                # emtgfile = "FAILURE_fake.emtg"

        # Initialize a new peatsa box object for this directory
        case = PEATSAbox.PEATSAbox(PEATSAorder,results_folder,emtgoptfile,emtgfile)

        # return the box
        return case
    
    def parse_fun_read_only(self,PEATSAorder,emtgoptfile,results_folder):
        """
        Description
        
        """
        import PEATSAbox
        import os
        
        # Get a listing of all the emtgopt files in this directory.
        #if PEATSAorder.parse_in_parallel:            
        # Check if this file is in the directory
        if emtgoptfile.rstrip("opt") in self.emtg_files:
            # It is! Hurray
            emtgfile = emtgoptfile.rstrip("opt")
        elif "FAILURE_" + emtgoptfile.rstrip("opt") in self.emtg_files:
            # This case didnt converge, but its in here
            emtgfile = "FAILURE_" + emtgoptfile.rstrip("opt")
        else:
            # Huh. Weird. Couldnt find the file. Lets just create fake one
            emtgfile = "FAILURE_fake.emtg"

        # else:
            # emtg_files = [file for file in os.listdir(results_folder) if file.endswith('.emtg')]
            
            # # Check if this file is in the directory
            # if emtgoptfile.rstrip("opt") in emtg_files:
                # # It is! Hurray
                # emtgfile = emtgoptfile.rstrip("opt")
            # elif "FAILURE_" + emtgoptfile.rstrip("opt") in emtg_files:
                # # This case didnt converge, but its in here
                # emtgfile = "FAILURE_" + emtgoptfile.rstrip("opt")
            # else:
                # # Huh. Weird. Couldnt find the file. Lets just create fake one
                # emtgfile = "FAILURE_fake.emtg"

        # Initialize a new peatsa box object for this directory
        # read only, don't parse!
        case = PEATSAbox.PEATSAbox(PEATSAorder,results_folder,emtgoptfile,emtgfile,initializationType = 1)

        # return the box
        return case
    
    def parse_fun_after_reading(self,PEATSAorder,PEATSAbox):
        """
        Call this only after parse_fun_read_only() has been called on a given PEATSAbox
        
        """
        import PEATSAbox
        import os
        
        case = PEATSAbox.initializeAfterReading(PEATSAorder)

        # return the box
        return case

    def sort_and_filter(self,PEATSAorder,PEATSAboxes,idx = 0):
        """
        Description.

        Parameters
        ----------
        PEATSAorder : TYPE
            DESCRIPTION.
        PEATSAboxes : TYPE
            DESCRIPTION.
        idx : TYPE, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        temp_cases : TYPE
            DESCRIPTION.

        """
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        import PEATSAbox
        logging.info("Sorting Results for Iteration " + str(PEATSAorder.iteration))
        
        # We want the best case to be first, so when we sort from min to max, we need to multiply
        # a maximization objective value by negative 1. This is not the objective function as evaluate
        # in EMTG, we are defining our own. 
        if PEATSAorder.max_or_min == "max":
            sign = -1
        elif PEATSAorder.max_or_min == "min":
            sign = 1
       
        # Get the seed criteria
        sort_criteria = PEATSAorder.seed_criteria[idx] 

        # Sort based on the fingerprint, and the objective value. if even the objective value is the same then put put the older case first (the one with the lesser value of first_iteration)
        PEATSAboxes = sorted(PEATSAboxes, key = lambda box: box.fingerprint[sort_criteria[0]] + (sign*box.objective_value,-box.iteration))
        
        logging.info("Filtering out sub-optimal cases with the same fingerprint")
        
        # Initialize the new filtered list
        temp_cases = [PEATSAboxes[0]]

        # get the last fingerprint to be added to the stack 
        last_fingerprint = PEATSAboxes[0].fingerprint[sort_criteria[0]]
        
        # Loop through all the cases
        import os
        for box in PEATSAboxes:
            
            # Get the fingerprint of this case
            current_fingerprint = box.fingerprint[sort_criteria[0]]
            
            # Check if it is a new fingerprint. If it is new, then
            # it will be added to the new stack. If it is not new, it
            # is going to get thrown away      
            if current_fingerprint != last_fingerprint:
                # It is new!

                # add this case to the stack
                temp_cases.append(box)
                
                # Update the last fingerprint 
                last_fingerprint = current_fingerprint
                
                #logging.info("busboy kept " + box.PEATSApath + box.PEATSAcrust_path)
                #logging.info("busboy kept " + box.PEATSApath + box.PEATSAdough_path)
            
            if box not in temp_cases:
                #logging.info("busboy did not keep " + box.PEATSApath + box.PEATSAcrust_path)
                #logging.info("busboy did not keep " + box.PEATSApath + box.PEATSAdough_path)
                
                if PEATSAorder.keep_only_current_and_previous:
                    # straight-up delete the files of ones that are not kept

                    # delete the .emtg, .emtgopt, and .emtg_spacecraftopt from results
                    if os.path.exists(PEATSAorder.working_directory + 'results/' + box.PEATSAcrust_path):
                        os.remove(PEATSAorder.working_directory + 'results/' + box.PEATSAcrust_path)
                        #logging.info("busboy removed " + box.PEATSApath + box.PEATSAcrust_path)
                    if os.path.exists(PEATSAorder.working_directory + 'results/' + box.PEATSAdough_path):
                        os.remove(PEATSAorder.working_directory + 'results/' + box.PEATSAdough_path)
                        #logging.info("busboy removed " + box.PEATSApath + box.PEATSAdough_path)

                    # use the emtgopt (dough) path instead of the emtg (crust) path
                    # because the emtg name might be FAILURE_fake, which is NOT the name of the spacecraft options file
                    spacecraftoptFile = PEATSAorder.working_directory + 'results/' + box.PEATSAdough_path.replace('.emtgopt', '.emtg_spacecraftopt')
                    #logging.info("want to remove " + spacecraftoptFile)
                    if os.path.exists(spacecraftoptFile):
                        os.remove(spacecraftoptFile)
                        #logging.info("busboy removed " + spacecraftoptFile)
                        
                    # delete the .emtgopt from cases directory
                    base_cases_dir = PEATSAorder.working_directory + 'cases/'
                    if "Previous" in box.PEATSApath:
                        emtgoptFile = base_cases_dir + "IterationPrevious/" + box.PEATSAdough_path
                    else:
                        emtgoptFile = base_cases_dir + "IterationCurrent/" + box.PEATSAdough_path
                    #logging.info("want to remove " + emtgoptFile)
                    if os.path.exists(emtgoptFile):
                        os.remove(emtgoptFile)
                        #logging.info("busboy removed " + emtgoptFile)
                
            
        # All done, return the result
        return temp_cases
        
    def write_to_csv(self,PEATSAorder,PEATSAboxes,history=None):
        """
        Writes the .csv summary file for each iteration.

        Parameters
        ----------
        PEATSAorder : TYPE
            DESCRIPTION.
        PEATSAboxes : TYPE
            DESCRIPTION.
        history : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        import PEATSAbox
        logging.info("Writing Results for Iteration " + str(PEATSAorder.iteration))
        
        # Open a handle to the file
        csv_file_handle = open(PEATSAorder.csv_file,'w')

        # Write the constant headers
        csv_file_handle.write("Folder,Filename,Mission Name,Instances Created,Instances Ran,Instances Finished,Instances Converged,First Feasible,Improved,Last Improved Iteration,Seed History,")
        for idx in range(len(PEATSAorder.seed_criteria)):
            csv_file_handle.write("SeedCriteria " + str(idx) + ",")
        csv_file_handle.write("Final Mass [kg],Flight Time [years],Feasibility,PEATSA objective,First NLP Solve Feasible,First NLP Solve Objective Function Value,Number of Solution Attempts,Best Solution Attempt,")
        
        # Loop through the user requested columns and add their headers
        for column in PEATSAorder.extra_csv_column_definitions:
            csv_file_handle.write(column[0] + ",")
        
        for box in PEATSAboxes:
            if box.converged:        
                # Always print the fingerprint last
                for idx,seed_criteria in enumerate(PEATSAorder.seed_criteria):
                    for jdx in range(len(box.fingerprint[seed_criteria[0]])):
                        csv_file_handle.write("S" + str(idx) + "F" + str(jdx) + ",")
                break
            
        # End the header line
        csv_file_handle.write("\n")

        # Loop through all the cases
        for box in PEATSAboxes:
            # Write the constant data
            csv_file_handle.write(box.PEATSApath + ",")
            csv_file_handle.write(box.PEATSAcrust_path + ",")
            csv_file_handle.write(box.mission_name + ",")
            if history==None:
                csv_file_handle.write("N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,")
            else:
                key = box.mission_name.split("_seeded_by")[0].split("_v")[0]
                if key in history["seed"][-1].keys():
                    csv_file_handle.write(str(len(history['seed'][-1][key])) + ",")
                else:
                    csv_file_handle.write("N/A,")
                if key in history["run"][-1].keys():
                    csv_file_handle.write(str(len(history['run'][-1][key])) + ",")
                else:
                    csv_file_handle.write("N/A,")
                numFinished = 0
                numConverged = 0
                if key in history["run"][-1].keys():
                    for instance in history['run'][PEATSAorder.iteration][key]:
                        if instance[0]:
                            numFinished += 1
                        if instance[1]:
                            numConverged += 1
                csv_file_handle.write(str(numFinished) + ",")
                csv_file_handle.write(str(numConverged) + ",")
                csv_file_handle.write(str(history['objective'][-1][key][2]) + ",")
                csv_file_handle.write(str(history['objective'][-1][key][3]) + ",")
                csv_file_handle.write(str(box.iteration) + ",")
                
                seed_code = ""
                if key in history["seed"][-1].keys():
                    for entry in history["seed"][-1][key]:
                        if entry == None:
                            seed_code += "0"
                        else:
                            times_using_seed = 1
                            if PEATSAorder.iteration > 0:
                                for iter_idx in range(2,PEATSAorder.iteration+1):
                                    if key in history["seed"][-iter_idx].keys():
                                        if entry in history["seed"][-iter_idx][key]:
                                            times_using_seed += 1
                            seed_code += str(times_using_seed)
                csv_file_handle.write(seed_code + ",")                     
                                
            for seed_criteria in PEATSAorder.seed_criteria:
                csv_file_handle.write(str(box.seed_criteria[seed_criteria[0]]) + ",")
            csv_file_handle.write(str(box.final_mass) + ",")  
            csv_file_handle.write(str(box.time_of_flight) + ",")  
            csv_file_handle.write(str(box.feasibility) + ",")
            csv_file_handle.write(str(box.objective_value) + ",")
            csv_file_handle.write(str(box.first_nlp_solve_feasible) + ",")
            csv_file_handle.write(str(box.objective_value_first_nlp_solve) + ",")
            csv_file_handle.write(str(box.number_of_solution_attempts) + ",")
            csv_file_handle.write(str(box.solution_attempt_index_that_produced_best_feasible_solution) + ",")
                
            # Print the data for each extra column
            for idx in range(len(PEATSAorder.extra_csv_column_definitions)):
                if len(box.extra_csv_column_data) > idx:
                    csv_file_handle.write(str(box.extra_csv_column_data[idx]).replace(",", "|"))
                csv_file_handle.write(",")
                
            # Write the fingerprint. 
            for seed_criteria in PEATSAorder.seed_criteria:
                for val in box.fingerprint[seed_criteria[0]][:-1]:
                    csv_file_handle.write(str(val).replace(",", "|") + ",")
            
            # End the line
            csv_file_handle.write("\n")
            
        # Close the file    
        csv_file_handle.close()
    
    def copy_results(self,PEATSAorder,boxes):
        """
        Description.
        
        """
        
        import shutil
                
        # Loop through all of the boxes
        for box in boxes:
            
            # copy the options file
            shutil.copyfile(box.PEATSApath + "/" + box.PEATSAdough_path,PEATSAorder.results_directory + box.PEATSAdough_path)
            
            # If the box has a mission object, copy the file
            if box.PEATSAcrust_path != None and box.PEATSAcrust_path != "FAILURE_fake.emtg":
                shutil.copyfile(box.PEATSApath + "/" + box.PEATSAcrust_path,PEATSAorder.results_directory + box.PEATSAcrust_path)
                
            # Update the locaton of the files    
            box.PEATSApath = PEATSAorder.results_directory
            
        # Send the results back    
        return boxes
            
    def write_history_file(self,PEATSAorder,history):
        """
        Writes/appends to the PEATSAhistory.csv file.

        Parameters
        ----------
        PEATSAorder : TYPE
            DESCRIPTION.
        history : TYPE
            DESCRIPTION.

        Returns
        -------
        PEATSAorder : TYPE
            DESCRIPTION.

        """
        import logging
            
        if PEATSAorder.iteration == 0:
            
            handle = open(PEATSAorder.csv_file.rstrip("Iteration1234567890.csv") + "PEATSAhistory.csv","w")
        
            handle.write("Iteration,\n")
            handle.write("Cases Created,\n")
            handle.write("Cases Run,\n")
            handle.write("Ran Unseeded,\n")
            handle.write("Cases Finished,\n")
            handle.write("Iteration Cases Converged,\n")
            handle.write("PEATSA Cases Converged,\n")
            handle.write("First Convergence this Iteration,\n")
            handle.write("Improvement this Iteration,\n")
            handle.write("Number of PEATSA Cases That Converged on First NLP Solve This Iteration,\n")
            handle.write("Average Number of Solution Attempts for All Cases This Iteration,\n")
            handle.write("Average Solution Attempt Number of Best Feasible Solution Among all Feasible Cases This Iteration,\n")
            handle.write("Times Using Seed Before Improvement,")
            
            handle.close()            
            
            PEATSAorder.improvement_tracker = [0] * len(PEATSAorder.constraint_walking)
            
        handle = open(PEATSAorder.csv_file.rstrip("Iteration1234567890.csv") + "PEATSAhistory.csv","r")
    
        current_lines = handle.readlines()
    
        handle.close()       

        old_lines = []
        improvement_count = []
    
        for line in current_lines:
            if line.startswith("Times Using Seed"):
                linesplit = line.split(",")
                for entry in linesplit[:-1]:
                    if "Times" not in entry:
                        improvement_count.append(int(entry))
            else:
                old_lines.append(line.rstrip("\n"))
                
        
        improvement_count.append(0)
        
        case_sum = 0
        unseeded_sum = 0
                
        for key in history['seed'][-1].keys():
            case_sum += len(history['seed'][-1][key])
            unseeded_sum += len([entry for entry in history['seed'][-1][key] if entry == None])
            
        run_sum = 0
        finished_sum = 0
        iteration_converged_sum = 0
        
        for key in history['run'][-1].keys():
            run_sum += len(history['run'][-1][key])
            finished_sum += len([entry for entry in history['run'][-1][key] if entry[0]])
            iteration_converged_sum += len([entry for entry in history['run'][-1][key] if entry[1]])
            
            for entry in history["run"][-1][key]:
                if entry[2]:
                    if entry[3] == None:
                        improvement_count[0] += 1
                    else:
                        times_using_seed = 1
                        # Check if this is the first iteration or not
                        if PEATSAorder.iteration > 0:
                            # It is not. Loop through the iterations
                            for iter_idx in range(2,PEATSAorder.iteration+1):
                                # Keep a flag if we need to break out of this loop
                                break_out = False
                                # Check if this case was run in the previous iteration
                                if key in history["seed"][-iter_idx].keys() and key in history["run"][-iter_idx].keys():
                                    # It was. Check if the seed was used in the previous iteration
                                    if entry[3] in history["seed"][-iter_idx][key]:
                                        # It was. Loop through the run history from that iteration
                                        for prev_entry in history["run"][-iter_idx][key]:
                                            # Check if this run is the one with the relevant seed
                                            if prev_entry[3] == entry[3]:
                                                # It is. Check if the run from the previous generation usign this seed improved
                                                if prev_entry[0] == 0:
                                                    # It didnt improve but it also didnt finish, so dont increment the counter, but dont stop searching backwards
                                                    pass
                                                elif prev_entry[2] == 0:
                                                    # IT did not. Increment the counter
                                                    times_using_seed += 1
                                                else:
                                                    # It did. Stop searching
                                                    break_out = True                                                    
                                                break
                                if break_out:
                                    break
                        improvement_count[times_using_seed] += 1
                
        peatsa_converged_sum = 0
        first_converged_sum = 0
        improvement_sum = 0
                
        for key in history['objective'][-1].keys():
            peatsa_converged_sum += history['objective'][-1][key][0]
            first_converged_sum += history['objective'][-1][key][2]
            improvement_sum += history['objective'][-1][key][3] 
            

        # number of cases that found a feasible solution on the first solution attempt
        firstSolveFeasibleSum = 0
        for key in history['firstSolveFeasible'][-1].keys():
            firstSolveFeasibleSum += history['firstSolveFeasible'][-1][key] # doesn't have [0] at the end because only one element in "tuple"
            
        # total number of solution attempts for all cases (converged or not)
        numSolutionAttemptsSum = 0
        for key in history['solutionAttempts'][-1].keys():
            numSolutionAttemptsSum += history['solutionAttempts'][-1][key] # doesn't have [0] at the end because only one element in "tuple"
        numSolutionAttemptsAvg = numSolutionAttemptsSum / run_sum
        
        # average index of best solution attempt amongst runs that converged
        bestSolutionIndexSum = 0
        for key in history['bestSolutionIndex'][-1].keys():
            bestSolutionIndexSum += history['bestSolutionIndex'][-1][key] # doesn't have [0] at the end because only one element in "tuple"
        
        # prevent divide-by-zero
        if iteration_converged_sum > 0:
            solutionIndexAvg = bestSolutionIndexSum / iteration_converged_sum
        else:
            solutionIndexAvg = -1
        
                        
        new_lines = []        
        
        for line in old_lines: 
            if line.startswith("Iteration,"):   
                new_lines.append(line + str(PEATSAorder.iteration) + ",\n")
            if line.startswith("Cases Created"):   
                new_lines.append(line + str(case_sum) + ",\n")
            if line.startswith("Cases Run"):   
                new_lines.append(line + str(run_sum) + ",\n")
            if line.startswith("Ran Unseeded"):   
                new_lines.append(line + str(unseeded_sum) + ",\n")
            if line.startswith("Cases Finished"):   
                new_lines.append(line + str(finished_sum) + ",\n") 
            if line.startswith("Iteration Cases Converged"):   
                new_lines.append(line + str(iteration_converged_sum) + ",\n")
            if line.startswith("PEATSA Cases Converged"):   
                new_lines.append(line + str(peatsa_converged_sum) + ",\n") 
            if line.startswith("First Convergence this Iteration"):   
                new_lines.append(line + str(first_converged_sum) + ",\n") 
            if line.startswith("Improvement this Iteration"):   
                new_lines.append(line + str(improvement_sum) + ",\n")
            if line.startswith("Number of PEATSA Cases That Converged on First NLP Solve"):
                new_lines.append(line + str(firstSolveFeasibleSum) + ",\n")
            if line.startswith("Average Number of Solution Attempts"):
                new_lines.append(line + str(round(numSolutionAttemptsAvg, 2)) + ",\n")
            if line.startswith("Average Solution Attempt Number"):
                new_lines.append(line + str(round(solutionIndexAvg, 2)) + ",\n")
        
        for idx in range(len(PEATSAorder.constraint_walking)):
            if improvement_sum > 0:
                PEATSAorder.improvement_tracker[idx] = 0
            else:
                if PEATSAorder.improvement_tracker[idx] == PEATSAorder.constraint_walking[idx][1]:
                    PEATSAorder.improvement_tracker[idx] = 1
                else:
                    PEATSAorder.improvement_tracker[idx] += 1                
        
        new_lines.append("Times Using Seed Before Improvement,")
        for entry in improvement_count:
            new_lines[-1] += str(entry) + ","

        handle = open(PEATSAorder.csv_file.rstrip("Iteration1234567890.csv") + "PEATSAhistory.csv","w")
        
        for line in new_lines:
            handle.write(line)
            
        handle.close()
            
        return PEATSAorder