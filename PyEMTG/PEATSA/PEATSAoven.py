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
PEATSAoven.py
========================
This file contains the PEATSAoven class, which is used to execute PEATSA cases.
"""


#PEATSAoven function
#function that runs all of the PEATSAdough (.emtgopt files) in the PEATSA meal
#and turns them into PEATSAcrusts (.emtg files)
#used to be step 2
#Jeremy Knittel 1/23/2018

import time

def executeSingleRunInParallel(PEATSAorder, box, killfileName, killTypeLineNumber, devnull):
    import logging
    f = open(killfileName, 'r')
    killfileLines = f.readlines()
    f.close()
    if "killiteration_nopostprocess" in killfileLines[killTypeLineNumber].lower() or "killrun_nopostprocess" in killfileLines[killTypeLineNumber].lower():
        return
        
    import subprocess

    sub = subprocess.Popen([PEATSAorder.emtg_root_directory + "/" + PEATSAorder.executable_name,PEATSAorder.cases_directory + box.PEATSAdough_path],stdout=devnull,stderr=devnull)
    
    tStart = time.time()
    logging.info('Started ' + PEATSAorder.cases_directory + box.PEATSAdough_path)
    while True:
        #logging.info('Checking whether to exit for ' + PEATSAorder.cases_directory + box.PEATSAdough_path)
        if sub.poll() is not None:
            #logging.info(PEATSAorder.cases_directory + box.PEATSAdough_path + ' ended naturally')
            break
        f = open(killfileName, 'r')
        killfileLines = f.readlines()
        f.close()
        if "killiteration_nopostprocess" in killfileLines[killTypeLineNumber].lower() or "killrun_nopostprocess" in killfileLines[killTypeLineNumber].lower():
            #logging.info(PEATSAorder.cases_directory + box.PEATSAdough_path + ' was killed by killfile')
            sub.kill()
            break
        try:
            sub.wait(5)
        except:
            t = time.time()
            if (t - tStart) > PEATSAorder.killtime:
                #logging.info(PEATSAorder.cases_directory + box.PEATSAdough_path + ' was killed by timeout')
                sub.kill()
                break

    return

def executeSingleRunInParellelWithPebble(PEATSAorder, box):
    """
    Routine to be run in parallel by a Pebble pool to execute an EMTG case

    Parameters
    ----------
    PEATSAorder : TYPE
        DESCRIPTION.
    box : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    import logging
    import os
    logging.info('Started ' + PEATSAorder.cases_directory + box.PEATSAdough_path)
    systemCall = PEATSAorder.emtg_root_directory + "/" + PEATSAorder.executable_name + ' ' + PEATSAorder.cases_directory + box.PEATSAdough_path
    runFlag = os.system(systemCall)
    
    return
    
def pebbleScheduledTaskDone(box, future):
    """
    Callback for pebble-based execution using schedule as opposed to map
    Taken from https://pypi.org/project/Pebble/
    """
    from concurrent.futures import TimeoutError
    import logging
    try:
        result = future.result() # blocks until results are ready
        logging.info(box.PEATSAdough_path + " ended naturally")
    except TimeoutError as error:
        logging.info("Killed " + box.PEATSAdough_path)
    except Exception as error:
        logging.info("PEATSAoven.executeWithPebble encountered an error")
        
    return

class PEATSAoven:
    """
    Used to execute PEATSA cases
    
    Parameters
    ----------
    None.

    Returns
    -------
    None.
    """
    
    def writeKillFile(self, killfileName):
        """
        Write the .DEATH file

        Parameters
        ----------
        killfileName : String
            Directory + name of file in which to write the .DEATH file

        Returns
        -------
        KillTypeLineNumber : Integer
            The line number whose contents are checked to see if we need to kill the iteration.

        """
        
        killfile = open(killfileName,'w')
        KillTypeLineNumber = 0 # increment after every \n. Do NOT increment after the final line where the user actually sets the input.
        killfile.write("# This file is used to end an iteration prematurely or an entire PEATSA run.\n")
        KillTypeLineNumber += 1
        killfile.write("# The contents of this file are checked every 5 seconds while PEATSA is running EMTG cases.\n")
        KillTypeLineNumber += 1
        killfile.write("# The contents of this file are NOT checked while PEATSA is doing python pre/post-processing.\n")
        KillTypeLineNumber += 1
        killfile.write("# To change PEATSA behavior, change the value assigned to KillType and save the file.\n")
        KillTypeLineNumber += 1
        killfile.write("# Valid options for KillType are given below and are case-insensitive.\n")
        KillTypeLineNumber += 1
        killfile.write("# None: (default) No effect. PEATSA continues normally.\n")
        KillTypeLineNumber += 1
        killfile.write("# KillRun_NoPostProcess: Ends the entire PEATSA run. Does NOT run post-processing routines.\n")
        KillTypeLineNumber += 1
        killfile.write("# KillIteration_NoPostProcess: End the current PEATSA iteration, then starts the next PEATSA iteration.\n")
        KillTypeLineNumber += 1
        killfile.write("# Does NOT run post-processing routines. This is useful if there is a midbake option you would like to have take effect immediately.\n")
        KillTypeLineNumber += 1
        killfile.write("KillType = None")
        killfile.close()
        
        return KillTypeLineNumber
    
    def executeWithPebble(self, PEATSAorder, PEATSAboxes):
        """
        Attempt to parallelize with Pebble rather than subprocess or multiprocessing

        Parameters
        ----------
        PEATSAorder : PEATSAmenu object
            The options file for the run as a whole.
        PEATSAboxes : List of PEATSAbox objects
            Options for each individual run.

        Returns
        -------
        None.

        """
        
        import logging
        logging.info("Entering the PEATSA oven to run cases for Iteration " + str(PEATSAorder.iteration))
        
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG/PEATSA")
        
        import pebble # for parallelization
        from pebble import ProcessPool, ProcessExpired, concurrent
        from concurrent.futures import TimeoutError
        import functools # for passing a single-arg function to pebble's pool.map()
        import PEATSAbox
        
        # Write out an iteration kill file
        killfileName = PEATSAorder.working_directory + "PEATSA_iteration_kill_file.DEATH"
        KillTypeLineNumber = self.writeKillFile(killfileName)
        
        devnull = open("/dev/null",'w')
            
        exe = PEATSAorder.emtg_root_directory + "/" + PEATSAorder.executable_name # the EMTG executable
        nCases = len(PEATSAboxes) # number of cases we need to run
        nCasesString = str(nCases)
        nProcessors = PEATSAorder.nCores[0][1] # number of parallel EMTGs to run
        
        # create a version of executeSingleRunInParellelWithPebble
        # that only requires a single input arg (the peatsabox)
        executeSingleRunInParellelWithPebbleSingleArg = functools.partial(executeSingleRunInParellelWithPebble, PEATSAorder)
        runKill = False
        tStart = time.time() # grab the current time
        if PEATSAorder.execution_type == 2:
            # using map
            with ProcessPool(max_workers = nProcessors) as pool:
                future = pool.map(executeSingleRunInParellelWithPebbleSingleArg, PEATSAboxes, timeout = PEATSAorder.killtime)
                
                iterator = future.result()
                
                while True:
                    try:
                        result = next(iterator)
                    except StopIteration:
                        break
                    except TimeoutError as error:
                        logging.info("Killed an EMTG")
                    except ProcessExpired as error:
                        logging.info("An EMTG exited")
                    except Exception as error:
                        logging.info("PEATSAoven.executeWithPebble encountered an error")
        elif PEATSAorder.execution_type == 3:
            # using schedule
            with ProcessPool(max_workers = nProcessors) as pool:
                for i in range(nCases):
                    tCurrent = time.time()
                    if (tCurrent - tStart > 5.0):
                        # check the .death file every 5 seconds
                        # does the break command get us out of this loop so that we don't kick off anymore emtgs?
                        f = open(killfileName, 'r')
                        contents = f.readlines()
                        f.close()
                        if "killiteration_nopostprocess" in contents[KillTypeLineNumber].lower():
                            logging.info("Killing iteration")
                            iterationKill = True
                            allCasesDone = True
                            break
                        elif "killrun_nopostprocess" in contents[KillTypeLineNumber].lower():
                            logging.info("Killing iteration")
                            runKill = True
                            allCasesDone = True
                            break
                        else:
                            tStart = time.time() # reset the timer
                        
                    # schedule a job within the pool
                    future = pool.schedule(executeSingleRunInParellelWithPebbleSingleArg, args=PEATSAboxes[i], timeout = PEATSAorder.killtime)
                    
                    # add a callback to be done when future is done
                    future.add_done_callback(functools.partial(pebbleScheduledTaskDone, PEATSAboxes[i]))
                    
        # what to do if iterationKill or runKill:
        if runKill == True:
            logging.info("Killing PEATSA")
            exit()
        
        return
    
    def executeV2(self, PEATSAorder, PEATSAboxes):
        import logging
        logging.info("Entering the PEATSA oven to run cases for Iteration " + str(PEATSAorder.iteration))
        
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG/PEATSA")
        
        import subprocess
        #import multiprocessing as mp
        #from multiprocessing.pool import ThreadPool
        #import MissionOptions
        import PEATSAbox
        
        # Write out an iteration kill file
        killfileName = PEATSAorder.working_directory + "PEATSA_iteration_kill_file.DEATH"
        KillTypeLineNumber = self.writeKillFile(killfileName)
        
        devnull = open("/dev/null",'w')
            
        exe = PEATSAorder.emtg_root_directory + "/" + PEATSAorder.executable_name # the EMTG executable
        nCases = len(PEATSAboxes) # number of cases we need to run
        nCasesString = str(nCases)
        nProcessors = PEATSAorder.nCores[0][1] # number of parallel EMTGs to run
        nCasesRunning = 0 # keep track of the number of cases currently executing
        mostRecentCase = 0
        iterationKill = False
        runKill = False
        allCasesDone = False
        subs = [] # list of current subprocesses. each element is also a list: [subprocess object, kickoff time, the options file that was executed]
        tStart = time.time() # grab the current time
        while mostRecentCase < nCases:
            # loop through all the cases that we need to run until they have all been kicked off
            fileName = PEATSAboxes[mostRecentCase].PEATSAdough_path
            optionsFileToExecute = PEATSAorder.cases_directory + fileName # grab the next case to kick off
            if nCasesRunning < nProcessors:
                # if we have a free processor, start a new case
                subs.append([subprocess.Popen([exe, optionsFileToExecute], stdout=devnull,stderr=devnull), time.time(), optionsFileToExecute])
                logging.info("Starting case " + str(mostRecentCase + 1) + " of " + nCasesString + " : " + fileName) # use +1 because of 0-indexing
                mostRecentCase += 1 # increment the most recent case
                nCasesRunning += 1 # increment the number of cases currently running
                
            subsToRemove = [] # list of subprocesses to remove because they have ended
            for i in range(len(subs)):
                # figure out which subprocesses to remove
                if subs[i][0].poll() is not None:
                    # this means that the subprocess ended naturally of its own accord
                    #logging.info('Naturally ended: ' + subs[i][2])
                    subsToRemove.append(i)
                elif (time.time() - subs[i][1]) > PEATSAorder.killtime:
                    # manually end process because it reached our killtime
                    subsToRemove.append(i)
                    
            subsToRemove.sort(reverse=True) # sort indices of subs to remove from greatest to least so we remove the highest ones first so that pop() doesn't mess with our ordering
            for subToRemoveIndex in subsToRemove:
                # remove subprocesses from the subs list
                # try removing the ones that ended naturally ... doesn't seem like they should need it, but do it anyway just in case
                try:
                    subs[subToRemoveIndex][0].kill()
                    logging.info('Killed ' + subs[subToRemoveIndex][2])
                except:
                    try:
                        logging.info('Ended without killing: ' + subs[subToRemoveIndex][2])
                    except:
                        logging.info('Ended without kill a process whose name we could not retrieve for some reason.')
                
                try:
                    r = subs.pop(subToRemoveIndex)
                    nCasesRunning -= 1 # decrement number of cases running
                except:
                    logging.info('Tried to pop case index ' + str(subToRemoveIndex) + ' but could not because it was out of range. Continuing.')
                
            if (time.time() - tStart) > 5.0:
                # every 5 seconds, check the killfile
                f = open(killfileName, 'r')
                contents = f.readlines()
                f.close()
                if "killiteration_nopostprocess" in contents[KillTypeLineNumber].lower():
                    logging.info("Killing iteration")
                    iterationKill = True
                    allCasesDone = True
                    # kill all processes
                    for i in range(len(subs)):
                        try:
                            subs[i][0].kill()
                            logging.info('Killed ' + subs[i][2])
                        except:
                            try:
                                logging.info('Ended without killing: ' + subs[i][2])
                            except:
                                logging.info('Ended without kill a process whose name we could not retrieve for some reason.')
                    break
                elif "killrun_nopostprocess" in contents[KillTypeLineNumber].lower():
                    logging.info("Killing iteration")
                    runKill = True
                    allCasesDone = True
                    # kill all processes
                    for i in range(len(subs)):
                        try:
                            subs[i][0].kill()
                            logging.info('Killed ' + subs[i][2])
                        except:
                            try:
                                logging.info('Ended without killing: ' + subs[i][2])
                            except:
                                logging.info('Ended without kill a process whose name we could not retrieve for some reason.')
                    break
                else:
                    tStart = time.time() # reset the timer
                    
                    
        # at this point, all cases have started
        tStart = time.time()
        logging.info("All cases have started")
        while not allCasesDone:
            
            # just keep looping through, figuring out who to remove and/or kill
            subsToRemove = []
            for i in range(len(subs)):
                # figure out which subprocesses to remove
                if subs[i][0].poll() is not None:
                    # this means that the subprocess ended naturally of its own accord
                    #logging.info('Naturally ended?: ' + subs[i][2])
                    subsToRemove.append(i)
                elif (time.time() - subs[i][1]) > PEATSAorder.killtime:
                    # manually end process because it reached our killtime
                    subsToRemove.append(i)
                    
            subsToRemove.sort(reverse=True) # sort indices of subs to remove from greatest to least so we remove the highest ones first so that pop() doesn't mess with our ordering
            for subToRemoveIndex in subsToRemove:
                # remove subprocesses from the subs list
                # try removing the ones that ended naturally ... doesn't seem like they should need it, but do it anyway just in case
                try:
                    subs[subToRemoveIndex][0].kill()
                    logging.info('Killed ' + subs[subToRemoveIndex][2])
                except:
                    try:
                        logging.info('Ended without killing: ' + subs[subToRemoveIndex][2])
                    except:
                        logging.info('Ended without kill a process whose name we could not retrieve for some reason.')
                
                try:
                    r = subs.pop(subToRemoveIndex)
                    nCasesRunning -= 1 # decrement number of cases running
                except:
                    logging.info('Tried to pop case index ' + str(subToRemoveIndex) + ' but could not because it was out of range. Continuing.')
                
            if (time.time() - tStart) > 5.0:
                # every 5 seconds, check the killfile
                f = open(killfileName, 'r')
                contents = f.readlines()
                f.close()
                if "killiteration_nopostprocess" in contents[KillTypeLineNumber].lower():
                    logging.info("Killing iteration")
                    iterationKill = True
                    # kill all processes
                    for i in range(len(subs)):
                        try:
                            subs[i][0].kill()
                            logging.info('Killed ' + subs[i][2])
                        except:
                            try:
                                logging.info('Ended without killing: ' + subs[i][2])
                            except:
                                logging.info('Ended without kill a process whose name we could not retrieve for some reason.')
                    break
                elif "killrun_nopostprocess" in contents[KillTypeLineNumber].lower():
                    logging.info("Killing iteration")
                    runKill = True
                    # kill all processes
                    for i in range(len(subs)):
                        try:
                            subs[i][0].kill()
                            logging.info('Killed ' + subs[i][2])
                        except:
                            try:
                                logging.info('Ended without killing: ' + subs[i][2])
                            except:
                                logging.info('Ended without kill a process whose name we could not retrieve for some reason.')
                    break
                else:
                    tStart = time.time() # reset the timer
            if subs == []:
                # everyone has finished
                allCasesDone = True
                
        # what to do if iterationKill or runKill:
        if runKill == True:
            logging.info("Killing PEATSA")
            exit()
        
        return
    
    def execute(self, PEATSAorder, PEATSAboxes):
        import logging
        logging.info("Entering the PEATSA oven to run cases for Iteration " + str(PEATSAorder.iteration))
        
        
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG/PEATSA")
        
        import subprocess
        import multiprocessing as mp
        from multiprocessing.pool import ThreadPool
        #import MissionOptions
        import PEATSAbox

        # Write out an iteration kill file
        killfileName = PEATSAorder.working_directory + "PEATSA_iteration_kill_file.DEATH"
        killfile = open(killfileName,'w')
        KillTypeLineNumber = 0 # increment after every \n. Do NOT increment after the final line where the user actually sets the input.
        killfile.write("# This file is used to end an iteration prematurely or an entire PEATSA run.\n")
        KillTypeLineNumber += 1
        killfile.write("# The contents of this file are checked every 5 seconds while PEATSA is running EMTG cases.\n")
        KillTypeLineNumber += 1
        killfile.write("# The contents of this file are NOT checked while PEATSA is doing python pre/post-processing.\n")
        KillTypeLineNumber += 1
        killfile.write("# To change PEATSA behavior, change the value assigned to KillType and save the file.\n")
        KillTypeLineNumber += 1
        killfile.write("# Valid options for KillType are given below and are case-insensitive.\n")
        KillTypeLineNumber += 1
        killfile.write("# None: (default) No effect. PEATSA continues normally.\n")
        KillTypeLineNumber += 1
        killfile.write("# KillRun_NoPostProcess: Ends the entire PEATSA run. Does NOT run post-processing routines.\n")
        KillTypeLineNumber += 1
        killfile.write("# KillIteration_NoPostProcess: End the current PEATSA iteration, then starts the next PEATSA iteration.\n")
        KillTypeLineNumber += 1
        killfile.write("# Does NOT run post-processing routines. This is useful if there is a midbake option you would like to have take effect immediately.\n")
        KillTypeLineNumber += 1
        killfile.write("KillType = None")
        killfile.close()
        
        devnull = open("/dev/null",'w')

        logging.info('nCores = ' + str(PEATSAorder.nCores[0][1]))
        threadPool = ThreadPool(PEATSAorder.nCores[0][1])
        #threadPool = ThreadPool(mp.cpu_count())
        for box in PEATSAboxes:
            x = threadPool.apply_async(executeSingleRunInParallel, (PEATSAorder, box, killfileName, KillTypeLineNumber, devnull))
        
        #logging.info('past loop')
        #while True:
        #    if x[-1].ready():
        #        logging.info('broke')
        #        break
        threadPool.close()
        logging.info('closed')
        threadPool.join()
        logging.info('joined')
        
        devnull.close()
        
        return
    
    def pop_and_kill(self,PEATSAorder,proclist,hostCt):
        """
        

        Parameters
        ----------
        PEATSAorder : TYPE
            DESCRIPTION.
        proclist : TYPE
            DESCRIPTION.
        hostCt : TYPE
            DESCRIPTION.

        Raises
        ------
        Exception
            DESCRIPTION.

        Returns
        -------
        proclist : TYPE
            DESCRIPTION.

        """
    
        import time
        import shutil
        import logging
        import os
    
        # Create a list of cases that need to be removed from the process list
        # ether because they finished successfully, or because they need to be killed
        poplist = []
    
        # Create a list of cases that need to be killed
        killlist = []
    
        # Loop through the process list
        for i in range(0,len(proclist)):
        
            # Grab the current process
            proc = proclist[i]
            # proc = [subprocess object, start time, case string, host string, host idx]
        
            if proc[0].poll() == None:
                # The case is still running
            
                # Check if it has been running too long
                if time.time() - proc[1] > PEATSAorder.killtime:
                    # It has. IT MUST DIE
                    logging.info("Killing: " + proc[2])
                
                    # The case needs to be added to the list of cases to remove from the process list and to the kill list
                    poplist.append(i)
                    killlist.append(i)
                else:
                    # It has been running for an appropriate period of time. Keep going
                    continue
            elif proc[0].poll() == 0:
                # The case has finished successfully
            
                logging.info("Case Finished: " + proc[2])
            
                # Add this to the list of cases to remove from the process list
                poplist.append(i)
            else:
                logging.info("Start time: " + str(proc[1]))
                logging.info("Emtgopt file: " + proc[2])
                logging.info("Host: " + proc[3])
                poplist.append(i)
                # raise Exception("Unexpected error. The process is not running, but it didnt end gracefully.")
    
        # Loop through all the cases that need to be killed
        # and MURDER THEM
        for popListIdx in range(len(killlist)):
            popIdx = killlist[popListIdx]
        
            # Check if the process is running locally or on another machine
            if proclist[popIdx][3] == "local":
            
                # Killing processes is a bit messy, so its safest to do so in a try block
                try:
                    # Python subprocess objects have a kill method. use it.
                    proclist[popIdx][0].kill()
                except:
                    # If we are here, also log the exception
                    logging.info(sys.exc_info()[0])
                    logging.info("Case appears to have finished now, so I cant murder it: " + proclist[popIdx][2])
            else:
                # IF we are not running locally, then the subprocess is ssh, not an emtg process. So we cant kill it directly, we need to 
                # load the process ID, which should have been written to file on the remote server, and then send a system command 
                # to kill it. NOTE: This is not OS safe, but you should not be here if you are just debugging on windows. Dont run
                # remote cases from a windows machine! 
            
                # Build the results folder location from 
                folder_name = PEATSAorder.results_dir + proclist[popIdx][2].split("/")[-1].rstrip("emtgopt").rstrip(".")
                # Open the pid file from the local machine. NOTE: THIS DOESNT WORK IF THERE ISNT SHARED MEMORY!
                pid_file = open(folder_name + "/PID.txt")
                # Read the file to get the process id number
                PIDno = pid_file.readline()
                # Send out a signal to the remote machine to kill the process
                os.system("ssh " + PEATSAorder.nCores[proclist[popIdx][3]][0] + " kill " + PIDno)
        
        # Loop backwards through the cases that we need to remove from the active process list. 
        # I go backwards because then the indexes that we added earlier stay valid as we go. Yes, I admit, 
        # there is almost certainly a better way to do this, but it does work just fine, and it isn't terribly
        # complicated
        for popListIdx in range(1,len(poplist)+1):
        
            # Get the index into the process list from the poplist
            popIdx = poplist[-popListIdx]
        
            # Reduce the number of cases running on this host
            hostCt[proclist[popIdx][3]] -= 1
        
            # Remove the process from the list
            proclist.pop(popIdx)
    
        # Give back the trimmed proclist
        return proclist

    def cook(self,PEATSAorder,PEATSAboxes):
        """
        

        Parameters
        ----------
        PEATSAorder : TYPE
            DESCRIPTION.
        PEATSAboxes : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG/PEATSA")
        import MissionOptions
        import PEATSAbox
        import logging
        import subprocess as sp
        import time
    
        logging.info("Entering the PEATSA oven to run cases for Iteration " + str(PEATSAorder.iteration))
    
        # Write out an iteration kill file
        killfile = open(PEATSAorder.working_directory + "PEATSA_iteration_kill_file.DEATH",'w')
        KillTypeLineNumber = 0 # increment after every \n. Do NOT increment after the final line where the user actually sets the input.
        killfile.write("# This file is used to end an iteration prematurely or an entire PEATSA run.\n")
        KillTypeLineNumber += 1
        killfile.write("# The contents of this file are checked every 5 seconds while PEATSA is running EMTG cases.\n")
        KillTypeLineNumber += 1
        killfile.write("# The contents of this file are NOT checked while PEATSA is doing python pre/post-processing.\n")
        KillTypeLineNumber += 1
        killfile.write("# To change PEATSA behavior, change the value assigned to KillType and save the file.\n")
        KillTypeLineNumber += 1
        killfile.write("# Valid options for KillType are given below and are case-insensitive.\n")
        KillTypeLineNumber += 1
        killfile.write("# None: (default) No effect. PEATSA continues normally.\n")
        KillTypeLineNumber += 1
        killfile.write("# KillRun_NoPostProcess: Ends the entire PEATSA run. Does NOT run post-processing routines.\n")
        KillTypeLineNumber += 1
        killfile.write("# KillIteration_NoPostProcess: End the current PEATSA iteration, then starts the next PEATSA iteration.\n")
        KillTypeLineNumber += 1
        killfile.write("# Does NOT run post-processing routines. This is useful if there is a midbake option you would like to have take effect immediately.\n")
        KillTypeLineNumber += 1
        killfile.write("KillType = None")
        killfile.close()
        
        ctr = 0
    
        # Create a list of all active processes
        proclist = []
    
        # Create a list to count how many cases a given host is running
        hostCt = {}
    
        # Count the total number of cores available
        totCores = 0
        for host in PEATSAorder.nCores:
            hostCt.update({host[0] : 0})
            totCores += host[1]
    
        # emtg output is going to go to devnull, so open a write on it
        devnull = open("/dev/null",'w')
        
        # Flag for killing prematurely
        kill_iteration = False
        
        # Flag for killing the entire peatsa run
        kill_peatsa = False
        
        nPeatsaBoxes = len(PEATSAboxes)
        nPeatsaBoxesString = str(nPeatsaBoxes)
    
        # Loop through each case
        for box in PEATSAboxes:
            PEATSAdough_path = box.PEATSAdough_path
            
            # increment the counter
            ctr += 1
           
            # Check if the number of cases running is greater than or equal to the number
            # of cores available     
            while len(proclist) >= totCores:
            
                # We do not currently have any cores open
            
                logging.info(str(len(proclist)) + " cases running")
            
                # Sleep for a few seconds so that we arent constantly actively looping
                time.sleep(5)
                
                # open the kill file
                killfile_handle = open(PEATSAorder.working_directory + "PEATSA_iteration_kill_file.DEATH",'r')
                killfileLines = killfile_handle.readlines()
                killfile_handle.close()
                                
                # get the lines from the kill file
                # Note: Since the default text contains the word 'kill' in it, the following check will not kill the iteration if the word 'Replace' is also contained in the file.
                if "killiteration_nopostprocess" in killfileLines[KillTypeLineNumber].lower():
                    logging.info("User typed 'killiteration_nopostprocess' in .DEATH file. Ending iteration.")
                    kill_iteration = True
                    break
                
                elif "killrun_nopostprocess" in killfileLines[KillTypeLineNumber].lower():
                    logging.info("User typed 'killrun_nopostprocess' in .DEATH file. Quitting PEATSA.")
                    kill_iteration = True
                    kill_peatsa = True
                    break
                
                            
                # Check if any cases finished and if any cases need to be killed
                proclist = self.pop_and_kill(PEATSAorder,proclist,hostCt)
            
            if kill_iteration:
                break
            
            # Loop through all of the hosts available
            for host in PEATSAorder.nCores:
            
                # Check if this host has idle cores
                if hostCt[host[0]] < host[1]:    
                    # It does!                    
                
                    logging.info("Starting case " + str(ctr) + " of " + nPeatsaBoxesString + " [" + host[0] +"]: " + PEATSAdough_path)
                
                    # Add one to the number of cases being run on this host
                    hostCt[host[0]] += 1
                
                    # Starting the case is different depending on if its local or not
                    if host[0] == "local":
                        # It is local. Just start the EMTG process, and add it along with its details to the process list
                        proclist.append([sp.Popen([PEATSAorder.emtg_root_directory + "/" + PEATSAorder.executable_name,PEATSAorder.cases_directory + PEATSAdough_path],stdout=devnull,stderr=devnull),time.time(),PEATSAdough_path, host[0]])
                    else:
                        # It is remote. Start the EMTG process on the remote machine by way of ssh, and add it along with its details to the process list
                        # This is not OS safe, but you shouldnt be kicking off runs remotely from windows.
                        proclist.append([sp.Popen(["ssh",PEATSAorder.nCores[idx][0],PEATSAorder.emtg_root_directory + "/" + PEATSAorder.executable_name,case],stdout=devnull,stderr=devnull),time.time(),PEATSAdough_path,host[0]])
                
                    # This case has been started, so no need to keep looping through hosts
                    break
        
        
        if kill_iteration:
            logging.info("Killing iteration")
            # Loop backwards through the cases that we need to remove from the active process list. 
            # I go backwards because then the indexes that we added earlier stay valid as we go. Yes, I admit, 
            # there is almost certainly a better way to do this, but it does work just fine, and it isn't terribly
            # complicated
            for proc in proclist:
                        
                # Check if the process is running locally or on another machine
                if proc[3] == "local":
            
                    # Killing processes is a bit messy, so its safest to do so in a try block
                    try:
                        # Python subprocess objects have a kill method. use it.
                        proc[0].kill()
                    except:
                        # If we are here, also log the exception
                        logging.info(sys.exc_info()[0])
                        logging.info("Case appears to have finished now, so I cant murder it: " + proc[2])
                else:
                    # IF we are not running locally, then the subprocess is ssh, not an emtg process. So we cant kill it directly, we need to 
                    # load the process ID, which should have been written to file on the remote server, and then send a system command 
                    # to kill it. NOTE: This is not OS safe, but you should not be here if you are just debugging on windows. Dont run
                    # remote cases from a windows machine! 
            
                    # Build the results folder location from 
                    folder_name = PEATSAorder.results_dir + proc[2].split("/")[-1].rstrip("emtgopt").rstrip(".")
                    # Open the pid file from the local machine. NOTE: THIS DOESNT WORK IF THERE ISNT SHARED MEMORY!
                    pid_file = open(folder_name + "/PID.txt")
                    # Read the file to get the process id number
                    PIDno = pid_file.readline()
                    # Send out a signal to the remote machine to kill the process
                    os.system("ssh " + PEATSAorder.nCores[proc[3]][0] + " kill " + PIDno)
            
            if kill_peatsa:
                exit()
        else:
                    
            logging.info("All cases started")
                        
            while len(proclist):    
                
                # We have kicked off all cases. Now we just need to wait
    
                logging.info(str(len(proclist)) + " cases running")
    
                # Sleep for a few seconds so that we arent constantly actively looping
                time.sleep(5)
                
                # open the kill file
                killfile_handle = open(PEATSAorder.working_directory + "PEATSA_iteration_kill_file.DEATH",'r')
                killfileLines = killfile_handle.readlines()
                killfile_handle.close()
                                
                # get the lines from the kill file
                # Note: Since the default text contains the word 'kill' in it, the following check will not kill the iteration if the word 'Replace' is also contained in the file.
                if "killiteration_nopostprocess" in killfileLines[KillTypeLineNumber].lower():
                    logging.info("User typed 'killiteration_nopostprocess' in .DEATH file. Ending iteration.")
                    kill_iteration = True
                    #break
                
                elif "killrun_nopostprocess" in killfileLines[KillTypeLineNumber].lower():
                    logging.info("User typed 'killrun_nopostprocess' in .DEATH file. Quitting PEATSA.")
                    kill_iteration = True
                    kill_peatsa = True
                    #break            
                    
                if kill_iteration:
                    logging.info("Killing iteration")
                    # Loop backwards through the cases that we need to remove from the active process list. 
                    # I go backwards because then the indexes that we added earlier stay valid as we go. Yes, I admit, 
                    # there is almost certainly a better way to do this, but it does work just fine, and it isn't terribly
                    # complicated
                    for proc in proclist:
                        
                        # Check if the process is running locally or on another machine
                        if proc[3] == "local":
            
                            # Killing processes is a bit messy, so its safest to do so in a try block
                            try:
                                # Python subprocess objects have a kill method. use it.
                                proc[0].kill()
                            except:
                                # If we are here, also log the exception
                                logging.info(sys.exc_info()[0])
                                logging.info("Case appears to have finished now, so I cant murder it: " + proc[2])
                        else:
                            # IF we are not running locally, then the subprocess is ssh, not an emtg process. So we cant kill it directly, we need to 
                            # load the process ID, which should have been written to file on the remote server, and then send a system command 
                            # to kill it. NOTE: This is not OS safe, but you should not be here if you are just debugging on windows. Dont run
                            # remote cases from a windows machine! 
            
                            # Build the results folder location from 
                            folder_name = PEATSAorder.results_dir + proc[2].split("/")[-1].rstrip("emtgopt").rstrip(".")
                            # Open the pid file from the local machine. NOTE: THIS DOESNT WORK IF THERE ISNT SHARED MEMORY!
                            pid_file = open(folder_name + "/PID.txt")
                            # Read the file to get the process id number
                            PIDno = pid_file.readline()
                            # Send out a signal to the remote machine to kill the process
                            os.system("ssh " + PEATSAorder.nCores[proc[3]][0] + " kill " + PIDno)
                    
                    if kill_peatsa:
                        exit()
                    break
                else:
                    # Check if any cases finished and if any cases need to be killed
                    proclist = self.pop_and_kill(PEATSAorder,proclist,hostCt)
