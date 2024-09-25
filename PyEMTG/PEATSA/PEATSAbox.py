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
PEATSAbox.py
==================

Class that is a sortable container for a mission and its options script

"""

class PEATSAbox(object):
    """
    A sortable container for a mission and its options script
    
    Parameters
    ----------
    PEATSAorder : PEATSAmenu object
        Holds options for the PEATSA run.
    PEATSAfolder : String
        Path to the output directory for this PEATSA run. Does NOT need to end with '/'.
    PEATSAdough_path : String, optional
        Path to .emtgopt file for this case relative to PEATSApath. The default is None.
    PEATSAcrust_path : String, optional
        Path to .emtg file for this case relative to PEATSApath. The default is None.

    Returns
    -------
    None.
    
    """
    def __init__(self, PEATSAorder, PEATSAfolder, PEATSAdough_path=None, PEATSAcrust_path=None, initializationType = 0):
        """
        Constructor for PEATSAbox

        Parameters
        ----------
        PEATSAorder : PEATSAmenu object
            Holds options for the PEATSA run.
        PEATSAfolder : String
            Path to the output directory for this PEATSA run. Does NOT need to end with '/'.
        PEATSAdough_path : String, optional
            Path to .emtgopt file for this case relative to PEATSApath. The default is None.
        PEATSAcrust_path : String, optional
            Path to .emtg file for this case relative to PEATSApath. The default is None.
        initializationType : Integer, optional
            0: Old way: all initialization at once
            1: New way (for parallelization). Only reading of files happens on construction. Rest of initialization requires calling initializeAfterReading()

        Returns
        -------
        None.

        """
        self.PEATSApath       = PEATSAfolder.rstrip("/") + "/"
        """Path to the output directory for this PEATSA run"""
        self.PEATSAcrust_path = PEATSAcrust_path
        """Path to .emtg file for this case relative to PEATSApath"""
        self.PEATSAdough_path = PEATSAdough_path
        self.PEATSAdough_fullpath = self.PEATSApath + self.PEATSAdough_path
        self.PEATSAdough_fullpath = self.PEATSAdough_fullpath.replace("results", "cases")
        """Path to .emtgopt file for this case relative to PEATSApath"""
        self.seed_boxes = []
        """Description."""
        self.fingerprint = {}
        """Fingerprint dictionary. Gets updated in initialize() based on PEATSAorder.seed_criteria."""
        self.converged = 0
        """If this cases produced a non-failed solution, gets set to 1 in initialize(). Otherwise, 0."""
        self.seed_criteria = {}
        self.mission_name = ""
        self.final_mass = -1e10
        """Final mass of spacecraft for this case (kg). Set in initialize() iff the case did not fail."""
        self.time_of_flight = 1e10
        """Time of flight for mission for this case (years). Set in initialize() iff the case did not fail."""
        self.extra_csv_column_data = []
        self.plot_fingerprint = []
        self.in_study = True
        self.feasibility = 1e10
        self.iteration = 0
        self.used_seed = None
        self.re_run_using_this_seed = {}
        self.objective_value = 1e100
        self.successfullyReadMissionAndMissionOptions = False
        self.first_nlp_solve_feasible = 0
        self.objective_value_first_nlp_solve = 1e100
        self.number_of_solution_attempts = 0
        self.solution_attempt_index_that_produced_best_feasible_solution = 0

        # Check if this is a minimzation or maximization problem and set the 
        # default obective value accordingly
        if PEATSAorder.max_or_min == "min":
            self.objective_value = 1e100
        elif PEATSAorder.max_or_min == "max":
            self.objective_value = -1e100
        
        if initializationType == 0:
            self.initialize(PEATSAorder)
        elif initializationType == 1:
            self.initializeReadOnly(PEATSAorder)
        return
    
    def initializeReadOnly(self, PEATSAorder):
        """
        A version of initialize() that only does the stuff that requires reading from file
        because that stuff is not a good candidate for parallelization, while the other
        stuff is.
        """
        
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        try:
            import PyHardware
            use_pyhardware = True
        except:
            use_pyhardware = False
        import MissionOptions
        import Mission
        import ThrottleTable
        import math
        import numpy as np
        import logging
        
        # Check if we have a mission file
        try:
            #MO = MissionOptions.MissionOptions(self.PEATSApath + self.PEATSAdough_path)
            MO = MissionOptions.MissionOptions(self.PEATSAdough_fullpath)
            # Write the mission options object to a class variable
            self.PEATSAdough = MO
        except:
            #logging.info('WARNING: Unable to load mission options file: ' + self.PEATSApath + self.PEATSAdough_path)
            logging.info('WARNING: Unable to load mission options file: ' + self.PEATSAdough_fullpath)
            return
        if self.PEATSAcrust_path != None and self.PEATSAcrust_path != "" and 'fake' not in self.PEATSAcrust_path:
            # Sometimes things go wrong loading the data, so put this in a try block
            try:
                # Load the mission object
                M = Mission.Mission(self.PEATSApath + self.PEATSAcrust_path)

                self.successfullyReadMissionAndMissionOptions = True

                # # Even if it didnt converge, we should be able to get feasiblity
                # self.feasibility = M.worst_violation
                
                # Write the Mission object to the box
                #if not PEATSAorder.thin_crust:
                self.PEATSAcrust = M
                
                # # If the emtg case converged, then it will not have "FAILURE" in the name
                # # and it will be safe to load the file into an EMTG mission object
                # # (Jacob) every once in a while Jacob writes MBH code wrong and "FAILURE" does not appear, so let's also check feasibility ourselves!
                # if "FAILURE" not in self.PEATSAcrust_path and self.feasibility < MO.snopt_feasibility_tolerance:
                #     # Get the final mass
                #     self.final_mass = M.spacecraft_dry_mass
        
                #     # Get the time of flight
                #     self.time_of_flight = M.total_flight_time_years
        
                #     # Update the flag to note that it converged
                #     self.converged = 1

            except:
                # Something is wrong with this file. Consider it a failure
                self.successfullyReadMissionAndMissionOptions = False
                self.converged = 0
                logging.info("WARNING: .emtg file exists but threw an error while being loaded: " + self.PEATSApath + self.PEATSAcrust_path)
        
        # Load the mission options object
        #MO = MissionOptions.MissionOptions(self.PEATSApath + self.PEATSAdough_path)
            
        
        
        # Load the spacecraft options if a spacecraft options file exists
        if MO.SpacecraftModelInput == 1 and use_pyhardware:
            self.PEATSAstone = PyHardware.SpacecraftOptions(MO.HardwarePath + "/" + MO.SpacecraftOptionsFile)
        
            # # Initialize an array to hold the throttle table files
            # self.PEATSAcooking_temperatures = []
            
            # # Loop through all of the stages to get all of the throttle table files
            # for stage_idx in range(self.PEATSAstone.getNumberOfStages()):
                
            #     # Grab the stage options
            #     stage = self.PEATSAstone.getStageOptions(stage_idx)
                
            #     # Grab the EP options
            #     EPopts = stage.getElectricPropulsionSystemOptions()
                
            #     # Append the throttle table files
            #     self.PEATSAcooking_temperatures.append(ThrottleTable.ThrottleTable(EPopts.getThrottleTableFile()))
        
        # Load the launch vehicle options, if a lv library exists
        if MO.SpacecraftModelInput != 2 and use_pyhardware:
            self.PEATSAovenmitt = PyHardware.LaunchVehicleOptionsLibrary(MO.HardwarePath + "/" + MO.LaunchVehicleLibraryFile)
                
        # # Get the mission name
        # self.mission_name = self.PEATSAdough.mission_name
        
        # # Note what iteration it was introducted
        # self.iteration = PEATSAorder.iteration
        
        # # Evaluate the objective.
        # if self.PEATSAcrust_path != None and self.PEATSAcrust_path != "" and "FAILURE" not in self.PEATSAcrust_path and self.feasibility < MO.snopt_feasibility_tolerance:
        #     try:
        #         M = Mission.Mission(self.PEATSApath + self.PEATSAcrust_path)
        #         # self.PEATSAcrust = M
        #         self.objective_value = eval(PEATSAorder.objective_formula)
        #     except:
        #         if 'M' in locals():
        #             raise Exception("The mission object exists, but the objective formula wasnt able to be evaluated")
        
        return

    def initializeAfterReading(self, PEATSAorder):
        """
        Do the parallelizable initialization stuff. Called only AFTER initializeReadOnly()
        """
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        try:
            import PyHardware
            use_pyhardware = True
        except:
            use_pyhardware = False
        import MissionOptions
        import Mission
        import ThrottleTable
        import math
        import numpy as np
        import logging
        
        if self.successfullyReadMissionAndMissionOptions:
            # Even if it didnt converge, we should be able to get feasiblity
            MO = self.PEATSAdough
            M = self.PEATSAcrust
            self.feasibility = M.worst_violation
            self.first_nlp_solve_feasible = M.first_nlp_solve_feasible
            self.objective_value_first_nlp_solve = M.objective_value_first_nlp_solve
            self.number_of_solution_attempts = M.number_of_solution_attempts
            self.solution_attempt_index_that_produced_best_feasible_solution = M.solution_attempt_index_that_produced_best_feasible_solution
            
            # If the emtg case converged, then it will not have "FAILURE" in the name
            # and it will be safe to load the file into an EMTG mission object
            # (Jacob) every once in a while Jacob writes MBH code wrong and "FAILURE" does not appear, so let's also check feasibility ourselves!
            if "FAILURE" not in self.PEATSAcrust_path and self.feasibility < MO.snopt_feasibility_tolerance:
                # Get the final mass
                self.final_mass = M.spacecraft_dry_mass
    
                # Get the time of flight
                self.time_of_flight = M.total_flight_time_years
    
                # Update the flag to note that it converged
                self.converged = 1
                
        # Load the spacecraft options if a spacecraft options file exists
        # Initialize an array to hold the throttle table files
        self.PEATSAcooking_temperatures = []
        
        # Loop through all of the stages to get all of the throttle table files
        for stage_idx in range(self.PEATSAstone.getNumberOfStages()):
            # Grab the stage options
            stage = self.PEATSAstone.getStageOptions(stage_idx)
            
            # Grab the EP options
            EPopts = stage.getElectricPropulsionSystemOptions()
            
            # Append the throttle table files
            self.PEATSAcooking_temperatures.append(ThrottleTable.ThrottleTable(EPopts.getThrottleTableFile()))
            
        # Get the mission name
        self.mission_name = self.PEATSAdough.mission_name
        
        # Note what iteration it was introducted
        self.iteration = PEATSAorder.iteration
        
        # Evaluate the objective.
        if self.PEATSAcrust_path != None and self.PEATSAcrust_path != "" and "FAILURE" not in self.PEATSAcrust_path and self.feasibility < MO.snopt_feasibility_tolerance:
            try:
                #M = Mission.Mission(self.PEATSApath + self.PEATSAcrust_path)
                # self.PEATSAcrust = M
                #M = self.PEATSAcrust
                self.objective_value = eval(PEATSAorder.objective_formula)
            except:
                # don't do anything ... if we made it 
                raise Exception("The mission object exists, but the objective formula wasnt able to be evaluated")
                
        # Loop through the different seed crieteria
        for seed_criteria in PEATSAorder.seed_criteria:
            
            # Add a new dictionary entry in the box's seed criteria by
            # evaluating the seed criteria. This is not in a try block
            # because it needs to work
            self.seed_criteria.update({seed_criteria[0] : eval(seed_criteria[0])})
            
            # Initialize a boolean for each seed
            self.re_run_using_this_seed.update({seed_criteria[0] : False})
            
            # Initialize the fingerprint for this seed criteria in the dictionary   
            # to an empty list
            self.fingerprint.update({seed_criteria[0]: () })    
                
            # Get the fingerprint
            for idx,formula in enumerate(PEATSAorder.fingerprint):
                # Add this formula. This is not in a try block, because
                # if PEATSA can't get the fingerprint for all of its cases
                # PEATSA cant function!
                if formula != seed_criteria[0]:
                    self.fingerprint[seed_criteria[0]] += (eval(formula),)
        
            # Add the seed criteria to the fingerprint
            self.fingerprint[seed_criteria[0]] += (self.seed_criteria[seed_criteria[0]],)
                
        # Get the extra user requested csv data
        for column in PEATSAorder.extra_csv_column_definitions:
            # This does not need to succeed, for PEATSA to work, so it will go in 
            # a try block
            try:
                # Try evaluating the formula
                self.extra_csv_column_data.append(eval(column[1]))
            except:
                # It didnt work. Oh well, just throw an empty string in the list
                self.extra_csv_column_data.append("")
                
                # from IPython import embed
                # embed()
                # stop
        
                # If we are debugging however, we want to know that the formula failed.
                if PEATSAorder.if_run_cases == 0:
                    if ('M.' in column[1] and "M" in locals()) or "M." not in column[1]:
                        raise Exception("Error in user csv column formula: " + column[1])
        
        return
    
    def initialize(self, PEATSAorder):
        """
        Description.

        Parameters
        ----------
        PEATSAorder : TYPE
            DESCRIPTION.

        Raises
        ------
        Exception
            DESCRIPTION.

        Returns
        -------
        None.

        """
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        try:
            import PyHardware
            use_pyhardware = True
        except:
            use_pyhardware = False
        import MissionOptions
        import Mission
        import ThrottleTable
        import math
        import numpy as np
        import logging
        
        # Check if we have a mission file
        if self.PEATSAcrust_path != None and self.PEATSAcrust_path != "" and 'fake' not in self.PEATSAcrust_path:
            # Sometimes things go wrong loading the data, so put this in a try block
            try:
                # Load the mission object
                M = Mission.Mission(self.PEATSApath + self.PEATSAcrust_path)
                #MO = MissionOptions.MissionOptions(self.PEATSApath + self.PEATSAdough_path)
                MO = MissionOptions.MissionOptions(self.PEATSAdough_fullpath)

                # Even if it didnt converge, we should be able to get hef easiblity
                self.feasibility = M.worst_violation
                self.first_nlp_solve_feasible = M.first_nlp_solve_feasible
                self.objective_value_first_nlp_solve = M.objective_value_first_nlp_solve
                self.number_of_solution_attempts = M.number_of_solution_attempts
                self.solution_attempt_index_that_produced_best_feasible_solution = M.solution_attempt_index_that_produced_best_feasible_solution
                
                # Write the Mission object to the box
                if not PEATSAorder.thin_crust:
                    self.PEATSAcrust = M
                
                # If the emtg case converged, then it will not have "FAILURE" in the name
                # and it will be safe to load the file into an EMTG mission object
                # (Jacob) every once in a while Jacob writes MBH code wrong and "FAILURE" does not appear, so let's also check feasibility ourselves!
                if "FAILURE" not in self.PEATSAcrust_path and self.feasibility < MO.snopt_feasibility_tolerance:
                    # Get the final mass
                    self.final_mass = M.spacecraft_dry_mass
        
                    # Get the time of flight
                    self.time_of_flight = M.total_flight_time_years
        
                    # Update the flag to note that it converged
                    self.converged = 1

            except:
                # Something is wrong with this file. Consider it a failure
                self.converged = 0
                logging.info("WARNING: .emtg file exists but threw an error while being loaded: " + self.PEATSApath + self.PEATSAcrust_path)
        
        # Load the mission options object
        #MO = MissionOptions.MissionOptions(self.PEATSApath + self.PEATSAdough_path)
        MO = MissionOptions.MissionOptions(self.PEATSAdough_fullpath)
            
        # Write the mission and mission object to a class variable
        self.PEATSAdough = MO
        
        # Load the spacecraft options if a spacecraft options file exists
        if MO.SpacecraftModelInput == 1 and use_pyhardware:
            self.PEATSAstone = PyHardware.SpacecraftOptions(MO.HardwarePath + "/" + MO.SpacecraftOptionsFile)
        
            # Initialize an array to hold the throttle table files
            self.PEATSAcooking_temperatures = []
            
            # Loop through all of the stages to get all of the throttle table files
            for stage_idx in range(self.PEATSAstone.getNumberOfStages()):
                
                # Grab the stage options
                stage = self.PEATSAstone.getStageOptions(stage_idx)
                
                # Grab the EP options
                EPopts = stage.getElectricPropulsionSystemOptions()
                
                # Append the throttle table files
                self.PEATSAcooking_temperatures.append(ThrottleTable.ThrottleTable(EPopts.getThrottleTableFile()))
        
        # Load the launch vehicle options, if a lv library exists
        if MO.SpacecraftModelInput != 2 and use_pyhardware:
            self.PEATSAovenmitt = PyHardware.LaunchVehicleOptionsLibrary(MO.HardwarePath + "/" + MO.LaunchVehicleLibraryFile)
                
        # Get the mission name
        self.mission_name = self.PEATSAdough.mission_name
        
        # Note what iteration it was introducted
        self.iteration = PEATSAorder.iteration
        
        # Evaluate the objective.
        if self.PEATSAcrust_path != None and self.PEATSAcrust_path != "" and "FAILURE" not in self.PEATSAcrust_path and self.feasibility < MO.snopt_feasibility_tolerance:
            try:
                M = Mission.Mission(self.PEATSApath + self.PEATSAcrust_path)
                # self.PEATSAcrust = M
                self.objective_value = eval(PEATSAorder.objective_formula)
            except:
                if 'M' in locals():
                    raise Exception("The mission object exists, but the objective formula wasnt able to be evaluated")

        # Loop through the different seed crieteria
        for seed_criteria in PEATSAorder.seed_criteria:
            
            # Add a new dictionary entry in the box's seed criteria by
            # evaluating the seed criteria. This is not in a try block
            # because it needs to work
            self.seed_criteria.update({seed_criteria[0] : eval(seed_criteria[0])})
            
            # Initialize a boolean for each seed
            self.re_run_using_this_seed.update({seed_criteria[0] : False})
            
            # Initialize the fingerprint for this seed criteria in the dictionary   
            # to an empty list
            self.fingerprint.update({seed_criteria[0]: () })    
                
            # Get the fingerprint
            for idx,formula in enumerate(PEATSAorder.fingerprint):
                # Add this formula. This is not in a try block, because
                # if PEATSA can't get the fingerprint for all of its cases
                # PEATSA cant function!
                if formula != seed_criteria[0]:
                    self.fingerprint[seed_criteria[0]] += (eval(formula),)
        
            # Add the seed criteria to the fingerprint
            self.fingerprint[seed_criteria[0]] += (self.seed_criteria[seed_criteria[0]],)
                
        # Get the extra user requested csv data
        for column in PEATSAorder.extra_csv_column_definitions:
            # This does not need to succeed, for PEATSA to work, so it will go in 
            # a try block
            try:
                # Try evaluating the formula
                self.extra_csv_column_data.append(eval(column[1]))
            except:
                # It didnt work. Oh well, just throw an empty string in the list
                self.extra_csv_column_data.append("")
                
                # from IPython import embed
                # embed()
                # stop
        
                # If we are debugging however, we want to know that the formula failed.
                if PEATSAorder.if_run_cases == 0:
                    if ('M.' in column[1] and "M" in locals()) or "M." not in column[1]:
                        raise Exception("Error in user csv column formula: " + column[1])
        
        
