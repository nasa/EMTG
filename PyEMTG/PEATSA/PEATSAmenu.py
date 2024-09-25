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
PEATSAmenu.py
==================

File holds the PEATSAmenu class AND the PEATSAinstagram class.

PEATSAmenu holds the options for the PEATSArun.

PEATSAinstagram holds plotting options.

"""

#PEATSAmenu class
#class that holds the PEATSAtoppings (options) and sets default values for those options
#for a PEATSAmeal. A PEATSA order overrides many of these default PEATSAtoppings.
#Jeremy Knittel 1-23-2018

# Define a class for plotting options
class PEATSAinstagram(object):
    """
    Holds plotting options
    
    Parameters
    ----------
    options_file : String, optional
        If provided, plotting options are loaded from the provided file. The default is None.

    Returns
    -------
    None.
    
    """
    def __init__(self,options_file = None):
        """
        Constructor
        """
        self.xlabel = ""
        self.ylabel = ""
        self.title = ""
        self.ylim = [0,0]
        self.xlim = [0,0]
        self.file_name_str = ""
        self.x_is_date = 0
        self.y_is_date = 0
        self.plot_commands = []
        self.lines_or_dots = "dots"
        self.x_variable = ""
        self.y_variable = ""
        self.groupings = []
        self.legend = False
        # A list of lists. Each sub-list must be teh same length as self.groupings. 
        # Each entry in the sub-list should be either a value of that grouping fingerprint that you want included in the plot
        # or "NA". Data must match only one of the sublists, so it uses or logic
        self.grouping_requirements = [] 
        
        if options_file != None:
            self.load_options_from_file(options_file)
        
    def write_options(self,filename):
        """
        Writes plotting options saved in object's attributes to a file.

        Parameters
        ----------
        filename : String
            Full path and file name to which plotting options are to be written.

        Returns
        -------
        None.

        """
        
        # Open the output file
        file_handle = open(filename,'w')
    
        # Loop through all members
        for member in self.__dict__.keys():
            
            # Check if the variable is a string or not and write it out
            if isinstance(eval("self." + member),str):
                file_handle.write(member + ' = "' + str(eval("self." + member)) + '"\n')
                file_handle.write("\n")
            else:
                file_handle.write(member + " = " + str(eval("self." + member)) + "\n")
                file_handle.write("\n")
                    
        # Close the file
        file_handle.close()
    
    def load_options_from_file(self,options_file):
        """
        Description.

        Parameters
        ----------
        options_file : String
            DESCRIPTION.

        Raises
        ------
        Exception
            DESCRIPTION.

        Returns
        -------
        None.

        """
        import logging
        
        # Try to read the options file        
        try:
            exec(open(options_file).read())
            logging.info("Successfully read plot options file")
        except:
            raise Exception('Not able to read/run that plot options script')
    
        # Get a dictionary of all local variables
        # The only local variables will be ones set from the options script
        local_var_dir = locals()
    
        # Loop through all local variables 
        for var in local_var_dir:
            if var != "self" and var != "options_file":
                if hasattr(self,var):
                    setattr(self,var,local_var_dir[var])
                else: 
                    logging.info("Unknown option: " + var)
                    
        logging.info("Successfully taken PEATSA delivery order")

class PEATSAmenu(object):
    def __init__(self,options_file = None):
        
        # Define the sections
        self.section_order = ['Path Options',"Start Options","Initial Case Creation Options","EMTG Control Options","Optimization Options","Sorting Options","Case Re-creation Options","Post-processing Options"]
        self.sections = {"Path Options" : ["run_name","working_directory","nCores","emtg_root_directory","executable_name","logfile","thin_crust"],
                         "Start Options" : ["execution_type", "max_iterations", "keep_only_current_and_previous", "start_type","PEATSA_type","if_run_cases","iteration","restart_run_root_directory","copy_previous_results"],
                         "Initial Case Creation Options" : ["custom_PEATSA_chef_path","nSamples","thrusting_standard_deviations","epoch_range","stop_after_journey","reference_files_location","reference_options_file_name","reference_mission_file_name","maneuver_file","reference_spk_data_file","trade_study_options_files","starting_case_index","missed_thrust_timestep","skip_forced_coast_periods","convert_3d_to_patched_conic","post_missed_thrust_electric_propellant_margin","pre_missed_thrust_electric_propellant_margin","vent_pre_missed_thrust_electric_propellant_margin","missed_thrust_recovery_forced_terminal_coast_reduction","override_propulated_duty_cycle"],
                         "EMTG Control Options" : ["killtime", "MBH_max_run_time"],
                         "Optimization Options" : ["objective_type","objective_formula","max_or_min","peatsa_goal_objective","polyfit_order","polyfit_margin","constraint_walking"],
                         "Sorting Options" : ["fingerprint"],
                         "Case Re-creation Options" : ["find_seeds_in_parallel", "find_seeds_in_parallel_nCores", "write_seeds_in_parallel", "write_seeds_in_parallel_nCores", "seed_folders","seed_from_seed_folders_on_fresh_start","allow_cases_to_seed_themselves","initial_guess_type","new_timestep","turn_off_seed_MBH_when_unseeded","seed_from_infeasible_self_when_unseeded","seed_from_infeasible_cases_if_no_feasible_seeds_found","seed_from_cases_that_havent_met_target","only_one_seed_in_all_seed_directions","seed_criteria","wait_until_seeded","only_rerun_if","override_options"],
                         "Post-processing Options" : ["parse_in_parallel","parse_in_parallel_nCores", "extra_csv_column_definitions","generate_default_plots","built_in_plotter_files","custom_post_processing_script","move_results"]}
        self.descriptions = {}
        self.requirements = {}
        self.available_midbake = {}
        
        # -------------- DEFINE ALL DEFAULT VALUES ----------------------------------                
        self.seed_from_seed_folders_on_fresh_start = False
        self.init_var(
            "seed_from_seed_folders_on_fresh_start",
            "self.start_type == 'Fresh'",
            "False",
       "# When fresh starting, if you have an external seed folder, would you like to use it before the first iteration runs?"
        )
        
        
        self.seed_from_infeasible_cases_if_no_feasible_seeds_found = False
        self.init_var(
            "seed_from_infeasible_cases_if_no_feasible_seeds_found",
            "True",
            "True",
       "# Would you like cases to be able to be seeded from infeasible seeds if no feasible seeds were found?"
        )
        
        self.override_propulated_duty_cycle = 0.0
        self.init_var(
            "override_propulated_duty_cycle",
            "self.PEATSA_type == 5",
            "True",
       "# Would you like to use the default duty cycle to propulate or a different value? If a different\n\
        # value, then enter the duty cycle (i.e. .95). Otherwise enter zero."
        )
        
        self.thin_crust = False
        self.init_var(
            "thin_crust",
            "True",
            "False",
       "# Would you like PEATSA to run in a slimmer memory mode? It will use less memory but will require more time between iterations."
        )
        
        self.keep_only_current_and_previous = False
        self.init_var(
            "keep_only_current_and_previous",
            "True",
            "False",
            "# Would you like PEATSA to keep only two Cases and Results subdirectories? (IterationPrevious and IterationCurrent).\n\
             # If true, you will use less disk space but lose the exhaustive history of PEATSA iterations.\n\
             # Currently can only be used with a Fresh start."
             )
        
        self.constraint_walking = []
        self.init_var(
            "constraint_walking",
            "True",
            "True",
       "# Would you like to release constraints if PEATSA stalls? Enter a list of tuples. Each tuple has the form:\n\
        # (string of constraint as a PyEMTG MissionOption,how many failed iterations to wait before loosening constraint, increment)\n\
        # The first entry should be a variable in the PyEMTG MissionOption class, such as: 'MO.Journeys[-1].arrival_date_bounds[1]'\n\
        # The second entry is how many times no PEATSA cases found improvement before the constraint will be loosened.\n\
        # Finally enter the increment to change the constraint each time this event is triggered."
        )
    
        self.move_results = ""
        self.init_var(
            "move_results",
            "True",
            "True",
       "# Would you like peatsa to move the latest results to another folder for easier ftp'ing?\n\
        # If so, provide a full path to the desired destination. If not, leave an empty string."
        )
        
        self.turn_off_seed_MBH_when_unseeded = 0
        self.init_var(
            "turn_off_seed_MBH_when_unseeded",
            "self.wait_until_seeded == 0",
            "True",
       "# If a case doesnt have a seed, rather than running unseeded, should it keep using whatever is in its trialX,\n\
        # presumably the seed from reference mission. Or do you want to turn off seed_MBH and actually start without an initial guess?"
        )

        self.seed_from_infeasible_self_when_unseeded = 0
        self.init_var(
            "seed_from_infeasible_self_when_unseeded",
            "self.wait_until_seeded == 0",
            "True",
       "# If a case doesnt have a seed, rather than running unseeded, should it seed from the end of its last run,\n\
        # even though it didnt find something feasible at the end of the last run?"
        )
        
        self.starting_case_index = 0
        self.init_var(
            "starting_case_index",
            "self.start_type == 'Fresh' and self.PEATSA_type == 2",
            "False",
       "# If you are just using this code to create new cases that will be added into another \n\
        # PEATSA run, then you can start the filename's X in 'TradeStudy_CaseX' at a nonzero value.\n\
        # What index should new cases start at?"
        )
        
        
        self.allow_cases_to_seed_themselves = 1
        self.init_var(
            "allow_cases_to_seed_themselves",
            "1",
            "True",
       "# Should cases be allowed to seed themselves? Or only their neighbors? Enter 1 or 0"
        )
        
        self.initial_guess_type = "copy"
        self.init_var(
            "initial_guess_type",
            "1",
            "True",
       "# How do you want to generate initial guesses?\n\
        # if initial_guess_type == 'x_days', then the number of timesteps in the mission itself\n\
        #                                    will be modified such that each timestep is ~ x days\n\
        #                                    (use new_timestep to define x)\n\
        # if initial_guess_type == 'n_steps', then the number of timesteps in the mission will\n\
        #                                     be kept the same or modified to be n steps\n\
        #                                     (use new_timestep to define x)\n\
        # if initial_guess_type == 'copy', then PEATSA will just copy the initial guess over without\n\
        #                                  attempting to interpolate"
        )
        
        self.new_timestep = -1
        self.init_var(
            "new_timestep",
            "self.initial_guess_type != 'copy'",
            "True",
       "# What do you want the new timestep to be?\n\
        # if initial_guess_type == 'x_days',\n\
        #          if new_timestep == -2, then the number of days in each journey control step will be copied from the seed mission\n\
        #          if new_timetsep > 0, then all journeys will be given this timestep in days\n\
        #          if new_timestep = [x,y,z], then the final 3 journeys will be given x,y,z days\n\
        #                                     timesteps respectively. if len(new_timestep) > len(new_timestep),\n\
        #                                     entries from new_timestep will be matched backwards. I.E:\n\
        #                                     MO.Journeys[-1].journey_number_of_timesteps = new_timestep[-1]\n\
        #                                     MO.Journeys[-2].journey_number_of_timesteps = new_timestep[-2], etc.\n\
        #                                     Note: len(new_timestep)  < len(Journeys) will throw an error\n\
        # if initial_guess_type == 'n_steps,\n\
        #          if new_timestep == -2, then the number of timesteps will be copied from seed mission\n\
        #          if new_timestep == -1, then the number of timesteps will not be modified from mission options\n\
        #          if new_timestep > 0, then all journeys will be given this number of timesteps\n\
        #          if new_timestep = [x,y,z], then the final 3 journeys will be given x,y,z\n\
        #                                     timesteps respectively. if len(new_timestep) > len(new_timestep),\n\
        #                                     entries from new_timestep will be matched backwards. I.E:\n\
        #                                     MO.Journeys[-1].number_of_timesteps = new_timestep[-1]\n\
        #                                     MO.Journeys[-2].number_of_timesteps = new_timestep[-2], etc.\n\
        #                                     Note: len(new_timestep)  < len(Journeys) will throw an error"
        )
        
        
        self.find_seeds_in_parallel = 0
        self.init_var(
            "find_seeds_in_parallel",
            "1",
            True,
            "# CURRENTLY BUGGY, DO NOT USE.\n\
            # Do you want to find seeds in parallel? (0 or 1).\n\
            # Uses find_seeds_in_parallel_nCores number of cores.")
            
        self.find_seeds_in_parallel_nCores = 1
        self.init_var(
            "find_seeds_in_parallel_nCores",
            "1",
            True,
            "# If find_seeds_in_parallel == 1, this is how many parallel processes are used to find seeds. Integer >= 1.")
            
        self.write_seeds_in_parallel = 0
        self.init_var(
            "write_seeds_in_parallel",
            "1",
            True,
            "# CURRENTLY BUGGY, DO NOT USE.\n\
            # Do you want to write seeds in parallel? (0 or 1).")
            
        self.write_seeds_in_parallel_nCores = 1
        self.init_var(
            "write_seeds_in_parallel_nCores",
            "1",
            True,
            "# If write_seeds_in_parallel == 1, this is how many parallel processes are used to write seeds. Integer >= 1.")
            
        self.seed_folders = []
        self.init_var(
            "seed_folders",
            "1",
            True,
       "# Are there additional folders that contain emtg cases that can be used as seeds\n\
        # for the current PEATSA run, but are not valid cases. For example, a group of results\n\
        # that seem feasible to PEATSA, but the EMTG feasibility criteria changed internally to\n\
        # EMTG, so all the old cases need to be tossed, but the initial guesses are still useful.\n\
        # Enter a list of tuples:\n\
        # (folder_type,full_folder_path)\n\
        # if folder_type == 0, then the folder is a peatsa results folder with many iteration\n\
        #                      subdirectories, all of which will be included\n\
        # if folder_type == 1, then the folder is an individual results folder and only its contents\n\
        #                      will be parsed, with all subdirectories ignored."
        )
        
        self.emtg_root_directory = "/Utilities/emtg/" 
        self.init_var(
            "emtg_root_directory", 
            "1",
            False,
       "# Where is EMTG located? The emtgv9 executable should be in this directory and \n\
        # the PyEMTG files should be located in: emtg_root_directory + '/PyEMTG'.\n\
        # Enter a string.")

        self.executable_name = 'EMTGv9'
        self.init_var(
            "executable_name",
            "1",
            False,
       "# What is the name of the executable? Enter as a string. Typically emtg or EMTGv9")

        self.skip_forced_coast_periods = 1
        self.init_var(
            "skip_forced_coast_periods",
            "self.start_type == 'Fresh' and (self.PEATSA_type == 1 or self.PEATSA_type == 3)",
            False,
       "# Should missed thrust recovery cases be started during forced coasting periods? Enter a 1 or 0.")

        self.convert_3d_to_patched_conic = 0
        self.init_var(
            "convert_3d_to_patched_conic",
            "self.start_type == 'Fresh' and (self.PEATSA_type == 1 or self.PEATSA_type == 3)",
            False,
       "# If the mission has coast phases, would you like to remove these from the missed thrust recovery cases?"
        )

        self.if_run_cases = 0
        self.init_var( 
            "if_run_cases",
            "1",
            False,
       "# Should cases be run or not? Should the PEATSA oven be called? Enter a 1 or 0.")

        self.logfile = "peatsa.out"
        self.init_var( 
            "logfile", 
            "1",
            False,
       "# Where would you like to log the results? Enter a string.")
        
        self.execution_type = 0
        self.init_var("execution_type",
                      "1",
                      True,
                      "#Choose which execution method is used.\n\
                        # 0: Traditional PEATSA\n\
                        # 1: Noble's new version\n\
                        # 2: Noble's new verion with pebble")

        self.start_type = "Fresh"
        self.init_var( 
            "start_type",
            "1",
            False,
       "# How should PEATSA start?\n\
        # start_type = 'Fresh', New PEATSA cases will be created\n\
        # start_type = 'Warm', PEATSA will start by running cases in an existing folder\n\
        # start_type = 'Hot', PEATSA will start by parsing a set of cases, re-seeding them\n\
        #                       and then running those new cases")

        self.run_name = "PEATSA_run"
        self.init_var( 
            "run_name",
            "1",
            False,
       "# Where should PEATSA work? Enter a string.")

        self.working_directory = "/Utilities/emtg/scratch/PEATSA"
        self.init_var( 
            "working_directory",
            "1",
            False,
       "# Where should PEATSA work? Enter a string.")

        self.copy_previous_results = 0
        self.init_var(
            "copy_previous_results",
            "self.start_type == 'Hot'",
            False,
       "# If start_type == 'Hot', do you want to copy the results from previous peatsa runs\n\
        # into the current run's iteration0 results folder? This will take a long time, but\n\
        # will save time in future hot starts. Enter a 1 or 0.")
        
        self.max_iterations = 1e10
        self.init_var(
            "max_iterations",
            "1",
            True,
            "# What is the maximum number of iterations PEATSA is allowed to run before exiting? (Integer >= 0)"
            )

        self.iteration = 0
        self.init_var( 
            "iteration",
            "self.start_type == 'Warm'",
            False,
       "# If start_type == 'Warm', what iteration to start from? Enter an integer.")

        self.restart_run_root_directory = [("/Utilities/emtg/scratch/PEATSA/PEATSA_run",'Full')]
        self.init_var(
            "restart_run_root_directory",
            "self.start_type == 'Warm' or self.start_type == 'Hot'",
            False,
       "# If start_type == 'Hot' or 'Warm', where are the results to start from? Enter a list of tuples:\n\
        # ('path to results','parse_type'\n\
        # parse_type = 'Full' : PEATSA will load every emtg and emtgopt in every subdirectory pointed to by the path/to/results\n\
        #                       It is assumed that the folder structure is: /path/to/results/Iteration#/EMTGandEMTGOPTfiles\n\
        # parse_type = 'Best' : PEATSA will only load the emtg and emtgopt files specified in the latest iteration csv in each /path/to/results/docs\n\
        #                       It is assumed that the path/to/results is actually a path to a peatsa root folder\n\
        # parse_type = 'Warm' : PEATSA will not load any results, it will just run the cases in /path/to/results/cases/"
        )

        self.PEATSA_type = 0
        self.init_var( 
            "PEATSA_type",
            "1",
            False,
       "# What type of PEATSA run is this? \n\
        # PEATSA_type = 0, custom \n\
        # PEATSA_type = 1, missed thrust\n\
        # PEATSA_type = 2, trade study\n\
        # PEATSA_type = 3, missed thrust trade study\n\
        # PEATSA_type = 4, batch run a single case\n\
        # PEATSA_type = 5, maneuver execution error monte carlo")
        
        self.maneuver_file = ""
        self.init_var(
            "maneuver_file",
            "1",
            False,
       "# Provide a path to the manuever_spec file for the reference mission"
        )
         
        self.nSamples = 100
        self.init_var(
            "nSamples",
            "self.PEATSA_type == 5",
            False,
       "# How many samples would you like to run in the monte carlo?"    
        )
        
        self.thrusting_standard_deviations = [0.0,0.0,0.0]
        self.init_var(
            "thrusting_standard_deviations",
            "self.PEATSA_type == 5 and self.start_type == 'Fresh'",
            False,
       "# What are the standard deviations of thrusting error? Enter a list:\n\
        # [% Magnitude error,Body Fixed Angle 1 [deg], Body Fixed Angle 2 [deg]]"    
        )
        
        self.epoch_range = [0,1e10]
        self.init_var(
            "epoch_range",
            "self.PEATSA_type == 5 and self.start_type == 'Fresh'",
            False,
       "# What is the range of epochs over which you want to monte carlo the maneuvers? Enter a list:\n\
        # [Earliest Julian Date, Latest Julian Date]"    
        )

        self.stop_after_journey = -1
        self.init_var(
            "stop_after_journey",
            "self.PEATSA_type == 5 and self.start_type == 'Fresh'",
            False,
       "# When optimizing monte carlo 'missions to go', do you want to optimize through the full reference mission,\n\
        # or stop after a specified journey? Enter the last journey you want in the 'mission to go' optimizations.\n\
        # Note: python negative indexing is fine"    
        )
        
        self.nCores = [("local",1)]
        self.init_var( 
            "nCores",
            "1",
            True,
       "# What machine(s) will run the emtg cases? \n\
        # This is a list of tuples. Each tuple is: \n\
        # (an IP address string or 'local' for the current machine, number of cores available on that host)\n\
        # NOTE! DO NOT USE ANYTHING BUT LOCAL UNLESS ALL MACHINES HAVE SHARED MEMORY AND\n\
        # working_directory is pointing to a shared folder on all machines")

        self.killtime = 1000000
        self.init_var( 
            "killtime",
            "1",
            True,
       "# How long (in seconds) before a case needs to be killed?\n\
       # Make this value at least the MBH max run time + the NLP max run time + some buffer time\n\
       # for writing to file. Buffer time depends on the complexity of the case, whether\n\
       # or not you are generating ephemeris files, etc. Recommend at least 30 seconds.")

        self.MBH_max_run_time = 1000000
        self.init_var( 
            "MBH_max_run_time",
            "self.start_type == 'Fresh'",
            False,
       "# How long (in seconds) should cases be run in EMTG?")

        self.objective_type = 0 
        self.init_var( 
            "objective_type",
            "self.PEATSA_type != 4",
            True,
       "# What is the objective type?\n\
        # objective_type == 0, objective_function better than peatsa_objective\n\
        # objective_type == 1, objective_function better than polyfit = objective_fun(seed_value)")

        self.objective_formula = ""
        self.init_var( 
            "objective_formula",
            "1",
            False,
       "# How can PEATSA evaluate the objective value?\n\
        # Input a string that can be evaluated. EMTG Mission data can be gathered using 'M.' and \n\
        # Mission option data is available using 'MO.'")

        self.max_or_min = "max"
        self.init_var( 
            "max_or_min",
            "1",
            False,
       "# Is this a maximization or minimization problem? Enter 'max' or 'min'.")

        self.peatsa_goal_objective = 0.0
        self.init_var( 
            "peatsa_goal_objective",
            "self.objective_type == 0",
            True,
       "# If PEATSA's objective_type == 0, what is the target objective value? I.E. for missed thrust, the \n\
        # target could be 60 days. Or for maximum mass, the target could be 100 \n\
        # kg. In some cases the goal objective will be the same as the EMTG\n\
        # objective_type, but it doesnt have to be. If a good target is not\n\
        # known, then simply set it to something impossible to reach. I.e. for \n\
        # missed thrust, 1e100, or for minumum time, -1 days. The value here should \n\
        # have the same units as whatever the 'objective_formula' will evaluate to.")

        self.polyfit_order = 1
        self.init_var( 
            "polyfit_order",
            "self.objective_type == 1 or 4 in [criteria[3] for criteria in self.seed_criteria]",
            True,
       "# If PEATSA's objective_type == 1, then an order for the polynomial fit is needed. Enter an integer.")

        self.polyfit_margin = .1
        self.init_var( 
            "polyfit_margin",
            "self.objective_type == 1",
            True,
       "# If PEATSA's objective_type == 1, then a margin is used so that the case value does not\n\
        # need to exceed the EXACT expected value from the polynomial fit. For example, if the polyfit\n\
        # predicts a value of 100, and the result for this case is 99, then we would need a margin of .01 \n\
        # or 1% to prevent the case from being re-run. If the result for the case is 95, then it still is likely \n\
        # 'good enough', and doesnt need to be re-run, so a margin of 5 to 10% is typically appropriate. Enter a float 0 to 1.")

        self.only_one_seed_in_all_seed_directions = 0
        self.init_var( 
            "only_one_seed_in_all_seed_directions",
            "self.PEATSA_type != 4",
            True,
       "# Flag to override the individual seed selection criterias specified in 'seed_criteria'.\n\
        # Instead, only one total seed will be selected regardless of which seed criteria was used. \n\
        # Note: This variable is a replacement for the old only_seed_from_best_in_all_seed_directions variable \n\
        # Valid choices: \n\
        # 0: Do not use \n\
        # 1: Only the best seed from all directions will be used (i.e., the old only_seed_from_best_in_all_seed_directions) \n\
        # 2: One randomly selected seed from amongst all seed directions that produced a valid seed will be used.")

        self.seed_from_cases_that_havent_met_target = 0
        self.init_var( 
            "seed_from_cases_that_havent_met_target",
            "1",
            True,
       "# Should cases that don't meet the goal criteria be considered as potential seeds?\n\
        # WARNING: potential seed cases will be evaluated against peatsa_goal_objective\n\
        # even if PEATSA_objective_type == 1. \n\
        # Enter a 1 or 0.")

        self.seed_criteria = []
        self.init_var( 
            "seed_criteria",
            "1",
            True,
       "# How can PEATSA evaluate the seed criteria?\n\
        # Input a list of tuples. Each tuple is: \n\
        # (a string that can be evaluated, max negative seed range, max positive seed range, seed selection criteria)\n\
        # For the evaluation string, EMTG Mission data can be gathered using\n\
        # 'MO.'. Don't use mission data, because even failed cases need to know how to\n\
        # pick a seed.\n\
        # For the max seed range, answer if there should be a maximum range for seeding, and if so, how long?\n\
        # max_negative_seed_range > 0, there is no maximum range in the negative direction\n\
        # max_negative_seed_range <= 0, this is the maximum range, in the units of whatever\n\
        #       seed_criteria evaluates in \n\
        # max_positive_seed_range < 0, there is no maximum range in the negative direction\n\
        # max_positive_seed_range >= 0, this is the maximum range, in the units of whatever\n\
        #       seed_criteria evaluates in \n\
        # For the seed selection criteria, specify how PEATSA should select new seeds?\n\
        # seed_selection_criteria == 0, seed from the closest potential seed\n\
        # seed_selection_criteria == 1, seed from all potential seeds\n\
        # seed_selection_criteria == 2, seed from only the best potential seed\n\
        # seed_selection_criteria == 3, seed from random potential seed\n\
        # seed_selection_criteria == 4, seed from a curve fit of all potential seeds.\n\
        #                               Each decision vector variable is a curve fit\n\
        #                               of a curve fit of the corresponding decision vector values\n\
        #                               from all of the seeds. Uses polyfit_order for the curve fit order.")

        self.extra_csv_column_definitions = []
        self.init_var( 
            "extra_csv_column_definitions",
            "1",
            True,
       "# Do you want any extra data printed to each iteration's csv summary?\n\
        # If not, leave this as an empty list\n\
        # If so, add tuples of:\n\
        # ('column header','string of formula to be evaluated')\n\
        # For data from an EMTG mission, use 'M.', and for an EMTG mission option, use 'MO.'")

        self.parse_in_parallel = 0
        self.init_var(
            "parse_in_parallel",
            "1",
            True,
       "# Should peatsa parse using multiple local cores? Enter a 1 or 0.\n\
        # A 1 means that peatsa will parse results in parallel using parse_in_parallel_nCores.\n\
        # Note that starting up parallelization is not free, which means that parallelizing will actually\n\
        # be SLOWER if the number of cases is not much greater than the number of cores being used.\n\
        # However, if you have thousands or tens of thousands of cases, you can definitely save time by parallelizing.")
        
        self.parse_in_parallel_nCores = 1
        self.init_var(
            "parse_in_parallel_nCores",
            "1",
            True,
       "# If parse_in_parallel == 1, how many cores should PEATSA use to parse in parallel?\n\
        # Enter an integer >= 1 and <= the number of available cores.")

        self.fingerprint = []
        self.init_var( 
            "fingerprint",  
            "self.PEATSA_type != 4",
            False,     
       "# What makes each case in this PEATSA run unique? Write a list of strings\n\
        # that will be evaluated for each case to give it a unique fingerprint. \n\
        # For data from an EMTG mission option, use 'MO.'. Do not use data\n\
        # from an EMTG mission. Cases need fingerprints even if they dont converge.\n\
        # Do not put any formulas that are already in the seed_criteria list.")

        self.reference_files_location = ""
        self.init_var( 
            "reference_files_location",
            "(self.start_type == 'Fresh' and (self.PEATSA_type == 1 or self.PEATSA_type == 4)) or self.PEATSA_type == 5",
            False,
       "# Where are the reference files? Enter a path as a string.")

        self.reference_options_file_name = ""
        self.init_var( 
            "reference_options_file_name",
            "(self.start_type == 'Fresh' and (self.PEATSA_type == 1 or self.PEATSA_type == 4)) or self.PEATSA_type == 5",
            False,  
       "# What is the name of the reference options file? Enter a string.")

        self.reference_mission_file_name = ""
        self.init_var( 
            "reference_mission_file_name",
            "self.start_type == 'Fresh' and (self.PEATSA_type == 1 or self.PEATSA_type == 3 or self.PEATSA_type == 5)",
            False,  
       "# What is the name of the reference mission file? Enter a string.")

        self.missed_thrust_timestep = 7
        self.init_var( 
            "missed_thrust_timestep",
            "self.start_type == 'Fresh' and (self.PEATSA_type == 1 or self.PEATSA_type == 3)",
            False,  
       "# If running a missed thrust case, how frequently, should cases be created?\n\
        # If missed_thrust_timestep <= 1 : PEATSA will create a timestep that is that percentage of the given phase's timestep\n\
        # If missed_thrust_timestep > 1 : PEATSA will use the timestep directly as number of days")

        self.reference_spk_data_file = ""
        self.init_var( 
            "reference_spk_data_file",
            "self.start_type == 'Fresh' and (self.PEATSA_type == 1 or self.PEATSA_type == 3 or self.PEATSA_type == 5)",
            False,  
       "# If running a missed thrust case, what is the path to the spk input data file\n\
        # that contains ephemeris for the reference trajectory. Enter a path as a string.\n\
        # PEATSA expects a .ephemeris file with a mass column.")

        self.post_missed_thrust_electric_propellant_margin = 0.0
        self.init_var( 
            "post_missed_thrust_electric_propellant_margin",
            "self.start_type == 'Fresh' and (self.PEATSA_type == 1 or self.PEATSA_type == 3)",
            False,  
       "# If running a missed thrust case, what propellant margin should the mission maintain?\n\
        # for the propellant remaining when the missed thrust event is over? Because Bruno Sarli\n\
        # is bad at math, you can specify units using either 'pct' or 'kg'. The whole thing\n\
        # must be a string. For example, '10 pct' or '100 kg'.\n\
        # If no units are specified, it will be assumed that the value is a percentage as a number from 0 to 1.")

        self.pre_missed_thrust_electric_propellant_margin = -1.0
        self.init_var( 
            "pre_missed_thrust_electric_propellant_margin",
            "self.start_type == 'Fresh' and (self.PEATSA_type == 1 or self.PEATSA_type == 3)",
            False,  
       "# If running a missed thrust case, what propellant margin should the mission consider\n\
        # itself to have used before the missed thrust event?\n\
        # For example, if the mission used 200 kg of propellant before the missed thrust\n\
        # event, and there are 100 kg of propellant remaining 'in the tank', we need to\n\
        # determine how much of the remaining 100 kg are actually usable. If the mission\n\
        # is to maintain a 10%% propellant margin over the entire mission, then only 80 kg\n\
        # of the remaning propellant are usable, because the margin on the 200 kg used\n\
        # so far is 20 kg.  \n\
        # pre_missed_thrust_electric_propellant_margin < 0, use the value from the reference mission\n\
        # pre_missed_thrust_electric_propellant_margin > 0 use the percentage specified (enter a number 0 to 1)")

        self.vent_pre_missed_thrust_electric_propellant_margin = 0
        self.init_var( 
            "vent_pre_missed_thrust_electric_propellant_margin",
            "self.start_type == 'Fresh' and (self.PEATSA_type == 1 or self.PEATSA_type == 3)",
            False,  
       "# If running a missed thrust case, should the pre-missed thrust event propellant\n\
        # be vented? Once the missed thrust event is over, the propellant that was book-kept\n\
        # as margin for the propellant used up until the event, is nominally converted to\n\
        # 'dry mass' in EMTG. However, this mass could be vented instead, which would reduce\n\
        # the mass on-board for missed thrust recovery\n\
        # Enter 1 or 0.")
        
        self.missed_thrust_recovery_forced_terminal_coast_reduction = 0.0
        self.init_var(
            "missed_thrust_recovery_forced_terminal_coast_reduction",
            "self.start_type == 'Fresh' and (self.PEATSA_type == 1 or self.PEATSA_type == 3)",
            False,  
       "# How many days should be removed from any terminal forced coasts for missed thrust recovery cases?\n\
        # If missed_thrust_recovery_forced_terminal_coast_reduction > 0 : PEATSA will reduce the coast by an amount that is that percentage of the given phase's timestep. Enter a float, but integers are probably smarter.\n\
        # If missed_thrust_recovery_forced_terminal_coast_reduction < 0 : PEATSA will use the reduction directly directly as a number of days.")
        
        self.override_options = [("1","MO.snopt_max_run_time = 300"),("1","MO.MBH_max_run_time = 600")]
        self.init_var(
            "override_options",
            "self.PEATSA_type != 4",
            True,  
       "# Are there any options that you want to override from the default options?\n\
        # Enter a tuple of: \n\
        # (string of condition of when to use the option, string of formula for option)\n\
        # Use MO for the mission options object.\n\
        # It is very likely that you want to at least override the directory paths for\n\
        # HardwareModels and Universe.\n\
        # Example: override_options = [('1','MO.launch_window_open_date = 23432'),('len(MO.Journeys) > 3','MO.Journeys[3].destination_list[1] = 2')]")

        self.wait_until_seeded = 0
        self.init_var(
            "wait_until_seeded",
            "self.PEATSA_type != 4",
            True,  
       "# NOTE: We are currently not sure whether or not this actually works.\n\
       # Should cases not be run until a seed can be found for them?\n\
        # Note: they will be run for the 0th iteration. But in the reseeding process,\n\
        # If a seed cant be found for it, then it will not be put back on the queue\n\
        # Enter 1 or 0")

        self.trade_study_options_files = []
        self.init_var(
            "trade_study_options_files",
            "(self.start_type == 'Fresh' or self.generate_default_plots) and (self.PEATSA_type == 2 or self.PEATSA_type == 3)",
            False,  
       "# If performing a trade study, where are the options being traded stored?\n\
        # Enter a list of tuples. Each tuple should have the format:\n\
        # (string of path to csv file,options file type)\n\
        # The path to the file should be full, not relative from anywhere.\n\
        # There are three formats that the file can be created. All three formats\n\
        # should be csv files though (It is assumed you will use excel to create the file\n\
        # and then save it as a csv). The three formats are:\n\
        # options file type = 0: Trade all values for a given variable individually against the default values of all other variables.\n\
        #                        First two rows specify the reference files location and emtgopt filename respectively\n\
        #                        Third row specifies the PyEMTG MissionOptions variables being traded. Use 'MO.' to access mission option data,\n\
        #                                                                                              Use 'SO.' to access spacecraft option data\n\
        #                                                                                              Use 'StgO[n].' to access stage n option data\n\
        #                        Fourth row lists the default value for the column it resides.\n\
        #                        Subsequent rows, as many as desired, indicate all other options for that variable.\n\
        #                        An example version of this file is in the PEATSA folder named 'opt_file_type_one.csv'\n\
        #                        The total number of cases generated will be equal to the sum of the number of all entries in the third row and below\n\
        # options file type = 1: Trade all values of each variable against all values of all other variables.\n\
        #                        First two rows specify the reference files location and emtgopt filename respectively\n\
        #                        Third row specifies the PyEMTG MissionOptions variables being traded. Use 'MO.' to access mission option data.\n\
        #                                                                                              Use 'SO.' to access spacecraft option data\n\
        #                                                                                              Use 'StgO[n].' to access stage n option data\n\
        #                        Subsequent rows, as many as desired, indicate all other options for that variable.\n\
        #                        An example version of this file is in the PEATSA folder named 'opt_file_type_two.csv'\n\
        #                        The total number of cases generated will be equal to the multiplicative of the number of entries in each column (not counting the top row)\n\
        # options file type = 2: List out the exact combinations of each variable.\n\
        #                        First two rows specify the reference files location and emtgopt filename respectively\n\
        #                        Third row specifies the PyEMTG MissionOptions variables being traded. Use 'MO.' to access mission option data.\n\
        #                                                                                              Use 'SO.' to access spacecraft option data\n\
        #                                                                                              Use 'StgO[n].' to access stage n option data\n\
        #                        Subsequent rows, as many as desired indicate the combinations of those variables to be considered.\n\
        #                        An example version of this file is in the PEATSA folder named 'opt_file_type_three.csv'\n\
        #                        The total number of cases generated will be equal to the number of rows (not counting the top row)"
        )

        self.custom_PEATSA_chef_path = ""
        self.init_var(
            "custom_PEATSA_chef_path",
            "self.start_type == 'Fresh' and self.PEATSA_type == 0",
            False,  
       "# If PEATSA_type == 0, enter the path to a python script which will generate\n\
        # the cases as part of the PEATSA run. This is meant to be used if the default\n\
        # PEATSA file generators for trade studies and missed thrust studies aren't\n\
        # sufficient for any reason. Enter a path as a string. The file must:\n\
        #      1)  Be valid python code, ending with .py\n\
        #      2)  Define a method called 'CustomPEATSAchef'\n\
        #      3)  CustomPEATSAchef must take a PEATSAmenu object as it's only input\n\
        #      4)  Write the mission objects into PEATSAorder.cases_directory\n\
        #      5)  Return a list of PEATSAbox objects with nothing but the PEATSAdough path set")

        self.only_rerun_if = []
        self.init_var(
            "only_rerun_if",
            "self.PEATSA_type != 4",
            True,  
       "# Enter a list of strings that will can be evaluated as a conditional statement. Each condition will filter out cases that should be ignored. These conditions\n\
        # are applied using 'or' logic, assuming that no cases will be re-run unless these conditions are met.\n\
        # If left empty, the normal conditions apply. If applied, the cases will still be compared to the PEATSA objective, so\n\
        # cases that meet the 'only_rerun_if', still might not be re-run if they also meet the stopping criteria.\n\
        # Note: this does not remove cases from the PEATSA run entirely. It only adds\n\
        # a condition, which if not met, will prevent cases from being re-run. But\n\
        # thoses filtered cases will still appear in csv summary files.\n\
        # Use 'M.' to access EMTG mission object data and 'MO.' to access EMTG mission options data.\n\
        # Example only_rerun_if = ['MO.launch_window_open_date > 64000','M.spacecraft_dry_mass < 100 and M.total_flight_time_years > 30']")

        self.generate_default_plots = 0
        self.init_var(
            "generate_default_plots",
            "self.PEATSA_type != 4",
            True,  
       "# NOTE: These plots are not necessarily all that useful, but see for yourself.\n\
       # Should the default plots be generated?\n\
        # For trade studies, each option being traded will be plotted \n\
        # on the x axis vs. the 'objective_value' on the y-axis. Note: this \n\
        # currently only works for trade study type 0 files.\n\
        # For missed thrust, missed thrust event date will be plotted on the x axis vs departure coast possible on the y-axis.")

        self.built_in_plotter_files = []
        self.init_var(
            "built_in_plotter_files",
            "self.PEATSA_type != 4",
            True,  
       "# Would you like to generate custom plots using the built-in plotter? If so, enter a list of paths as strings to csv files which\n\
        # contain options for a plot")

        self.custom_post_processing_script = ""
        self.init_var(
            "custom_post_processing_script",
            "1",
            True,  
       "# Would you like to call a custom post processing script? Enter a path as a string.\n\
        # This script will be called after each iteration. It can generate extra summary files, in addition to\n\
        # the .csv file generated already, or extra plots, in addition to the defaults or any made with the built-in plotting routine.\n\
        # Enter the path of a python module. The file must:\n\
        #      1)  Be valid python code, ending with .py\n\
        #      2)  Define a method called 'PEATSAPostProcess'\n\
        #      3)  PEATSAPostProcess must take a PEATSAmenu object as its first input\n\
        #      4)  PEATSAPostProcess must take a python list of PEATSAbox objects as its second input\n\
        # It is not necessary, but data files can be saved to 'PEATSAorder.docs_dir' and plots can be saved to 'PEATSAorder.images_dir'")
                   
        # --------------- LOAD USER OPTIONS ---------------------------
        if options_file != None:
            self.load_options_from_file(options_file)

    def load_options_from_file(self,options_file,if_mid_bake = False):
        # Read the options file        
        exec(open(options_file).read())
        if if_mid_bake == True:
            import logging
            logging.info("Successfully read options file")
        else:
            print("Successfully read options file")
    
        # Get a dictionary of all local variables
        # The only local variables will be ones set from the options script
        local_var_dir = locals()
    
        # Loop through all local variables 
        for var in local_var_dir:
            if var != "self" and var != "options_file" and var != "if_mid_bake":
                if hasattr(self,var):
                    setattr(self,var,local_var_dir[var])
                else: 
                    if if_mid_bake == True:
                        logging.info("Unknown option: " + var)
                    else:
                        print("Unknown option: " + var)
        
        needed_vars = self.get_needed_options()
        
        any_missing = False
        
        for var in needed_vars:
            if var not in local_var_dir:
                print("WARNING: Missing variable: " + var)
                any_missing = True
                
        if any_missing:
            print("WARNING: Maybe you need to run /PyEMTG/PEATSA/update_existing_PEATSA_script.py")
        
        # for idx,seed_crit in enumerate(self.seed_criteria):
        #     if len(seed_crit) == 3:
        #         seed_crit += ([],)
        #         self.seed_criteria[idx] = seed_crit      
                    
        if if_mid_bake == True:
            import logging
            logging.info("Successfully taken PEATSA order")
        else:
            print("Successfully taken PEATSA order")
    
    def get_needed_options(self):
        
        needed_vars = []
        
        allvars = dir(self)
        
        ignore_list = ["sections","section_order","descriptions","requirements","available_midbake","load_options_from_file","get_needed_options","init_var","write_to_file","conditionally_write_to_file"]
        
        for var in allvars:
            if var.startswith("__"):
                continue
            if var in ignore_list:
                continue
            
            if var not in self.requirements.keys():
                continue
            
            if (eval(self.requirements[var])):
                needed_vars.append(var)
        
        return needed_vars
            
    def init_var(self,var_name,requirements_string,midbake,description_string):
        self.descriptions.update({var_name : description_string})
        self.requirements.update({var_name : requirements_string})
        self.available_midbake.update({var_name : midbake})
        
    def write_to_file(self,filename,requirement_check = False,mid_bake_file = False):
        
        from math import floor, ceil
        
        # Get all of the member variables from the class, excpet for the descriptions dictionary
        # members = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__") and not isinstance(eval("self." + attr),dict)]
    
        # Open the output file
        file_handle = open(filename,'w')
    
        # Write the header for the sample options file
        file_handle.write("# PEATSA ")
        if mid_bake_file == True: 
            file_handle.write("mid-bake ")
        file_handle.write("options script" + "\n")
        file_handle.write("# written by PEATSA menu class" + "\n")
        file_handle.write("\n")
        file_handle.write("\n")

        # If this is a midbake file, let the user know these options can be changed while peatsa is running
        if mid_bake_file == True:
            file_handle.write("# These options can be modified while a PEATSA run is baking, and they will be updated as soon as the current iteration cases are done cooking" + "\n")
            file_handle.write("\n")
            file_handle.write("\n")
    
        # Define the length of the section headers
        section_header_length = 100
    
        # Loop through all the sections
        for section in self.section_order:
            
            # Get the length of the section title itself
            sec_len = len(section)
            
            file_handle.write("# ")
            # Loop through the front dashes
            for dash_idx in range(int(floor((section_header_length-sec_len)/2.0))):
                file_handle.write("-")
            
            # Write the section title
            file_handle.write(" " + section + " ")
            
            # Loop through the end dashes
            for dash_idx in range(int(ceil((section_header_length-sec_len)/2.0))):
                file_handle.write("-")
            
            file_handle.write("\n\n")
    
            # Loop through all member variables
            for member in self.sections[section]:
            
                # Check if this is a midbake file
                if mid_bake_file == True:
                    # It is.
                    
                    # Check if this variable is available mid-bake
                    if self.available_midbake[member] == False:
                        # It is not
                        continue
                        
                # Check if we are checking conditions for the option to be present
                if requirement_check:
                    # We are
                
                    # Check if the requirement is set
                    if not (eval(self.requirements[member])):
                        # Condition is not met, go to the next member variable
                        continue
            
                # Write the description comment to file
                file_handle.write(self.descriptions[member])
                file_handle.write("\n")

                # Check if the member variable is a string, and then write it out with or without quotes as appropriate
                if isinstance(eval("self." + member),str):
                    file_handle.write(member + " = '" + str(eval("self." + member)) + "'\n")
                    file_handle.write("\n")
                elif isinstance(eval("self." + member),list):
                    listvals = eval("self." + member)
                    file_handle.write(member + " = [" + "\n")
                    for eidx,entry in enumerate(listvals):
                        if isinstance(entry,str):
                            if '"' in entry:
                                file_handle.write(" " * len(member) + "    '" + str(entry) + "'")
                            else:
                                file_handle.write(" " * len(member) + '    "' + str(entry) + '"')
                        else:
                            file_handle.write(" " * len(member) + "    " + str(entry))
                        if eidx + 1 < len(listvals):
                            file_handle.write(",")
                        file_handle.write("\n")
                    file_handle.write(" " * len(member) + "   ]" + "\n")
                    file_handle.write("\n")
                else:
                    file_handle.write(member + " = " + str(eval("self." + member)) + "\n")
                    file_handle.write("\n")
    
            file_handle.write("\n")
    
        # Close the file
        file_handle.close()
        
    def conditionally_write_to_file(self,filename):
        # Nothing to do but call the normal file writer, but with the requirements flag set
        self.write_to_file(filename,True)