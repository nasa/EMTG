# PEATSA options script
# written by PEATSA menu class


# -------------------------------------------- Path Options --------------------------------------------

# Where should PEATSA work? Enter a string.
run_name = 'peatsa_test'

# Where should PEATSA work? Enter a string.
working_directory = '/Utilities/emtg/projects/peatsa_tests/'

# What machine(s) will run the emtg cases? 
        # This is a list of tuples. Each tuple is: 
        # (an IP address string or 'local' for the current machine, number of cores available on that host)
        # NOTE! DO NOT USE ANYTHING BUT LOCAL UNLESS ALL MACHINES HAVE SHARED MEMORY AND
        # working_directory is pointing to a shared folder on all machines
nCores = [
          ('local', 4)
         ]

# Where is EMTG located? The emtgv9 executable should be in this directory and 
        # the PyEMTG files should be located in: emtg_root_directory + '/PyEMTG'.
        # Enter a string.
emtg_root_directory = '/Utilities/emtg/bin/'

# What is the name of the executable? Enter as a string. Typically emtg or EMTGv9
executable_name = 'EMTGv9'

# Where would you like to log the results? Enter a string.
logfile = 'peatsa.out'

# Would you like PEATSA to run in a slimmer memory mode? It will use less memory but will require more time between iterations.
thin_crust = False


# ------------------------------------------- Start Options --------------------------------------------

#Choose which execution method is used.
                        # 0: Traditional PEATSA
                        # 1: Noble's new version
execution_type = 0

# How should PEATSA start?
        # start_type = 'Fresh', New PEATSA cases will be created
        # start_type = 'Warm', PEATSA will start by running cases in an existing folder
        # start_type = 'Hot', PEATSA will start by parsing a set of cases, re-seeding them
        #                       and then running those new cases
start_type = 'Fresh'

# What type of PEATSA run is this? 
        # PEATSA_type = 0, custom 
        # PEATSA_type = 1, missed thrust
        # PEATSA_type = 2, trade study
        # PEATSA_type = 3, missed thrust trade study
        # PEATSA_type = 4, batch run a single case
        # PEATSA_type = 5, maneuver execution error monte carlo
PEATSA_type = 2

# Should cases be run or not? Should the PEATSA oven be called? Enter a 1 or 0.
if_run_cases = 1

# What is the maximum number of iterations PEATSA is allowed to run before exiting? (Integer >= 0)
max_iterations = 3


# ----------------------------------- Initial Case Creation Options ------------------------------------

# Provide a path to the manuever_spec file for the reference mission
maneuver_file = ''

# If performing a trade study, where are the options being traded stored?
        # Enter a list of tuples. Each tuple should have the format:
        # (string of path to csv file,options file type)
        # The path to the file should be full, not relative from anywhere.
        # There are three formats that the file can be created. All three formats
        # should be csv files though (It is assumed you will use excel to create the file
        # and then save it as a csv). The three formats are:
        # options file type = 0: Trade all values for a given variable individually against the default values of all other variables.
        #                        First two rows specify the reference files location and emtgopt filename respectively
        #                        Third row specifies the PyEMTG MissionOptions variables being traded. Use 'MO.' to access mission option data,
        #                                                                                              Use 'SO.' to access spacecraft option data
        #                                                                                              Use 'StgO[n].' to access stage n option data
        #                        Fourth row lists the default value for the column it resides.
        #                        Subsequent rows, as many as desired, indicate all other options for that variable.
        #                        An example version of this file is in the PEATSA folder named 'opt_file_type_one.csv'
        #                        The total number of cases generated will be equal to the sum of the number of all entries in the third row and below
        # options file type = 1: Trade all values of each variable against all values of all other variables.
        #                        First two rows specify the reference files location and emtgopt filename respectively
        #                        Third row specifies the PyEMTG MissionOptions variables being traded. Use 'MO.' to access mission option data.
        #                                                                                              Use 'SO.' to access spacecraft option data
        #                                                                                              Use 'StgO[n].' to access stage n option data
        #                        Subsequent rows, as many as desired, indicate all other options for that variable.
        #                        An example version of this file is in the PEATSA folder named 'opt_file_type_two.csv'
        #                        The total number of cases generated will be equal to the multiplicative of the number of entries in each column (not counting the top row)
        # options file type = 2: List out the exact combinations of each variable.
        #                        First two rows specify the reference files location and emtgopt filename respectively
        #                        Third row specifies the PyEMTG MissionOptions variables being traded. Use 'MO.' to access mission option data.
        #                                                                                              Use 'SO.' to access spacecraft option data
        #                                                                                              Use 'StgO[n].' to access stage n option data
        #                        Subsequent rows, as many as desired indicate the combinations of those variables to be considered.
        #                        An example version of this file is in the PEATSA folder named 'opt_file_type_three.csv'
        #                        The total number of cases generated will be equal to the number of rows (not counting the top row)
trade_study_options_files = [
                             ('/Utilities/emtg/projects/peatsa_tests/peatsa_test.csv', 1)
                            ]

# If you are just using this code to create new cases that will be added into another 
        # PEATSA run, then you can start the filename's X in 'TradeStudy_CaseX' at a nonzero value.
        # What index should new cases start at?
starting_case_index = 0


# ---------------------------------------- EMTG Control Options ----------------------------------------

# How long (in seconds) before a case needs to be killed?
killtime = 20

# How long (in seconds) should cases be run in EMTG?
MBH_max_run_time = 11


# ---------------------------------------- Optimization Options ----------------------------------------

# What is the objective type?
        # objective_type == 0, objective_function better than peatsa_objective
        # objective_type == 1, objective_function better than polyfit = objective_fun(seed_value)
objective_type = 0

# How can PEATSA evaluate the objective value?
        # Input a string that can be evaluated. EMTG Mission data can be gathered using 'M.' and 
        # Mission option data is available using 'MO.'
objective_formula = 'M.total_deterministic_deltav'

# Is this a maximization or minimization problem? Enter 'max' or 'min'.
max_or_min = 'min'

# If PEATSA's objective_type == 0, what is the target objective value? I.E. for missed thrust, the 
        # target could be 60 days. Or for maximum mass, the target could be 100 
        # kg. In some cases the goal objective will be the same as the EMTG
        # objective_type, but it doesnt have to be. If a good target is not
        # known, then simply set it to something impossible to reach. I.e. for 
        # missed thrust, 1e100, or for minumum time, -1 days. The value here should 
        # have the same units as whatever the 'objective_formula' will evaluate to.
peatsa_goal_objective = 0

# Would you like to release constraints if PEATSA stalls? Enter a list of tuples. Each tuple has the form:
        # (string of constraint as a PyEMTG MissionOption,how many failed iterations to wait before loosening constraint, increment)
        # The first entry should be a variable in the PyEMTG MissionOption class, such as: 'MO.Journeys[-1].arrival_date_bounds[1]'
        # The second entry is how many times no PEATSA cases found improvement before the constraint will be loosened.
        # Finally enter the increment to change the constraint each time this event is triggered.
constraint_walking = [
                     ]


# ------------------------------------------ Sorting Options -------------------------------------------

# What makes each case in this PEATSA run unique? Write a list of strings
        # that will be evaluated for each case to give it a unique fingerprint. 
        # For data from an EMTG mission option, use 'MO.'. Do not use data
        # from an EMTG mission. Cases need fingerprints even if they dont converge.
        # Do not put any formulas that are already in the seed_criteria list.
fingerprint = [
               "MO.launch_window_open_date"
              ]


# -------------------------------------- Case Re-creation Options --------------------------------------

# CURRENTLY BUGGY, DO NOT USE.
            # Do you want to find seeds in parallel? (0 or 1).
            # Uses find_seeds_in_parallel_nCores number of cores.
find_seeds_in_parallel = 0

# If find_seeds_in_parallel == 1, this is how many parallel processes are used to find seeds. Integer >= 1.
find_seeds_in_parallel_nCores = 1

# CURRENTLY BUGGY, DO NOT USE.
            # Do you want to write seeds in parallel? (0 or 1).
write_seeds_in_parallel = 0

# If write_seeds_in_parallel == 1, this is how many parallel processes are used to write seeds. Integer >= 1.
write_seeds_in_parallel_nCores = 1

# Are there additional folders that contain emtg cases that can be used as seeds
        # for the current PEATSA run, but are not valid cases. For example, a group of results
        # that seem feasible to PEATSA, but the EMTG feasibility criteria changed internally to
        # EMTG, so all the old cases need to be tossed, but the initial guesses are still useful.
        # Enter a list of tuples:
        # (folder_type,full_folder_path)
        # if folder_type == 0, then the folder is a peatsa results folder with many iteration
        #                      subdirectories, all of which will be included
        # if folder_type == 1, then the folder is an individual results folder and only its contents
        #                      will be parsed, with all subdirectories ignored.
seed_folders = [
               ]

# When fresh starting, if you have an external seed folder, would you like to use it before the first iteration runs?
seed_from_seed_folders_on_fresh_start = False

# Should cases be allowed to seed themselves? Or only their neighbors? Enter 1 or 0
allow_cases_to_seed_themselves = 1

# How do you want to generate initial guesses?
        # if initial_guess_type == 'x_days', then the number of timesteps in the mission itself
        #                                    will be modified such that each timestep is ~ x days
        #                                    (use new_timestep to define x)
        # if initial_guess_type == 'n_steps', then the number of timesteps in the mission will
        #                                     be kept the same or modified to be n steps
        #                                     (use new_timestep to define x)
        # if initial_guess_type == 'copy', then PEATSA will just copy the initial guess over without
        #                                  attempting to interpolate
initial_guess_type = 'copy'

# If a case doesnt have a seed, rather than running unseeded, should it keep using whatever is in its trialX,
        # presumably the seed from reference mission. Or do you want to turn off seed_MBH and actually start without an initial guess?
turn_off_seed_MBH_when_unseeded = 0

# If a case doesnt have a seed, rather than running unseeded, should it seed from the end of its last run,
        # even though it didnt find something feasible at the end of the last run?
seed_from_infeasible_self_when_unseeded = 0

# Would you like cases to be able to be seeded from infeasible seeds if no feasible seeds were found?
seed_from_infeasible_cases_if_no_feasible_seeds_found = False

# Should cases that don't meet the goal criteria be considered as potential seeds?
        # WARNING: potential seed cases will be evaluated against peatsa_goal_objective
        # even if PEATSA_objective_type == 1. 
        # Enter a 1 or 0.
seed_from_cases_that_havent_met_target = 1

# Flag to override the individual seed selection criterias specified in 'seed_criteria'.
        # Instead, only one total seed will be selected regardless of which seed criteria was used. 
        # Enter a 1 or 0.
only_one_seed_from_best_in_all_seed_directions = 1

# How can PEATSA evaluate the seed criteria?
        # Input a list of tuples. Each tuple is: 
        # (a string that can be evaluated, max negative seed range, max positive seed range, seed selection criteria)
        # For the evaluation string, EMTG Mission data can be gathered using
        # 'MO.'. Don't use mission data, because even failed cases need to know how to
        # pick a seed.
        # For the max seed range, answer if there should be a maximum range for seeding, and if so, how long?
        # max_negative_seed_range > 0, there is no maximum range in the negative direction
        # max_negative_seed_range <= 0, this is the maximum range, in the units of whatever
        #       seed_criteria evaluates in 
        # max_positive_seed_range < 0, there is no maximum range in the negative direction
        # max_positive_seed_range >= 0, this is the maximum range, in the units of whatever
        #       seed_criteria evaluates in 
        # For the seed selection criteria, specify how PEATSA should select new seeds?
        # seed_selection_criteria == 0, seed from the closest potential seed
        # seed_selection_criteria == 1, seed from all potential seeds
        # seed_selection_criteria == 2, seed from only the best potential seed
        # seed_selection_criteria == 3, seed from random potential seed
        # seed_selection_criteria == 4, seed from a curve fit of all potential seeds.
        #                               Each decision vector variable is a curve fit
        #                               of a curve fit of the corresponding decision vector values
        #                               from all of the seeds. Uses polyfit_order for the curve fit order.
seed_criteria = [
                 ('MO.launch_window_open_date', -1.1, 1.1, 2)
                ]

# Should cases not be run until a seed can be found for them?
        # Note: they will be run for the 0th iteration. But in the reseeding process,
        # If a seed cant be found for it, then it will not be put back on the queue
        # Enter 1 or 0
wait_until_seeded = 0

# Enter a list of strings that will can be evaluated as a conditional statement. Each condition will filter out cases that should be ignored. These conditions
        # are applied using 'or' logic, assuming that no cases will be re-run unless these conditions are met.
        # If left empty, the normal conditions apply. If applied, the cases will still be compared to the PEATSA objective, so
        # cases that meet the 'only_rerun_if', still might not be re-run if they also meet the stopping criteria.
        # Note: this does not remove cases from the PEATSA run entirely. It only adds
        # a condition, which if not met, will prevent cases from being re-run. But
        # thoses filtered cases will still appear in csv summary files.
        # Use 'M.' to access EMTG mission object data and 'MO.' to access EMTG mission options data.
        # Example only_rerun_if = ['MO.launch_window_open_date > 64000','M.spacecraft_dry_mass < 100 and M.total_flight_time_years > 30']
only_rerun_if = [
                ]

# Are there any options that you want to override from the default options?
        # Enter a tuple of: 
        # (string of condition of when to use the option, string of formula for option)
        # Use MO for the mission options object.
        # Example: override_options = [('1','MO.launch_window_open_date = 23432'),('len(MO.Journeys) > 3','MO.Journeys[3].destination_list[1] = 2')]
override_options = [
                    ('1', 'MO.objective_type=0'),
                    ('1', 'MO.run_inner_loop=3'),
                    ('1', 'MO.snopt_feasibility_tolerance = 1e-5'),
                    ('1', 'MO.snopt_optimality_tolerance = 1e-4'),
                    ('1', 'MO.quiet_basinhopping=1'),
                    ('1', 'MO.quiet_NLP=1'),
                    ('1', 'MO.snopt_max_run_time=10'),
                    ('1', 'MO.MBH_max_run_time=11'),
                    ('1', 'MO.short_output_file_names=1'),
                    ('1', 'MO.print_only_non_default_options=0'),
                    ('1', 'MO.MBH_RNG_seed=-1'),
                    ('1', 'MO.skip_first_nlp_run=0'),
                    ('1', 'MO.generate_forward_integrated_ephemeris=0'),
                    ('1', 'MO.output_maneuver_and_target_spec_files=0'),
                    ('1', 'MO.universe_folder = "/Utilities/emtg/OSIRIS-REx_Tutorial/universe/"'),
                    ('1', 'MO.HardwarePath = "/Utilities/emtg/OSIRIS-REx_Tutorial/HardwareModels/"')
                   ]


# -------------------------------------- Post-processing Options ---------------------------------------

# By default, PEATSA will parse previous results using just one core, in parallel.
        # Should peatsa instead parse using all allocated local cores? Enter a 1 or 0.
parse_in_parallel = 0

# Do you want any extra data printed to each iteration's csv summary?
        # If not, leave this as an empty list
        # If so, add tuples of:
        # ('column header','string of formula to be evaluated')
        # For data from an EMTG mission, use 'M.', and for an EMTG mission option, use 'MO.'
extra_csv_column_definitions = [
                                ('Gregorian_Launch_Date', 'M.Journeys[0].missionevents[0].GregorianDate'),
                                ('Julian_Launch_Date', 'M.Journeys[0].missionevents[0].JulianDate'),
                                ('Launch_C3_[km^2/s^2]', 'M.Journeys[0].missionevents[0].C3'),
                                ('Delivered_Dry_mass[kg]', 'M.final_mass_including_propellant_margin'),
                                ('Total_Deterministic_DV[km/s]', 'M.total_deterministic_deltav')
                               ]

# Should the default plots be generated?
        # For trade studies, each option being traded will be plotted 
        # on the x axis vs. the 'objective_value' on the y-axis. Note: this 
        # currently only works for trade study type 0 files.
        # For missed thrust, missed thrust event date will be plotted on the x axis vs departure coast possible on the y-axis.
generate_default_plots = 0

# Would you like to generate custom plots using the built-in plotter? If so, enter a list of paths as strings to csv files which
        # contain options for a plot
built_in_plotter_files = [
                         ]

# Would you like to call a custom post processing script? Enter a path as a string.
        # This script will be called after each iteration. It can generate extra summary files, in addition to
        # the .csv file generated already, or extra plots, in addition to the defaults or any made with the built-in plotting routine.
        # Enter the path of a python module. The file must:
        #      1)  Be valid python code, ending with .py
        #      2)  Define a method called 'PEATSAPostProcess'
        #      3)  PEATSAPostProcess must take a PEATSAmenu object as its first input
        #      4)  PEATSAPostProcess must take a python list of PEATSAbox objects as its second input
        # It is not necessary, but data files can be saved to 'PEATSAorder.docs_dir' and plots can be saved to 'PEATSAorder.images_dir'
custom_post_processing_script = ''

# Would you like peatsa to move the latest results to another folder for easier ftp'ing?
        # If so, provide a full path to the desired destination. If not, leave an empty string.
move_results = ''


