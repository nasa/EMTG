import os
import JourneyOptionsV1 as JO
import platform

class MissionOptions(object):

    #************************************************************************************constructor
    def __init__(self, input_file_name='default.emtgopt'):
        self.filename = ""
        self.success = []
    
        #problem type    
        self.problem_type = 0
    
        #outer loop solver settings
        self.run_outerloop = 0 #whether or not to run the outer loop if false
        self.outerloop_popsize = 200 #population size
        self.outerloop_genmax = 40 #number of generations
        self.outerloop_tournamentsize = 4 #tournament size for selection
        self.outerloop_CR = 0.3 #math::crossover ratio
        self.outerloop_mu = 0.5 #mutation rate
        self.outerloop_stallmax = 20 #maximum number of stall generations
        self.outerloop_tolfit = 1.0e-4 #fitness tolerance
        self.outerloop_ntrials = 1 #how many times to run the outer loop
        self.outerloop_elitecount = 1 #how many elite individuals to retain
        self.outerloop_useparallel = 0 #whether or not to use the parallel outer-loop
        self.outerloop_warmstart = 0 #this will be the number of generations that have elapsed in the prior run of the GA
        self.outerloop_warm_population = "none"
        self.outerloop_warm_archive = "none"
        self.outerloop_reevaluate_full_population = 0
        self.quiet_outerloop = 1#if true, suppress all text outputs except error catches
    
        #outer loop selectable options settings
        self.outerloop_vary_spacecraft = 0
        self.outerloop_vary_power_system = 0
        self.outerloop_vary_chemical_fuel_tank_capacity = 0
        self.outerloop_vary_chemical_oxidizer_tank_capacity = 0
        self.outerloop_vary_electric_propellant_tank_capacity = 0
        self.outerloop_vary_electric_propulsion_system = 0
        self.outerloop_vary_number_of_electric_propulsion_systems = 0
        self.outerloop_vary_launch_vehicle = 0
        self.outerloop_vary_duty_cycle = 0
        self.outerloop_vary_launch_epoch = 0
        self.outerloop_vary_flight_time_upper_bound = 0
        self.outerloop_vary_departure_C3 = 0
        self.outerloop_vary_arrival_C3 = 0
        self.outerloop_vary_arrival_declination = 0
        self.outerloop_vary_final_journey_entry_interface_velocity = 0
        self.outerloop_vary_number_of_DSMs = 0
        self.outerloop_vary_phase_type = 0
        self.outerloop_vary_stop_after_journey = 0
        self.outerloop_spacecraft_choices = ['default.emtg_spacecraftopt']
        self.outerloop_power_system_choices = ['DefaultPowerSystem']
        self.outerloop_chemical_fuel_tank_capacity_choices = [1000.0]
        self.outerloop_chemical_oxidizer_tank_capacity_choices = [1000.0]
        self.outerloop_electric_propellant_tank_capacity_choices = [1000.0]
        self.outerloop_electric_propulsion_system_choices = ['NSTAR']
        self.outerloop_number_of_electric_propulsion_systems_choices = [1]
        self.outerloop_launch_vehicle_choices = ['Atlas_V_401']
        self.outerloop_duty_cycle_choices = [1.0]
        self.outerloop_launch_epoch_choices = [51545.0]
        self.outerloop_flight_time_upper_bound_choices = [365.25]
        self.outerloop_departure_C3_choices = [0]
        self.outerloop_arrival_C3_choices = [0]
        self.outerloop_arrival_declination_choices = [0]
        self.outerloop_final_journey_entry_interface_velocity_choices = [0]
        self.outerloop_phase_type_choices = [2, 6]
        self.outerloop_restrict_flight_time_lower_bound = 0
    
        #outer-loop variable number of journeys
        self.outerloop_vary_number_of_journeys = 0
        self.outerloop_maximum_number_of_variable_journeys = 8
        self.outerloop_insert_variable_journeys_after_fixed_journey_index = 0
        self.outerloop_inherit_variable_journey_settings_from_fixed_journey_index = 0
    
        #outer-loop point group settings
        self.outerloop_point_groups_number_to_score = [1]
        self.outerloop_point_groups_values = [1]
        self.outerloop_point_groups_members = [[1, 2, 3]]
    
        #outer-loop filter settings
        self.outerloop_filter_successive_destinations_on_inclination = 0
        self.outerloop_destination_inclination_filter_bandpass = 10
        self.outerloop_filter_by_groups = 0
        self.outerloop_group_filter_number_of_groups = 1
        self.outerloop_group_filter_lists = [[1, 2]]
    
        #outerloop objective settings
        self.outerloop_objective_function_choices = [2, 6]
    
        #inner loop solver settings
        self.NLP_solver_type = 0
        self.NLP_solver_mode = 1
        self.quiet_NLP = 1
        self.print_NLP_movie_frames = 0
        self.enable_NLP_chaperone = 1
        self.quiet_basinhopping = 0
        self.ACE_feasible_point_finder = 1
        self.MBH_max_not_improve = 100000
        self.MBH_max_trials = 100000
        self.MBH_max_run_time = 600
        self.MBH_max_step_size = 1.0
        self.MBH_hop_distribution = 1
        self.MBH_time_hop_probability = 0.05
        self.MBH_Pareto_alpha = 1.4
        self.MBH_always_write_archive = 0
        self.MBH_archive_state_vector = 0
        self.MBH_write_every_improvement = 0
        self.MBH_RNG_seed = 0
        self.snopt_feasibility_tolerance = 1.0e-5
        self.snopt_optimality_tolerance = 1.0e-6
        self.NLP_max_step = 1.0
        self.snopt_major_iterations = 8000
        self.snopt_minor_iterations = 500
        self.snopt_max_run_time = 15
        self.seed_MBH = 0
        self.skip_first_nlp_run = 0
        self.NLP_stop_on_goal_attain = 0
        self.NLP_objective_goal = 0.0
    
        self.objective_journey = 0
    
        #problem settings set by the user
    
        #ephemeris data
        self.ephemeris_source = 1
        self.SPICE_leap_seconds_kernel = "naif0009.tls"
        self.SPICE_reference_frame_kernel = "pck00010.tpc"
        self.universe_folder = "../Universe/"
    
        #integrator options
        self.integrator_tolerance = 1.0e-8
        self.propagatorType = 0
        self.integratorType = 0
        self.integration_time_step_size = 10.0 * 86400.0

        #state representation options       
        self.PSFBstateRepresentation = 2#Cartesian, SphericalRADEC, SphericalAZFPA
        self.PeriapseBoundaryStateRepresentation = 1#Cartesian, SphericalRADEC, SphericalAZFPA
    
        #low thrust solver parameters
        self.num_timesteps = 10 #number of timesteps per phase
        self.spiral_segments = 10
    
        #impulsive thrust solver parameters
        self.maximum_number_of_impulses_per_phase = 1
    
        #vehicle parameters
        self.maximum_mass = 1000 #the maximum possible mass of the spacecraft (negative number means use LV max)
        self.allow_initial_mass_to_vary = 0 #flag on whether or not the solver can choose the initial mass (make the spacecraft wet mass lighter)
        self.LV_margin = 0.0 #launch vehicle margin
        self.LV_adapter_mass = 0.0 #launch vehicle adapter mass (kg)
        self.IspLT = 1 #specific impulse of the engine used for low-thrust maneuvers
        self.IspLT_minimum = 1 #minimum Isp for VSI systems
        self.IspChem = 10000 #specific impulse of the engine used for impulsive maneuvers
        self.Thrust = 0.01 #thrust of the spacecraft, in Newtons
    
        self.engine_type = 0
        self.number_of_electric_propulsion_systems = 1
        self.throttle_logic_mode = 0
        self.throttle_sharpness = 100.0
        self.engine_duty_cycle = 1.0 #percentage of time that engine can operate
        self.duty_cycle_type = 0
        self.thrust_scale_factor = 1.0#percentage of thrust that goes through the CG
        self.power_at_1_AU = 10.0 #in kW
        self.power_source_type = 1
        self.solar_power_model_type = 0#0: classic Sauer, 1: polynomial
        self.solar_power_gamma = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #coefficients for solar panels
        self.power_margin = 0.0 #for propulsion, as a fraction
        self.engine_input_thrust_coefficients = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.engine_input_mass_flow_rate_coefficients = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.engine_input_power_bounds = [1.0, 5.0]
        self.user_defined_engine_efficiency = 0.7
        self.spacecraft_power_coefficients = [0.0, 0.0, 0.0]
        self.spacecraft_power_model_type = 0
        self.power_decay_rate = 0.0 #percent per year
    
        #TCMs
        self.TCM_Isp = 200.0
        self.TCM_post_launch = 0.0
        self.TCM_pre_flyby = 0.0
        self.TCM_maneuver_fraction = 0.0
    
        #ACS
        self.trackACS = 0
        self.ACS_kg_per_day = 0.0
    
        #spacecraft object fields
        self.SpacecraftModelInput = 2 #SpacecraftModelInputType::AssembleFromMissionOptions
        self.HardwarePath = "C:/Utilities/HardwareModels/"
        self.ThrottleTableFile = "NEXT_TT11_NewFrontiers_EOL_1_3_2017.ThrottleTable"
        self.LaunchVehicleLibraryFile = "NLSII_April2017.emtg_launchvehicleopt"
        self.PowerSystemsLibraryFile = "default.emtg_powersystemsopt"
        self.PropulsionSystemsLibraryFile = "4_18_2017.emtg_propulsionsystemopt"
        self.SpacecraftOptionsFile = "default.emtg_spacecraftopt"
        self.LaunchVehicleKey = "Atlas_V_401"
        self.PowerSystemKey = "5kW_basic"
        self.ElectricPropulsionSystemKey = "NSTAR"
        self.ChemicalPropulsionSystemKey = "DefaultChemicalPropulsionSystem"
        
        #minimum dry mass constraint and related parameters
        self.final_mass_constraint_bounds = [0.0, 0.0] #in kg
        self.constrain_final_mass = 0
        self.constrain_dry_mass = 0
    
        #propellant tank constraints
        self.enable_electric_propellant_tank_constraint = False
        self.maximum_electric_propellant = 1000.0
        self.electric_propellant_margin = 0.0
        self.enable_chemical_propellant_tank_constraint = False
        self.maximum_chemical_fuel = 1000.0
        self.maximum_chemical_oxidizer = 1000.0
        self.bipropellant_mixture_ratio = 0.925
        self.chemical_propellant_margin = 0.0
        
        
        #perturbation settings
        self.perturb_SRP = 0
        self.perturb_thirdbody = 0
        self.perturb_J2 = 0
        self.spacecraft_area = 1 #in m^2
        self.coefficient_of_reflectivity = 1
    
        #global problem settings
        self.number_of_journeys = 1
        self.max_phases_per_journey = 8
        self.include_initial_impulse_in_cost = 0
        self.global_timebounded = 1#0: unbounded, 1: bounded total time (note that the global arrival date bound is by definition the same as the last journey"s arrival date bound and is not duplicated
        self.launch_window_open_date = 60000.0#MJD
        self.total_flight_time_bounds = [0, 1000.0]#[2] days
        self.objective_type = 2 #0: minimum deltaV, 1: minimum time, #2: maximum final mass
        self.DLA_bounds = [-28.5, 28.5] #DLA in degrees
        self.RLA_bounds = [-2880.0, 2880.0]
        self.mission_name = "default"
        self.mission_type = 2
        self.forced_post_launch_coast = 0.0
        self.forced_pre_flyby_coast = 0.0
        self.forced_post_flyby_coast = 0.0
        self.waypoint_file_path = './banana.ephemeris'
        self.covariance_file_path = './dbanana_dbanana.covariance'
        
        #array of JourneyOptions objects
        self.Journeys = []
        self.ActiveJourney = 0
        
        #output format settings
        self.output_file_frame = 1
        self.post_mission_wait_time = 0.0
        self.override_working_directory = 0;
        self.forced_working_directory = "..//EMTG_v9_Results"
        self.override_mission_subfolder = 0;
        self.forced_mission_subfolder = "mission_subfolder"
        self.short_output_file_names = 0
        self.generate_forward_integrated_ephemeris = 0#0 :no, 1: yes
        self.add_control_switch_line_to_ephemeris = 0
        self.append_mass_to_ephemeris_output = 0
        self.append_control_to_ephemeris_output = 0
        self.append_thrust_to_ephemeris_output = 0
        self.append_mdot_to_ephemeris_output = 0
        self.append_Isp_to_ephemeris_output = 0
        self.append_number_of_active_engines_to_ephemeris_output = 0
        self.append_active_power_to_ephemeris_output = 0
        self.append_throttle_level_to_ephemeris_output = 0
        self.background_mode = 0 #0: no, 1: yes
        self.output_dormant_journeys = 0
        self.output_SPK = 0
        self.output_STMs = 0
        self.output_maneuver_and_target_spec_files = 0
    
        #manual outer-loop control code
        self.stop_after_journey = 32767
    
        #manual inner-loop control code
        self.check_derivatives = 0
        self.run_inner_loop = 1
        self.trialX = []
    
        #maneuver constraint definitions
        self.ManeuverConstraintDefinitions = []
    
        #boundary constraint definitions
        self.BoundaryConstraintDefinitions = []

        #phase distance constraint definitions
        self.PhaseDistanceConstraintDefinitions = []
        
        self.user_data = {}
        
        #now run some stuff
        self.parse_options_file(input_file_name)
        

    #************************************************************************************parse options file
    def parse_options_file(self, input_file_name):
        #Step 1: open the file

        self.filename = input_file_name

        if os.path.isfile(self.filename):
            inputfile = open(input_file_name, "r")
            self.success = 1
        else:
            print("Unable to open", input_file_name, "EMTG Error")
            self.success = 0
            return
        
        maneuver_constraint_line_flag = 0
        boundary_constraint_line_flag = 0
        distance_constraint_line_flag = 0
        perturb_line_flag = 0
        point_group_members_flag = 0
        group_filter_flag = 0
        flyby_choice_line_flag = 0
        destination_choice_line_flag = 0
        arrival_type_choice_line_flag = 0
        enable_periapse_burn_line_flag = 0
        trialX_line_flag = 0
        bad_option_flag = 0


        #Step 2: scan through the file
        linenumber = 0
        for line in inputfile:
            #strip off the newline character
            line = line.rstrip("\n\r ")
            linenumber = linenumber + 1

            if line.strip('\r') != "":
                if line[0] != "#":
                    #this is an active line, so it is space delimited
                    linecell = [entry.rstrip(" \r\n") for entry in line.split(" ")]
                                                        
                    choice = linecell[0].strip(' ')

                    if bad_option_flag > 0:
                        print('continuing to strip bad option')

                    elif perturb_line_flag > 0:
                        self.Journeys[perturb_line_flag - 1].journey_perturbation_bodies = []
                        for x in linecell:
                            self.Journeys[perturb_line_flag - 1].journey_perturbation_bodies.append(int(x))
                        perturb_line_flag += 1

                    elif group_filter_flag > 0:
                        group_filter_flag += 1
                        temp_filter_group = []
                        for entry in linecell[1:]:
                            temp_filter_group.append(int(entry))
                        self.outerloop_group_filter_lists.append(temp_point_group)

                    elif point_group_members_flag > 0:
                        point_group_members_flag += 1
                        temp_point_group = []
                        for entry in linecell[1:]:
                            temp_point_group.append(int(entry))
                        self.outerloop_point_groups_members.append(temp_point_group)

                    elif flyby_choice_line_flag > 0:
                        self.Journeys[flyby_choice_line_flag - 1].outerloop_journey_flyby_sequence_choices = []
                        for entry in linecell[1:]:
                            self.Journeys[flyby_choice_line_flag - 1].outerloop_journey_flyby_sequence_choices.append(int(entry))
                        flyby_choice_line_flag += 1

                    elif destination_choice_line_flag > 0:
                        self.Journeys[destination_choice_line_flag - 1].outerloop_journey_destination_choices = []
                        for entry in linecell[1:]:
                            self.Journeys[destination_choice_line_flag - 1].outerloop_journey_destination_choices.append(int(entry))
                        destination_choice_line_flag += 1

                    elif arrival_type_choice_line_flag > 0:
                        self.Journeys[arrival_type_choice_line_flag - 1].outerloop_journey_arrival_type_choices = []
                        for entry in linecell[1:]:
                            self.Journeys[arrival_type_choice_line_flag - 1].outerloop_journey_arrival_type_choices.append(int(entry))
                        arrival_type_choice_line_flag += 1

                    elif maneuver_constraint_line_flag > 0:
                        if linecell[0] == "END_MANEUVER_CONSTRAINT_BLOCK":
                             maneuver_constraint_line_flag = 0
                        else:   
                            self.ManeuverConstraintDefinitions.append(linecell[0])
                            maneuver_constraint_line_flag += 1
                            
                    elif boundary_constraint_line_flag > 0:
                        if linecell[0] == "END_BOUNDARY_CONSTRAINT_BLOCK":
                             boundary_constraint_line_flag = 0
                        else:   
                            self.BoundaryConstraintDefinitions.append(linecell[0])
                            boundary_constraint_line_flag += 1

                    elif distance_constraint_line_flag > 0:
                        if linecell[0] == "END_PHASE_DISTANCE_CONSTRAINT_BLOCK":
                             distance_constraint_line_flag = 0
                        else:   
                            self.PhaseDistanceConstraintDefinitions.append(linecell[0])
                            distance_constraint_line_flag += 1

                    elif trialX_line_flag > 0:
                        if "END_TRIALX" in linecell[0]:
                            trialX_line_flag = 0
                            self.ConvertDecisionVector()
                        else:
                            commalinecell = line.split(',')
                            self.trialX.append(commalinecell)
                            trialX_line_flag += 1

                    elif choice == "problem_type":
                        self.problem_type = int(linecell[1])

                    #outer-loop solver settings
                    elif choice == "run_outerloop":
                        self.run_outerloop = int(linecell[1])
                    elif choice == "outerloop_popsize":
                        self.outerloop_popsize = int(linecell[1])
                    elif choice == "outerloop_genmax":
                        self.outerloop_genmax = int(linecell[1])
                    elif choice == "outerloop_tournamentsize":
                        self.outerloop_tournamentsize = int(linecell[1])
                    elif choice == "outerloop_CR":
                        self.outerloop_CR = float(linecell[1])
                    elif choice == "outerloop_mu":
                        self.outerloop_mu = float(linecell[1])
                    elif choice == "outerloop_stallmax":
                        self.outerloop_stallmax = int(linecell[1])
                    elif choice == "outerloop_tolfit":
                        self.outerloop_tolfit = float(linecell[1])
                    elif choice == "outerloop_ntrials":
                        self.outerloop_ntrials = int(linecell[1])
                    elif choice == "outerloop_elitecount":
                        self.outerloop_elitecount = int(linecell[1])
                    elif choice == "outerloop_useparallel":
                        self.outerloop_useparallel = int(linecell[1])
                    elif choice == "outerloop_warmstart":
                        self.outerloop_warmstart = int(linecell[1])
                    elif choice == "outerloop_warm_population":
                        self.outerloop_warm_population = linecell[1]
                    elif choice == "outerloop_warm_archive":
                        self.outerloop_warm_archive = linecell[1]
                    elif choice == "outerloop_reevaluate_full_population":
                        self.outerloop_reevaluate_full_population = int(linecell[1])
                    elif choice == "quiet_outerloop":
                        self.quiet_outerloop = int(linecell[1])

                    #outer loop selectable options settings
                    elif choice == "outerloop_vary_spacecraft":
                        self.outerloop_vary_spacecraft = int(linecell[1])
                    elif choice == "outerloop_vary_power_system":
                        self.outerloop_vary_power_system = int(linecell[1])
                    elif choice == "outerloop_vary_chemical_fuel_tank_capacity":
                        self.outerloop_vary_chemical_fuel_tank_capacity = int(linecell[1])    
                    elif choice == "outerloop_vary_chemical_oxidizer_tank_capacity":
                        self.outerloop_vary_chemical_oxidizer_tank_capacity = int(linecell[1])          
                    elif choice == "outerloop_vary_electric_propellant_tank_capacity":
                        self.outerloop_vary_electric_propellant_tank_capacity = int(linecell[1])                 
                    elif choice == "outerloop_vary_launch_epoch":
                        self.outerloop_vary_launch_epoch = int(linecell[1])
                    elif choice == "outerloop_vary_flight_time_upper_bound":
                        self.outerloop_vary_flight_time_upper_bound = int(linecell[1])
                    elif choice == "outerloop_restrict_flight_time_lower_bound":
                        self.outerloop_restrict_flight_time_lower_bound = int(linecell[1])
                    elif choice == "outerloop_vary_electric_propulsion_system":
                        self.outerloop_vary_electric_propulsion_system = int(linecell[1])
                    elif choice == "outerloop_vary_number_of_electric_propulsion_systems":
                        self.outerloop_vary_number_of_electric_propulsion_systems = int(linecell[1])
                    elif choice == "outerloop_vary_duty_cycle":
                        self.outerloop_vary_duty_cycle = int(linecell[1])
                    elif choice == "outerloop_vary_launch_vehicle":
                        self.outerloop_vary_launch_vehicle = int(linecell[1])
                    elif choice == "outerloop_vary_departure_C3":
                        self.outerloop_vary_departure_C3 = int(linecell[1])
                    elif choice == "outerloop_vary_arrival_C3":
                        self.outerloop_vary_arrival_C3 = int(linecell[1])
                    elif choice == "outerloop_vary_arrival_declination":
                        self.outerloop_vary_arrival_declination = int(linecell[1])
                    elif choice == "outerloop_vary_final_journey_entry_interface_velocity":
                        self.outerloop_vary_final_journey_entry_interface_velocity = int(linecell[1])
                    elif choice == "outerloop_vary_number_of_DSMs":
                        self.outerloop_vary_number_of_DSMs = int(linecell[1])
                    elif choice == "outerloop_vary_phase_type":
                        self.outerloop_vary_phase_type = int(linecell[1])
                    elif choice == "outerloop_vary_journey_destination":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].outerloop_vary_journey_destination = int(linecell[j+1])
                    elif choice == "outerloop_vary_journey_flyby_sequence":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].outerloop_vary_journey_flyby_sequence = int(linecell[j+1])
                    elif choice == "outerloop_vary_journey_arrival_type":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].outerloop_vary_journey_arrival_type = int(linecell[j+1])
                    elif choice == "outerloop_vary_stop_after_journey":
                        self.outerloop_vary_stop_after_journey = int(linecell[1])
                    elif choice == "outerloop_spacecraft_choices":
                        self.outerloop_spacecraft_choices = []
                        for x in linecell[1:]:
                            self.outerloop_spacecraft_choices.append(x)
                    elif choice == "outerloop_power_system_choices":
                        self.outerloop_power_system_choices = []
                        for x in linecell[1:]:
                            self.outerloop_power_system_choices.append(x)
                    elif choice == "outerloop_chemical_fuel_tank_capacity_choices":  
                        self.outerloop_chemical_fuel_tank_capacity_choices = []
                        for x in linecell[1:]:
                            self.outerloop_chemical_fuel_tank_capacity_choices.append(float(x))
                    elif choice == "outerloop_chemical_oxidizer_tank_capacity_choices":  
                        self.outerloop_chemical_oxidizer_tank_capacity_choices = []
                        for x in linecell[1:]:
                            self.outerloop_chemical_oxidizer_tank_capacity_choices.append(float(x))
                    elif choice == "outerloop_electric_propellant_tank_capacity_choices":  
                        self.outerloop_electric_propellant_tank_capacity_choices = []
                        for x in linecell[1:]:
                            self.outerloop_electric_propellant_tank_capacity_choices.append(float(x))
                    elif choice == "outerloop_duty_cycle_choices":
                        self.outerloop_duty_cycle_choices = []
                        for x in linecell[1:]:
                            self.outerloop_duty_cycle_choices.append(float(x))
                    elif choice == "outerloop_launch_epoch_choices":
                        self.outerloop_launch_epoch_choices = []
                        for x in linecell[1:]:
                            self.outerloop_launch_epoch_choices.append(float(x))
                    elif choice == "outerloop_flight_time_upper_bound_choices":
                        self.outerloop_flight_time_upper_bound_choices = []
                        for x in linecell[1:]:
                            self.outerloop_flight_time_upper_bound_choices.append(float(x))
                    elif choice == "outerloop_electric_propulsion_system_choices":
                        self.outerloop_electric_propulsion_system_choices = []
                        for x in linecell[1:]:
                            self.outerloop_electric_propulsion_system_choices.append(x)
                    elif choice == "outerloop_number_of_electric_propulsion_systems_choices":
                        self.outerloop_number_of_electric_propulsion_systems_choices = []
                        for x in linecell[1:]:
                            self.outerloop_number_of_electric_propulsion_systems_choices.append(int(float(x)))
                    elif choice == "outerloop_launch_vehicle_choices":
                        self.outerloop_launch_vehicle_choices = []
                        for x in linecell[1:]:
                            self.outerloop_launch_vehicle_choices.append(x)
                    elif choice == "outerloop_departure_C3_choices":
                        self.outerloop_departure_C3_choices = []
                        for x in linecell[1:]:
                            self.outerloop_departure_C3_choices.append(float(x))
                    elif choice == "outerloop_arrival_C3_choices":
                        self.outerloop_arrival_C3_choices = []
                        for x in linecell[1:]:
                            self.outerloop_arrival_C3_choices.append(float(x))
                    elif choice == "outerloop_arrival_declination_choices":
                        self.outerloop_arrival_declination_choices = []
                        for x in linecell[1:]:
                            self.outerloop_arrival_declination_choices.append(float(x))
                    elif choice == "outerloop_final_journey_entry_interface_velocity_choices":
                        self.outerloop_final_journey_entry_interface_velocity_choices = []
                        for x in linecell[1:]:
                            self.outerloop_final_journey_entry_interface_velocity_choices.append(float(x))
                    elif choice == "outerloop_phase_type_choices":
                        self.outerloop_phase_type_choices = []
                        for x in linecell[1:]:
                            self.outerloop_phase_type_choices.append(int(float(x)))
                    elif choice == "outerloop_journey_flyby_sequence_choices":
                        flyby_choice_line_flag = 1
                    elif choice == "outerloop_journey_destination_choices":
                        destination_choice_line_flag = 1
                    elif choice == "outerloop_journey_maximum_number_of_flybys":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].outerloop_journey_maximum_number_of_flybys = int(float(linecell[j+1]))
                    elif choice == "outerloop_journey_arrival_type_choices":
                        arrival_type_choice_line_flag = 1

                    #outer-loop variable number of journeys
                    elif choice == "outerloop_vary_number_of_journeys":
                        self.outerloop_vary_number_of_journeys = int(float(linecell[1]))
                    elif choice == "outerloop_maximum_number_of_variable_journeys":
                        self.outerloop_maximum_number_of_variable_journeys = int(float(linecell[1]))
                    elif choice == "outerloop_insert_variable_journeys_after_fixed_journey_index":
                        self.outerloop_insert_variable_journeys_after_fixed_journey_index = int(float(linecell[1]))
                    elif choice == "outerloop_inherit_variable_journey_settings_from_fixed_journey_index":
                            self.outerloop_inherit_variable_journey_settings_from_fixed_journey_index = int(float(linecell[1]))

                    #outer-loop point groups settings
                    elif choice == "outerloop_point_groups_values":
                        self.outerloop_point_groups_values = []
                        for x in linecell[1:]:
                            self.outerloop_point_groups_values.append(int(x))

                        #start reading point group members
                        self.outerloop_point_groups_members = []
                        point_group_members_flag = 1

                    elif choice == "outerloop_point_groups_number_to_score":
                        self.outerloop_point_groups_number_to_score = []
                        for x in linecell[1:]:
                            self.outerloop_point_groups_number_to_score.append(int(x))

                    #outer-loop filter settings
                    elif choice == "outerloop_filter_successive_destinations_on_inclination":
                        self.outerloop_filter_successive_destinations_on_inclination = int(linecell[1])
                    elif choice == "outerloop_destination_inclination_filter_bandpass":
                        self.outerloop_destination_inclination_filter_bandpass = float(linecell[1])
                    elif choice == "outerloop_filter_by_groups":
                        self.outerloop_filter_by_groups = int(linecell[1])
                    elif choice == "outerloop_group_filter_lists":
                        self.outerloop_group_filter_number_of_groups = int(linecell[1])

                        #start reading point group members
                        self.outerloop_group_filter_lists = []
                        group_filter_flag = 1

                    #outerloop objective settings
                    elif choice == "outerloop_objective_function_choices":
                        self.outerloop_objective_function_choices = []
                        for x in linecell[1:]:
                            self.outerloop_objective_function_choices.append(int(x))

                    
                    #inner loop solver settings
                    elif choice == "NLP_solver_type":
                        self.NLP_solver_type = int(linecell[1])
                    elif choice == "NLP_solver_mode":
                        self.NLP_solver_mode = int(linecell[1])
                    elif choice == "ACE_feasible_point_finder":
                        self.ACE_feasible_point_finder = int(linecell[1])
                    elif choice == "quiet_NLP":
                        self.quiet_NLP = int(linecell[1])
                    elif choice == "print_NLP_movie_frames":
                        self.print_NLP_movie_frames = int(linecell[1])
                    elif choice == "enable_NLP_chaperone":
                        self.enable_NLP_chaperone = int(linecell[1])
                    elif choice == "NLP_stop_on_goal_attain":
                        self.NLP_stop_on_goal_attain = int(linecell[1])
                    elif choice == "NLP_objective_goal":
                        self.NLP_objective_goal = float(linecell[1])
                    elif choice == "quiet_basinhopping":
                        self.quiet_basinhopping = int(linecell[1])
                    elif choice ==  "MBH_max_not_improve":
                        self.MBH_max_not_improve = int(linecell[1])
                    elif choice ==  "MBH_max_trials":
                        self.MBH_max_trials = int(linecell[1])
                    elif choice ==  "MBH_max_run_time":
                        self.MBH_max_run_time = int(linecell[1])
                    elif choice ==  "MBH_max_step_size":
                        self.MBH_max_step_size = float(linecell[1])
                    elif choice == "MBH_hop_distribution":
                        self.MBH_hop_distribution = int(linecell[1])
                    elif choice == "MBH_Pareto_alpha":
                        self.MBH_Pareto_alpha = float(linecell[1])
                    elif choice == "MBH_time_hop_probability":
                        self.MBH_time_hop_probability = float(linecell[1])
                    elif choice == "MBH_always_write_archive":
                        self.MBH_always_write_archive = int(float(linecell[1]))
                    elif choice == "MBH_archive_state_vector":
                        self.MBH_archive_state_vector = int(float(linecell[1]))
                    elif choice == "MBH_write_every_improvement":
                        self.MBH_write_every_improvement = int(float(linecell[1]))
                    elif choice == "MBH_RNG_seed":
                        self.MBH_RNG_seed = int(float(linecell[1]))
                    elif choice ==  "snopt_feasibility_tolerance":
                        self.snopt_feasibility_tolerance = float(linecell[1])
                    elif choice ==  "snopt_optimality_tolerance":
                        self.snopt_optimality_tolerance = float(linecell[1])
                    elif choice == "NLP_max_step":
                        self.NLP_max_step = float(linecell[1])
                    elif choice ==  "snopt_major_iterations":
                        self.snopt_major_iterations = int(linecell[1])
                    elif choice ==  "snopt_minor_iterations":
                        self.snopt_minor_iterations = int(linecell[1])
                    elif choice ==  "snopt_max_run_time":
                        self.snopt_max_run_time = int(linecell[1])
                    elif choice ==  "seed_MBH":
                        self.seed_MBH = int(linecell[1])
                    elif choice == "skip_first_nlp_run":
                        self.skip_first_nlp_run = int(linecell[1])

                    elif choice == "objective_journey":
                        self.objective_journey = int(linecell[1])

                    #problem settings set by the user
                    elif choice ==  "ephemeris_source":
                        self.ephemeris_source = int(linecell[1])
                    elif choice ==  "SPICE_leap_seconds_kernel":
                        self.SPICE_leap_seconds_kernel = linecell[1]
                    elif choice ==  "SPICE_reference_frame_kernel":
                        self.SPICE_reference_frame_kernel = linecell[1]
                    elif choice ==  "universe_folder":
                        self.universe_folder = linecell[1]

                    elif choice == "integrator_tolerance":
                        self.integrator_tolerance = float(linecell[1])
                    elif choice == "propagatorType":
                        self.propagatorType = int(linecell[1])
                    elif choice == "integratorType":
                        self.integratorType = int(linecell[1])
                    elif choice == "integration_time_step_size":
                        self.integration_time_step_size = float(linecell[1])

                    #state representation settings
                    elif choice == "PSFBstateRepresentation":
                        self.PSFBstateRepresentation = int(linecell[1])
                    elif choice == "PeriapseBoundaryStateRepresentation":
                        self.PeriapseBoundaryStateRepresentation = int(linecell[1])

                    #low thrust solver parameters
                    elif choice ==  "num_timesteps":
                        self.num_timesteps = int(linecell[1])
                    elif choice == "spiral_segments":
                        self.spiral_segments = int(linecell[1])

                    #impulsive thrust solver parameters
                    elif choice == "maximum_number_of_impulses_per_phase":
                        self.maximum_number_of_impulses_per_phase = int(linecell[1])
                    elif choice == "enable_periapse_burns":
                        print("Deprecated option: enable_periapse_burns. Use journey_enable_periapse_burns instead. If journey_enable_periapse_burns was not set also, periapse burns will default to off for all flybys.")

                    #vehicle parameters
                    elif choice == "maximum_mass":
                        self.maximum_mass = float(linecell[1])
                    elif choice == "enable_maximum_propellant_mass_constraint":
                        print("Deprecated option: enable_maximum_propellant_mass_constraint. Use individual electric and chemical constraints instead.")
                    elif choice == "maximum_propellant_mass":
                        print("Deprecated option: maximum_propellant_mass. Use individual electric and chemical constraints instead.")
                    elif choice == "LV_margin":
                        self.LV_margin = float(linecell[1])
                    elif choice == "LV_adapter_mass":
                        self.LV_adapter_mass = float(linecell[1])
                    elif choice == "IspLT":
                        self.IspLT = float(linecell[1])
                    elif choice == "IspLT_minimum":
                        self.IspLT_minimum = float(linecell[1])
                    elif choice == "IspChem":
                        self.IspChem = float(linecell[1])
                    elif choice == "Thrust":
                        self.Thrust = float(linecell[1])
                    elif choice == "engine_type":
                        self.engine_type = int(linecell[1])
                    elif choice == "number_of_electric_propulsion_systems":
                        self.number_of_electric_propulsion_systems = int(linecell[1])
                    elif choice == "throttle_logic_mode":
                        self.throttle_logic_mode = int(linecell[1])
                        if self.throttle_logic_mode > 1:
                            self.throttle_logic_mode = 1
                    elif choice == "throttle_sharpness":
                        self.throttle_sharpness = float(linecell[1])
                    elif choice == "engine_duty_cycle":
                        self.engine_duty_cycle = float(linecell[1])
                    elif choice == "duty_cycle_type":
                        self.duty_cycle_type = int(float(linecell[1]))
                    elif choice == "thrust_scale_factor":
                        self.thrust_scale_factor = float(linecell[1])
                    elif choice == "power_at_1_AU":
                        self.power_at_1_AU = float(linecell[1])
                    elif choice == "power_source_type":
                        self.power_source_type = int(linecell[1])
                    elif choice == "solar_power_model_type":
                        self.solar_power_model_type = int(linecell[1])
                    elif choice == "solar_power_gamma":
                        self.solar_power_gamma = []
                        for x in linecell[1:]:
                            self.solar_power_gamma.append(float(x))
                        for k in range(len(self.solar_power_gamma),7):
                            self.solar_power_gamma.append(0.0)
                    elif choice == "power_margin":
                        self.power_margin = float(linecell[1])
                    elif choice == "power_decay_rate":
                        self.power_decay_rate = float(linecell[1])
                    elif choice == "spacecraft_power_coefficients":
                        self.spacecraft_power_coefficients = []
                        for x in linecell[1:]:
                            self.spacecraft_power_coefficients.append(float(x))
                    elif choice == "engine_input_thrust_coefficients":
                        self.engine_input_thrust_coefficients = []
                        for x in linecell[1:]:
                            self.engine_input_thrust_coefficients.append(float(x))
                        while len(self.engine_input_thrust_coefficients) < 7:
                            self.engine_input_thrust_coefficients.append(0.0)
                    elif choice == "engine_input_mass_flow_rate_coefficients":
                        self.engine_input_mass_flow_rate_coefficients = []
                        for x in linecell[1:]:
                            self.engine_input_mass_flow_rate_coefficients.append(float(x))
                        while len(self.engine_input_mass_flow_rate_coefficients) < 7:
                            self.engine_input_mass_flow_rate_coefficients.append(0.0)
                    elif choice == "engine_input_power_bounds":
                        self.engine_input_power_bounds = []
                        for x in linecell[1:]:
                            self.engine_input_power_bounds.append(float(x))
                    elif choice == "user_defined_engine_efficiency":
                        self.user_defined_engine_efficiency = float(linecell[1])
                    elif choice == "spacecraft_power_model_type":
                        self.spacecraft_power_model_type = int(linecell[1])
                    elif choice == "allow_initial_mass_to_vary":
                        self.allow_initial_mass_to_vary = int(linecell[1])

                    elif choice == "TCM_Isp":
                        self.TCM_Isp = float(linecell[1])
                    elif choice == "TCM_post_launch":
                        self.TCM_post_launch = float(linecell[1])
                    elif choice == "TCM_pre_flyby":
                        self.TCM_pre_flyby = float(linecell[1])
                    elif choice == "TCM_maneuver_fraction":
                        self.TCM_maneuver_fraction = float(linecell[1])

                    elif choice == "trackACS":
                        self.trackACS = int(float(linecell[1]))
                    elif choice == "ACS_kg_per_day":
                        self.ACS_kg_per_day = float(linecell[1])

                    elif choice == "constrain_dry_mass":
                        self.constrain_dry_mass = int(float(linecell[1]))
                    elif choice == "constrain_final_mass":
                        self.constrain_final_mass = int(float(linecell[1]))
                    elif choice == "constrain_dry_mass":
                        self.constrain_dry_mass = int(float(linecell[1]))
                    elif choice == "final_mass_constraint_bounds":
                        self.final_mass_constraint_bounds = [float(linecell[1]), float(linecell[2])]
                    elif choice == "propellant_margin":
                        print("Deprecated option: propellant_margin. Use individual electric and chemical constraints instead.")

                    #propellant tank constraints
                    elif choice == "enable_electric_propellant_tank_constraint":
                        self.enable_electric_propellant_tank_constraint = int(linecell[1])
                    elif choice == "maximum_electric_propellant":
                        self.maximum_electric_propellant = float(linecell[1])
                    elif choice == "electric_propellant_margin":
                        self.electric_propellant_margin = float(linecell[1])
                    elif choice == "enable_chemical_propellant_tank_constraint":
                        self.enable_chemical_propellant_tank_constraint = int(linecell[1])
                    elif choice == "maximum_chemical_fuel":
                        self.maximum_chemical_fuel = float(linecell[1])
                    elif choice == "maximum_chemical_oxidizer":
                        self.maximum_chemical_oxidizer = float(linecell[1])
                    elif choice == "bipropellant_mixture_ratio":
                        self.bipropellant_mixture_ratio = float(linecell[1])
                    elif choice == "chemical_propellant_margin":
                        self.chemical_propellant_margin = float(linecell[1])

                    #spacecraft object-related fields
                    elif choice == "SpacecraftModelInput":
                        self.SpacecraftModelInput = int(linecell[1])
                    elif choice == "HardwarePath":
                        self.HardwarePath = linecell[1].rstrip(" \r\n")
                    elif choice == "ThrottleTableFile":
                        self.ThrottleTableFile = linecell[1].rstrip(" \r\n")
                    elif choice == "LaunchVehicleLibraryFile":
                        self.LaunchVehicleLibraryFile = linecell[1].rstrip(" \r\n")
                    elif choice == "PowerSystemsLibraryFile":
                        self.PowerSystemsLibraryFile = linecell[1].rstrip(" \r\n")
                    elif choice == "PropulsionSystemsLibraryFile":
                        self.PropulsionSystemsLibraryFile = linecell[1].rstrip(" \r\n")
                    elif choice == "SpacecraftOptionsFile":
                        self.SpacecraftOptionsFile = linecell[1].rstrip(" \r\n")
                    elif choice == "LaunchVehicleKey":
                        self.LaunchVehicleKey = linecell[1].rstrip(" \r\n")
                    elif choice == "PowerSystemKey":
                        self.PowerSystemKey = linecell[1].rstrip(" \r\n")
                    elif choice == "ElectricPropulsionSystemKey":
                        self.ElectricPropulsionSystemKey = linecell[1].rstrip(" \r\n")
                    elif choice == "ChemicalPropulsionSystemKey":
                        self.ChemicalPropulsionSystemKey = linecell[1].rstrip(" \r\n")

                    elif choice == "number_of_journeys":
                        self.number_of_journeys = int(linecell[1])
                        for j in range(0, self.number_of_journeys):
                            temp_JourneyOptions = JO.JourneyOptions(self.mission_type)
                            self.Journeys.append(temp_JourneyOptions)

                    elif choice == "max_phases_per_journey":
                        self.max_phases_per_journey = int(linecell[1])
                    elif choice == "destination_list":
                        for j in range(0,self.number_of_journeys):
                            self.Journeys[j].destination_list = [int(linecell[2*j+1]), int(linecell[2*j+2])]
                    elif choice == "include_initial_impulse_in_cost":
                        self.include_initial_impulse_in_cost = int(float(linecell[1]))
                    elif choice == "global_timebounded":
                        self.global_timebounded = int(linecell[1])
                    elif choice == "launch_window_open_date":
                        self.launch_window_open_date = float(linecell[1])
                    elif choice == "total_flight_time_bounds":
                        self.total_flight_time_bounds = [float(linecell[1]), float(linecell[2])]
                    elif choice == "objective_type":
                        self.objective_type = int(linecell[1])
                    elif choice == "waypoint_file_path":
                        self.waypoint_file_path = linecell[1]
                    elif choice == "covariance_file_path":
                        self.covariance_file_path = linecell[1]
                    elif choice == "DLA_bounds":
                        self.DLA_bounds = [float(linecell[1]), float(linecell[2])]
                    elif choice == "RLA_bounds":
                        self.RLA_bounds = [float(linecell[1]), float(linecell[2])]
                    elif choice == "mission_name":
                        self.mission_name = linecell[1]
                    elif choice == "mission_type":
                        self.mission_type = int(linecell[1])
                    elif choice == "forced_post_launch_coast":
                        self.forced_post_launch_coast = float(linecell[1])
                    elif choice == "forced_pre_flyby_coast":
                        self.forced_pre_flyby_coast = float(linecell[1])
                    elif choice == "forced_post_flyby_coast":
                        self.forced_post_flyby_coast = float(linecell[1])

                    #parse all of the journey options and load them
                    #into the journey objects
                    elif choice == "journey_override_num_steps":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_override_num_steps = int(float(linecell[j+1]))

                    elif choice == "journey_number_of_steps":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_number_of_steps = int(float(linecell[j+1]))

                    elif choice == "journey_num_interior_control_points":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_num_interior_control_points = int(float(linecell[j+1]))

                    elif choice == "journey_freeze_decision_variables":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_freeze_decision_variables = int(float(linecell[j+1]))

                    elif choice == "journey_names":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_names = linecell[j+1]

                    elif choice == "journey_starting_mass_increment":
                        print("deprecated option: journey_starting_mass_increment. Use \'journey_fixed_starting_mass_increment\'")
                        
                    elif choice == "journey_fixed_ending_mass_increment":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_fixed_ending_mass_increment = float(linecell[j+1])
                            
                    elif choice == "journey_fixed_starting_mass_increment":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_fixed_starting_mass_increment = float(linecell[j+1])
                        
                    elif choice == "journey_minimum_starting_mass_increment":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_minimum_starting_mass_increment = float(linecell[j+1])
                            
                    elif choice == "journey_maximum_starting_mass_increment":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_maximum_starting_mass_increment = float(linecell[j+1])
                            
                    elif choice == "journey_variable_mass_increment":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_variable_mass_increment = int(float(linecell[j+1]))
                            
                    elif choice == "journey_constrain_initial_mass":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_constrain_initial_mass = int(linecell[j+1])
                    
                    elif choice == "journey_maximum_initial_mass":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_maximum_initial_mass = float(linecell[j+1])
                                    
                    elif choice == "journey_timebounded":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_timebounded = int(float(linecell[j+1]))

                    elif choice == "journey_bounded_departure_date":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_bounded_departure_date = int(float(linecell[j+1]))

                    elif choice == "journey_departure_date_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_date_bounds = [float(linecell[2*j+1]), float(linecell[2*j+2])]
                    
                    elif choice == "journey_wait_time_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_wait_time_bounds = [float(linecell[2*j+1]), float(linecell[2*j+2])]
                    
                    elif choice == "journey_flight_time_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_flight_time_bounds = [float(linecell[2*j+1]), float(linecell[2*j+2])]
                    
                    elif choice == "journey_arrival_date_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_date_bounds = [float(linecell[2*j+1]), float(linecell[2*j+2])]
                    
                    elif choice == "journey_initial_impulse_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_initial_impulse_bounds = [float(linecell[2*j+1]), float(linecell[2*j+2])]
                    
                    elif choice == "journey_departure_type":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_type = int(float(linecell[j+1]))
                    
                    elif choice == "journey_departure_class":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_class = int(float(linecell[j+1]))

                    elif choice == "journey_departure_ellipsoid_axes":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_ellipsoid_axes = [float(linecell[3*j+1]), float(linecell[3*j+2]), float(linecell[3*j+3])]
                    
                    elif choice == "journey_departure_elements_type":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_elements_type = int(float(linecell[j+1]))

                    elif choice == "journey_departure_elements_frame":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_elements_frame = int(float(linecell[j+1]))
                    
                    elif choice == "journey_departure_elements":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_elements = [float(linecell[6*j+1]), float(linecell[6*j+2]), float(linecell[6*j+3]), float(linecell[6*j+4]), float(linecell[6*j+5]), float(linecell[6*j+6])]
                    
                    elif choice == "journey_departure_elements_vary_flag":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_elements_vary_flag = [int(float(linecell[6*j+1])), int(float(linecell[6*j+2])), int(float(linecell[6*j+3])), int(float(linecell[6*j+4])), int(float(linecell[6*j+5])), int(float(linecell[6*j+6]))]
                    
                    elif choice == "journey_departure_elements_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_elements_bounds = [float(linecell[12*j+1]), float(linecell[12*j+2]), float(linecell[12*j+3]), float(linecell[12*j+4]), float(linecell[12*j+5]), float(linecell[12*j+6]), float(linecell[12*j+7]), float(linecell[12*j+8]), float(linecell[12*j+9]), float(linecell[12*j+10]), float(linecell[12*j+11]), float(linecell[12*j+12])]
                            
                    elif choice == "AllowJourneyFreePointDepartureToPropagate":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].AllowJourneyFreePointDepartureToPropagate = int(float(linecell[j+1]))
                            
                    elif choice == "journey_departure_elements_reference_epoch":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_departure_elements_reference_epoch = float(linecell[j+1])

                    elif choice == "journey_override_flyby_altitude_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_override_flyby_altitude_bounds = int(float(linecell[j + 1]))

                    elif choice == "journey_flyby_altitude_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_flyby_altitude_bounds = [float(linecell[2 * j + 1]), float(linecell[2 * j + 2])]
                                                
                    elif choice == "journey_override_PropagatorType":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_override_PropagatorType = int(float(linecell[j + 1]))

                    elif choice == "journey_PropagatorType":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_PropagatorType = int(float(linecell[j + 1]))

                    elif choice == "journey_override_integration_step_size":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_override_integration_step_size = int(float(linecell[j + 1]))

                    elif choice == "journey_integration_step_size":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_integration_step_size = float(linecell[j + 1])

                    elif choice == "journey_CoastPhaseMatchPointFraction":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_CoastPhaseMatchPointFraction = float(linecell[j + 1])

                    elif choice == "journey_CoastPhaseForwardIntegrationStepLength":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_CoastPhaseForwardIntegrationStepLength = float(linecell[j + 1])

                    elif choice == "journey_CoastPhaseBackwardIntegrationStepLength":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_CoastPhaseBackwardIntegrationStepLength = float(linecell[j + 1])
                    
                    elif choice == "journey_arrival_type":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_type = int(float(linecell[j+1]))
                    
                    elif choice == "journey_arrival_class":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_class = int(float(linecell[j+1]))

                    elif choice == "journey_arrival_ellipsoid_axes":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_ellipsoid_axes = [float(linecell[3*j+1]), float(linecell[3*j+2]), float(linecell[3*j+3])]

                    elif choice == "journey_PeriapseArrival_override_altitude":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_PeriapseArrival_override_altitude = int(float(linecell[j+1]))

                    elif choice == "journey_PeriapseArrival_altitude_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_PeriapseArrival_altitude_bounds = [float(linecell[2*j+1]), float(linecell[2*j+2])]
                    
                    elif choice == "journey_PeriapseDeparture_altitude_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_PeriapseDeparture_altitude_bounds = [float(linecell[2*j+1]), float(linecell[2*j+2])]
                    
                    elif choice == "journey_arrival_elements_type":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_elements_type = int(float(linecell[j+1]))

                    elif choice == "journey_arrival_elements_frame":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_elements_frame = int(float(linecell[j+1]))
                    
                    elif choice == "journey_arrival_elements":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_elements = [float(linecell[6*j+1]), float(linecell[6*j+2]), float(linecell[6*j+3]), float(linecell[6*j+4]), float(linecell[6*j+5]), float(linecell[6*j+6])]
                    
                    elif choice == "journey_arrival_elements_vary_flag":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_elements_vary_flag = [int(float(linecell[6*j+1])), int(float(linecell[6*j+2])), int(float(linecell[6*j+3])), int(float(linecell[6*j+4])), int(float(linecell[6*j+5])), int(float(linecell[6*j+6]))]
                    
                    elif choice == "journey_arrival_elements_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_elements_bounds = [float(linecell[12*j+1]), float(linecell[12*j+2]), float(linecell[12*j+3]), float(linecell[12*j+4]), float(linecell[12*j+5]), float(linecell[12*j+6]), float(linecell[12*j+7]), float(linecell[12*j+8]), float(linecell[12*j+9]), float(linecell[12*j+10]), float(linecell[12*j+11]), float(linecell[12*j+12])]
                    
                    elif choice == "AllowJourneyFreePointArrivalToPropagate":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].AllowJourneyFreePointArrivalToPropagate = int(float(linecell[j+1]))
                            
                    elif choice == "journey_arrival_elements_reference_epoch":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_elements_reference_epoch = float(linecell[j+1])

                    elif choice == "journey_central_body":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_central_body = linecell[j+1]

                    elif choice == "journey_final_velocity":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_final_velocity = [float(linecell[j*3+1]), float(linecell[j*3+2]), float(linecell[j*3+3])]

                    elif choice =="journey_arrival_declination_constraint_flag":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_declination_constraint_flag = int(float(linecell[j+1]))

                    elif choice == "journey_arrival_declination_bounds":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_arrival_declination_bounds = [float(linecell[j*2+1]), float(linecell[j*2+2])]

                    elif choice == "journey_escape_spiral_starting_radius":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_escape_spiral_starting_radius = float(linecell[j+1])
                    elif choice == "journey_capture_spiral_final_radius":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_capture_spiral_final_radius = float(linecell[j+1])
                    elif choice == "journey_forced_terminal_coast":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_forced_terminal_coast = float(linecell[j+1])
                    elif choice == "journey_forced_initial_coast":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_forced_initial_coast = float(linecell[j+1])

                    elif choice == "journey_end_deltav":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_end_deltav = float(linecell[j+1])
                    elif choice == "journey_end_propulsion_system":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_end_propulsion_system = int(linecell[j+1])
                    elif choice == "journey_end_TCM":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_end_TCM = float(linecell[j+1])

                    elif choice == "journey_override_duty_cycle":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_override_duty_cycle = int(linecell[j+1])
                    elif choice == "journey_duty_cycle":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_duty_cycle = float(linecell[j+1])
                                
                    #staging-related quantities
                    elif choice == "journey_stage_after_departure":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_stage_after_departure = int(linecell[j+1])
                    elif choice == "journey_stage_before_arrival":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_stage_before_arrival = int(linecell[j+1])
                    elif choice == "journey_stage_after_arrival":
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_stage_after_arrival = int(linecell[j+1])

                    #specialized maneuver constraints (epoch and magnitude)      
                    elif choice == "BEGIN_MANEUVER_CONSTRAINT_BLOCK":
                        self.ManeuverConstraintDefinitions = []
                        maneuver_constraint_line_flag = 1

                    
                    #specialized boundary constraints     
                    elif choice == "BEGIN_BOUNDARY_CONSTRAINT_BLOCK":
                        self.BoundaryConstraintDefinitions = []
                        boundary_constraint_line_flag = 1

                        
                    #phase distance constraints    
                    elif choice == "BEGIN_PHASE_DISTANCE_CONSTRAINT_BLOCK":
                        self.PhaseDistanceConstraintDefinitions = []
                        distance_constraint_line_flag = 1
                                                 
                    #perturbation-related quantities    
                    elif choice == "perturb_SRP":
                        self.perturb_SRP = int(linecell[1])
                    elif choice == "perturb_thirdbody":
                        self.perturb_thirdbody = int(linecell[1])
                    elif choice == "perturb_J2":
                        self.perturb_J2 = int(linecell[1])
                    elif choice == "journey_perturbation_bodies":
                        #first get the number of perturbation bodies for each journey
                        for j in range(0, self.number_of_journeys):
                            self.Journeys[j].journey_number_of_perturbation_bodies = int(linecell[j+1])

                        #next read the bodies line for each journey
                        perturb_line_flag = 1
                    
                    elif choice == "spacecraft_area":
                        self.spacecraft_area = float(linecell[1])
                    elif choice == "coefficient_of_reflectivity":
                        self.coefficient_of_reflectivity = float(linecell[1])
                            
                    #output format settings
                    elif choice == "output_file_frame":
                        self.output_file_frame = int(linecell[1])
                    elif choice == "post_mission_wait_time":
                        self.post_mission_wait_time = float(linecell[1])
                    elif choice == "override_working_directory":
                        self.override_working_directory = int(linecell[1])
                    elif choice == "forced_working_directory":
                        self.forced_working_directory = linecell[1]
                    elif choice == "override_mission_subfolder":
                        self.override_mission_subfolder = int(linecell[1])
                    elif choice == "forced_mission_subfolder":
                        self.forced_mission_subfolder = linecell[1]
                    elif choice == "generate_forward_integrated_ephemeris":
                        self.generate_forward_integrated_ephemeris = int(linecell[1])
                    elif choice == "add_control_switch_line_to_ephemeris":
                        self.add_control_switch_line_to_ephemeris = int(linecell[1])
                    elif choice == "append_mass_to_ephemeris_output":
                        self.append_mass_to_ephemeris_output = int(linecell[1])
                    elif choice == "append_control_to_ephemeris_output":
                        self.append_control_to_ephemeris_output = int(linecell[1])
                    elif choice == "append_thrust_to_ephemeris_output":
                        self.append_thrust_to_ephemeris_output = int(linecell[1])
                    elif choice == "append_mdot_to_ephemeris_output":
                        self.append_mdot_to_ephemeris_output = int(linecell[1])
                    elif choice == "append_Isp_to_ephemeris_output":
                        self.append_Isp_to_ephemeris_output = int(linecell[1])
                    elif choice == "append_number_of_active_engines_to_ephemeris_output":
                        self.append_number_of_active_engines_to_ephemeris_output = int(linecell[1])
                    elif choice == "append_active_power_to_ephemeris_output":
                        self.append_active_power_to_ephemeris_output = int(linecell[1])
                    elif choice == "append_throttle_level_to_ephemeris_output":
                        self.append_throttle_level_to_ephemeris_output = int(linecell[1])
                    elif choice == "background_mode":
                        self.background_mode = int(linecell[1])
                    elif choice == "output_dormant_journeys":
                        self.output_dormant_journeys = int(linecell[1])
                    elif choice == "short_output_file_names":
                        self.short_output_file_names = int(linecell[1])
                    elif choice == "output_STMs":
                        self.output_STMs = int(linecell[1])
                    elif choice == "output_maneuver_and_target_spec_files":
                        self.output_maneuver_and_target_spec_files = int(linecell[1])

                    #trialX, sequence input, etc
                    elif choice == "check_derivatives":
                        self.check_derivatives = int(linecell[1])
                    elif choice == "run_inner_loop":
                        self.run_inner_loop = int(linecell[1])
                    elif choice == "sequence":
                        for j in range(0, self.number_of_journeys):
                            seq = [0] * self.max_phases_per_journey
                            for p in range (0, self.max_phases_per_journey):
                                seq[p] = int(linecell[self.max_phases_per_journey*j+p + 1])
                            
                            self.Journeys[j].sequence = seq
                            self.Journeys[j].phase_type = [0] * (self.max_phases_per_journey + 1)
                            self.Journeys[j].impulses_per_phase = [self.maximum_number_of_impulses_per_phase] * (self.max_phases_per_journey + 1)
                            self.Journeys[j].journey_enable_periapse_burns = [0] * (self.max_phases_per_journey)
                            self.Journeys[j].number_of_phases = sum(1 for x in seq if x > 0)

                    elif choice == "stop_after_journey":
                        self.stop_after_journey = int(linecell[1])
                       
                    elif choice == "phase_type":                        
                        for j in range(0, self.number_of_journeys):
                            iphase_type = [0] * (self.max_phases_per_journey + 1)
                            for p in range (0, self.max_phases_per_journey + 1):
                                iphase_type[p] = int(linecell[(self.max_phases_per_journey + 1) * j + p + 1])

                            self.Journeys[j].phase_type = iphase_type
                            
                    elif choice == "impulses_per_phase":                        
                        for j in range(0, self.number_of_journeys):
                            iimpulses = [0] * (self.max_phases_per_journey + 1)
                            for p in range (0, self.max_phases_per_journey + 1):
                                iimpulses[p] = int(linecell[(self.max_phases_per_journey + 1)*j+p + 1])

                            self.Journeys[j].impulses_per_phase = iimpulses
                            
                    elif choice == "journey_enable_periapse_burns":                        
                        for j in range(0, self.number_of_journeys):
                            iseq = [0] * (self.max_phases_per_journey)
                            for p in range (0, self.max_phases_per_journey):
                                iseq[p] = int(linecell[self.max_phases_per_journey*j+p + 1])

                            self.Journeys[j].journey_enable_periapse_burns = iseq
                            
                    elif choice == "BEGIN_TRIALX":
                        self.trialX = []
                        trialX_line_flag = 1
                    
                    elif choice == "user_data":
                        self.user_data = dict()
                        full_notes = line.lstrip("user_data").lstrip(" ").rstrip(" \r\n")
                        if ":" in full_notes or full_notes.replace(" ","") != "":
                            full_notes = full_notes.split(":")
                                                                                  
                            for note in full_notes:
                                var = note.lstrip('("').split(",")[0].rstrip('"')
                                val = eval(note.lstrip("(").lstrip(var + '"').lstrip(", ").rstrip(") "))
                                self.user_data.update({var:val})

                    #if option is not recognized
                    else:
                        errorstring = "Option not recognized: " + str(linecell[0]) + " on line " + str(linenumber) + ". Option will be stripped so that parsing can continue."
                        #print(errorstring)
                        bad_option_flag = 1
                else:
                    maneuver_constraint_line_flag = 0
                    boundary_constraint_line_flag = 0
                    distance_constraint_line_flag = 0
                    perturb_line_flag = 0
                    point_group_members_flag = 0
                    group_filter_flag = 0
                    sequence_line_flag = 0
                    phase_type_flag = 0
                    number_of_impulses_line_flag = 0
                    trialX_line_flag = 0
                    flyby_choice_line_flag = 0
                    destination_choice_line_flag = 0
                    arrival_type_choice_line_flag = 0
                    enable_periapse_burn_line_flag = 0
                    trialX_line_flag = 0
                    bad_option_flag = 0
            else:
                maneuver_constraint_line_flag = 0
                boundary_constraint_line_flag = 0
                distance_constraint_line_flag = 0
                perturb_line_flag = 0
                point_group_members_flag = 0
                group_filter_flag = 0
                sequence_line_flag = 0
                phase_type_flag = 0
                number_of_impulses_line_flag = 0
                trialX_line_flag = 0
                flyby_choice_line_flag = 0
                destination_choice_line_flag = 0
                arrival_type_choice_line_flag = 0
                enable_periapse_burn_line_flag = 0
                bad_option_flag = 0
        inputfile.close()

    def write_options_file(self, output_file_name):
        #first make some error-preventing correction
        if (self.IspChem < 1.0):
            self.IspChem = 1.0
            
        #first open the file for writing
        outputfile = open(output_file_name, "w")
        
        outputfile.write("##Options file for EMTG_v9\n")
        outputfile.write("\n")
            
        outputfile.write("##problem type\n")
        outputfile.write("#0: standard EMTG mission\n")
        outputfile.write("problem_type " + str(self.problem_type) + "\n")
        outputfile.write("\n")

        outputfile.write("##outer-loop solver settings\n")
        outputfile.write("#Do you want to run an outer-loop?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: Genetic algorithm (number of objective functions determines which GA to run)\n")
        outputfile.write("run_outerloop " + str(self.run_outerloop) + "\n")
        outputfile.write("#outer-loop population size\n")   
        outputfile.write("outerloop_popsize " + str(self.outerloop_popsize) + "\n")
        outputfile.write("#maximum number of outer-loop generations\n") 
        outputfile.write("outerloop_genmax " + str(self.outerloop_genmax) + "\n")
        outputfile.write("#tournament size for selection\n")    
        outputfile.write("outerloop_tournamentsize " + str(self.outerloop_tournamentsize) + "\n")
        outputfile.write("#crossover ratio\n")  
        outputfile.write("outerloop_CR " + str(self.outerloop_CR) + "\n")
        outputfile.write("#mutation rate\n")    
        outputfile.write("outerloop_mu " + str(self.outerloop_mu) + "\n")
        outputfile.write("#maximum number of stall generations\n")  
        outputfile.write("outerloop_stallmax " + str(self.outerloop_stallmax) + "\n")
        outputfile.write("#fitness tolerance for the outer-loop\n") 
        outputfile.write("outerloop_tolfit " + str(self.outerloop_tolfit) + "\n")
        outputfile.write("#how many elite individuals to retain\n") 
        outputfile.write("outerloop_elitecount " + str(self.outerloop_elitecount) + "\n")
        outputfile.write("#how many times to run the outer-loop\n") 
        outputfile.write("outerloop_ntrials " + str(self.outerloop_ntrials) + "\n")
        outputfile.write("#whether or not to use the parallel outer-loop\n")    
        outputfile.write("outerloop_useparallel " + str(self.outerloop_useparallel) + "\n")
        outputfile.write("#whether or not to perform an outer loop warm start\n")
        outputfile.write("outerloop_warmstart " + str(self.outerloop_warmstart) + "\n")
        outputfile.write("#Population file for outerloop warm start (set to none if not warm starting)\n")
        outputfile.write("outerloop_warm_population " + self.outerloop_warm_population + "\n")
        outputfile.write("#Archive file for outerloop warm start (set to none if not warm starting)\n")
        outputfile.write("outerloop_warm_archive " + self.outerloop_warm_archive + "\n")
        outputfile.write("#Re-evaluate the entire outerloop each generation? Otherwise read from the archive.\n")
        outputfile.write("outerloop_reevaluate_full_population " + str(self.outerloop_reevaluate_full_population) + "\n")
        outputfile.write("#Quiet outer-loop?\n")
        outputfile.write("quiet_outerloop " + str(self.quiet_outerloop) + "\n")
        outputfile.write("\n")

        outputfile.write("##inner-loop solver settings\n")
        outputfile.write("#NLP solver type\n")
        outputfile.write("#0: SNOPT\n")
        outputfile.write("#1: WORHP\n")
        outputfile.write("NLP_solver_type " + str(self.NLP_solver_type) + "\n")
        outputfile.write("#NLP solver mode\n")
        outputfile.write("#0: find feasible point only\n")
        outputfile.write("#1: find optimal solution\n")
        outputfile.write("#2: satisfy equality constraints\n")
        outputfile.write("NLP_solver_mode " + str(self.NLP_solver_mode) + "\n")
        outputfile.write("#Quiet NLP solver?\n")
        outputfile.write("quiet_NLP " + str(self.quiet_NLP) + "\n")
        outputfile.write("#Print NLP movie frames?\n")
        outputfile.write("print_NLP_movie_frames " + str(self.print_NLP_movie_frames) + "\n")
        outputfile.write("#Enable NLP chaperone?\n")
        outputfile.write("enable_NLP_chaperone " + str(self.enable_NLP_chaperone) + "\n")
        outputfile.write("#Stop NLP upon attaining objective goal?\n")
        outputfile.write("NLP_stop_on_goal_attain " + str(self.NLP_stop_on_goal_attain) + "\n")
        outputfile.write("#NLP objective goal\n")
        outputfile.write("NLP_objective_goal " + str(self.NLP_objective_goal) + "\n")
        outputfile.write("#Quiet MBH?\n")
        outputfile.write("quiet_basinhopping " + str(self.quiet_basinhopping) + "\n")
        outputfile.write("#Enable ACE feasible point finder?\n")
        outputfile.write("ACE_feasible_point_finder " + str(self.ACE_feasible_point_finder) + "\n")
        outputfile.write("#quantity Max_not_improve for MBH\n")
        outputfile.write("MBH_max_not_improve " + str(self.MBH_max_not_improve) + "\n")
        outputfile.write("#maximum number of trials for MBH\n")
        outputfile.write("MBH_max_trials " + str(self.MBH_max_trials) + "\n")
        outputfile.write("#maximum run time for MBH, in seconds\n")
        outputfile.write("MBH_max_run_time " + str(self.MBH_max_run_time) + "\n")
        outputfile.write("#maximum step size for uniform MBH, or scaling factor for Cauchy MBH\n")
        outputfile.write("MBH_max_step_size " + str(self.MBH_max_step_size) + "\n")
        outputfile.write("#MBH hop probabilty distribution\n")
        outputfile.write("#0: uniform\n")
        outputfile.write("#1: Cauchy\n")
        outputfile.write("#2: Pareto\n")
        outputfile.write("#3: Gaussian\n")
        outputfile.write("MBH_hop_distribution " +str(self.MBH_hop_distribution) + "\n")
        outputfile.write("#Pareto distribution alpha\n")
        outputfile.write("MBH_Pareto_alpha " + str(self.MBH_Pareto_alpha) + "\n");
        outputfile.write("#probability of MBH time hop operation\n")
        outputfile.write("MBH_time_hop_probability " + str(self.MBH_time_hop_probability) + "\n")
        outputfile.write("#should the MBH archive file always be written?\n")
        outputfile.write("MBH_always_write_archive " + str(self.MBH_always_write_archive) + "\n")
        outputfile.write("#include state vector in MBH archive file?\n")
        outputfile.write("MBH_archive_state_vector " + str(self.MBH_archive_state_vector) + "\n")
        outputfile.write("#Write every MBH improvement for later animation?\n")
        outputfile.write("MBH_write_every_improvement " + str(self.MBH_write_every_improvement) + "\n")
        outputfile.write("#MBH RNG seed (negative number means system clock)\n")
        outputfile.write("MBH_RNG_seed " + str(self.MBH_RNG_seed) + "\n")
        outputfile.write("#feasibility tolerance\n")
        outputfile.write("snopt_feasibility_tolerance " + str(self.snopt_feasibility_tolerance) + "\n")
        outputfile.write("#optimality tolerance\n")
        outputfile.write("snopt_optimality_tolerance " + str(self.snopt_optimality_tolerance) + "\n")
        outputfile.write("#NLP max step\n")
        outputfile.write("NLP_max_step " + str(self.NLP_max_step) + "\n")
        outputfile.write("#maximum number of major iterations for SNOPT\n")
        outputfile.write("snopt_major_iterations " + str(self.snopt_major_iterations) + "\n")
        outputfile.write("#maximum number of minor iterations for SNOPT\n")
        outputfile.write("snopt_minor_iterations " + str(self.snopt_minor_iterations) + "\n")
        outputfile.write("#Maximum run time, in seconds, for a single call to SNOPT\n")
        outputfile.write("snopt_max_run_time " + str(self.snopt_max_run_time) + "\n")
        outputfile.write("#Will MBH be seeded with an initial point? Otherwise MBH starts from a completely random point.\n")
        outputfile.write("seed_MBH " + str(self.seed_MBH) + "\n")
        outputfile.write("#Which journey to use for objective types that reference a specific journey?\n")
        outputfile.write("objective_journey " + str(self.objective_journey) + "\n")
        outputfile.write("#If seed_MBH is on, would you like to skip the first nlp run and hop right away?\n")
        outputfile.write("skip_first_nlp_run " + str(self.skip_first_nlp_run) + "\n")
        outputfile.write("\n")

        outputfile.write("##low-thrust solver parameters\n")    
        outputfile.write("#number of time steps per phase\n")   
        outputfile.write("num_timesteps " + str(self.num_timesteps) + "\n")
        outputfile.write("#Number of spiral time-segments\n")
        outputfile.write("spiral_segments " + str(self.spiral_segments) + "\n")
        outputfile.write("\n")

        outputfile.write("##impulsive-thrust solver parameters\n")
        outputfile.write("#maximum number of impulses for MGAnDSM\n")
        outputfile.write("maximum_number_of_impulses_per_phase " + str(self.maximum_number_of_impulses_per_phase) + "\n")
        outputfile.write("\n")

        outputfile.write("##ephemeris data\n")  
        outputfile.write("#ephemeris source\n") 
        outputfile.write("#0: static\n")    
        outputfile.write("#1: SPICE (default to static if no SPICE file supplied for a body)\n")    
        outputfile.write("#2: SplineEphem (default to static if no SPICE file supplied for a body)\n")    
        outputfile.write("ephemeris_source " + str(self.ephemeris_source) + "\n")
        outputfile.write("#Universe folder\n")
        outputfile.write("universe_folder " + str(self.universe_folder) + "\n")
        outputfile.write("#SPICE leap seconds kernel - required for SPICE to work\n")
        outputfile.write("SPICE_leap_seconds_kernel " + str(self.SPICE_leap_seconds_kernel) + "\n")
        outputfile.write("#SPICE_reference_frame_kernel\n")
        outputfile.write("SPICE_reference_frame_kernel " + str(self.SPICE_reference_frame_kernel) + "\n")
        outputfile.write("\n")

        outputfile.write("##integrator options\n")
        outputfile.write("#integration tolerance\n")
        outputfile.write("integrator_tolerance " + str(self.integrator_tolerance) + "\n")
        outputfile.write("#Propagator type\n")
        outputfile.write("#0: Keplerian propagator\n")
        outputfile.write("#1: Integrated propagator\n")
        outputfile.write("propagatorType " + str(int(self.propagatorType)) + "\n")
        outputfile.write("#Integrator type\n")
        outputfile.write("#0: rk7813M adaptive step\n")
        outputfile.write("#1: rk8 fixed step\n")
        outputfile.write("integratorType " + str(int(self.integratorType)) + "\n")
        outputfile.write("#Fixed-step integrator time step size (seconds)\n")
        outputfile.write("integration_time_step_size " + str(self.integration_time_step_size) + "\n")
        outputfile.write("\n")
        
        outputfile.write("##state representation options\n")
        outputfile.write("#state representation for periapse boundary conditions\n")
        outputfile.write("PeriapseBoundaryStateRepresentation " + str(self.PeriapseBoundaryStateRepresentation) + "\n")
        outputfile.write("#state representation for PSFB\n")
        outputfile.write("PSFBstateRepresentation " + str(self.PSFBstateRepresentation) + "\n")
        outputfile.write("\n")
            
        outputfile.write("##vehicle parameters\n")  
        outputfile.write("#the maximum possible mass in kg of the spacecraft (negative number means use LV max)\n") 
        outputfile.write("maximum_mass " + str(self.maximum_mass) + "\n")
        outputfile.write("#Launch vehicle margin (0.0 - 1.0)\n")
        outputfile.write("LV_margin " + str(self.LV_margin) + "\n")
        outputfile.write("#Launch vehicle adapter mass (kg)\n")
        outputfile.write("LV_adapter_mass " + str(self.LV_adapter_mass) + "\n")
        outputfile.write("\n")

        outputfile.write("##parameters that are only relevant for missions that use chemical propulsion\n") 
        outputfile.write("##dummy values should be used if the mission does not use chemical propulsion but are not strictly necessary\n")  
        outputfile.write("#specific impulse in seconds of the engine used for impulsive maneuvers\n")   
        outputfile.write("IspChem " + str(self.IspChem) + "\n")
        outputfile.write("\n")
            
        outputfile.write("##parameters that are only relevant for missions that use low-thrust\n")  
        outputfile.write("##dummy values should be used if the mission does not use low-thrust but are not strictly necessary\n")   
        outputfile.write("#specific impulse in seconds of the engine used for low-thrust maneuvers\n")  
        outputfile.write("#for VSI systems, this represents maximum Isp\n")
        outputfile.write("IspLT " + str(self.IspLT) + "\n")
        outputfile.write("#minimum Isp for VSI systems\n")
        outputfile.write("IspLT_minimum " + str(self.IspLT_minimum) + "\n")
        outputfile.write("#thrust of the spacecraft low-thrust motor, in Newtons\n")    
        outputfile.write("Thrust " + str(self.Thrust) + "\n")
        outputfile.write("#low-thrust engine type\n")
        outputfile.write("#0: fixed thrust/Isp\n")
        outputfile.write("#1: constant Isp, efficiency, EMTG computes input power\n")
        outputfile.write("#2: choice of power model, constant efficiency, EMTG chooses Isp\n")
        outputfile.write("#3: choice of power model, constant efficiency and Isp\n")
        outputfile.write("#4: continuously-varying specific impulse\n")
        outputfile.write("#5: custom thrust and mass flow rate polynomial\n")
        outputfile.write("#6: NSTAR\n")
        outputfile.write("#7: XIPS-25\n")
        outputfile.write("#8: BPT-4000 High-Isp\n")
        outputfile.write("#9: BPT-4000 High-Thrust\n")
        outputfile.write("#10: BPT-4000 Ex-High-Isp\n")
        outputfile.write("#11: NEXT high-Isp v9\n")
        outputfile.write("#12: VASIMR (argon, using analytical model)\n")
        outputfile.write("#13: Hall Thruster (Xenon, using analytical model)\n")
        outputfile.write("#14: NEXT high-ISP v10\n")
        outputfile.write("#15: NEXT high-thrust v10\n")
        outputfile.write("#16: BPT-4000 MALTO\n")
        outputfile.write("#17: NEXIS Cardiff 8-15-201\n")
        outputfile.write("#18: H6MS Cardiff 8-15-2013\n")
        outputfile.write("#19: BHT20K Cardiff 8-16-2013\n")
        outputfile.write("#20: Aerojet HiVHAC EM\n")
        outputfile.write("#21: 13 kW STMD Hall high-Isp (not available in open-source)\n")
        outputfile.write("#22: 13 kW STMD Hall high-thrust (not available in open-source)\n")
        outputfile.write("#23: NEXT TT11 High-Thrust\n")
        outputfile.write("#24: NEXT TT11 High-Isp\n")
        outputfile.write("#25: NEXT TT11 Expanded Throttle Table\n")
        outputfile.write("#26: 13 kW STMD Hall high-Isp 10-1-2014 (not available in open-source)\n")
        outputfile.write("#27: 13 kW STMD Hall medium-thrust 10-1-2014 (not available in open-source)\n")
        outputfile.write("#28: 13 kW STMD Hall high-thrust 10-1-2014 (not available in open-source)\n")
        outputfile.write("#29: 2D Throttle table\n")
        outputfile.write("#30: 1D Throttle table high-thrust\n")
        outputfile.write("#31: 1D Throttle table high-Isp\n")
        outputfile.write("#32: 2D polynomial fit\n")
        outputfile.write("engine_type " + str(self.engine_type) + "\n")
        outputfile.write("#Custom engine thrust coefficients (T = A + BP + C*P^2 + D*P^3 + E*P^4 + G*P^5 + H*P^6)\n")
        outputfile.write("engine_input_thrust_coefficients")
        for k in range(0,7):
            outputfile.write(" " + str(self.engine_input_thrust_coefficients[k]))
        outputfile.write("\n")
        outputfile.write("#Custom engine mass flow rate coefficients (mdot = A + BP + C*P^2 + D*P^3 + E*P^4 + G*P^5 + H*P^6)\n")
        outputfile.write("engine_input_mass_flow_rate_coefficients")
        for k in range(0,7):
            outputfile.write(" " + str(self.engine_input_mass_flow_rate_coefficients[k]))
        outputfile.write("\n")
        outputfile.write("#Custom engine lower and upper bounds on input power (per engine, in kW)\n")
        outputfile.write("engine_input_power_bounds " + str(self.engine_input_power_bounds[0]) + " " + str(self.engine_input_power_bounds[1]) + "\n")
        outputfile.write("#Custom engine input efficiency\n")
        outputfile.write("user_defined_engine_efficiency " + str(self.user_defined_engine_efficiency) + "\n")
        outputfile.write("#number of low-thrust propulsion systems\n")
        outputfile.write("number_of_electric_propulsion_systems " + str(self.number_of_electric_propulsion_systems) + "\n")
        outputfile.write("#Throttle logic mode\n")
        outputfile.write("#0: maximum number of thrusters\n")
        outputfile.write("#1: minimum number of thrusters\n")
        outputfile.write("throttle_logic_mode " + str(self.throttle_logic_mode) + "\n")
        outputfile.write("#Throttle sharpness (higher means more precise, lower means smoother)\n")
        outputfile.write("throttle_sharpness " + str(self.throttle_sharpness) + "\n")
        outputfile.write("#engine duty cycle [0,1]\n")
        outputfile.write("engine_duty_cycle " + str(self.engine_duty_cycle) + "\n")
        outputfile.write("#duty cycle type\n")
        outputfile.write("#0: averaged\n")
        outputfile.write("#1: realistic\n")
        outputfile.write("duty_cycle_type " + str(self.duty_cycle_type) + "\n")
        outputfile.write("#Thrust scale factor, for simulating cosine loss [0, 1]\n")
        outputfile.write("thrust_scale_factor " + str(self.thrust_scale_factor) + "\n")
        outputfile.write("#electrical power available at 1 AU (kW)\n")
        outputfile.write("power_at_1_AU " + str(self.power_at_1_AU) + "\n")
        outputfile.write("#power source type, 0: fixed, 1: solar\n")
        outputfile.write("power_source_type " + str(self.power_source_type) + "\n")
        outputfile.write("#Solar power model type\n")
        outputfile.write("#0: classic Sauer model\n")
        outputfile.write("#1: polynomial (0th order on the left)\n")
        outputfile.write("solar_power_model_type " + str(self.solar_power_model_type) + "\n")
        outputfile.write("#solar power coefficients gamma_1 through gamma_5\n")
        outputfile.write("#if all gamma = [1 0 0 0 0], then solar power is a simple 1/r^2\n")
        outputfile.write("solar_power_gamma")
        for k in range(0,7):
            outputfile.write(" " + str(self.solar_power_gamma[k]))
        outputfile.write("\n")
        outputfile.write("#Power margin (for thrusters, as a fraction)\n")
        outputfile.write("power_margin " + str(self.power_margin) + "\n")
        outputfile.write("#Power system decay rate (per year)\n")
        outputfile.write("power_decay_rate " + str(self.power_decay_rate) + "\n")
        outputfile.write("#spacecraft power coefficients A, B, and C\n")
        outputfile.write("#represent the power requirements of the spacecraft at a distance r from the sun\n")
        outputfile.write("#i.e. heaters, communications, etc\n")
        outputfile.write("spacecraft_power_coefficients " + str(self.spacecraft_power_coefficients[0]) + " " + str(self.spacecraft_power_coefficients[1]) + " " + str(self.spacecraft_power_coefficients[2]) + "\n")
        outputfile.write("#spacecraft power model type\n")
        outputfile.write("#0: P_sc = A + B/r + C/r^2\n")
        outputfile.write("#1: P_sc = A if P > A, A + B(C - P) otherwise\n")
        outputfile.write("spacecraft_power_model_type " + str(self.spacecraft_power_model_type) + "\n")
        outputfile.write("#Allow initial mass to vary, up to maximum possible mass? (only relevant for MGALT and FBLT)\n")
        outputfile.write("allow_initial_mass_to_vary " + str(self.allow_initial_mass_to_vary) + "\n")

        outputfile.write("#Isp for TCMs (s)\n")
        outputfile.write("TCM_Isp " + str(self.TCM_Isp) + "\n")
        outputfile.write("#Magnitude of post-launch TCM (km/s)\n")
        outputfile.write("TCM_post_launch " + str(self.TCM_post_launch) + "\n")
        outputfile.write("#Magnitude of pre-encounter TCMs (km/s)\n")
        outputfile.write("TCM_pre_flyby " + str(self.TCM_pre_flyby) + "\n")
        outputfile.write("#Magnitude of post-DSM TCMs as a fraction of DSM magnitude\n")
        outputfile.write("TCM_maneuver_fraction " + str(self.TCM_maneuver_fraction) + "\n")

        
        outputfile.write("#Track ACS propellant consumption?\n")
        outputfile.write("trackACS " + str(int(self.trackACS)) + "\n")
        outputfile.write("#ACS propellant use per day (kg)\n")
        outputfile.write("ACS_kg_per_day " + str(self.ACS_kg_per_day) + "\n")

        outputfile.write("#Final mass constraint bounds\n")
        outputfile.write("final_mass_constraint_bounds " + str(self.final_mass_constraint_bounds[0]) + " " + str(self.final_mass_constraint_bounds[1]) + "\n")
        outputfile.write("#Constrain final mass?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("constrain_final_mass " + str(int(self.constrain_final_mass)) + "\n")
        outputfile.write("#Constrain dry mass?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("constrain_dry_mass " + str(int(self.constrain_dry_mass)) + "\n")
        outputfile.write("#Enable maximum electric propulsion propellant constraint?\n")
        outputfile.write("enable_electric_propellant_tank_constraint " + str(int(self.enable_electric_propellant_tank_constraint)) + "\n")
        outputfile.write("#Maximum electric propulsion propellant (kg)\n")
        outputfile.write("maximum_electric_propellant " + str(self.maximum_electric_propellant) + "\n")
        outputfile.write("#Electric propulsion propellant margin\n")
        outputfile.write("electric_propellant_margin " + str(self.electric_propellant_margin) + "\n")
        outputfile.write("#Enable maximum chemical propulsion propellant constraint?\n")
        outputfile.write("enable_chemical_propellant_tank_constraint " + str(int(self.enable_chemical_propellant_tank_constraint)) + "\n")
        outputfile.write("#Maximum chemical fuel (kg)\n")
        outputfile.write("maximum_chemical_fuel " + str(self.maximum_chemical_fuel) + "\n")
        outputfile.write("#Maximum chemical oxidizer (kg)\n")
        outputfile.write("maximum_chemical_oxidizer " + str(self.maximum_chemical_oxidizer) + "\n")
        outputfile.write("#Bipropellant mixture ratio\n")
        outputfile.write("bipropellant_mixture_ratio " + str(self.bipropellant_mixture_ratio) + "\n")
        outputfile.write("#Chemical propulsion propellant margin\n")
        outputfile.write("chemical_propellant_margin " + str(self.chemical_propellant_margin) + "\n")
        outputfile.write("\n")

        outputfile.write("##Spacecraft Object Settings\n")
        outputfile.write("#Spacecraft object input type\n")
        outputfile.write("#0: Assemble from libraries\n")
        outputfile.write("#1: Read .emtg_spacecraftoptions file\n")
        outputfile.write("#2: Assemble from missionoptions object\n")
        outputfile.write("SpacecraftModelInput " + str(self.SpacecraftModelInput) + "\n")
        outputfile.write("HardwarePath " + self.HardwarePath + "\n")
        outputfile.write("ThrottleTableFile " + self.ThrottleTableFile + "\n")
        outputfile.write("LaunchVehicleLibraryFile " + self.LaunchVehicleLibraryFile + "\n")
        outputfile.write("PowerSystemsLibraryFile " + self.PowerSystemsLibraryFile + "\n")
        outputfile.write("PropulsionSystemsLibraryFile " + self.PropulsionSystemsLibraryFile + "\n")
        outputfile.write("SpacecraftOptionsFile " + self.SpacecraftOptionsFile + "\n")
        outputfile.write("LaunchVehicleKey " + self.LaunchVehicleKey + "\n")
        outputfile.write("PowerSystemKey " + self.PowerSystemKey + "\n")
        outputfile.write("ElectricPropulsionSystemKey " + self.ElectricPropulsionSystemKey + "\n")
        outputfile.write("ChemicalPropulsionSystemKey " + self.ChemicalPropulsionSystemKey + "\n")
        outputfile.write("\n")

        outputfile.write("##Global problem settings\n")
        outputfile.write("#mission name\n")
        outputfile.write("mission_name " + str(self.mission_name) + "\n")
        outputfile.write("#mission type\n")
        outputfile.write("#0: MGALTS\n")
        outputfile.write("#1: FBLTS\n")
        outputfile.write("#2: MGALT\n")
        outputfile.write("#3: FBLT\n")
        outputfile.write("#4: PSBI\n")
        outputfile.write("#5: PSFB\n")
        outputfile.write("#6: MGAnDSMs\n")
        outputfile.write("#7: CoastPhase\n")
        outputfile.write("#8: SundmanCoastPhase\n")
        outputfile.write("#9: Variable\n")
        outputfile.write("mission_type " + str(self.mission_type) + "\n")
        outputfile.write("#number of journeys (user-defined endpoints)\n")
        outputfile.write("#Each journey has a central body and two boundary points\n")
        outputfile.write("#Each central body has a menu of destinations which is used to choose the boundary points. Every menu is structured:\n")
        outputfile.write("#-1: Boundary at a point in space, either fixed or free\n")
        outputfile.write("#0: Nothing happens. This code is only used to signify ""no flyby"" and should NEVER be coded as a destination.\n")
        outputfile.write("#1: Body 1 (i.e. Mercury, Io, etc)\n")
        outputfile.write("#2: Body 2 (i.e. Venus, Europa, etc)\n")
        outputfile.write("#...\n")
        outputfile.write("#N: Body N\n")
        outputfile.write("number_of_journeys " + str(self.number_of_journeys) + "\n")
        outputfile.write("#maximum number of phases allowed per journey\n")
        outputfile.write("max_phases_per_journey " + str(self.max_phases_per_journey) + "\n")
        outputfile.write("#destination list \n")
        outputfile.write("destination_list")
        for j in range(0, self.number_of_journeys):
            # if j > 0:
                # if self.Journeys[j].destination_list[0] != self.Journeys[j-1].destination_list[1]:
                    # print("Second entry in destination list for Journey ", j, " does not match first entry in destination list for Journey ", j-1, ". Are you sure that you want to do that?")
            outputfile.write(" " + str(self.Journeys[j].destination_list[0]) + " " + str(self.Journeys[j].destination_list[1]))
        outputfile.write("\n")
        outputfile.write("#the following option is relevant only if optimizing over total deltaV, should the initial impulse be included in the cost?\n")
        outputfile.write("include_initial_impulse_in_cost " + str(self.include_initial_impulse_in_cost) + "\n")
        outputfile.write("#global time bounds\n")
        outputfile.write("#0: unbounded\n")
        outputfile.write("#1: bounded total time (note that the global arrival date bound is by definition the same as the last journey arrival date bound and is not duplicated\n")
        outputfile.write("global_timebounded " + str(self.global_timebounded) + "\n")
        outputfile.write("#MJD of the opening of the launch window\n")
        outputfile.write("launch_window_open_date " + str(self.launch_window_open_date) + "\n")
        outputfile.write("#total flight time bounds, in days\n")
        outputfile.write("total_flight_time_bounds " + str(self.total_flight_time_bounds[0]) + " " + str(self.total_flight_time_bounds[1]) + "\n")
        outputfile.write("#objective function type\n")
        outputfile.write("#0: minimum deltaV\n")
        outputfile.write("#1: minimum time\n")
        outputfile.write("#2: maximum final mass\n")    
        outputfile.write("#3: maximize initial mass\n")
        outputfile.write("#4: depart as late as possible in the window\n")
        outputfile.write("#5: depart as early as possible in the window\n")
        outputfile.write("#6: maximize orbit energy\n")
        outputfile.write("#7: minimize launch mass\n")
        outputfile.write("#8: arrive as early as possible\n")
        outputfile.write("#9: arrive as late as possible\n")
        outputfile.write("#10: minimum propellant (not the same as #2)\n")
        outputfile.write("#11: maximum dry/wet ratio\n")
        outputfile.write("#12: maximum arrival kinetic energy\n")
        outputfile.write("#13: minimum BOL power\n")
        outputfile.write("#14: maximize log_10(final mass)\n")
        outputfile.write("#15: maximum log_e(final mass)\n")
        outputfile.write("#16: maximum dry mass margin\n")
        outputfile.write("#17: maximum dry mass\n")
        outputfile.write("#18: maximum log_10(dry mass)\n")
        outputfile.write("#19: maximum log_e(dry mass)\n")
        outputfile.write("#20: minimize chemical fuel\n")
        outputfile.write("#21: minimize chemical oxidizer\n")
        outputfile.write("#22: minimize electric propellant\n")
        outputfile.write("#23: minimize total propellant\n")
        outputfile.write("#24: minimize waypoint tracking error\n")
        outputfile.write("objective_type " + str(self.objective_type) + "\n")
        outputfile.write("#bounds on the DLA, in degrees (typically set to declination of your launch site)\n") 
        outputfile.write("DLA_bounds " + str(self.DLA_bounds[0]) + " " + str(self.DLA_bounds[1]) + "\n")
        outputfile.write("#bounds on the RLA, in degrees (typically set to a greater than 360 degree range)\n")
        outputfile.write("RLA_bounds " + str(self.RLA_bounds[0]) + " " + str(self.RLA_bounds[1]) + "\n")
        outputfile.write("#Forced post-launch coast (in days, to be enforced after launch)\n")
        outputfile.write("forced_post_launch_coast " + str(self.forced_post_launch_coast) + "\n")
        outputfile.write("#Forced pre flyby coast (in days, to be enforced before each flyby)\n")
        outputfile.write("forced_pre_flyby_coast " + str(self.forced_pre_flyby_coast) + "\n")
        outputfile.write("#Forced post flyby coast (in days, to be enforced after each flyby)\n")
        outputfile.write("forced_post_flyby_coast " + str(self.forced_post_flyby_coast) + "\n")
        outputfile.write("#Waypoint file path\n")
        outputfile.write("waypoint_file_path " + self.waypoint_file_path + "\n")
        outputfile.write("#Covariance file path\n")
        outputfile.write("covariance_file_path " + self.waypoint_file_path + "\n")
        outputfile.write("\n")

        outputfile.write("##Settings for each journey\n")   
        outputfile.write("##dummy values should be used - they should not be necessary but testing was not exhaustive so please use them\n")
        outputfile.write("#journey names\n")
        outputfile.write("journey_names")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_names))
        outputfile.write("\n")
        outputfile.write("#Override this journey's number of steps?\n")
        outputfile.write("journey_override_num_steps")
        for journey in self.Journeys:
            outputfile.write(" " + str(journey.journey_override_num_steps))
        outputfile.write("\n")
        outputfile.write("#Number of time steps for this journey, if overriden\n")
        outputfile.write("journey_number_of_steps")
        for journey in self.Journeys:
            outputfile.write(" " + str(journey.journey_number_of_steps))
        outputfile.write("\n")
        outputfile.write( "#Number of interior control points for parallel shooting phase types\n")
        outputfile.write( "journey_num_interior_control_points")
        for journey in self.Journeys:
            outputfile.write(" " + str(journey.journey_num_interior_control_points))
        outputfile.write("\n")
        outputfile.write("#Freeze this journey's decision variables?\n")
        outputfile.write("journey_freeze_decision_variables")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(int(self.Journeys[j].journey_freeze_decision_variables)))
        outputfile.write("\n")
        outputfile.write("#How much mass to add to the spacecraft at the end of the journey (a negative number indicates a mass drop)\n")
        outputfile.write("journey_fixed_ending_mass_increment")
        for j in range (0,self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_fixed_ending_mass_increment))
        outputfile.write('\n')
        outputfile.write("#How much mass to add to the spacecraft at the beginning of the journey (a negative number indicates a mass drop)\n")
        outputfile.write("journey_fixed_starting_mass_increment")
        for j in range (0,self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_fixed_starting_mass_increment))
        outputfile.write('\n')
        outputfile.write("#Is the mass increment variable (i.e. can the optimizer choose how much mass to add)\n")
        outputfile.write("#This option is ignored for journeys with zero or negative mass increment\n")
        outputfile.write("journey_variable_mass_increment")
        for j in range (0,self.number_of_journeys):
            outputfile.write(" " + str(int(self.Journeys[j].journey_variable_mass_increment)))
        outputfile.write('\n')
        outputfile.write("#What is the minimum amount of mass to drop if you have a variable mass increment.\n")
        outputfile.write("#Put the value here with the smallest absolute value\n")
        outputfile.write("journey_minimum_starting_mass_increment")
        for j in range (0,self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_minimum_starting_mass_increment));
        outputfile.write('\n')
        outputfile.write("#What is the maximum amount of mass to drop if you have a variable mass increment.\n")
        outputfile.write("#Put the value here with the largest absolute value\n")
        outputfile.write("journey_maximum_starting_mass_increment")
        for j in range (0,self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_maximum_starting_mass_increment))
        outputfile.write('\n')
        outputfile.write("#Do you want to constrain the maximum initial mass for this journey. Meant to be used\n")
        outputfile.write("#in concert with journey_variable_mass_increment.\n")
        outputfile.write("journey_constrain_initial_mass")
        for j in range (0,self.number_of_journeys):
            outputfile.write(" " + str(int(self.Journeys[j].journey_constrain_initial_mass)))
        outputfile.write('\n')
        outputfile.write("#If initial mass for this journey is constrained, enter the constraint value here.\n")
        outputfile.write("journey_maximum_initial_mass")
        for j in range (0,self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_maximum_initial_mass))
        outputfile.write('\n')
        outputfile.write("#is each journey time bounded (one value per journey)\n")
        outputfile.write("#0: unbounded\n")
        outputfile.write("#1: bounded flight time\n")
        outputfile.write("#2: bounded arrival date\n")
        outputfile.write("#3: bounded aggregate flight time\n")
        outputfile.write("journey_timebounded")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_timebounded))
        outputfile.write("\n")
        outputfile.write("#Does each journey have a bounded departure date? (note this is redundant for the first journey)\n")
        outputfile.write("journey_bounded_departure_date")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(int(self.Journeys[j].journey_bounded_departure_date)))
        outputfile.write("\n")
        outputfile.write("#Lower and upper bounds on journey departure date, in MJD. Two numbers per journey\n")
        outputfile.write("journey_departure_date_bounds")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_departure_date_bounds[0]) + " " + str(self.Journeys[j].journey_departure_date_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#what are the wait time lower and upper bounds, in days, for each journey (two numbers per journey)\n")   
        outputfile.write("journey_wait_time_bounds")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_wait_time_bounds[0]) + " " + str(self.Journeys[j].journey_wait_time_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#what are the flight time bounds for each journey (two numbers per journey, use dummy values if no flight time bounds)\n")    
        outputfile.write("journey_flight_time_bounds")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_flight_time_bounds[0]) + " " + str(self.Journeys[j].journey_flight_time_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#what are the arrival date bounds for each journey (two numbers per journey, use dummy values if no flight time bounds)\n")   
        outputfile.write("journey_arrival_date_bounds")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_arrival_date_bounds[0]) + " " + str(self.Journeys[j].journey_arrival_date_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#what are the bounds on the initial impulse for each journey in km/s (two numbers per journey)\n")    
        outputfile.write("#you can set a very high upper bound if you are using a launchy vehicle model - the optimizer will find the correct value\n") 
        outputfile.write("journey_initial_impulse_bounds")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_initial_impulse_bounds[0]) + " " + str(self.Journeys[j].journey_initial_impulse_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#journey departure type (one value per journey)\n")   
        outputfile.write("#0: launch or direct insertion\n")
        outputfile.write("#1: depart from parking orbit (you can use this one in place of a launch vehicle model, and the departure burn will be done with the EDS motor)\n")
        outputfile.write("#2: free direct departure, i.e. do not burn to get the departure v_infinity (used for when operations about a small body are not modeled but the departure velocity is known)\n")
        outputfile.write("#3: flyby (only valid for successive journeys)\n")
        outputfile.write("#4: flyby with fixed v-infinity-out (only valid for successive journeys)\n")
        outputfile.write("#5: spiral-out from circular orbit (low-thrust missions only)\n")
        outputfile.write("#6: zero-turn flyby (for small bodies)\n")
        outputfile.write("journey_departure_type")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_departure_type))
        outputfile.write("\n")
        outputfile.write("#journey departure boundary class (one value per journey)\n")
        outputfile.write("#0: Ephemeris-pegged (default EMTG)\n")
        outputfile.write("#1: Free point\n")
        outputfile.write("#2: Ephemeris-referenced\n")
        outputfile.write("#3: Periapse\n")
        outputfile.write("journey_departure_class")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_departure_class))
        outputfile.write("\n")
        outputfile.write("#journey departure ellipsoid axes (3 per journey)\n")
        outputfile.write("journey_departure_ellipsoid_axes")
        for thisJourney in self.Journeys:
            for entry in thisJourney.journey_departure_ellipsoid_axes:
                outputfile.write(' ' + str(entry))
        outputfile.write("\n")
        outputfile.write("#PeriapseDeparture altitude bounds (two per journey, in km)\n")
        outputfile.write("journey_PeriapseDeparture_altitude_bounds")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_PeriapseDeparture_altitude_bounds[0]) + " " + str(self.Journeys[j].journey_PeriapseDeparture_altitude_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#journey arrival type (one value per journey)\n")
        outputfile.write("#0: insertion into parking orbit (use chemical Isp)\n")
        outputfile.write("#1: rendezvous (use chemical Isp)\n")
        outputfile.write("#2: intercept with bounded V_infinity\n")
        outputfile.write("#3: low-thrust rendezvous (does not work if terminal phase is not low-thrust)\n")
        outputfile.write("#4: match final v-infinity vector\n")
        outputfile.write("#5: match final v-infinity vector (low-thrust)\n")
        outputfile.write("#6: capture spiral\n")
        outputfile.write("journey_arrival_type")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_arrival_type))
        outputfile.write("\n")
        outputfile.write("#journey arrival boundary class (one value per journey)\n")
        outputfile.write("#0: Ephemeris-pegged (default EMTG)\n")
        outputfile.write("#1: Free point\n")
        outputfile.write("#2: Ephemeris-referenced\n")
        outputfile.write("#3: Periapse\n")
        outputfile.write("journey_arrival_class")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_arrival_class))
        outputfile.write("\n")
        outputfile.write("#journey arrival ellipsoid axes (3 per journey)\n")
        outputfile.write("journey_arrival_ellipsoid_axes")
        for thisJourney in self.Journeys:
            for entry in thisJourney.journey_arrival_ellipsoid_axes:
                outputfile.write(' ' + str(entry))
        outputfile.write("\n")
        outputfile.write("#Journey forced initial coast (in days)\n")
        outputfile.write("journey_forced_initial_coast")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_forced_initial_coast))
        outputfile.write("\n")
        outputfile.write("#override PeriapseArrival altitude bounds?\n")
        outputfile.write("journey_PeriapseArrival_override_altitude")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(int(self.Journeys[j].journey_PeriapseArrival_override_altitude)))
        outputfile.write("\n")
        outputfile.write("#PeriapseArrival altitude bounds (two per journey, in km)\n")
        outputfile.write("journey_PeriapseArrival_altitude_bounds")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_PeriapseArrival_altitude_bounds[0]) + " " + str(self.Journeys[j].journey_PeriapseArrival_altitude_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#Override journey flyby altitude? (one per journey, affects only first phase)\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("journey_override_flyby_altitude_bounds")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(int(self.Journeys[j].journey_override_flyby_altitude_bounds)))
        outputfile.write("\n")
        outputfile.write("#Lower and upper bound on journey flyby altitude (two per journey, in km, only affects first phase)\n")
        outputfile.write("journey_flyby_altitude_bounds")
        for journey in self.Journeys:
            outputfile.write(" " + str(journey.journey_flyby_altitude_bounds[0]) + " " + str(journey.journey_flyby_altitude_bounds[1]))
        outputfile.write("\n")
        outputfile.write("#type of orbit elements specified at beginning of journey(0: inertial, 1: COE)\n")
        outputfile.write("journey_departure_elements_type")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_departure_elements_type))
        outputfile.write("\n")
        outputfile.write("#reference frame for journey departure elements (0: J2000_ICRF, 1: J2000_BCI, 2: J2000_BCF, 3: TrueOfDate_BCI, 4: TrueOfDate_BCF, 5: Principle Axes, 6: Topocentric, 7: Polar)\n")
        outputfile.write("journey_departure_elements_frame")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_departure_elements_frame))
        outputfile.write("\n")
        outputfile.write("#orbit elements at beginning of journey (a(km), e, i, RAAN, AOP, MA) supply angles in degrees\n")
        outputfile.write("journey_departure_elements")
        for j in range (0, self.number_of_journeys):
            for k in range(0, 6):
                outputfile.write(" " + str(self.Journeys[j].journey_departure_elements[k]))
        outputfile.write("\n")
        outputfile.write("#Vary journey departure elements? (one entry per element per journey: 0 means no, 1 means yes)\n")
        outputfile.write("journey_departure_elements_vary_flag")
        for j in range (0, self.number_of_journeys):
            for k in range(0, 6):
                outputfile.write(" " + str(self.Journeys[j].journey_departure_elements_vary_flag[k]))
        outputfile.write("\n")
        outputfile.write("#Lower and upper bounds on journey departure elements (two per element per journey, ignored if vary flag is off for that element)\n")
        outputfile.write("journey_departure_elements_bounds")
        for j in range (0, self.number_of_journeys):
            for k in range(0, 12):
                outputfile.write(" " + str(self.Journeys[j].journey_departure_elements_bounds[k]))
        outputfile.write("\n")    
        outputfile.write("#Reference epoch (MJD) for journey departure elements\n")
        outputfile.write("journey_departure_elements_reference_epoch")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_departure_elements_reference_epoch))
        outputfile.write("\n")
        outputfile.write("#Allow journey departure free point boundary to propagate (otherwise it is a fixed waypoint)\n")
        outputfile.write("AllowJourneyFreePointDepartureToPropagate")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].AllowJourneyFreePointDepartureToPropagate))
        outputfile.write("\n")

        outputfile.write("#type of orbit elements specified at end of journey(0: inertial, 1: COE)\n")
        outputfile.write("journey_arrival_elements_type")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_arrival_elements_type))
        outputfile.write("\n")
        outputfile.write("#reference frame for journey arrival elements (0: J2000_ICRF, 1: J2000_BCI, 2: J2000_BCF, 3: TrueOfDate_BCI, 4: TrueOfDate_BCF, 5: Principle Axes, 6: Topocentric, 7: Polar)\n")
        outputfile.write("journey_arrival_elements_frame")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_arrival_elements_frame))
        outputfile.write("\n")
        outputfile.write("#orbit elements at end of journey (a(km), e, i, RAAN, AOP, MA) supply angles in degrees\n")
        outputfile.write("journey_arrival_elements")
        for j in range (0, self.number_of_journeys):
            for k in range(0, 6):
                outputfile.write(" " + str(self.Journeys[j].journey_arrival_elements[k]))
        outputfile.write("\n")
        outputfile.write("#Vary journey arrival elements? (one entry per element per journey: 0 means no, 1 means yes)\n")
        outputfile.write("journey_arrival_elements_vary_flag")
        for j in range (0, self.number_of_journeys):
            for k in range(0, 6):
                outputfile.write(" " + str(self.Journeys[j].journey_arrival_elements_vary_flag[k]))
        outputfile.write("\n")
        outputfile.write("#Lower and upper bounds on journey arrival elements (two per element per journey, ignored if vary flag is off for that element)\n")
        outputfile.write("journey_arrival_elements_bounds")
        for j in range (0, self.number_of_journeys):
            for k in range(0, 12):
                outputfile.write(" " + str(self.Journeys[j].journey_arrival_elements_bounds[k]))
        outputfile.write("\n")
        outputfile.write("#Reference epoch (MJD) for journey arrival elements\n")
        outputfile.write("journey_arrival_elements_reference_epoch")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_arrival_elements_reference_epoch))
        outputfile.write("\n")
        outputfile.write("#Allow journey arrival free point boundary to propagate (otherwise it is a fixed waypoint)\n")
        outputfile.write("AllowJourneyFreePointArrivalToPropagate")
        for j in range (0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].AllowJourneyFreePointArrivalToPropagate))
        outputfile.write("\n")

        outputfile.write("#journey central body\n")
        outputfile.write("#Must match the name of a Universe file in the Universe folder\n")
        outputfile.write("journey_central_body")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_central_body))
        outputfile.write("\n")
        outputfile.write("#final VHP for journeys that end in intercepts, in km/s (three numbers per journey)\n")
        outputfile.write("journey_final_velocity")
        for j in range(0, self.number_of_journeys):
            for k in range(0,3):
                outputfile.write(" " + str(self.Journeys[j].journey_final_velocity[k]))
        outputfile.write("\n")
        outputfile.write("#Impose arrival declination constraint on each journey?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("journey_arrival_declination_constraint_flag")
        for journey in self.Journeys:
            outputfile.write(" "  + str(journey.journey_arrival_declination_constraint_flag))
        outputfile.write("\n")
        outputfile.write("#Arrival declination bounds for each journey\n")
        outputfile.write("#Two numbers per journey, in degrees\n")
        outputfile.write("journey_arrival_declination_bounds")
        for journey in self.Journeys:
            outputfile.write(" " + str(journey.journey_arrival_declination_bounds[0]) + " " + str(journey.journey_arrival_declination_bounds[1]))
        outputfile.write("\n")

        outputfile.write("#Starting orbital radius for an escape spiral at the beginning of the journey\n")
        outputfile.write("journey_escape_spiral_starting_radius")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_escape_spiral_starting_radius))
        outputfile.write("\n")
        outputfile.write("#Final orbital radius for a capture spiral at the end of the journey\n")
        outputfile.write("journey_capture_spiral_final_radius")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_capture_spiral_final_radius))
        outputfile.write("\n")
        outputfile.write("#Journey forced terminal coast (in days)\n")
        outputfile.write("journey_forced_terminal_coast")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_forced_terminal_coast))
        outputfile.write("\n")

        outputfile.write("#Journey-end delta-v (km/s)\n")
        outputfile.write("journey_end_deltav")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_end_deltav))
        outputfile.write("\n")
        outputfile.write("#Propulsion system for journey-end maneuver\n")
        outputfile.write("#0: Monoprop chemical\n")
        outputfile.write("#1: Biprop chemical\n")
        outputfile.write("#2: Electric\n")
        outputfile.write("journey_end_propulsion_system")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_end_propulsion_system))
        outputfile.write("\n")
        outputfile.write("#Journey override global duty cycle\n")
        outputfile.write("journey_override_duty_cycle")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(int(self.Journeys[j].journey_override_duty_cycle)))
        outputfile.write("\n")
        outputfile.write("#Journey duty cycle\n")
        outputfile.write("journey_duty_cycle")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_duty_cycle))
        outputfile.write("\n")
        outputfile.write("#Journey-end TCM magnitude (km/s)\n");
        outputfile.write("journey_end_TCM");
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_end_TCM))
        outputfile.write("\n")
        outputfile.write("#Journey override propagator type?\n")
        outputfile.write("journey_override_PropagatorType")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_override_PropagatorType))
        outputfile.write("\n")
        outputfile.write("#Journey propagator type?\n")
        outputfile.write("journey_PropagatorType")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_PropagatorType))
        outputfile.write("\n")
        outputfile.write("#Journey override integration step size\n")
        outputfile.write("journey_override_integration_step_size")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(int(self.Journeys[j].journey_override_integration_step_size)))
        outputfile.write("\n")
        outputfile.write("#Journey integration step size (days)\n")
        outputfile.write("journey_integration_step_size")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_integration_step_size))
        outputfile.write("\n")
        outputfile.write("#Journey coast phase match point fraction\n")
        outputfile.write("journey_CoastPhaseMatchPointFraction")        
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_CoastPhaseMatchPointFraction))
        outputfile.write("\n")
        outputfile.write("#Journey coast phase forward integration step length (seconds)\n")
        outputfile.write("journey_CoastPhaseForwardIntegrationStepLength")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_CoastPhaseForwardIntegrationStepLength))
        outputfile.write("\n")
        outputfile.write("#Journey coast phase backward integration step length (seconds)\n")
        outputfile.write("journey_CoastPhaseBackwardIntegrationStepLength")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_CoastPhaseBackwardIntegrationStepLength))
        outputfile.write("\n")

        outputfile.write("\n")
        outputfile.write("##Staging options\n")
        outputfile.write("#Stage after departure? (one entry per journey)\n")
        outputfile.write("journey_stage_after_departure")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(int(self.Journeys[j].journey_stage_after_departure)))
        outputfile.write("\n")
        outputfile.write("#Stage before arrival? (one entry per journey)\n")
        outputfile.write("journey_stage_before_arrival")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(int(self.Journeys[j].journey_stage_before_arrival)))
        outputfile.write("\n")
        outputfile.write("#Stage after arrival? (one entry per journey)\n")
        outputfile.write("journey_stage_after_arrival")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(int(self.Journeys[j].journey_stage_after_arrival)))
        outputfile.write("\n")
            
        outputfile.write("##Perturbation settings\n")
        outputfile.write("#Enable solar radiation pressure?\n")
        outputfile.write("perturb_SRP " + str(self.perturb_SRP) + "\n")
        outputfile.write("#Enable third-body perturbations?\n")
        outputfile.write("perturb_thirdbody " + str(self.perturb_thirdbody) + "\n")
        outputfile.write("#Enable J2 perturbations?\n")
        outputfile.write("perturb_J2 " + str(self.perturb_J2) + "\n")
        outputfile.write("#Journey perturbation bodies. One line per journey. The numbers in the line correspond to\n")
        outputfile.write("#bodies in the journey""s Universe file. If perturbations are off, each line should just have a zero\n")
        outputfile.write("#the numbers in the first line are the number of perturbation bodies for each journey\n")
        outputfile.write("journey_perturbation_bodies")
        for j in range(0, self.number_of_journeys):
            outputfile.write(" " + str(self.Journeys[j].journey_number_of_perturbation_bodies))
        outputfile.write("\n")
        for j in range(0, self.number_of_journeys):
            if self.Journeys[j].journey_number_of_perturbation_bodies > 0:
                outputfile.write(str(self.Journeys[j].journey_perturbation_bodies[0]))
                for b in range(1, self.Journeys[j].journey_number_of_perturbation_bodies):
                    outputfile.write(" " + str(self.Journeys[j].journey_perturbation_bodies[b]))
            else:
                outputfile.write("0")
            outputfile.write("\n")
        outputfile.write("#end_journey_perturbation_bodies\n")
        outputfile.write("#Spacecraft area (in m^2)\n")
        outputfile.write("spacecraft_area " + str(self.spacecraft_area) + "\n")
        outputfile.write("#Coefficient of reflectivity\n")
        outputfile.write("#0.0: perfectly translucent\n")
        outputfile.write("#1.0: perfectly absorbing\n")
        outputfile.write("#2.0: perfectly reflecting\n")
        outputfile.write("coefficient_of_reflectivity " + str(self.coefficient_of_reflectivity) + "\n")
        outputfile.write("\n")

        outputfile.write("##Outer-loop selectable options settings\n")
        outputfile.write("#Allow outer-loop to vary spacecraft file (only relevant if SpacecraftModelInput is set to ReadSpacecraftFile)\n")
        outputfile.write("outerloop_vary_spacecraft " + str(self.outerloop_vary_spacecraft) + "\n")
        outputfile.write("#Allow outer-loop to vary power system?\n")
        outputfile.write("outerloop_vary_power_system " + str(self.outerloop_vary_power_system) + "\n")
        outputfile.write("#Allow outer-loop to vary chemical fuel tank capacity?\n")
        outputfile.write("outerloop_vary_chemical_fuel_tank_capacity " + str(self.outerloop_vary_chemical_fuel_tank_capacity) + "\n")
        outputfile.write("#Allow outer-loop to vary chemical oxidizer tank capacity?\n")
        outputfile.write("outerloop_vary_chemical_oxidizer_tank_capacity " + str(self.outerloop_vary_chemical_oxidizer_tank_capacity) + "\n")
        outputfile.write("#Allow outer-loop to vary electric propellant tank capacity?\n")
        outputfile.write("outerloop_vary_electric_propellant_tank_capacity " + str(self.outerloop_vary_electric_propellant_tank_capacity) + "\n")
        outputfile.write("#Allow outer-loop to vary electric propulsion system?\n")
        outputfile.write("outerloop_vary_electric_propulsion_system " + str(self.outerloop_vary_electric_propulsion_system) + "\n")
        outputfile.write("#Allow outer-loop to vary number of electric propulsion systems?\n")
        outputfile.write("outerloop_vary_number_of_electric_propulsion_systems " + str(self.outerloop_vary_number_of_electric_propulsion_systems) + "\n")
        outputfile.write("#Allow outer-loop to vary thruster duty cycle?\n")
        outputfile.write("outerloop_vary_duty_cycle " + str(self.outerloop_vary_duty_cycle) + "\n")
        outputfile.write("#Allow outer-loop to vary launch vehicle?\n")
        outputfile.write("outerloop_vary_launch_vehicle " + str(self.outerloop_vary_launch_vehicle) + "\n")
        outputfile.write("#Allow outer-loop to vary launch epoch?\n")
        outputfile.write("outerloop_vary_launch_epoch " + str(self.outerloop_vary_launch_epoch) + "\n")
        outputfile.write("#Allow outer-loop to vary flight time upper bound?\n")
        outputfile.write("outerloop_vary_flight_time_upper_bound " + str(self.outerloop_vary_flight_time_upper_bound) + "\n")
        outputfile.write("#Restrict flight-time lower bound when running outer-loop?\n")
        outputfile.write("outerloop_restrict_flight_time_lower_bound " + str(self.outerloop_restrict_flight_time_lower_bound) + "\n")
        outputfile.write("#Allow outer-loop to vary first journey departure C3?\n")
        outputfile.write("outerloop_vary_departure_C3 " + str(self.outerloop_vary_departure_C3) + "\n")
        outputfile.write("#Allow outer-loop to vary last journey arrival C3?\n")
        outputfile.write("outerloop_vary_arrival_C3 " + str(self.outerloop_vary_arrival_C3) + "\n")
        outputfile.write("#Allow outer-loop to vary last journey arrival declination?\n")
        outputfile.write("outerloop_vary_arrival_declination " + str(self.outerloop_vary_arrival_declination) + "\n")
        outputfile.write("#Allow outer-loop to vary last journey entry interface velocity?\n")
        outputfile.write("outerloop_vary_final_journey_entry_interface_velocity " + str(self.outerloop_vary_final_journey_entry_interface_velocity) + "\n")
        outputfile.write("#Allow outer-loop to vary number of DSMs in each phase?\n")
        outputfile.write("outerloop_vary_number_of_DSMs " + str(self.outerloop_vary_number_of_DSMs) + "\n")
        outputfile.write("#Allow outer-loop to vary phase type\n")
        outputfile.write("outerloop_vary_phase_type " + str(self.outerloop_vary_phase_type) + "\n")
        outputfile.write("#Allow outer-loop to vary journey destination? (one value per journey)\n")
        outputfile.write("outerloop_vary_journey_destination")
        for journey in self.Journeys:
            outputfile.write(" " + str(journey.outerloop_vary_journey_destination))
        outputfile.write("\n")
        outputfile.write("#Allow outer-loop to vary journey flyby sequence? (one value per journey)\n")
        outputfile.write("outerloop_vary_journey_flyby_sequence")
        for journey in self.Journeys:
            outputfile.write(" " + str(journey.outerloop_vary_journey_flyby_sequence));
        outputfile.write("\n")
        outputfile.write("#Allow outer-loop to vary journey arrival type? (one value per journey)\n")
        outputfile.write("outerloop_vary_journey_arrival_type")
        for journey in self.Journeys:
            outputfile.write(" " + str(journey.outerloop_vary_journey_arrival_type));
        outputfile.write("\n")
        outputfile.write("#Allow outer-loop to choose which journey to stop after?\n")
        outputfile.write("outerloop_vary_stop_after_journey " + str(self.outerloop_vary_stop_after_journey) + "\n")
        outputfile.write("\n")
        outputfile.write("#Outer-loop spacecraft file choices\n")
        outputfile.write("outerloop_spacecraft_choices")
        for entry in self.outerloop_spacecraft_choices:
            outputfile.write(" " + entry)
        outputfile.write("\n")
        outputfile.write("#Outer-loop power system choices\n")
        outputfile.write("outerloop_power_system_choices")
        for entry in self.outerloop_power_system_choices:
            outputfile.write(" " + entry)
        outputfile.write("\n")
        outputfile.write("#Outer-loop chemical fuel tank capacity choices (in kg)\n")
        outputfile.write("outerloop_chemical_fuel_tank_capacity_choices")
        for entry in self.outerloop_chemical_fuel_tank_capacity_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop chemical oxidizer tank capacity choices (in kg)\n")
        outputfile.write("outerloop_chemical_oxidizer_tank_capacity_choices")
        for entry in self.outerloop_chemical_oxidizer_tank_capacity_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop electric propulsion system choices (in order of most to least preferable)\n")
        outputfile.write("outerloop_electric_propulsion_system_choices")
        for entry in self.outerloop_electric_propulsion_system_choices:
            outputfile.write(" " + entry)
        outputfile.write("\n")
        outputfile.write("#Outer-loop number of electric propulsion systems choices\n")
        outputfile.write("outerloop_number_of_electric_propulsion_systems_choices")
        for entry in self.outerloop_number_of_electric_propulsion_systems_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop duty cycle choices\n")
        outputfile.write("outerloop_duty_cycle_choices")
        for entry in self.outerloop_duty_cycle_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop launch vehicle choices (in order of most to least preferable)\n")
        outputfile.write("outerloop_launch_vehicle_choices")
        for entry in self.outerloop_launch_vehicle_choices:
            outputfile.write(" " + entry)
        outputfile.write("\n")
        outputfile.write("#Outer-loop launch window open epoch choices (in MJD)\n")
        outputfile.write("outerloop_launch_epoch_choices")
        for entry in self.outerloop_launch_epoch_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop flight time upper bound choices (in days)\n")
        outputfile.write("outerloop_flight_time_upper_bound_choices")
        for entry in self.outerloop_flight_time_upper_bound_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop first journey departure C3 choices\n")
        outputfile.write("outerloop_departure_C3_choices")
        for entry in self.outerloop_departure_C3_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop last journey arrival C3 choices\n")
        outputfile.write("outerloop_arrival_C3_choices")
        for entry in self.outerloop_arrival_C3_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop last journey arrival declination choices\n")
        outputfile.write("outerloop_arrival_declination_choices")
        for entry in self.outerloop_arrival_declination_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("##Outer-loop last journey entry interface velocity choices (km/s)\n")
        outputfile.write("outerloop_final_journey_entry_interface_velocity_choices")
        for entry in self.outerloop_final_journey_entry_interface_velocity_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop phase type choices\n")
        outputfile.write("outerloop_phase_type_choices")
        for entry in self.outerloop_phase_type_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("#Outer-loop maximum number of flybys (one value for each journey)\n")
        outputfile.write("outerloop_journey_maximum_number_of_flybys")
        for journey in self.Journeys:
            outputfile.write(" " + str(journey.outerloop_journey_maximum_number_of_flybys))
        outputfile.write("\n")
        outputfile.write("#Outer-loop journey destination choices (one line for each journey)\n")
        outputfile.write("outerloop_journey_destination_choices\n")
        for journey in self.Journeys:
            for entry in journey.outerloop_journey_destination_choices:
                outputfile.write(" " + str(entry))
            outputfile.write("\n")
        outputfile.write("#Outer-loop flyby sequence choices (one line for each journey)\n")
        outputfile.write("outerloop_journey_flyby_sequence_choices\n")
        for journey in self.Journeys:
            for entry in journey.outerloop_journey_flyby_sequence_choices:
                outputfile.write(" " + str(entry))
            outputfile.write("\n")
        outputfile.write("#Outer-loop journey arrival type choices (one line for each journey)\n")
        outputfile.write("outerloop_journey_arrival_type_choices\n")
        for journey in self.Journeys:
            for entry in journey.outerloop_journey_arrival_type_choices:
                outputfile.write(" " + str(entry))
            outputfile.write("\n")
        outputfile.write("\n")

        outputfile.write("##Settings for outer-loop variable number of journeys\n")
        outputfile.write("#Allow the outer-loop to vary the number of journeys?\n")
        outputfile.write("outerloop_vary_number_of_journeys " + str(self.outerloop_vary_number_of_journeys) + "\n")
        outputfile.write("#Maximum number of variable journeys that the outer-loop may add\n")
        outputfile.write("outerloop_maximum_number_of_variable_journeys " + str(self.outerloop_maximum_number_of_variable_journeys) + "\n")
        outputfile.write("#After which fixed journey should the variable journeys be inserted (numbers greater than the number of fixed journeys will place all variable journeys at the end\n")
        outputfile.write("outerloop_insert_variable_journeys_after_fixed_journey_index " + str(self.outerloop_insert_variable_journeys_after_fixed_journey_index) + "\n")
        outputfile.write("#Which fixed journey should be the template for all journey-specific options in the variable journeys?\n")
        outputfile.write("outerloop_inherit_variable_journey_settings_from_fixed_journey_index " + str(self.outerloop_inherit_variable_journey_settings_from_fixed_journey_index) + "\n")
        outputfile.write("\n")

        outputfile.write("##Outer-loop objective function settings\n")
        outputfile.write("#Pick as many as you want. The Pareto surface will be generated in these dimensions\n")
        outputfile.write("#0: BOL power at 1 AU (kW)\n")
        outputfile.write("#1: Launch epoch (MJD)\n")
        outputfile.write("#2: Flight time (days)\n")
        outputfile.write("#3: Thruster preference\n")
        outputfile.write("#4: Number of thrusters\n")
        outputfile.write("#5: Launch vehicle preference\n")
        outputfile.write("#6: Delivered mass to final target (kg)\n")
        outputfile.write("#7: Final journey mass increment (for maximizing sample return)\n")
        outputfile.write("#8: First journey departure C3 (km^2/s^2)\n")
        outputfile.write("#9: Final journey arrival C3 (km^2/s^2)\n")
        outputfile.write("#10: Total deterministic delta-v (km/s)\n")
        outputfile.write("#11: Inner-loop objective (whatever it was)\n")
        outputfile.write("#12: Point-group value\n")
        outputfile.write("#13: Total propellant mass including margin (kg)\n")
        outputfile.write("#14: Number of journeys\n")
        outputfile.write("#15: Dry mass margin\n")
        outputfile.write("#16: bus power (kW)\n")
        outputfile.write("#17: Final journey interface velocity (km/s)\n")
        outputfile.write("#18: Thruster duty cycle\n")
        outputfile.write("#19: Normalized aggregate control\n")
        outputfile.write("#20: final journey arrival declination\n")
        outputfile.write("outerloop_objective_function_choices")
        for entry in self.outerloop_objective_function_choices:
            outputfile.write(" " + str(entry))
        outputfile.write("\n")
        outputfile.write("\n")
        
        outputfile.write("##Outer-loop point group settings\n")
        outputfile.write("#Point group values and members\n")
        outputfile.write("outerloop_point_groups_values")
        for g in range(0, len(self.outerloop_point_groups_values)):
            outputfile.write(" " + str(self.outerloop_point_groups_values[g]))
        outputfile.write("\n")
        for g in range(0, len(self.outerloop_point_groups_values)):
            for m in range(0, len(self.outerloop_point_groups_members[g])):
                outputfile.write(" " + str(self.outerloop_point_groups_members[g][m]))
            outputfile.write("\n")
        outputfile.write("#How many members to score from each point group (additional members add no more points)\n")
        outputfile.write("outerloop_point_groups_number_to_score")
        for g in range(0, len(self.outerloop_point_groups_values)):
            outputfile.write(" " + str(self.outerloop_point_groups_number_to_score[g]))
        outputfile.write("\n")
        outputfile.write("\n")

        outputfile.write("##Outer-loop filter settings\n")
        outputfile.write("#Relative inclination filter\n")
        outputfile.write("outerloop_filter_successive_destinations_on_inclination " + str(self.outerloop_filter_successive_destinations_on_inclination) + "\n")
        outputfile.write("#Bandpass for outer-loop destination inclination filter (degrees)\n")
        outputfile.write("outerloop_destination_inclination_filter_bandpass " + str(self.outerloop_destination_inclination_filter_bandpass) + "\n")
        outputfile.write("#Group filter\n")
        outputfile.write("outerloop_filter_by_groups " + str(self.outerloop_filter_by_groups) + "\n")
        outputfile.write("#Group filter lists (first number is how many lists, successive lines are what is in the list)" + "\n")
        outputfile.write("outerloop_group_filter_lists "  + str(self.outerloop_group_filter_number_of_groups) + "\n")
        for g in range(0, self.outerloop_group_filter_number_of_groups):
            for m in range(0, len(self.outerloop_group_filter_lists[g])):
                outputfile.write(" " + str(self.outerloop_group_filter_lists[g][m]))
            outputfile.write("\n")
        outputfile.write("\n")
        outputfile.write("\n")
            
        outputfile.write("##output format settings\n")
        outputfile.write("#Output file frame\n")
        outputfile.write("#0: ICRF, 1: J2000_BCI, 2: J2000_BCF, 3: TrueOfDate_BCI, 4: TrueOfDate_BCF, 5: Principle Axes, 6: Topocentric, 7: Polar\n")
        outputfile.write("output_file_frame " + str(self.output_file_frame) + "\n")
        outputfile.write("#Output journey entries for wait times at intermediate and final target?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("output_dormant_journeys " + str(self.output_dormant_journeys) + "\n")
        outputfile.write("#Post-mission wait time at the final target (if zero, no post-mission ephemeris will be printed)\n")
        outputfile.write("post_mission_wait_time " + str(self.post_mission_wait_time) + "\n")
        outputfile.write("#Override working directory?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("override_working_directory " + str(self.override_working_directory) + "\n")
        outputfile.write("#Custom working directory\n")
        outputfile.write("forced_working_directory " + self.forced_working_directory + "\n")
        outputfile.write("#Override working directory?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("override_mission_subfolder " + str(self.override_mission_subfolder) + "\n")
        outputfile.write("#Custom working directory\n")
        outputfile.write("forced_mission_subfolder " + self.forced_mission_subfolder + "\n")
        outputfile.write("#Short output file names?\n")
        outputfile.write("short_output_file_names " + str(self.short_output_file_names) + "\n")
        outputfile.write("#Generate forward integrated ephemeris (STK compatible)?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("generate_forward_integrated_ephemeris " + str(self.generate_forward_integrated_ephemeris) + "\n")
        outputfile.write("#Add an extra line to the forward integrated ephemeris output to facilitate control switches?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("add_control_switch_line_to_ephemeris " + str(self.add_control_switch_line_to_ephemeris) + "\n")
        outputfile.write("#Append mass to forward integrated ephemeris output?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("append_mass_to_ephemeris_output " + str(self.append_mass_to_ephemeris_output) + "\n")
        outputfile.write("#Append control to forward integrated ephemeris output?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("append_control_to_ephemeris_output " + str(self.append_control_to_ephemeris_output) + "\n")
        outputfile.write("#Append thrust to forward integrated ephemeris output?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("append_thrust_to_ephemeris_output " + str(self.append_thrust_to_ephemeris_output) + "\n")
        outputfile.write("#Append mdot to forward integrated ephemeris output?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("append_mdot_to_ephemeris_output " + str(self.append_mdot_to_ephemeris_output) + "\n")
        outputfile.write("#Append Isp to forward integrated ephemeris output?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("append_Isp_to_ephemeris_output " + str(self.append_Isp_to_ephemeris_output) + "\n")
        outputfile.write("#Append number of active engines to forward integrated ephemeris output?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("append_number_of_active_engines_to_ephemeris_output " + str(self.append_number_of_active_engines_to_ephemeris_output) + "\n")
        outputfile.write("#Append active power to forward integrated ephemeris output?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("append_active_power_to_ephemeris_output " + str(self.append_active_power_to_ephemeris_output) + "\n")
        outputfile.write("#Append throttle level to forward integrated ephemeris output?\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("append_throttle_level_to_ephemeris_output " + str(self.append_throttle_level_to_ephemeris_output) + "\n")

        
        outputfile.write("#Enable background mode (do not ask for key press on exit)\n")
        outputfile.write("#0: no\n")
        outputfile.write("#1: yes\n")
        outputfile.write("background_mode " + str(self.background_mode) + "\n")
        outputfile.write("#Output STMs?\n")
        outputfile.write("output_STMs " + str(self.output_STMs) + "\n")
        outputfile.write("#Output maneuver and target spec files?\n")
        outputfile.write("output_maneuver_and_target_spec_files " + str(self.output_maneuver_and_target_spec_files) + "\n")
        outputfile.write("\n")

        outputfile.write("##manual outer-loop control code\n")
        outputfile.write("#stop after i-th journey? (indexed from 0)\n")
        outputfile.write("stop_after_journey " + str(self.stop_after_journey) + "\n")
        outputfile.write("#sequence, must have (max_phases_per_journey) entries for each journey. Use 0 to encode no flyby\n")
        outputfile.write("#integer codes represent planets\n")
        outputfile.write("#this option is NOT used if the outer-loop is turned on\n")
        outputfile.write("sequence ")
        for j in range(0, self.number_of_journeys):
            #first check to make sure that this journey has a sequence
            #line to print
            if j == 0:
                outputfile.write(str(self.Journeys[j].sequence[0]))
            else:
                outputfile.write(" " + str(self.Journeys[j].sequence[0]))

            for p in range(1, self.max_phases_per_journey):
                if p < len(self.Journeys[j].sequence):
                    outputfile.write(" " + str(self.Journeys[j].sequence[p]))
                else:
                    outputfile.write(" 0")
        outputfile.write("\n")
        outputfile.write("#phase type, must have one entry for each phase in the mission\n")
        outputfile.write("#this option allows you to have different phases use different propulsion systems\n")
        outputfile.write("#0: MGALTS\n")
        outputfile.write("#1: FBLTS\n")
        outputfile.write("#2: MGALT\n")
        outputfile.write("#3: FBLT\n")
        outputfile.write("#4: PSBI\n")
        outputfile.write("#5: PSFB\n")
        outputfile.write("#6: MGAnDSMs\n")
        outputfile.write("#7: CoastPhase\n")
        outputfile.write("#8: SundmanCoastPhase\n")
        outputfile.write("#this option is only read if mission_type is set to Variable\n")
        outputfile.write("#if mission_type == Variable and the outer-loop is ON, then the following option is ignored\n")
        outputfile.write("phase_type")
        for j in range(0, self.number_of_journeys):
            for p in range(0, self.max_phases_per_journey + 1):
                outputfile.write(" " + str(self.Journeys[j].phase_type[p]))
        outputfile.write("\n")
        outputfile.write("#Number of impulses per phase\n")
        outputfile.write("impulses_per_phase ")
        for j in range(0, self.number_of_journeys):
            #first check to make sure that this journey has a sequence
            #line to print
            if j == 0:
                outputfile.write(str(self.Journeys[j].impulses_per_phase[0]))
            else:
                outputfile.write(" " + str(self.Journeys[j].impulses_per_phase[0]))

            for p in range(1, self.max_phases_per_journey + 1):
                if p < len(self.Journeys[j].impulses_per_phase):
                    outputfile.write(" " + str(self.Journeys[j].impulses_per_phase[p]))
                else:
                    outputfile.write(" 0")
        outputfile.write("\n")
        outputfile.write("#powered flyby flag for each phase. The first value for each journey is ignored. The first flyby corresponds to the second flag\n")
        outputfile.write("#with one entry for every phase. As of now, this is not part of the outerloop, so only one line is allowed\n")
        outputfile.write("journey_enable_periapse_burns ")
        for j in range(0, self.number_of_journeys):
            #first check to make sure that this journey has a sequence
            #line to print
            if j == 0:
                outputfile.write(str(self.Journeys[j].journey_enable_periapse_burns[0]))
            else:
                outputfile.write(" " + str(self.Journeys[j].journey_enable_periapse_burns[0]))

            for p in range(1, self.max_phases_per_journey):
                if p < len(self.Journeys[j].journey_enable_periapse_burns):
                    outputfile.write(" " + str(self.Journeys[j].journey_enable_periapse_burns[p]))
                else:
                    outputfile.write(" 0")
        outputfile.write("\n")
        outputfile.write("\n")        
                
        outputfile.write("\n")
        outputfile.write("#Maneuver constraint code\n")
        outputfile.write("#Works for absolute and relative epochs and also magnitudes\n")
        outputfile.write("BEGIN_MANEUVER_CONSTRAINT_BLOCK\n")
        for ManeuverConstraintDefinition in self.ManeuverConstraintDefinitions:
            outputfile.write(ManeuverConstraintDefinition + "\n")
        outputfile.write("END_MANEUVER_CONSTRAINT_BLOCK\n")
        outputfile.write("\n")     
                
        outputfile.write("\n")
        outputfile.write("#Boundary constraint code\n")
        outputfile.write("BEGIN_BOUNDARY_CONSTRAINT_BLOCK\n")
        for BoundaryConstraintDefinition in self.BoundaryConstraintDefinitions:
            outputfile.write(BoundaryConstraintDefinition + "\n")
        outputfile.write("END_BOUNDARY_CONSTRAINT_BLOCK\n")
        outputfile.write("\n")
        
        outputfile.write("\n")
        outputfile.write("#Phase distance constraint code\n")
        outputfile.write("BEGIN_PHASE_DISTANCE_CONSTRAINT_BLOCK\n")
        for PhaseDistanceConstraintDefinition in self.PhaseDistanceConstraintDefinitions:
            outputfile.write(PhaseDistanceConstraintDefinition + "\n")
        outputfile.write("END_PHASE_DISTANCE_CONSTRAINT_BLOCK\n")
        outputfile.write("\n")

        outputfile.write("##manual inner-loop control code\n")
        outputfile.write("#Check derivatives against finite differencing?\n")
        outputfile.write("check_derivatives " + str(self.check_derivatives) + "\n")
        outputfile.write("#which inner loop solver to run?\n")  
        outputfile.write("#0: none, evaluate trialX\n")
        outputfile.write("#1: run MBH\n")
        outputfile.write("#2: run constrained DE\n")
        outputfile.write("#3: run NLP using trialX as initial guess\n")        
        outputfile.write("#4: filament walker\n")
        outputfile.write("run_inner_loop " + str(self.run_inner_loop) + "\n")
            
        if len(self.trialX) > 0:
            outputfile.write("#trial decision vector\n")
            outputfile.write("BEGIN_TRIALX\n")
            for entry in self.trialX:
                outputfile.write(entry[0] + "," + '%17.20f' % float(entry[1]) + "\n")
            outputfile.write("END_TRIALX\n")
        outputfile.write("\n")
            
        outputfile.write("#Enter any user data that should be appended to the .emtg file.\n")
        outputfile.write("#This is typically used in python wrappers\n")
        outputfile.write("user_data ")
        first_entry = True
        for entry in self.user_data.keys():
            if first_entry == True:
                first_entry = False
            else:
                outputfile.write(":")
            outputfile.write('("' + entry + '",')
            if isinstance(self.user_data[entry],str):
                outputfile.write("'" + str(self.user_data[entry]) + "')") 
            else:
                outputfile.write(str(self.user_data[entry]) + ")")    
        outputfile.write("\n")            
        outputfile.write("\n")
            
        outputfile.write("#end options file\n")

        outputfile.close()

    

    def ConvertDecisionVector(self):
        if self.PSFBstateRepresentation == 1:#SphericalRADEC, so we need to convert any SphericalAZFPA segments to SphericalRADEC
            from math import sin, cos, atan2, asin
            for Xindex in range(0, len(self.trialX)):
                description = self.trialX[Xindex][0]
                prefix = description.split(':')[0]

                if 'PSFB_Step' in description and 'left state AZ' in description:
                    #extract the SphericalAZFPA state
                    r = float(self.trialX[Xindex - 4][1])
                    RA = float(self.trialX[Xindex - 3][1])
                    DEC = float(self.trialX[Xindex - 2][1])
                    v = float(self.trialX[Xindex - 1][1])
                    AZ = float(self.trialX[Xindex][1])
                    FPA = float(self.trialX[Xindex + 1][1])
                    #convert to cartesian
                    cosRA = cos(RA)
                    sinRA = sin(RA)
                    cosDEC = cos(DEC)
                    sinDEC = sin(DEC)
                    cosAZ = cos(AZ)
                    sinAZ = sin(AZ)
                    cosFPA = cos(FPA)
                    sinFPA = sin(FPA)
                
                    xdot = -v * (sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA); 
                    ydot =  v * (sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA); 
                    zdot =  v * (cosFPA*sinDEC + cosDEC * cosAZ*sinFPA); 

                    #compute vRA and vDEC
                    vRA = atan2(ydot, xdot)
                    vDEC = asin(zdot / v)

                    #write vRA and vDEC into the decision vector
                    self.trialX[Xindex][0] = self.trialX[Xindex][0].replace('AZ','vRA')
                    self.trialX[Xindex][1] = str(vRA)
                    self.trialX[Xindex + 1][0] = self.trialX[Xindex + 1][0].replace('FPA','vDEC')
                    self.trialX[Xindex + 1][1] = str(vDEC)
        elif self.PSFBstateRepresentation == 2:#SphericalAZFPA - convert all SphericalRADEC to this format
            from math import sin, cos, atan2, asin, acos, pi
            import numpy
            for Xindex in range(0, len(self.trialX)):
                description = self.trialX[Xindex][0]
                prefix = description.split(':')[0]

                if 'PSFB_Step' in description and 'left state vRA' in description:
                    #extract the SphericalAZFPA state
                    r = float(self.trialX[Xindex - 4][1])
                    RA = float(self.trialX[Xindex - 3][1])
                    DEC = float(self.trialX[Xindex - 2][1])
                    v = float(self.trialX[Xindex - 1][1])
                    vRA = float(self.trialX[Xindex][1])
                    vDEC = float(self.trialX[Xindex + 1][1])
                    #convert to cartesian
                    cosRA = cos(RA)
                    sinRA = sin(RA)
                    cosDEC = cos(DEC)
                    sinDEC = sin(DEC)
                    cosvRA = cos(vRA)
                    sinvRA = sin(vRA)
                    cosvDEC = cos(vDEC)
                    sinvDEC = sin(vDEC)
                
                    x = r * cosRA * cosDEC
                    y = r * sinRA * cosDEC
                    z = r * sinDEC;
                    xdot = v * cosvRA * cosvDEC
                    ydot = v * sinvRA * cosvDEC
                    zdot = v * sinvDEC;

                    #compute AZ and FPA
                    FPA = acos( (x*xdot + y*ydot + z*zdot) / r / v )
        
                    #azimuth is complicated
                    xhat = numpy.matrix([cos(RA)*cos(DEC), sin(RA)*cos(DEC), sin(DEC)]).T
                    yhat = numpy.matrix([cos(RA + pi / 2.0), sin(RA + pi / 2.0), 0.0]).T
                    zhat = numpy.matrix([-cos(RA)*sin(DEC), -sin(RA)*sin(DEC), cos(DEC)]).T
                    R = numpy.hstack([xhat, yhat, zhat]).T
                    V = numpy.matrix([xdot, ydot, zdot]).T
                    Vprime = R * V
                    AZ = atan2(Vprime[1], Vprime[2])

                    #write vRA and vDEC into the decision vector
                    self.trialX[Xindex][0] = self.trialX[Xindex][0].replace('vRA','AZ')
                    self.trialX[Xindex][1] = str(AZ)
                    self.trialX[Xindex + 1][0] = self.trialX[Xindex + 1][0].replace('vDEC','FPA')
                    self.trialX[Xindex + 1][1] = str(FPA)

        if self.PeriapseBoundaryStateRepresentation == 1:#SphericalRADEC - convert all SphericalAZFPA to SphericalRADEC
            from math import sin, cos, atan2, asin
            for Xindex in range(0, len(self.trialX)):
                description = self.trialX[Xindex][0]
                prefix = description.split(':')[0]

                if 'Periapse' in description and 'event left state AZ' in description:
                    #extract the SphericalAZFPA state
                    r = float(self.trialX[Xindex - 4][1])
                    RA = float(self.trialX[Xindex - 3][1])
                    DEC = float(self.trialX[Xindex - 2][1])
                    v = float(self.trialX[Xindex - 1][1])
                    AZ = float(self.trialX[Xindex][1])
                    FPA = float(self.trialX[Xindex + 1][1])
                    #convert to cartesian
                    cosRA = cos(RA)
                    sinRA = sin(RA)
                    cosDEC = cos(DEC)
                    sinDEC = sin(DEC)
                    cosAZ = cos(AZ)
                    sinAZ = sin(AZ)
                    cosFPA = cos(FPA)
                    sinFPA = sin(FPA)
                
                    xdot = -v * (sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA); 
                    ydot =  v * (sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA); 
                    zdot =  v * (cosFPA*sinDEC + cosDEC * cosAZ*sinFPA); 

                    #compute vRA and vDEC
                    vRA = atan2(ydot, xdot)
                    vDEC = asin(zdot / v)

                    #write vRA and vDEC into the decision vector
                    self.trialX[Xindex][0] = self.trialX[Xindex][0].replace('AZ','vRA')
                    self.trialX[Xindex][1] = str(vRA)
                    self.trialX[Xindex + 1][0] = self.trialX[Xindex + 1][0].replace('FPA','vDEC')
                    self.trialX[Xindex + 1][1] = str(vDEC)
        elif self.PeriapseBoundaryStateRepresentation == 2:#SphericalAZFPA - convert all SphericalRADEC to this format
            from math import sin, cos, atan2, asin, acos, pi
            import numpy
            for Xindex in range(0, len(self.trialX)):
                description = self.trialX[Xindex][0]
                prefix = description.split(':')[0]

                if 'Periapse' in description and 'event left state vRA' in description:
                    #extract the SphericalAZFPA state
                    r = float(self.trialX[Xindex - 4][1])
                    RA = float(self.trialX[Xindex - 3][1])
                    DEC = float(self.trialX[Xindex - 2][1])
                    v = float(self.trialX[Xindex - 1][1])
                    vRA = float(self.trialX[Xindex][1])
                    vDEC = float(self.trialX[Xindex + 1][1])
                    #convert to cartesian
                    cosRA = cos(RA)
                    sinRA = sin(RA)
                    cosDEC = cos(DEC)
                    sinDEC = sin(DEC)
                    cosvRA = cos(vRA)
                    sinvRA = sin(vRA)
                    cosvDEC = cos(vDEC)
                    sinvDEC = sin(vDEC)
                
                    x = r * cosRA * cosDEC
                    y = r * sinRA * cosDEC
                    z = r * sinDEC;
                    xdot = v * cosvRA * cosvDEC
                    ydot = v * sinvRA * cosvDEC
                    zdot = v * sinvDEC;

                    #compute AZ and FPA
                    FPA = acos( (x*xdot + y*ydot + z*zdot) / r / v )
        
                    #azimuth is complicated
                    xhat = numpy.matrix([cos(RA)*cos(DEC), sin(RA)*cos(DEC), sin(DEC)]).T
                    yhat = numpy.matrix([cos(RA + pi / 2.0), sin(RA + pi / 2.0), 0.0]).T
                    zhat = numpy.matrix([-cos(RA)*sin(DEC), -sin(RA)*sin(DEC), cos(DEC)]).T
                    R = numpy.hstack([xhat, yhat, zhat]).T
                    V = numpy.matrix([xdot, ydot, zdot]).T
                    Vprime = R * V
                    AZ = atan2(Vprime[1], Vprime[2])

                    #write vRA and vDEC into the decision vector
                    self.trialX[Xindex][0] = self.trialX[Xindex][0].replace('vRA','AZ')
                    self.trialX[Xindex][1] = str(AZ)
                    self.trialX[Xindex + 1][0] = self.trialX[Xindex + 1][0].replace('vDEC','FPA')
                    self.trialX[Xindex + 1][1] = str(FPA)
        return