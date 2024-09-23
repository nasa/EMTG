#conversion script from the .emtgopt v1 format to v2
#ingests a filename for the v1 file and returns a v2 MissionOptions object

import MissionOptions as MissionOptionsV2
import JourneyOptions as JourneyOptionsV2
import MissionOptionsV1

def Convert_emtgopt_v1_to_v2(fileName):
    #Step 1: determine the format of the input file
    input_file_version = 1
    with open(fileName, "r") as inputFile:
        if '#EMTGv9 .emtgopt file version 2' in inputFile.readline():
            input_file_version = 2

    #Step 2: if the file is already v2, create and return a MissionOptions object. Otherwise perform the conversion.
    if input_file_version == 2:
        print('Input file ' + fileName + ' is already in .emtgopt v2 format. No conversion is necessary. Returning v2 MissionOptions object')

        return MissionOptionsV2.MissionOptions(fileName)
    elif input_file_version == 1:
        originalOptions = MissionOptionsV1.MissionOptions(fileName)

        if originalOptions.success != 1:
            print('Input file ' + fileName + ' is not in a recognized format! Aborting.')
            return

        print('Input file ' + fileName + ' is in .emtgopt v1 format. Converting...')

        #instantiate a v2 MissionOptions object
        newOptions = MissionOptionsV2.MissionOptions()

        #convert global mission fields
        newOptions.mission_name  = originalOptions.mission_name 
        newOptions.objective_type  = originalOptions.objective_type 
        newOptions.objective_journey  = originalOptions.objective_journey 
        newOptions.include_initial_impulse_in_cost  = originalOptions.include_initial_impulse_in_cost 
        newOptions.global_timebounded  = originalOptions.global_timebounded 
        newOptions.launch_window_open_date  = originalOptions.launch_window_open_date 
        newOptions.total_flight_time_bounds  = originalOptions.total_flight_time_bounds 
        newOptions.DLA_bounds  = originalOptions.DLA_bounds 
        newOptions.RLA_bounds  = originalOptions.RLA_bounds 
        newOptions.mission_type  = originalOptions.mission_type 
        newOptions.NLP_solver_type  = originalOptions.NLP_solver_type 
        newOptions.NLP_solver_mode  = originalOptions.NLP_solver_mode 
        newOptions.quiet_NLP  = originalOptions.quiet_NLP 
        newOptions.ACE_feasible_point_finder  = originalOptions.ACE_feasible_point_finder 
        newOptions.MBH_always_write_archive  = originalOptions.MBH_always_write_archive 
        newOptions.MBH_archive_state_vector  = originalOptions.MBH_archive_state_vector 
        newOptions.MBH_max_not_improve  = originalOptions.MBH_max_not_improve 
        newOptions.MBH_max_trials  = originalOptions.MBH_max_trials 
        newOptions.MBH_max_run_time  = originalOptions.MBH_max_run_time 
        newOptions.MBH_max_step_size  = originalOptions.MBH_max_step_size 
        newOptions.MBH_hop_distribution  = originalOptions.MBH_hop_distribution 
        newOptions.MBH_Pareto_alpha  = originalOptions.MBH_Pareto_alpha 
        newOptions.MBH_write_every_improvement  = originalOptions.MBH_write_every_improvement 
        newOptions.MBH_time_hop_probability  = originalOptions.MBH_time_hop_probability 
        newOptions.snopt_feasibility_tolerance  = originalOptions.snopt_feasibility_tolerance 
        newOptions.snopt_optimality_tolerance  = originalOptions.snopt_optimality_tolerance 
        newOptions.NLP_max_step  = originalOptions.NLP_max_step 
        newOptions.snopt_major_iterations  = originalOptions.snopt_major_iterations 
        newOptions.snopt_minor_iterations  = originalOptions.snopt_minor_iterations 
        newOptions.snopt_max_run_time  = originalOptions.snopt_max_run_time
        newOptions.enable_NLP_chaperone  = originalOptions.enable_NLP_chaperone 
        newOptions.seed_MBH  = originalOptions.seed_MBH 
        newOptions.skip_first_nlp_run  = originalOptions.skip_first_nlp_run 
        newOptions.NLP_stop_on_goal_attain  = originalOptions.NLP_stop_on_goal_attain 
        newOptions.NLP_objective_goal  = originalOptions.NLP_objective_goal 
        newOptions.MBH_RNG_seed  = originalOptions.MBH_RNG_seed 
        newOptions.print_NLP_movie_frames  = originalOptions.print_NLP_movie_frames 
        newOptions.quiet_basinhopping  = originalOptions.quiet_basinhopping 
        newOptions.SPICE_leap_seconds_kernel  = originalOptions.SPICE_leap_seconds_kernel 
        newOptions.SPICE_reference_frame_kernel  = originalOptions.SPICE_reference_frame_kernel 
        newOptions.universe_folder  = originalOptions.universe_folder 
        newOptions.ephemeris_source  = originalOptions.ephemeris_source
        newOptions.propagatorType  = originalOptions.propagatorType 
        newOptions.integratorType  = originalOptions.integratorType 
        newOptions.integrator_tolerance  = originalOptions.integrator_tolerance 
        newOptions.integration_time_step_size  = originalOptions.integration_time_step_size 
        newOptions.num_timesteps  = originalOptions.num_timesteps 
        newOptions.spiral_segments  = originalOptions.spiral_segments 
        newOptions.allow_initial_mass_to_vary  = originalOptions.allow_initial_mass_to_vary 
        newOptions.maximum_mass  = originalOptions.maximum_mass 
        newOptions.IspLT  = originalOptions.IspLT 
        newOptions.IspLT_minimum  = originalOptions.IspLT_minimum 
        newOptions.IspChem  = originalOptions.IspChem 
        newOptions.Thrust  = originalOptions.Thrust 
        newOptions.LV_margin  = originalOptions.LV_margin 
        newOptions.LV_adapter_mass  = originalOptions.LV_adapter_mass 
        newOptions.engine_type  = originalOptions.engine_type 
        newOptions.number_of_electric_propulsion_systems  = originalOptions.number_of_electric_propulsion_systems 
        newOptions.engine_duty_cycle  = originalOptions.engine_duty_cycle 
        newOptions.duty_cycle_type  = originalOptions.duty_cycle_type 
        newOptions.thrust_scale_factor  = originalOptions.thrust_scale_factor 
        newOptions.power_at_1_AU  = originalOptions.power_at_1_AU 
        newOptions.power_source_type  = originalOptions.power_source_type 
        newOptions.solar_power_model_type  = originalOptions.solar_power_model_type 
        newOptions.solar_power_gamma  = originalOptions.solar_power_gamma 
        newOptions.power_margin  = originalOptions.power_margin 
        newOptions.power_decay_rate  = originalOptions.power_decay_rate 
        newOptions.throttle_sharpness  = originalOptions.throttle_sharpness 
        newOptions.throttle_logic_mode  = originalOptions.throttle_logic_mode 
        newOptions.spacecraft_power_coefficients  = originalOptions.spacecraft_power_coefficients 
        newOptions.engine_input_thrust_coefficients  = originalOptions.engine_input_thrust_coefficients 
        newOptions.engine_input_mass_flow_rate_coefficients  = originalOptions.engine_input_mass_flow_rate_coefficients 
        newOptions.engine_input_power_bounds  = originalOptions.engine_input_power_bounds 
        newOptions.user_defined_engine_efficiency  = originalOptions.user_defined_engine_efficiency 
        newOptions.spacecraft_power_model_type  = originalOptions.spacecraft_power_model_type 
        newOptions.TCM_Isp  = originalOptions.TCM_Isp 
        newOptions.TCM_post_launch  = originalOptions.TCM_post_launch 
        newOptions.TCM_pre_flyby  = originalOptions.TCM_pre_flyby 
        newOptions.TCM_maneuver_fraction  = originalOptions.TCM_maneuver_fraction 
        newOptions.trackACS  = originalOptions.trackACS 
        newOptions.ACS_kg_per_day  = originalOptions.ACS_kg_per_day 
        newOptions.final_mass_constraint_bounds  = originalOptions.final_mass_constraint_bounds 
        newOptions.constrain_final_mass  = originalOptions.constrain_final_mass 
        newOptions.constrain_dry_mass  = originalOptions.constrain_dry_mass 
        newOptions.enable_electric_propellant_tank_constraint  = originalOptions.enable_electric_propellant_tank_constraint 
        newOptions.maximum_electric_propellant  = originalOptions.maximum_electric_propellant 
        newOptions.electric_propellant_margin  = originalOptions.electric_propellant_margin 
        newOptions.enable_chemical_propellant_tank_constraint  = originalOptions.enable_chemical_propellant_tank_constraint 
        newOptions.maximum_chemical_fuel  = originalOptions.maximum_chemical_fuel 
        newOptions.maximum_chemical_oxidizer  = originalOptions.maximum_chemical_oxidizer 
        newOptions.bipropellant_mixture_ratio  = originalOptions.bipropellant_mixture_ratio 
        newOptions.chemical_propellant_margin  = originalOptions.chemical_propellant_margin 
        newOptions.SpacecraftModelInput  = originalOptions.SpacecraftModelInput 
        newOptions.HardwarePath  = originalOptions.HardwarePath 
        newOptions.ThrottleTableFile  = originalOptions.ThrottleTableFile 
        newOptions.LaunchVehicleLibraryFile  = originalOptions.LaunchVehicleLibraryFile 
        newOptions.PowerSystemsLibraryFile  = originalOptions.PowerSystemsLibraryFile 
        newOptions.PropulsionSystemsLibraryFile  = originalOptions.PropulsionSystemsLibraryFile 
        newOptions.SpacecraftOptionsFile  = originalOptions.SpacecraftOptionsFile 
        newOptions.LaunchVehicleKey  = originalOptions.LaunchVehicleKey 
        newOptions.PowerSystemKey  = originalOptions.PowerSystemKey 
        newOptions.ElectricPropulsionSystemKey  = originalOptions.ElectricPropulsionSystemKey 
        newOptions.ChemicalPropulsionSystemKey  = originalOptions.ChemicalPropulsionSystemKey 
        newOptions.perturb_SRP  = originalOptions.perturb_SRP 
        newOptions.perturb_thirdbody  = originalOptions.perturb_thirdbody 
        newOptions.perturb_J2  = originalOptions.perturb_J2 
        newOptions.spacecraft_area  = originalOptions.spacecraft_area 
        newOptions.coefficient_of_reflectivity  = originalOptions.coefficient_of_reflectivity 
        newOptions.forced_post_launch_coast  = originalOptions.forced_post_launch_coast 
        newOptions.forced_pre_flyby_coast  = originalOptions.forced_pre_flyby_coast 
        newOptions.forced_post_flyby_coast  = originalOptions.forced_post_flyby_coast 
        newOptions.waypoint_file_path  = originalOptions.waypoint_file_path 
        newOptions.covariance_file_path  = originalOptions.covariance_file_path 
        newOptions.ParallelShootingStateRepresentation  = originalOptions.PSFBstateRepresentation 
        newOptions.PeriapseBoundaryStateRepresentation  = originalOptions.PeriapseBoundaryStateRepresentation 
        newOptions.output_file_frame  = originalOptions.output_file_frame
        newOptions.output_dormant_journeys  = originalOptions.output_dormant_journeys 
        newOptions.post_mission_wait_time  = originalOptions.post_mission_wait_time 
        newOptions.override_working_directory  = originalOptions.override_working_directory 
        newOptions.forced_working_directory  = originalOptions.forced_working_directory 
        newOptions.override_mission_subfolder  = originalOptions.override_mission_subfolder 
        newOptions.forced_mission_subfolder  = originalOptions.forced_mission_subfolder 
        newOptions.short_output_file_names  = originalOptions.short_output_file_names 
        newOptions.generate_forward_integrated_ephemeris  = originalOptions.generate_forward_integrated_ephemeris 
        newOptions.add_control_switch_line_to_ephemeris  = originalOptions.add_control_switch_line_to_ephemeris 
        newOptions.append_mass_to_ephemeris_output  = originalOptions.append_mass_to_ephemeris_output 
        newOptions.append_control_to_ephemeris_output  = originalOptions.append_control_to_ephemeris_output 
        newOptions.append_thrust_to_ephemeris_output  = originalOptions.append_thrust_to_ephemeris_output 
        newOptions.append_mdot_to_ephemeris_output  = originalOptions.append_mdot_to_ephemeris_output 
        newOptions.append_Isp_to_ephemeris_output  = originalOptions.append_Isp_to_ephemeris_output 
        newOptions.append_active_power_to_ephemeris_output  = originalOptions.append_active_power_to_ephemeris_output 
        newOptions.append_number_of_active_engines_to_ephemeris_output  = originalOptions.append_number_of_active_engines_to_ephemeris_output 
        newOptions.append_throttle_level_to_ephemeris_output  = originalOptions.append_throttle_level_to_ephemeris_output 
        newOptions.background_mode  = originalOptions.background_mode 
        newOptions.output_STMs  = originalOptions.output_STMs 
        newOptions.output_maneuver_and_target_spec_files  = originalOptions.output_maneuver_and_target_spec_files 
        newOptions.stop_after_journey  = originalOptions.stop_after_journey 
        newOptions.run_inner_loop  = originalOptions.run_inner_loop 
        newOptions.check_derivatives  = originalOptions.check_derivatives 
        newOptions.user_data  = originalOptions.user_data 

        #convert journey-specific fields
        newOptions.Journeys = []
        for originalJourneyOptions in originalOptions.Journeys:            
            newOptions.Journeys.append(JourneyOptionsV2.JourneyOptions())

            #options that don't change
            newOptions.Journeys[-1].journey_central_body = originalJourneyOptions.journey_central_body
            newOptions.Journeys[-1].journey_end_TCM = originalJourneyOptions.journey_end_TCM
            newOptions.Journeys[-1].journey_end_deltav = originalJourneyOptions.journey_end_deltav
            newOptions.Journeys[-1].journey_end_propulsion_system = originalJourneyOptions.journey_end_propulsion_system
            newOptions.Journeys[-1].destination_list = originalJourneyOptions.destination_list
            newOptions.Journeys[-1].AllowJourneyFreePointDepartureToPropagate = originalJourneyOptions.AllowJourneyFreePointDepartureToPropagate
            newOptions.Journeys[-1].AllowJourneyFreePointArrivalToPropagate = originalJourneyOptions.AllowJourneyFreePointArrivalToPropagate

            #options that change format
            newOptions.Journeys[-1].sequence = [b for b in originalJourneyOptions.sequence if b > 0]
            newOptions.Journeys[-1].phase_type = originalJourneyOptions.phase_type[0]
            newOptions.Journeys[-1].impulses_per_phase = originalJourneyOptions.impulses_per_phase[0]
            newOptions.Journeys[-1].enable_periapse_burns = originalJourneyOptions.journey_enable_periapse_burns[0]

            #options that change name
            newOptions.Journeys[-1].journey_name = originalJourneyOptions.journey_names
            newOptions.Journeys[-1].override_num_steps = originalJourneyOptions.journey_override_num_steps
            newOptions.Journeys[-1].number_of_steps = originalJourneyOptions.journey_number_of_steps
            newOptions.Journeys[-1].override_duty_cycle = originalJourneyOptions.journey_override_duty_cycle
            newOptions.Journeys[-1].duty_cycle = originalJourneyOptions.journey_duty_cycle
            newOptions.Journeys[-1].override_PropagatorType = originalJourneyOptions.journey_override_PropagatorType
            newOptions.Journeys[-1].propagatorType = originalJourneyOptions.journey_PropagatorType
            newOptions.Journeys[-1].override_integration_step_size = originalJourneyOptions.journey_override_integration_step_size
            newOptions.Journeys[-1].integration_step_size = originalJourneyOptions.journey_integration_step_size
            newOptions.Journeys[-1].override_flyby_altitude_bounds = originalJourneyOptions.journey_override_flyby_altitude_bounds
            newOptions.Journeys[-1].flyby_altitude_bounds = originalJourneyOptions.journey_flyby_altitude_bounds
            newOptions.Journeys[-1].PeriapseArrival_override_altitude = originalJourneyOptions.journey_PeriapseArrival_override_altitude
            newOptions.Journeys[-1].PeriapseArrival_altitude_bounds = originalJourneyOptions.journey_PeriapseArrival_altitude_bounds
            newOptions.Journeys[-1].PeriapseDeparture_altitude_bounds = originalJourneyOptions.journey_PeriapseDeparture_altitude_bounds
            newOptions.Journeys[-1].num_interior_control_points = originalJourneyOptions.journey_num_interior_control_points
            newOptions.Journeys[-1].CoastPhaseMatchPointFraction = originalJourneyOptions.journey_CoastPhaseMatchPointFraction
            newOptions.Journeys[-1].CoastPhaseForwardIntegrationStepLength = originalJourneyOptions.journey_CoastPhaseForwardIntegrationStepLength
            newOptions.Journeys[-1].CoastPhaseBackwardIntegrationStepLength = originalJourneyOptions.journey_CoastPhaseBackwardIntegrationStepLength
            newOptions.Journeys[-1].bounded_departure_date = originalJourneyOptions.journey_bounded_departure_date
            newOptions.Journeys[-1].timebounded = originalJourneyOptions.journey_timebounded
            newOptions.Journeys[-1].departure_date_bounds = originalJourneyOptions.journey_departure_date_bounds
            newOptions.Journeys[-1].wait_time_bounds = originalJourneyOptions.journey_wait_time_bounds
            newOptions.Journeys[-1].flight_time_bounds = originalJourneyOptions.journey_flight_time_bounds
            newOptions.Journeys[-1].arrival_date_bounds = originalJourneyOptions.journey_arrival_date_bounds
            newOptions.Journeys[-1].departure_type = originalJourneyOptions.journey_departure_type
            newOptions.Journeys[-1].initial_impulse_bounds = originalJourneyOptions.journey_initial_impulse_bounds
            newOptions.Journeys[-1].departure_elements_vary_flag = originalJourneyOptions.journey_departure_elements_vary_flag
            newOptions.Journeys[-1].departure_elements = originalJourneyOptions.journey_departure_elements
            newOptions.Journeys[-1].departure_elements_bounds = originalJourneyOptions.journey_departure_elements_bounds
            newOptions.Journeys[-1].departure_elements_reference_epoch = originalJourneyOptions.journey_departure_elements_reference_epoch
            newOptions.Journeys[-1].departure_elements_frame = originalJourneyOptions.journey_departure_elements_frame
            newOptions.Journeys[-1].maximum_starting_mass_increment = originalJourneyOptions.journey_maximum_starting_mass_increment
            newOptions.Journeys[-1].minimum_starting_mass_increment = originalJourneyOptions.journey_minimum_starting_mass_increment
            newOptions.Journeys[-1].fixed_starting_mass_increment = originalJourneyOptions.journey_fixed_starting_mass_increment
            newOptions.Journeys[-1].fixed_ending_mass_increment = originalJourneyOptions.journey_fixed_ending_mass_increment
            newOptions.Journeys[-1].variable_mass_increment = originalJourneyOptions.journey_variable_mass_increment
            newOptions.Journeys[-1].constrain_initial_mass = originalJourneyOptions.journey_constrain_initial_mass
            newOptions.Journeys[-1].maximum_initial_mass = originalJourneyOptions.journey_maximum_initial_mass
            newOptions.Journeys[-1].departure_class = originalJourneyOptions.journey_departure_class
            newOptions.Journeys[-1].departure_ellipsoid_axes = originalJourneyOptions.journey_departure_ellipsoid_axes
            newOptions.Journeys[-1].arrival_type = originalJourneyOptions.journey_arrival_type
            newOptions.Journeys[-1].arrival_elements_vary_flag = originalJourneyOptions.journey_arrival_elements_vary_flag
            newOptions.Journeys[-1].arrival_elements = originalJourneyOptions.journey_arrival_elements
            newOptions.Journeys[-1].arrival_elements_bounds = originalJourneyOptions.journey_arrival_elements_bounds
            newOptions.Journeys[-1].arrival_elements_reference_epoch = originalJourneyOptions.journey_arrival_elements_reference_epoch
            newOptions.Journeys[-1].arrival_elements_frame = originalJourneyOptions.journey_arrival_elements_frame
            newOptions.Journeys[-1].final_velocity = originalJourneyOptions.journey_final_velocity
            newOptions.Journeys[-1].forced_terminal_coast = originalJourneyOptions.journey_forced_terminal_coast
            newOptions.Journeys[-1].forced_initial_coast = originalJourneyOptions.journey_forced_initial_coast
            newOptions.Journeys[-1].arrival_class = originalJourneyOptions.journey_arrival_class
            newOptions.Journeys[-1].arrival_ellipsoid_axes = originalJourneyOptions.journey_arrival_ellipsoid_axes
            newOptions.Journeys[-1].escape_spiral_starting_radius = originalJourneyOptions.journey_escape_spiral_starting_radius
            newOptions.Journeys[-1].capture_spiral_final_radius = originalJourneyOptions.journey_capture_spiral_final_radius
            newOptions.Journeys[-1].perturbation_bodies = list(filter(lambda b: b != 0, originalJourneyOptions.journey_perturbation_bodies))
            newOptions.Journeys[-1].stage_after_departure = originalJourneyOptions.journey_stage_after_departure
            newOptions.Journeys[-1].stage_before_arrival = originalJourneyOptions.journey_stage_before_arrival
            newOptions.Journeys[-1].stage_after_arrival = originalJourneyOptions.journey_stage_after_arrival
            newOptions.Journeys[-1].freeze_decision_variables = originalJourneyOptions.journey_freeze_decision_variables

        newOptions.number_of_journeys = len(newOptions.Journeys)

        #initial guess - have to distribute to the journeys
        newOptions.trialX = originalOptions.trialX
        newOptions.DisassembleMasterDecisionVector()

        #script constraints - have to distribute to the journeys
        newOptions.ManeuverConstraintDefinitions = originalOptions.ManeuverConstraintDefinitions
        newOptions.BoundaryConstraintDefinitions = originalOptions.BoundaryConstraintDefinitions
        newOptions.PhaseDistanceConstraintDefinitions = originalOptions.PhaseDistanceConstraintDefinitions
        newOptions.DisassembleMasterConstraintVectors()

        print('Finished converting. Returning v2-formatted options object.')

        return newOptions

    else:
        print('Input file ' + fileName + ' is not in a recognized format! Aborting.')
        return