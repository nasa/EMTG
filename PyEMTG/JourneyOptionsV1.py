#Python class file for EMTG JourneyOptions
class JourneyOptions(object):

    #************************************************************************************constructor
    def __init__(self, phase_type=6):
        self.journey_names = 'New_Journey'
        self.journey_override_num_steps = 0
        self.journey_number_of_steps = 20
        self.journey_bounded_departure_date = 0
        self.journey_timebounded = 0#0: unbounded, 1: bounded flight time, 2: bounded arrival date
        self.journey_departure_date_bounds = [0.0]*2
        self.journey_wait_time_bounds = [0.0]*2#days
        self.journey_flight_time_bounds = [0.0]*2
        self.journey_arrival_date_bounds = [0.0]*2
        self.journey_initial_impulse_bounds = [0.0]*2 #in km/s
        self.journey_arrival_type = 1 #0: orbit insertion (use chemical Isp), 1: rendezvous (use chemical Isp), 2: flyby with bounded VHP, 3: low-thrust rendezvous (does not work if terminal phase is not low-thrust), 4: match v-infinity, 5: match v-infinity low-thrust, 6: E=0, 7: capture spiral
        self.journey_arrival_class = 0
        self.journey_arrival_ellipsoid_axes = [1.0e6, 1.0e6, 1.0e6]
        self.journey_departure_type = 0 #0: launch or direct insertion, 1: depart from parking orbit (you can use this one in place of a launch vehicle model, and the departure burn will be done with the EDS motor), 2: 'free' direct departure, i.e. do not burn to get the departure v_infinity, 3: Start from Sphere of Influence (use SOI angles chosen by previous journey's endpoint, i.e. after a spiral-out or fully modeled departure from parking orbit), 3: flyby (only valid for successive journeys), 4: flyby with fixed VHP, 5: escape spiral
        self.journey_departure_class = 0
        self.journey_departure_ellipsoid_axes = [1.0e6, 1.0e6, 1.0e6]
        self.journey_arrival_elements_type = 1 #0: cartesian, 1: COE
        self.journey_arrival_elements_frame = 1 #0: ICRF, 1: J2000_BCI, 2: J2000_BCF, 3: TrueOfDate_BCI, 4: TrueOfDate_BCF
        self.journey_arrival_elements = [0.0]*6 #a(km), e, i, RAAN, AOP, TA
        self.journey_arrival_elements_bounds = [0.0]*12
        self.journey_arrival_elements_vary_flag = [0]*6    
        self.journey_arrival_elements_reference_epoch = 51544.5
        self.AllowJourneyFreePointArrivalToPropagate = 0
        self.journey_departure_elements_type = 1 #0: cartesian, 1: COE
        self.journey_departure_elements_frame = 1 #0: ICRF, 1: J2000_BCI, 2: J2000_BCF, 3: TrueOfDate_BCI, 4: TrueOfDate_BCF
        self.journey_departure_elements = [0.0]*6 #a(km), e, i, RAAN, AOP, TA
        self.journey_departure_elements_bounds = [0.0]*12
        self.journey_departure_elements_vary_flag = [0]*6
        self.AllowJourneyFreePointDepartureToPropagate = 1
        self.journey_departure_elements_reference_epoch = 51544.5
        self.journey_override_flyby_altitude_bounds = False
        self.journey_flyby_altitude_bounds = [1000.0, 10000.0]
        self.journey_central_body = 'Sun' #spice names
        self.journey_final_velocity = [0.0, 20.0, 0.0] #in km/s
        self.journey_fixed_ending_mass_increment = 0 #in kg
        self.journey_fixed_starting_mass_increment = 0 #in kg
        self.journey_minimum_starting_mass_increment = 0 #in kg
        self.journey_maximum_starting_mass_increment = 0 #in kg
        self.journey_variable_mass_increment = 0 #whether or not the optimizer can choose the mass increment (ignored for non-positive mass increments)
        self.journey_maximum_initial_mass = 0
        self.journey_constrain_initial_mass = 0
        self.journey_arrival_declination_constraint_flag = 0
        self.journey_arrival_declination_bounds = [0.0]*2#in degrees
        self.journey_number_of_perturbation_bodies = 1
        self.journey_perturbation_bodies = [0]
        self.journey_escape_spiral_starting_radius = 6678#in km
        self.journey_capture_spiral_final_radius = 6678#in km
        self.journey_forced_terminal_coast = 0.0 #in days
        self.journey_forced_initial_coast = 0.0 #in days
    
        self.journey_end_deltav = 0 #in km/s, delta-v to be performed after the journey is over (i.e. prox ops, divert maneuvers, etc)
        self.journey_end_propulsion_system = 0 #monoprop

        self.journey_end_TCM = 0.0# in km/s

        #staging information
        self.journey_stage_after_departure = False
        self.journey_stage_before_arrival = False
        self.journey_stage_after_arrival = False

        
        self.journey_PeriapseArrival_override_altitude = False
        self.journey_PeriapseArrival_altitude_bounds = [300.0, 100000.0]
        self.journey_PeriapseDeparture_altitude_bounds = [300.0, 300.0]
        
        #sequence information
        self.number_of_phases = 1
        self.sequence = [0]*8
        self.destination_list = [1, 1]
        self.impulses_per_phase = [1]*9
        self.journey_enable_periapse_burns = [0]*8

        #outer loop selectable options settings
        self.outerloop_vary_journey_destination = 0
        self.outerloop_vary_journey_flyby_sequence = 0
        self.outerloop_vary_journey_arrival_type = 0
        self.outerloop_journey_destination_choices = [1]
        self.outerloop_journey_flyby_sequence_choices = [1]
        self.outerloop_journey_maximum_number_of_flybys = 8
        self.outerloop_journey_arrival_type_choices = [1]

        #solver information
        self.journey_freeze_decision_variables = 0
    
        self.journey_override_duty_cycle = False
        self.journey_duty_cycle = 0.0

        self.journey_override_PropagatorType = 0
        self.journey_PropagatorType = 0#Keplerian
        self.journey_override_integration_step_size = False
        self.journey_integration_step_size = 1.0

        self.journey_number_of_perturbation_bodies = 1
        self.journey_perturbation_bodies = [0]

        if phase_type in [2, 3, 4]:
            self.journey_arrival_type = 3
        self.phase_type = [phase_type] * 9

        #coast phase propagation
        self.journey_CoastPhaseMatchPointFraction = 0.5
        self.journey_CoastPhaseForwardIntegrationStepLength = 86400.0
        self.journey_CoastPhaseBackwardIntegrationStepLength = 86400.0

        #controls        
        self.journey_num_interior_control_points = 1



