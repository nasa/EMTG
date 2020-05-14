# EMTG: Evolutionary Mission Trajectory Generator
# An open-source global optimization tool for preliminary mission design
# Provided by NASA Goddard Space Flight Center
#
# Copyright (c) 2013 - 2020 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Other Rights Reserved.
#
# Licensed under the NASA Open Source License (the "License"); 
# You may not use this file except in compliance with the License. 
# You may obtain a copy of the License at:
# https://opensource.org/licenses/NASA-1.3
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
# express or implied.   See the License for the specific language
# governing permissions and limitations under the License.

import wx
import wx.adv
import wx.lib.scrolledpanel
import platform

import JourneyOptions as JO
import Universe
import BodyPicker
import DepartureElementsPanel
import ArrivalElementsPanel
import JourneyStagingPanel

import os

class JourneyOptionsPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, missionoptions):
        self.missionoptions = missionoptions
        self.parent = parent
        
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)

        self.JourneyList = []

        self.JourneySelectBox = wx.ListBox(self, -1, choices=self.JourneyList, size=(300,200), style=wx.LB_SINGLE)
        self.btnAddNewJourney = wx.Button(self, -1, "New Journey", size=(200,-1))
        self.btnDeleteJourney = wx.Button(self, -1, "Delete Journey", size=(200,-1))
        self.btnMoveJourneyUp = wx.Button(self, -1, "Move Journey Up", size=(200,-1))
        self.btnMoveJourneyDown = wx.Button(self, -1, "Move Journey Down", size=(200,-1))

        buttonstacksizer = wx.BoxSizer(wx.VERTICAL)
        buttonstacksizer.AddMany([self.btnAddNewJourney, self.btnDeleteJourney, self.btnMoveJourneyUp, self.btnMoveJourneyDown])

        JourneySelectionSizer = wx.BoxSizer(wx.HORIZONTAL)
        JourneySelectionSizer.Add(self.JourneySelectBox)
        JourneySelectionSizer.AddSpacer(5)
        JourneySelectionSizer.Add(buttonstacksizer)

        self.lbljourney_name = wx.StaticText(self, -1, "Journey name")
        self.txtjourney_name = wx.TextCtrl(self, -1, "journey_name", size=(300,-1))

        self.lblPhaseType = wx.StaticText(self, -1, "Phase type")
        phasetypes = ['MGALTS (not yet implemented)', 'FBLTS (not yet implemented)', 'MGALT', 'FBLT', 'PSBI', 'PSFB', 'MGAnDSMs', 'CoastPhase', 'SundmanCoastPhase']
        self.cmbPhaseType = wx.ComboBox(self, -1, choices=phasetypes, style=wx.CB_READONLY)

        self.lbloverride_num_steps = wx.StaticText(self, -1, "Override journey number of time steps")
        self.chkoverride_num_steps = wx.CheckBox(self, -1)
        self.lblnumber_of_steps = wx.StaticText(self, -1, "Number of time steps")
        self.txtnumber_of_steps = wx.TextCtrl(self, -1, "number_of_steps")

        self.lblnum_interior_control_points = wx.StaticText(self, -1, "Number of interior control points")
        self.txtnum_interior_control_points = wx.TextCtrl(self, -1, "num_interior_control_points")
        
        self.lblimpulses_per_phase = wx.StaticText(self, -1, "Impulses per phase")
        self.txtimpulses_per_phase = wx.TextCtrl(self, -1, "impulses_per_phase")

        self.lblforce_unit_magnitude_control = wx.StaticText(self, -1, "Force 100% control in this journey?")
        self.cmbforce_unit_magnitude_control = wx.ComboBox(self, -1, choices=['up to unit magnitude', 'unit magnitude', 'zero magnitude'], style=wx.CB_READONLY)
        
        self.lblforce_fixed_inertial_control = wx.StaticText(self, -1, "Force single intertial control vector across each phase?")
        self.chkforce_fixed_inertial_control = wx.CheckBox(self, -1)

        self.lbljourney_central_body = wx.StaticText(self, -1, "Central body")
        self.txtjourney_central_body = wx.TextCtrl(self, -1, "journey_central_body", size=(300,-1))
        self.btnjourney_central_body = wx.Button(self, -1, "...")
        journey_central_body_box = wx.BoxSizer(wx.HORIZONTAL)
        journey_central_body_box.Add(self.txtjourney_central_body)
        journey_central_body_box.AddSpacer(5)
        journey_central_body_box.Add(self.btnjourney_central_body)

        self.lbldestination_list = wx.StaticText(self, -1, "Destination list")
        self.txtdestination_list = wx.TextCtrl(self, -1, "destination_list")
        self.btndestination_list = wx.Button(self, -1, "...")
        destination_list_box = wx.BoxSizer(wx.HORIZONTAL)
        destination_list_box.Add(self.txtdestination_list)
        destination_list_box.AddSpacer(5)
        destination_list_box.Add(self.btndestination_list)
        
        self.lblfixed_ending_mass_increment = wx.StaticText(self, -1, "Fixed ending mass increment (kg)")
        self.txtfixed_ending_mass_increment = wx.TextCtrl(self, -1, "fixed_ending_mass_increment")
        
        self.lblfixed_starting_mass_increment = wx.StaticText(self, -1, "Fixed starting mass increment (kg)")
        self.txtfixed_starting_mass_increment = wx.TextCtrl(self, -1, "fixed_starting_mass_increment")

        self.lblminimum_starting_mass_increment = wx.StaticText(self, -1, "Minimum starting mass increment (kg)")
        self.txtminimum_starting_mass_increment = wx.TextCtrl(self, -1, "minimum_starting_mass_increment")
        
        self.lblmaximum_starting_mass_increment = wx.StaticText(self, -1, "Maximum starting mass increment (kg)")
        self.txtmaximum_starting_mass_increment = wx.TextCtrl(self, -1, "maximum_starting_mass_increment")
        
        self.lblvariable_mass_increment = wx.StaticText(self, -1, "Variable mass increment")
        self.chkvariable_mass_increment = wx.CheckBox(self, -1)
        
        self.lblconstrain_initial_mass = wx.StaticText(self, -1, "Constrain initial mass")
        self.chkconstrain_initial_mass = wx.CheckBox(self, -1)
        
        self.lblmaximum_initial_mass = wx.StaticText(self, -1, "Maximum initial mass (kg)")
        self.txtmaximum_initial_mass = wx.TextCtrl(self, -1, "maximum_initial_mass")

        self.lblwait_time_bounds = wx.StaticText(self, -1, "Wait time bounds (days)")
        self.txtwait_time_bounds_lower = wx.TextCtrl(self, -1, "wait_time_bounds[0]")
        self.txtwait_time_bounds_upper = wx.TextCtrl(self, -1, "wait_time_bounds[1]")
        wait_time_sizer = wx.BoxSizer(wx.HORIZONTAL)
        wait_time_sizer.AddMany([self.txtwait_time_bounds_lower, self.txtwait_time_bounds_upper])

        self.lbltimebounded = wx.StaticText(self, -1, "Journey time bounds")
        time_bounds_choices = ['unbounded','bounded flight time','bounded arrival date','bounded aggregate flight time']
        self.cmbtimebounded = wx.ComboBox(self, -1, choices=time_bounds_choices, style=wx.CB_READONLY)

        self.lblflight_time_bounds = wx.StaticText(self, -1, "Journey flight time bounds (days)")
        self.txtflight_time_bounds_lower = wx.TextCtrl(self, -1, "flight_time_bounds[0]")
        self.txtflight_time_bounds_upper = wx.TextCtrl(self, -1, "flight_time_bounds[1]")
        flight_time_sizer = wx.BoxSizer(wx.HORIZONTAL)
        flight_time_sizer.AddMany([self.txtflight_time_bounds_lower, self.txtflight_time_bounds_upper])

        self.lblarrival_date_bounds = wx.StaticText(self, -1, "Journey arrival date bounds")
        self.txtarrival_date_bounds_lower = wx.TextCtrl(self, -1, "arrival_date_bounds[0]")
        self.txtarrival_date_bounds_upper = wx.TextCtrl(self, -1, "arrival_date_bounds[1]")
        self.ArrivalDateLowerCalendar = wx.adv.CalendarCtrl(self, -1)
        self.ArrivalDateUpperCalendar = wx.adv.CalendarCtrl(self, -1)
        arrival_date_sizer = wx.BoxSizer(wx.HORIZONTAL)
        arrival_date_sizer.AddMany([self.txtarrival_date_bounds_lower, self.ArrivalDateLowerCalendar, self.txtarrival_date_bounds_upper, self.ArrivalDateUpperCalendar])

        self.lblbounded_departure_date = wx.StaticText(self, -1, "Bounded journey departure date?")
        self.chkbounded_departure_date = wx.CheckBox(self, -1)

        self.lbldeparture_date_bounds = wx.StaticText(self, -1, "Journey departure date bounds")
        self.txtdeparture_date_bounds_lower = wx.TextCtrl(self, -1, "departure_date_bounds[0]")
        self.txtdeparture_date_bounds_upper = wx.TextCtrl(self, -1, "departure_date_bounds[1]")
        self.departureDateLowerCalendar = wx.adv.CalendarCtrl(self, -1)
        self.departureDateUpperCalendar = wx.adv.CalendarCtrl(self, -1)
        departure_date_sizer = wx.BoxSizer(wx.HORIZONTAL)
        departure_date_sizer.AddMany([self.txtdeparture_date_bounds_lower, self.departureDateLowerCalendar, self.txtdeparture_date_bounds_upper, self.departureDateUpperCalendar])

        self.lblinitial_impulse_bounds = wx.StaticText(self, -1, "Journey initial impulse bounds (km/s)")
        self.txtinitial_impulse_bounds_lower = wx.TextCtrl(self, -1, "initial_impulse_bounds[0]")
        self.txtinitial_impulse_bounds_upper = wx.TextCtrl(self, -1, "initial_impulse_bounds[1]")
        initial_impulse_sizer = wx.BoxSizer(wx.HORIZONTAL)
        initial_impulse_sizer.AddMany([self.txtinitial_impulse_bounds_lower, self.txtinitial_impulse_bounds_upper])

        self.lblforce_free_point_direct_insertion_along_velocity_vector = wx.StaticText(self, -1, "Force departure impulse along velocity vector?")
        self.chkforce_free_point_direct_insertion_along_velocity_vector = wx.CheckBox(self, -1)

        self.lbldeparture_type = wx.StaticText(self, -1, "Journey departure type")
        departure_type_choices = ['0: launch or direct insertion','1: depart from parking orbit','2: free direct departure',
                                        '3: flyby','4: flyby with fixed v-infinity-out','5: Spiral-out from circular orbit','6: zero-turn flyby (for small bodies)']
        self.cmbdeparture_type = wx.ComboBox(self, -1, choices=departure_type_choices, style=wx.CB_READONLY)

        self.lbldeparture_class = wx.StaticText(self, -1, "Journey departure class")
        departure_class_choices = ['0: Ephemeris-pegged', '1: Free point', '2: Ephemeris-referenced','3: Periapse']
        self.cmbdeparture_class = wx.ComboBox(self, -1, choices=departure_class_choices, style=wx.CB_READONLY)

        self.lbldeparture_ellipsoid_axes = wx.StaticText(self, -1, "Departure reference ellipsoid axes")
        self.txtdeparture_ellipsoid_axes = wx.TextCtrl(self, -1, "departure_ellipsoid_axes", size=(300, -1))
        self.btndeparture_ellipsoid_axes = wx.Button(self, -1, "Default")
        departure_ellipsoid_sizer = wx.BoxSizer(wx.HORIZONTAL)
        departure_ellipsoid_sizer.AddMany([self.txtdeparture_ellipsoid_axes, self.btndeparture_ellipsoid_axes])

        self.lblforced_initial_coast = wx.StaticText(self, -1, "Forced initial coast (days)")
        self.txtforced_initial_coast = wx.TextCtrl(self, -1, "forced_initial_coast")

        self.lbloverride_flyby_altitude_bounds = wx.StaticText(self, -1, "Override flyby altitude bounds for first phase")
        self.chkoverride_flyby_altitude_bounds = wx.CheckBox(self, -1)
        self.lblflyby_altitude_bounds = wx.StaticText(self, -1, "First phase flyby altitude bounds (km)")
        self.txtflyby_altitude_boundsLower = wx.TextCtrl(self, -1, "flyby_altitude_bounds[0]")
        self.txtflyby_altitude_boundsUpper = wx.TextCtrl(self, -1, "flyby_altitude_bounds[1]")
        flyby_altitude_box = wx.BoxSizer(wx.HORIZONTAL)
        flyby_altitude_box.AddMany([self.txtflyby_altitude_boundsLower, self.txtflyby_altitude_boundsUpper])

        self.lblescape_spiral_starting_radius = wx.StaticText(self, -1, "Orbital radius for beginning of escape spiral (km)")
        self.txtescape_spiral_starting_radius = wx.TextCtrl(self, -1, "escape_spiral_starting_radius")
        self.lblescape_spiral_final_radius = wx.StaticText(self, -1, "Orbital radius for end of escape spiral (km)")
        self.txtescape_spiral_final_radius = wx.TextCtrl(self, -1, "escape_spiral_final_radius")

        self.lblzero_turn_flyby_distance = wx.StaticText(self, -1, "Flyby distance-from-center (km)")
        self.txtzero_turn_flyby_distance = wx.TextCtrl(self, -1, "zero_turn_flyby_distance")

        self.lblarrival_type = wx.StaticText(self, -1, "Journey arrival type")
        arrival_type_choices = ['0: insertion into parking orbit (use chemical Isp)','1: rendezvous (with chemical maneuver)','2: intercept with bounded V_infinity',
                                        '3: rendezvous (no maneuver)','4: match final v-infinity vector (with chemical maneuver)',
                                        '5: match final v-infinity vector (no maneuver)','6: capture spiral', '7: momentum transfer']
        self.cmbarrival_type = wx.ComboBox(self, -1, choices=arrival_type_choices, style=wx.CB_READONLY)

        self.lblimpact_momentum_enhancement_factor = wx.StaticText(self, -1, "Impact momentum enhancement factor (beta)")
        self.txtimpact_momentum_enhancement_factor = wx.TextCtrl(self, -1, "impact_momentum_enhancement_factor")

        self.lblarrival_class = wx.StaticText(self, -1, "Journey arrival class")
        arrival_class_choices = ['0: Ephemeris-pegged', '1: Free point', '2: Ephemeris-referenced','3: Periapse']
        self.cmbarrival_class = wx.ComboBox(self, -1, choices=arrival_class_choices, style=wx.CB_READONLY)

        self.lblephemeris_pegged_orbit_insertion_SMA = wx.StaticText(self, -1, "Insertion orbit SMA (km)")
        self.txtephemeris_pegged_orbit_insertion_SMA = wx.TextCtrl(self, -1, "Insertion orbit SMA (km)")
        self.lblephemeris_pegged_orbit_insertion_ECC = wx.StaticText(self, -1, "Insertion orbit ECC")
        self.txtephemeris_pegged_orbit_insertion_ECC = wx.TextCtrl(self, -1, "Insertion orbit ECC")

        self.lblarrival_ellipsoid_axes = wx.StaticText(self, -1, "arrival reference ellipsoid axes")
        self.txtarrival_ellipsoid_axes = wx.TextCtrl(self, -1, "arrival_ellipsoid_axes", size=(300, -1))
        self.btnarrival_ellipsoid_axes = wx.Button(self, -1, "Default")
        arrival_ellipsoid_sizer = wx.BoxSizer(wx.HORIZONTAL)
        arrival_ellipsoid_sizer.AddMany([self.txtarrival_ellipsoid_axes, self.btnarrival_ellipsoid_axes])
        
        self.lblcapture_spiral_starting_radius = wx.StaticText(self, -1, "Orbital radius for beginning of capture spiral (km)")
        self.txtcapture_spiral_starting_radius = wx.TextCtrl(self, -1, "capture_spiral_starting_radius")
        self.lblcapture_spiral_final_radius = wx.StaticText(self, -1, "Orbital radius for end of capture spiral (km)")
        self.txtcapture_spiral_final_radius = wx.TextCtrl(self, -1, "capture_spiral_final_radius")
        
        self.lblforced_terminal_coast = wx.StaticText(self, -1, "Forced terminal coast (days)")
        self.txtforced_terminal_coast = wx.TextCtrl(self, -1, "forced_terminal_coast")

        self.lblfinal_velocity = wx.StaticText(self, -1, "Journey final velocity vector (km/s)")
        self.txtfinal_velocity0 = wx.TextCtrl(self, -1, "final_velocity[0]")
        self.txtfinal_velocity1 = wx.TextCtrl(self, -1, "final_velocity[1]")
        self.txtfinal_velocity2 = wx.TextCtrl(self, -1, "final_velocity[2]")
        final_velocity_box = wx.BoxSizer(wx.HORIZONTAL)
        final_velocity_box.AddMany([self.txtfinal_velocity0, self.txtfinal_velocity1, self.txtfinal_velocity2])

        self.lblFreePointArrival_print_target_spec = wx.StaticText(self, -1, "Print a target spec line for free point arrival?")
        self.chkFreePointArrival_print_target_spec = wx.CheckBox(self, -1)
 
        
        self.lbljourney_end_deltav = wx.StaticText(self, -1, "Journey-end delta-v (km/s)")
        self.txtjourney_end_deltav = wx.TextCtrl(self, -1, "journey_end_deltav")
        self.lbljourney_end_TCM = wx.StaticText(self, -1, "Journey-end TCM magnitude (km/s)")
        self.txtjourney_end_TCM = wx.TextCtrl(self, -1, "journey_end_TCM")

        self.lbloverride_duty_cycle = wx.StaticText(self,-1,"Override this journey's duty cycle")
        self.chkoverride_duty_cycle = wx.CheckBox(self,-1)
        self.lblduty_cycle = wx.StaticText(self,-1,"Journey duty cycle")
        self.txtduty_cycle = wx.TextCtrl(self,-1,"txtduty_cycle")

        self.lblsequence = wx.StaticText(self, -1, "Flyby sequence")
        self.txtsequence = wx.TextCtrl(self, -1, "sequence", size=(300,60), style=wx.TE_MULTILINE)
        self.btnsequence = wx.Button(self, -1, "...")
        sequence_box = wx.BoxSizer(wx.HORIZONTAL)
        sequence_box.Add(self.txtsequence)
        sequence_box.AddSpacer(5)
        sequence_box.Add(self.btnsequence)

        self.lblperiapse_burns = wx.StaticText(self, -1, "Enable powered flybys?")
        self.chkperiapse_burns = wx.CheckBox(self, -1)

        self.lblperturbation_bodies = wx.StaticText(self, -1, "Perturbation_bodies")
        self.txtperturbation_bodies = wx.TextCtrl(self, -1, "perturbation_bodies", size=(300,-1))
        self.btnperturbation_bodies = wx.Button(self, -1, "...")
        perturbation_bodies_box = wx.BoxSizer(wx.HORIZONTAL)
        perturbation_bodies_box.Add(self.txtperturbation_bodies)
        perturbation_bodies_box.AddSpacer(5)
        perturbation_bodies_box.Add(self.btnperturbation_bodies)

        #self.btnEditJourneyDistanceConstraints = wx.Button(self, -1, "Edit journey distance constraints")

        self.lblfreeze_decision_variables = wx.StaticText(self, -1, "Freeze this journey's decision variables?")
        self.chkfreeze_decision_variables = wx.CheckBox(self, -1)

        self.lbloverride_integration_step_size = wx.StaticText(self, -1, "Override journey's integration step size?")
        self.chkoverride_integration_step_size = wx.CheckBox(self, -1)
        self.lblintegration_step_size = wx.StaticText(self, -1, "Integration step size (seconds)")
        self.txtintegration_step_size = wx.TextCtrl(self, -1, "integration_step_size")
        
        self.lblCoastPhaseMatchPointFraction = wx.StaticText(self, -1, "Where do you want to place the match point? (fraction)")
        self.txtCoastPhaseMatchPointFraction = wx.TextCtrl(self, -1, "CoastPhaseMatchPointFraction")
        self.lblCoastPhaseForwardIntegrationStepLength = wx.StaticText(self, -1, "Integration step size for the forward half-phase (seconds)")
        self.txtCoastPhaseForwardIntegrationStepLength = wx.TextCtrl(self, -1, "CoastPhaseForwardIntegrationStepLength ")
        self.lblCoastPhaseBackwardIntegrationStepLength = wx.StaticText(self, -1, "Integration step size for the backward half-phase (seconds)")
        self.txtCoastPhaseBackwardIntegrationStepLength = wx.TextCtrl(self, -1, "CoastPhaseBackwardIntegrationStepLength")

        
        self.lbloverride_PropagatorType = wx.StaticText(self, -1, "Override journey's propagator type?")
        self.chkoverride_PropagatorType = wx.CheckBox(self, -1)
        self.lblpropagatorType = wx.StaticText(self, -1, "Propagator type")
        propagatorType_choices = ["Keplerian", "Integrator"]
        self.cmbpropagatorType = wx.ComboBox(self, -1, choices=propagatorType_choices, style=wx.CB_READONLY)

        ##############
        # drag options
        ##############

        # enable drag
        self.lblEnableDrag = wx.StaticText(self, -1, "Enable aerodynamic drag?")
        self.chkEnableDrag = wx.CheckBox(self, -1)

        # drag area
        self.lblSpacecraftDragArea = wx.StaticText(self, -1, "Spacecraft drag area (m^2)")
        self.txtSpacecraftDragArea = wx.TextCtrl(self, -1, "spacecraftDragArea")

        # drag coefficient
        self.lblSpacecraftDragCoefficient = wx.StaticText(self, -1, "Spacecraft drag coefficient")
        self.txtSpacecraftDragCoefficient = wx.TextCtrl(self, -1, "spacecraftDragCoefficient")

        # density model
        self.lblAtmosphericDensityModel = wx.StaticText(self, -1, "Atmospheric density model")
        self.atmosphericDensityModelChoices = ["Exponential"]
        self.cmbAtmosphericDensityModel = wx.ComboBox(self, -1, choices=self.atmosphericDensityModelChoices, style=wx.CB_READONLY)

        # density model data file
        self.lblAtmosphericDensityModelDataFile = wx.StaticText(self, -1, "Atmospheric density model data file")
        self.txtAtmosphericDensityModelDataFile = wx.TextCtrl(self, -1, "atmosphericDensityModelDataFile", size=(300,-1))
        self.btnAtmosphericDensityModelDataFile = wx.Button(self, -1, "...")
        self.atmosphericDensityModelDataFileBox = wx.BoxSizer(wx.HORIZONTAL)
        self.atmosphericDensityModelDataFileBox.Add(self.txtAtmosphericDensityModelDataFile)
        self.atmosphericDensityModelDataFileBox.AddSpacer(5)
        self.atmosphericDensityModelDataFileBox.Add(self.btnAtmosphericDensityModelDataFile)

        JourneyInformationGrid = wx.FlexGridSizer(100,2,5,5)
        JourneyInformationGrid.AddMany([self.lbljourney_name, self.txtjourney_name,
                                        self.lblPhaseType, self.cmbPhaseType,
                                        self.lblfreeze_decision_variables, self.chkfreeze_decision_variables,
                                        self.lbloverride_num_steps, self.chkoverride_num_steps,
                                        self.lblnumber_of_steps, self.txtnumber_of_steps,
                                        self.lblnum_interior_control_points, self.txtnum_interior_control_points,
                                        self.lblimpulses_per_phase, self.txtimpulses_per_phase,
                                        self.lblforce_unit_magnitude_control, self.cmbforce_unit_magnitude_control,
                                        self.lblforce_fixed_inertial_control, self.chkforce_fixed_inertial_control,
                                        self.lbljourney_central_body, journey_central_body_box,
                                        self.lbldestination_list, destination_list_box,
                                        self.lblfixed_starting_mass_increment, self.txtfixed_starting_mass_increment,
                                        self.lblminimum_starting_mass_increment, self.txtminimum_starting_mass_increment,
                                        self.lblmaximum_starting_mass_increment, self.txtmaximum_starting_mass_increment,
                                        self.lblvariable_mass_increment, self.chkvariable_mass_increment,
                                        self.lblfixed_ending_mass_increment, self.txtfixed_ending_mass_increment,
                                        self.lblconstrain_initial_mass, self.chkconstrain_initial_mass,
                                        self.lblmaximum_initial_mass, self.txtmaximum_initial_mass,
                                        self.lblwait_time_bounds, wait_time_sizer,
                                        self.lblbounded_departure_date, self.chkbounded_departure_date,
                                        self.lbldeparture_date_bounds, departure_date_sizer,
                                        self.lbltimebounded, self.cmbtimebounded,
                                        self.lblflight_time_bounds, flight_time_sizer,
                                        self.lblarrival_date_bounds, arrival_date_sizer,
                                        self.lblinitial_impulse_bounds, initial_impulse_sizer,
                                        self.lblforce_free_point_direct_insertion_along_velocity_vector, self.chkforce_free_point_direct_insertion_along_velocity_vector,
                                        self.lbldeparture_type, self.cmbdeparture_type,
                                        self.lbldeparture_class, self.cmbdeparture_class,
                                        self.lbldeparture_ellipsoid_axes, departure_ellipsoid_sizer,
                                        self.lblzero_turn_flyby_distance, self.txtzero_turn_flyby_distance,
                                        self.lblforced_initial_coast, self.txtforced_initial_coast,
                                        self.lbloverride_flyby_altitude_bounds, self.chkoverride_flyby_altitude_bounds,
                                        self.lblflyby_altitude_bounds, flyby_altitude_box,
                                        self.lblescape_spiral_starting_radius, self.txtescape_spiral_starting_radius,
                                        self.lblescape_spiral_final_radius, self.txtescape_spiral_final_radius,
                                        self.lblarrival_type, self.cmbarrival_type,
                                        self.lblarrival_class, self.cmbarrival_class,
                                        self.lblimpact_momentum_enhancement_factor, self.txtimpact_momentum_enhancement_factor,
                                        self.lblephemeris_pegged_orbit_insertion_SMA, self.txtephemeris_pegged_orbit_insertion_SMA,
                                        self.lblephemeris_pegged_orbit_insertion_ECC, self.txtephemeris_pegged_orbit_insertion_ECC,
                                        self.lblarrival_ellipsoid_axes, arrival_ellipsoid_sizer,
                                        self.lblforced_terminal_coast, self.txtforced_terminal_coast,
                                        self.lblcapture_spiral_starting_radius, self.txtcapture_spiral_starting_radius,
                                        self.lblcapture_spiral_final_radius, self.txtcapture_spiral_final_radius,
                                        self.lblfinal_velocity, final_velocity_box,  
                                        self.lblFreePointArrival_print_target_spec, self.chkFreePointArrival_print_target_spec,
                                        self.lbljourney_end_deltav, self.txtjourney_end_deltav,
                                        self.lbljourney_end_TCM, self.txtjourney_end_TCM,
                                        self.lbloverride_duty_cycle, self.chkoverride_duty_cycle,
                                        self.lblduty_cycle, self.txtduty_cycle,
                                        self.lblsequence, sequence_box,
                                        self.lblperiapse_burns, self.chkperiapse_burns,
                                        self.lblperturbation_bodies, perturbation_bodies_box,
                                        self.lbloverride_PropagatorType, self.chkoverride_PropagatorType,
                                        self.lblpropagatorType, self.cmbpropagatorType,
                                        self.lbloverride_integration_step_size, self.chkoverride_integration_step_size,
                                        self.lblintegration_step_size, self.txtintegration_step_size,
                                        self.lblCoastPhaseMatchPointFraction, self.txtCoastPhaseMatchPointFraction,
                                        self.lblCoastPhaseForwardIntegrationStepLength, self.txtCoastPhaseForwardIntegrationStepLength,
                                        self.lblCoastPhaseBackwardIntegrationStepLength, self.txtCoastPhaseBackwardIntegrationStepLength,
                                        self.lblEnableDrag, self.chkEnableDrag,
                                        self.lblSpacecraftDragArea, self.txtSpacecraftDragArea,
                                        self.lblSpacecraftDragCoefficient, self.txtSpacecraftDragCoefficient,
                                        self.lblAtmosphericDensityModel, self.cmbAtmosphericDensityModel,
                                        self.lblAtmosphericDensityModelDataFile, self.atmosphericDensityModelDataFileBox])

        JourneyInformationStacker = wx.BoxSizer(wx.VERTICAL)
        JourneyInformationStacker.AddMany([JourneySelectionSizer, JourneyInformationGrid])


        #custom orbit elements
        self.DepartureElementsPanel = DepartureElementsPanel.DepartureElementsPanel(self, self.missionoptions)
        self.ArrivalElementsPanel = ArrivalElementsPanel.ArrivalElementsPanel(self, self.missionoptions)
        
        #staging
        self.StagingPanel = JourneyStagingPanel.JourneyStagingPanel(self, self.missionoptions)

        RightStacker = wx.BoxSizer(wx.VERTICAL)
        RightStacker.AddMany([self.DepartureElementsPanel, self.ArrivalElementsPanel, self.StagingPanel])

        self.mainbox = wx.BoxSizer(wx.HORIZONTAL)
        self.mainbox.AddMany([JourneyInformationStacker, RightStacker])
        self.SetSizer(self.mainbox)
        self.SetupScrolling()

        #bindings
        self.JourneySelectBox.Bind(wx.EVT_LISTBOX, self.ChangeJourneySelectBoxChoice)
        self.btnAddNewJourney.Bind(wx.EVT_BUTTON, self.ClickAddNewJourney)
        self.btnDeleteJourney.Bind(wx.EVT_BUTTON, self.ClickDeleteJourney)
        self.btnMoveJourneyUp.Bind(wx.EVT_BUTTON, self.ClickMoveJourneyUp)
        self.btnMoveJourneyDown.Bind(wx.EVT_BUTTON, self.ClickMoveJourneyDown)

        self.txtjourney_name.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_name)
        self.cmbPhaseType.Bind(wx.EVT_COMBOBOX, self.ChangePhaseType)
        self.chkoverride_num_steps.Bind(wx.EVT_CHECKBOX, self.Changeoverride_num_steps)
        self.txtnumber_of_steps.Bind(wx.EVT_KILL_FOCUS, self.Changenumber_of_steps)
        self.txtnum_interior_control_points.Bind(wx.EVT_KILL_FOCUS, self.Changenum_interior_control_points)
        self.txtimpulses_per_phase.Bind(wx.EVT_KILL_FOCUS, self.Changeimpulses_per_phase)
        self.cmbforce_unit_magnitude_control.Bind(wx.EVT_COMBOBOX, self.Changeforce_unit_magnitude_control)
        self.chkforce_fixed_inertial_control.Bind(wx.EVT_CHECKBOX, self.Changeforce_fixed_inertial_control)
        self.txtjourney_central_body.Bind(wx.EVT_KILL_FOCUS,self.Changejourney_central_body)
        self.txtdestination_list.Bind(wx.EVT_KILL_FOCUS,self.Changedestination_list)
        self.txtfixed_ending_mass_increment.Bind(wx.EVT_KILL_FOCUS,self.Changefixed_ending_mass_increment)
        self.txtfixed_starting_mass_increment.Bind(wx.EVT_KILL_FOCUS,self.Changefixed_starting_mass_increment)
        self.txtminimum_starting_mass_increment.Bind(wx.EVT_KILL_FOCUS,self.Changeminimum_starting_mass_increment)
        self.txtmaximum_starting_mass_increment.Bind(wx.EVT_KILL_FOCUS,self.Changemaximum_starting_mass_increment)
        self.chkvariable_mass_increment.Bind(wx.EVT_CHECKBOX,self.Changevariable_mass_increment)
        self.chkconstrain_initial_mass.Bind(wx.EVT_CHECKBOX,self.Changeconstrain_initial_mass)
        self.txtmaximum_initial_mass.Bind(wx.EVT_KILL_FOCUS,self.Changemaximum_initial_mass)
        self.txtwait_time_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.Changewait_time_bounds_lower)
        self.txtwait_time_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.Changewait_time_bounds_upper)
        self.cmbtimebounded.Bind(wx.EVT_COMBOBOX,self.Changetimebounded)
        self.txtflight_time_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.Changeflight_time_bounds_lower)
        self.txtflight_time_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.Changeflight_time_bounds_upper)
        self.txtarrival_date_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.Changearrival_date_bounds_lower)
        self.ArrivalDateLowerCalendar.Bind(wx.adv.EVT_CALENDAR_SEL_CHANGED, self.ChangeArrivalDateLowerCalendar)
        self.ArrivalDateUpperCalendar.Bind(wx.adv.EVT_CALENDAR_SEL_CHANGED, self.ChangeArrivalDateUpperCalendar)
        self.txtarrival_date_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.Changearrival_date_bounds_upper)
        self.chkbounded_departure_date.Bind(wx.EVT_CHECKBOX, self.Changebounded_departure_date)
        self.txtdeparture_date_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.Changedeparture_date_bounds_lower)
        self.departureDateLowerCalendar.Bind(wx.adv.EVT_CALENDAR_SEL_CHANGED, self.ChangedepartureDateLowerCalendar)
        self.departureDateUpperCalendar.Bind(wx.adv.EVT_CALENDAR_SEL_CHANGED, self.ChangedepartureDateUpperCalendar)
        self.txtdeparture_date_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.Changedeparture_date_bounds_upper)
        self.txtinitial_impulse_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.Changeinitial_impulse_bounds_lower)
        self.txtinitial_impulse_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.Changeinitial_impulse_bounds_upper)
        self.chkforce_free_point_direct_insertion_along_velocity_vector.Bind(wx.EVT_CHECKBOX, self.Changeforce_free_point_direct_insertion_along_velocity_vector)
        self.cmbdeparture_type.Bind(wx.EVT_COMBOBOX,self.Changedeparture_type)
        self.cmbdeparture_class.Bind(wx.EVT_COMBOBOX,self.Changedeparture_class)
        self.txtdeparture_ellipsoid_axes.Bind(wx.EVT_KILL_FOCUS, self.Changedeparture_ellipsoid_axes)
        self.btndeparture_ellipsoid_axes.Bind(wx.EVT_BUTTON, self.ClickButtondeparture_ellipsoid_axes)
        self.txtzero_turn_flyby_distance.Bind(wx.EVT_KILL_FOCUS, self.Changezero_turn_flyby_distance)
        self.txtforced_initial_coast.Bind(wx.EVT_KILL_FOCUS, self.Changeforced_initial_coast)

        self.chkoverride_flyby_altitude_bounds.Bind(wx.EVT_CHECKBOX,self.Changeoverride_flyby_altitude_bounds)
        self.txtflyby_altitude_boundsLower.Bind(wx.EVT_KILL_FOCUS,self.Changeflyby_altitude_boundsLower)
        self.txtflyby_altitude_boundsUpper.Bind(wx.EVT_KILL_FOCUS,self.Changeflyby_altitude_boundsUpper)

        self.txtescape_spiral_starting_radius.Bind(wx.EVT_KILL_FOCUS, self.Changeescape_spiral_starting_radius)
        self.txtescape_spiral_final_radius.Bind(wx.EVT_KILL_FOCUS, self.Changeescape_spiral_final_radius)
        self.cmbarrival_type.Bind(wx.EVT_COMBOBOX,self.Changearrival_type)
        self.cmbarrival_class.Bind(wx.EVT_COMBOBOX,self.Changearrival_class)
        self.txtarrival_ellipsoid_axes.Bind(wx.EVT_KILL_FOCUS, self.Changearrival_ellipsoid_axes)
        self.btnarrival_ellipsoid_axes.Bind(wx.EVT_BUTTON, self.ClickButtonarrival_ellipsoid_axes)
        self.txtephemeris_pegged_orbit_insertion_SMA.Bind(wx.EVT_KILL_FOCUS, self.Changeephemeris_pegged_orbit_insertion_SMA)
        self.txtephemeris_pegged_orbit_insertion_ECC.Bind(wx.EVT_KILL_FOCUS, self.Changeephemeris_pegged_orbit_insertion_ECC)
        self.txtcapture_spiral_starting_radius.Bind(wx.EVT_KILL_FOCUS, self.Changecapture_spiral_starting_radius)
        self.txtcapture_spiral_final_radius.Bind(wx.EVT_KILL_FOCUS, self.Changecapture_spiral_final_radius)
        self.txtforced_terminal_coast.Bind(wx.EVT_KILL_FOCUS, self.Changeforced_terminal_coast)
        self.txtfinal_velocity0.Bind(wx.EVT_KILL_FOCUS,self.Changefinal_velocity0)
        self.txtfinal_velocity1.Bind(wx.EVT_KILL_FOCUS,self.Changefinal_velocity1)
        self.txtfinal_velocity2.Bind(wx.EVT_KILL_FOCUS,self.Changefinal_velocity2)
        self.chkFreePointArrival_print_target_spec.Bind(wx.EVT_KILL_FOCUS, self.ChangeFreePointArrival_print_target_spec)

        self.txtimpact_momentum_enhancement_factor.Bind(wx.EVT_KILL_FOCUS, self.Changeimpact_momentum_enhancement_factor)

                
        self.chkoverride_duty_cycle.Bind(wx.EVT_CHECKBOX, self.Changeoverride_duty_cycle)
        self.txtduty_cycle.Bind(wx.EVT_KILL_FOCUS, self.Changeduty_cycle)

        self.txtsequence.Bind(wx.EVT_KILL_FOCUS,self.Changesequence)
        self.chkperiapse_burns.Bind(wx.EVT_CHECKBOX,self.Changeperiapse_burns)
        self.txtperturbation_bodies.Bind(wx.EVT_KILL_FOCUS,self.Changeperturbation_bodies)
        
        self.btndestination_list.Bind(wx.EVT_BUTTON,self.Clickdestination_list)
        self.btnjourney_central_body.Bind(wx.EVT_BUTTON,self.Clickjourney_central_body)
        self.btnsequence.Bind(wx.EVT_BUTTON,self.Clicksequence)
        self.btnperturbation_bodies.Bind(wx.EVT_BUTTON, self.Clickperturbation_bodies)
        self.chkfreeze_decision_variables.Bind(wx.EVT_CHECKBOX, self.Changefreeze_decision_variables)

        self.txtjourney_end_deltav.Bind(wx.EVT_KILL_FOCUS, self.Changejourney_end_deltav)
        self.txtjourney_end_TCM.Bind(wx.EVT_KILL_FOCUS, self.Changejourney_end_TCM)

        self.chkoverride_PropagatorType.Bind(wx.EVT_CHECKBOX, self.Changeoverride_PropagatorType)
        self.cmbpropagatorType.Bind(wx.EVT_COMBOBOX, self.ChangepropagatorType)

        self.chkoverride_integration_step_size.Bind(wx.EVT_CHECKBOX, self.Changeoverride_integration_step_size)
        self.txtintegration_step_size.Bind(wx.EVT_KILL_FOCUS, self.Changeintegration_step_size)

        
        self.txtCoastPhaseMatchPointFraction.Bind(wx.EVT_KILL_FOCUS, self.ChangeCoastPhaseMatchPointFraction)
        self.txtCoastPhaseForwardIntegrationStepLength.Bind(wx.EVT_KILL_FOCUS, self.ChangeCoastPhaseForwardIntegrationStepLength)
        self.txtCoastPhaseBackwardIntegrationStepLength.Bind(wx.EVT_KILL_FOCUS, self.ChangeCoastPhaseBackwardIntegrationStepLength)

        ###############
        # drag bindings
        ###############

        self.chkEnableDrag.Bind(wx.EVT_CHECKBOX, self.ChangeEnableDrag)
        self.txtSpacecraftDragArea.Bind(wx.EVT_KILL_FOCUS, self.ChangeSpacecraftDragArea)
        self.txtSpacecraftDragCoefficient.Bind(wx.EVT_KILL_FOCUS, self.ChangeSpacecraftDragCoefficient)
        self.cmbAtmosphericDensityModel.Bind(wx.EVT_COMBOBOX, self.ChangeAtmosphericDensityModel)
        self.btnAtmosphericDensityModelDataFile.Bind(wx.EVT_BUTTON, self.ClickAtmosphericDensityModelDataFile)
        self.txtAtmosphericDensityModelDataFile.Bind(wx.EVT_KILL_FOCUS,self.ChangeAtmosphericDensityModelDataFile)

    def update(self):
        self.StagingPanel.update()

        self.Journeylist = []
        for j in range(0, self.missionoptions.number_of_journeys):
            self.Journeylist.append(self.missionoptions.Journeys[j].journey_name)

        self.JourneySelectBox.SetItems(self.Journeylist)
        self.JourneySelectBox.SetSelection(self.missionoptions.ActiveJourney)

        self.txtjourney_name.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_name))        
        self.cmbPhaseType.SetSelection(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].phase_type)
        self.chkoverride_num_steps.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_num_steps)
        self.txtnumber_of_steps.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].number_of_steps))
        self.txtnum_interior_control_points.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].num_interior_control_points))
        self.txtimpulses_per_phase.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].impulses_per_phase))
        self.cmbforce_unit_magnitude_control.SetSelection(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].force_unit_magnitude_control)
        self.chkforce_fixed_inertial_control.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].force_fixed_inertial_control)
        self.txtjourney_central_body.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body))
        self.txtdestination_list.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list))
        self.txtfixed_starting_mass_increment.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].fixed_starting_mass_increment))
        self.txtminimum_starting_mass_increment.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].minimum_starting_mass_increment))
        self.txtmaximum_starting_mass_increment.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].maximum_starting_mass_increment))
        self.chkvariable_mass_increment.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].variable_mass_increment)    
        self.txtfixed_ending_mass_increment.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].fixed_ending_mass_increment))        
        self.chkconstrain_initial_mass.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].constrain_initial_mass)
        self.txtmaximum_initial_mass.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].maximum_initial_mass))
        self.chkconstrain_initial_mass.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].constrain_initial_mass)
        self.txtmaximum_initial_mass.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].maximum_initial_mass))
        self.txtwait_time_bounds_lower.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].wait_time_bounds[0]))
        self.txtwait_time_bounds_upper.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].wait_time_bounds[1]))
        self.cmbtimebounded.SetSelection(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].timebounded)
        self.txtflight_time_bounds_lower.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].flight_time_bounds[0]))
        self.txtflight_time_bounds_upper.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].flight_time_bounds[1]))
        self.txtarrival_date_bounds_lower.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_date_bounds[0]))
        self.chkbounded_departure_date.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].bounded_departure_date)
        date = wx.DateTime.FromJDN(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_date_bounds[0] + 2400000.5)
        self.ArrivalDateLowerCalendar.SetDate(date.MakeUTC())
        self.txtarrival_date_bounds_upper.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_date_bounds[1]))
        date = wx.DateTime.FromJDN(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_date_bounds[1] + 2400000.5)
        self.ArrivalDateUpperCalendar.SetDate(date.MakeUTC())
        self.txtdeparture_date_bounds_lower.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_date_bounds[0]))
        date = wx.DateTime.FromJDN(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_date_bounds[0] + 2400000.5)
        self.departureDateLowerCalendar.SetDate(date.MakeUTC())
        self.txtdeparture_date_bounds_upper.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_date_bounds[1]))
        date = wx.DateTime.FromJDN(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_date_bounds[1] + 2400000.5)
        self.departureDateUpperCalendar.SetDate(date.MakeUTC())
        self.txtinitial_impulse_bounds_lower.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].initial_impulse_bounds[0]))
        self.txtinitial_impulse_bounds_upper.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].initial_impulse_bounds[1]))
        self.chkforce_free_point_direct_insertion_along_velocity_vector.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].force_free_point_direct_insertion_along_velocity_vector)
        self.cmbdeparture_type.SetSelection(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_type)
        self.cmbdeparture_class.SetSelection(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_class)
        self.txtdeparture_ellipsoid_axes.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_ellipsoid_axes))
        self.txtzero_turn_flyby_distance.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].zero_turn_flyby_distance))
        self.txtforced_initial_coast.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].forced_initial_coast))
        self.chkoverride_flyby_altitude_bounds.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_flyby_altitude_bounds)
        self.txtflyby_altitude_boundsLower.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].flyby_altitude_bounds[0]))
        self.txtflyby_altitude_boundsUpper.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].flyby_altitude_bounds[1]))
        self.txtescape_spiral_starting_radius.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].escape_spiral_starting_radius))
        self.txtescape_spiral_final_radius.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].escape_spiral_final_radius))
        self.cmbarrival_type.SetSelection(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_type)
        self.cmbarrival_class.SetSelection(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_class)
        self.txtephemeris_pegged_orbit_insertion_SMA.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].ephemeris_pegged_orbit_insertion_SMA))
        self.txtephemeris_pegged_orbit_insertion_ECC.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].ephemeris_pegged_orbit_insertion_ECC))
        self.txtarrival_ellipsoid_axes.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_ellipsoid_axes))
        self.txtcapture_spiral_starting_radius.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].capture_spiral_starting_radius))
        self.txtcapture_spiral_final_radius.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].capture_spiral_final_radius))
        self.txtforced_terminal_coast.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].forced_terminal_coast))
        self.txtfinal_velocity0.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].final_velocity[0]))
        self.txtfinal_velocity1.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].final_velocity[1]))
        self.txtfinal_velocity2.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].final_velocity[2]))
        self.chkFreePointArrival_print_target_spec.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].FreePointArrival_print_target_spec)
 
        self.txtimpact_momentum_enhancement_factor.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].impact_momentum_enhancement_factor))
 
        
        self.txtsequence.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].sequence))
        self.chkperiapse_burns.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].enable_periapse_burns)
        self.txtperturbation_bodies.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].perturbation_bodies))
        self.chkfreeze_decision_variables.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].freeze_decision_variables)
        
        self.txtjourney_end_deltav.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_end_deltav))
        self.txtjourney_end_TCM.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_end_TCM))
        
        self.chkoverride_duty_cycle.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_duty_cycle)
        self.txtduty_cycle.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].duty_cycle))
        self.chkoverride_PropagatorType.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_PropagatorType)
        self.cmbpropagatorType.SetSelection(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].propagatorType)
        self.chkoverride_integration_step_size.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_integration_step_size)
        self.txtintegration_step_size.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].integration_step_size))      
        self.txtCoastPhaseMatchPointFraction.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].CoastPhaseMatchPointFraction))
        self.txtCoastPhaseForwardIntegrationStepLength.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].CoastPhaseForwardIntegrationStepLength))
        self.txtCoastPhaseBackwardIntegrationStepLength.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].CoastPhaseBackwardIntegrationStepLength))

        ##############
        # drag updates
        ##############
        self.chkEnableDrag.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].perturb_drag)
        self.txtSpacecraftDragArea.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].spacecraft_drag_area))
        self.txtSpacecraftDragCoefficient.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].coefficient_of_drag))
        self.cmbAtmosphericDensityModel.SetSelection(
            self.atmosphericDensityModelChoices.index(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].AtmosphericDensityModelKey))
        #self.cmbAtmosphericDensityModel.SetSelection(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].AtmosphericDensityModelKey)
        self.txtAtmosphericDensityModelDataFile.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].AtmosphericDensityModelDataFile))

        #if there is only one journey in the list then disable delete, up, and down
        if self.missionoptions.number_of_journeys == 1:
            self.btnDeleteJourney.Disable()
            self.btnMoveJourneyUp.Disable()
            self.btnMoveJourneyDown.Disable()
        else:
            self.btnDeleteJourney.Enable()

            #if the first journey in the list is active then you cannot move up
            if self.missionoptions.ActiveJourney == 0:
                self.btnMoveJourneyUp.Disable()
            else:
                self.btnMoveJourneyUp.Enable()

            #if the last journey in the list is active then you cannot move down
            if self.missionoptions.ActiveJourney == self.missionoptions.number_of_journeys - 1:
                self.btnMoveJourneyDown.Disable()
            else:
                self.btnMoveJourneyDown.Enable()

        #only show phase type if mission type is variable
        if self.missionoptions.mission_type == 9:
            self.lblPhaseType.Show(True)
            self.cmbPhaseType.Show(True)
        else:
            self.lblPhaseType.Show(False)
            self.cmbPhaseType.Show(False)

        #only show number of steps if override is on
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_num_steps:
            self.lblnumber_of_steps.Show(True)
            self.txtnumber_of_steps.Show(True)
        else:
            self.lblnumber_of_steps.Show(False)
            self.txtnumber_of_steps.Show(False)

        #do we want to show the number of interior control points?
        if self.missionoptions.mission_type in [4, 5] or (self.missionoptions.mission_type == 9 and self.missionoptions.Journeys[self.missionoptions.ActiveJourney].phase_type in [4, 5]):
            self.lblnum_interior_control_points.Show(True)
            self.txtnum_interior_control_points.Show(True)
        else:
            self.lblnum_interior_control_points.Show(False)
            self.txtnum_interior_control_points.Show(False)

        #do we want to show the option to force unit magnitude control and forced fixed inertial control?
        if self.missionoptions.mission_type in [0, 1, 2, 3, 4, 5] or (self.missionoptions.mission_type == 9 and self.missionoptions.Journeys[self.missionoptions.ActiveJourney].phase_type in [0, 1, 2, 3, 4, 5]):
            self.lblforce_unit_magnitude_control.Show(True)
            self.cmbforce_unit_magnitude_control.Show(True)
            self.lblforce_fixed_inertial_control.Show(True)
            self.chkforce_fixed_inertial_control.Show(True)
        else:
            self.lblforce_unit_magnitude_control.Show(False)
            self.cmbforce_unit_magnitude_control.Show(False)
            self.lblforce_fixed_inertial_control.Show(False)
            self.chkforce_fixed_inertial_control.Show(False)

        #do we want to show number of impulses?
        if self.missionoptions.mission_type == 6 or (self.missionoptions.mission_type == 9 and self.missionoptions.Journeys[self.missionoptions.ActiveJourney].phase_type == 6):
            self.lblimpulses_per_phase.Show(True)
            self.txtimpulses_per_phase.Show(True)
        else:
            self.lblimpulses_per_phase.Show(False)
            self.txtimpulses_per_phase.Show(False)

        #hide or show flight time and arrival date bounds
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].timebounded == 0:
            self.lblflight_time_bounds.Show(False)
            self.txtflight_time_bounds_lower.Show(False)
            self.txtflight_time_bounds_upper.Show(False)
            self.lblarrival_date_bounds.Show(False)
            self.txtarrival_date_bounds_lower.Show(False)
            self.txtarrival_date_bounds_upper.Show(False)
            self.ArrivalDateLowerCalendar.Show(False)
            self.ArrivalDateUpperCalendar.Show(False)
        elif self.missionoptions.Journeys[self.missionoptions.ActiveJourney].timebounded == 1 or self.missionoptions.Journeys[self.missionoptions.ActiveJourney].timebounded == 3:
            self.lblflight_time_bounds.Show(True)
            self.txtflight_time_bounds_lower.Show(True)
            self.txtflight_time_bounds_upper.Show(True)
            self.lblarrival_date_bounds.Show(False)
            self.txtarrival_date_bounds_lower.Show(False)
            self.txtarrival_date_bounds_upper.Show(False)
            self.ArrivalDateLowerCalendar.Show(False)
            self.ArrivalDateUpperCalendar.Show(False)
        elif self.missionoptions.Journeys[self.missionoptions.ActiveJourney].timebounded == 2:
            self.lblflight_time_bounds.Show(False)
            self.txtflight_time_bounds_lower.Show(False)
            self.txtflight_time_bounds_upper.Show(False)
            self.lblarrival_date_bounds.Show(True)
            self.txtarrival_date_bounds_lower.Show(True)
            self.txtarrival_date_bounds_upper.Show(True)
            self.ArrivalDateLowerCalendar.Show(True)
            self.ArrivalDateUpperCalendar.Show(True)

        #hide or show departure date bounds
        if self.missionoptions.ActiveJourney == 0 \
                or self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_type == 3 \
                or self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_type == 4 \
                or self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_type == 6:
            self.lblbounded_departure_date.Show(False)
            self.chkbounded_departure_date.Show(False)
            self.lbldeparture_date_bounds.Show(False)
            self.txtdeparture_date_bounds_lower.Show(False)
            self.txtdeparture_date_bounds_upper.Show(False)
            self.departureDateLowerCalendar.Show(False)
            self.departureDateUpperCalendar.Show(False)
        elif self.missionoptions.Journeys[self.missionoptions.ActiveJourney].bounded_departure_date == 0:
            self.lblbounded_departure_date.Show(True)
            self.chkbounded_departure_date.Show(True)
            self.lbldeparture_date_bounds.Show(False)
            self.txtdeparture_date_bounds_lower.Show(False)
            self.txtdeparture_date_bounds_upper.Show(False)
            self.departureDateLowerCalendar.Show(False)
            self.departureDateUpperCalendar.Show(False)
        else:
            self.lblbounded_departure_date.Show(True)
            self.chkbounded_departure_date.Show(True)
            self.lbldeparture_date_bounds.Show(True)
            self.txtdeparture_date_bounds_lower.Show(True)
            self.txtdeparture_date_bounds_upper.Show(True)
            self.departureDateLowerCalendar.Show(True)
            self.departureDateUpperCalendar.Show(True)

        #enable or disable the orbit elements selection boxes as appropriate
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_class == 1:
            #enable departure orbit elements box
            self.DepartureElementsPanel.Show(True)
            self.DepartureElementsPanel.update()
        else:
            self.DepartureElementsPanel.Show(False)
        
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_class == 1:
            #enable arrival orbit elements box
            self.ArrivalElementsPanel.Show(True)
            self.lblFreePointArrival_print_target_spec.Show(True)
            self.chkFreePointArrival_print_target_spec.Show(True)
            self.ArrivalElementsPanel.update()
        else:
            self.ArrivalElementsPanel.Show(False)
            self.lblFreePointArrival_print_target_spec.Show(False)
            self.chkFreePointArrival_print_target_spec.Show(False)

        #impact momentum enhancement factor
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_class == 0 and self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_type == 7:
            self.lblimpact_momentum_enhancement_factor.Show(True)
            self.txtimpact_momentum_enhancement_factor.Show(True)
        else:
            self.lblimpact_momentum_enhancement_factor.Show(False)
            self.txtimpact_momentum_enhancement_factor.Show(False)

        #ephemeris pegged orbit insertion
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_class == 0 and self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_type == 0:
            self.txtephemeris_pegged_orbit_insertion_SMA.Show(True)
            self.lblephemeris_pegged_orbit_insertion_SMA.Show(True)
            self.txtephemeris_pegged_orbit_insertion_ECC.Show(True)
            self.lblephemeris_pegged_orbit_insertion_ECC.Show(True)
        else:
            self.txtephemeris_pegged_orbit_insertion_SMA.Show(False)
            self.lblephemeris_pegged_orbit_insertion_SMA.Show(False)
            self.txtephemeris_pegged_orbit_insertion_ECC.Show(False)
            self.lblephemeris_pegged_orbit_insertion_ECC.Show(False)

        #only show the flyby altitude override options if we are a flyby or a zero-turn flyby
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_type in [3, 4, 6]:
            self.lbloverride_flyby_altitude_bounds.Show(True)
            self.chkoverride_flyby_altitude_bounds.Show(True)
            if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_flyby_altitude_bounds:                    
                self.lblflyby_altitude_bounds.Show(True)
                self.txtflyby_altitude_boundsLower.Show(True)
                self.txtflyby_altitude_boundsUpper.Show(True)
            else:                    
                self.lblflyby_altitude_bounds.Show(False)
                self.txtflyby_altitude_boundsLower.Show(False)
                self.txtflyby_altitude_boundsUpper.Show(False)
        else:
            self.lbloverride_flyby_altitude_bounds.Show(False)
            self.chkoverride_flyby_altitude_bounds.Show(False)
            self.lblflyby_altitude_bounds.Show(False)
            self.txtflyby_altitude_boundsLower.Show(False)
            self.txtflyby_altitude_boundsUpper.Show(False)

        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_type in [3, 4, 6]:
            self.lblwait_time_bounds.Show(False)
            self.txtwait_time_bounds_lower.Show(False)
            self.txtwait_time_bounds_upper.Show(False)
            self.lblinitial_impulse_bounds.Show(False)
            self.txtinitial_impulse_bounds_lower.Show(False)
            self.txtinitial_impulse_bounds_upper.Show(False)
        else:
            self.lblwait_time_bounds.Show(True)
            self.txtwait_time_bounds_lower.Show(True)
            self.txtwait_time_bounds_upper.Show(True)
            self.lblinitial_impulse_bounds.Show(True)
            self.txtinitial_impulse_bounds_lower.Show(True)
            self.txtinitial_impulse_bounds_upper.Show(True)

        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_type in [0] \
            and self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_class == 1:

            self.lblforce_free_point_direct_insertion_along_velocity_vector.Show(True)
            self.chkforce_free_point_direct_insertion_along_velocity_vector.Show(True)
        else:
            self.lblforce_free_point_direct_insertion_along_velocity_vector.Show(False)
            self.chkforce_free_point_direct_insertion_along_velocity_vector.Show(False)


        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_type in [4, 5]:
            self.lblfinal_velocity.Show(True)
            self.lblfinal_velocity.SetLabel("Journey final velocity vector")
            self.txtfinal_velocity0.Show(True)
            self.txtfinal_velocity1.Show(True)
            self.txtfinal_velocity2.Show(True)
        elif self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_type in [1]:
            self.lblfinal_velocity.Show(True)
            self.lblfinal_velocity.SetLabel("Journey final impulse bounds (km/s)")
            self.txtfinal_velocity0.Show(True)
            self.txtfinal_velocity1.Show(True)
            self.txtfinal_velocity2.Show(False)
        elif self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_type in [2, 7]:
            self.lblfinal_velocity.Show(True)
            self.lblfinal_velocity.SetLabel("Journey final velocity bounds (km/s)")
            self.txtfinal_velocity0.Show(True)
            self.txtfinal_velocity1.Show(True)
            self.txtfinal_velocity2.Show(False)
        else:
            self.lblfinal_velocity.Show(False)
            self.txtfinal_velocity0.Show(False)
            self.txtfinal_velocity1.Show(False)
            self.txtfinal_velocity2.Show(False)

 
        #only show the pre-intercept coast  flag if this is a bounded v-infinity intercept or orbit insertion 
        #same rule applies for pre-intercept coast 
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_type in [0, 1, 2]:  
            self.lblforced_terminal_coast.Show(True) 
            self.txtforced_terminal_coast.Show(True) 
        else:
            self.lblforced_terminal_coast.Show(False) 
            self.txtforced_terminal_coast.Show(False)
 

        #options for an escape spiral
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_type == 5:
            self.lblinitial_impulse_bounds.Show(False)
            self.txtinitial_impulse_bounds_lower.Show(False)
            self.txtinitial_impulse_bounds_upper.Show(False)
            self.lblescape_spiral_starting_radius.Show(True)
            self.txtescape_spiral_starting_radius.Show(True)
            self.lblescape_spiral_final_radius.Show(True)
            self.txtescape_spiral_final_radius.Show(True)
        #free direct departure
        elif self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_type == 2:
            self.lblinitial_impulse_bounds.Show(False)
            self.txtinitial_impulse_bounds_lower.Show(False)
            self.txtinitial_impulse_bounds_upper.Show(False)
            self.lblescape_spiral_starting_radius.Show(False)
            self.txtescape_spiral_starting_radius.Show(False)
            self.lblescape_spiral_final_radius.Show(False)
            self.txtescape_spiral_final_radius.Show(False)
        else:
            #hide the initial v-infinity options
            self.lblescape_spiral_starting_radius.Show(False)
            self.txtescape_spiral_starting_radius.Show(False)
            self.lblescape_spiral_final_radius.Show(False)
            self.txtescape_spiral_final_radius.Show(False)

        #options for a capture spiral
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_type == 6:
            self.lblcapture_spiral_starting_radius.Show(True)
            self.txtcapture_spiral_starting_radius.Show(True)
            self.lblcapture_spiral_final_radius.Show(True)
            self.txtcapture_spiral_final_radius.Show(True)
        else:
            self.lblcapture_spiral_starting_radius.Show(False)
            self.txtcapture_spiral_starting_radius.Show(False)
            self.lblcapture_spiral_final_radius.Show(False)
            self.txtcapture_spiral_final_radius.Show(False)


        #options for third body perturbations
        if self.missionoptions.perturb_thirdbody == 1:
            self.lblperturbation_bodies.Show(True)
            self.txtperturbation_bodies.Show(True)
            self.btnperturbation_bodies.Show(True)
        else:
            self.lblperturbation_bodies.Show(False)
            self.txtperturbation_bodies.Show(False)
            self.btnperturbation_bodies.Show(False)

        if (self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_duty_cycle):
            self.lblduty_cycle.Show(True)
            self.txtduty_cycle.Show(True)
        else:
            self.lblduty_cycle.Show(False)
            self.txtduty_cycle.Show(False)

        #show journey-end maneuver components if it is needed
        #should point to journey-end propulsion system here
        #if (self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_end_deltav > 0.0):
         
        #else:
          
        #journey-end TCM, only needed for intercepts, orbit insertions, and chemical rendezvous
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_type in [0, 1, 2]:
            self.lbljourney_end_TCM.Show(True)
            self.txtjourney_end_TCM.Show(True)
        else:
            self.lbljourney_end_TCM.Show(False)
            self.txtjourney_end_TCM.Show(False)

        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].variable_mass_increment:
            self.txtmaximum_starting_mass_increment.Show(True)
            self.txtminimum_starting_mass_increment.Show(True)
            self.txtfixed_starting_mass_increment.Show(False)
            self.lblmaximum_starting_mass_increment.Show(True)
            self.lblminimum_starting_mass_increment.Show(True)
            self.lblfixed_starting_mass_increment.Show(False)
        else:
            self.txtmaximum_starting_mass_increment.Show(False)
            self.txtminimum_starting_mass_increment.Show(False)
            self.txtfixed_starting_mass_increment.Show(True)
            self.lblmaximum_starting_mass_increment.Show(False)
            self.lblminimum_starting_mass_increment.Show(False)
            self.lblfixed_starting_mass_increment.Show(True)
        
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].constrain_initial_mass:
            self.txtmaximum_initial_mass.Show(True)
            self.lblmaximum_initial_mass.Show(True)
        else:
            self.txtmaximum_initial_mass.Show(False)
            self.lblmaximum_initial_mass.Show(False)

        if self.missionoptions.mission_type in [1, 3, 5, 6, 7] or (self.missionoptions.mission_type == 9 and self.missionoptions.Journeys[self.missionoptions.ActiveJourney].phase_type in [1, 3, 5, 6, 7]):
            self.lbloverride_PropagatorType.Show(True)
            self.chkoverride_PropagatorType.Show(True)

            if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_PropagatorType:
                self.lblpropagatorType.Show(True)
                self.cmbpropagatorType.Show(True)
            else:
                self.lblpropagatorType.Show(False)
                self.cmbpropagatorType.Show(False)
        else:
            self.lbloverride_PropagatorType.Show(False)
            self.chkoverride_PropagatorType.Show(False)
            self.lblpropagatorType.Show(False)
            self.cmbpropagatorType.Show(False)

        if self.missionoptions.mission_type in [7] or (self.missionoptions.mission_type == 9 and self.missionoptions.Journeys[self.missionoptions.ActiveJourney].phase_type in [7]):
            self.lblCoastPhaseMatchPointFraction.Show(True)
            self.txtCoastPhaseMatchPointFraction.Show(True)
            self.lblCoastPhaseForwardIntegrationStepLength.Show(True)
            self.txtCoastPhaseForwardIntegrationStepLength.Show(True)
            self.lblCoastPhaseBackwardIntegrationStepLength.Show(True)
            self.txtCoastPhaseBackwardIntegrationStepLength.Show(True)
            self.lbloverride_integration_step_size.Show(False)
            self.chkoverride_integration_step_size.Show(False)
            self.lblintegration_step_size.Show(False)
            self.txtintegration_step_size.Show(False)
        elif self.missionoptions.mission_type in [1, 3, 5, 8] \
            or (self.missionoptions.mission_type in [6] and (self.missionoptions.propagatorType == 1 or self.missionoptions.Journeys[self.missionoptions.ActiveJourney].propagatorType == 1)) \
            or (self.missionoptions.mission_type == 9 and (self.missionoptions.Journeys[self.missionoptions.ActiveJourney].phase_type in [1, 3, 5, 8] or (self.missionoptions.Journeys[self.missionoptions.ActiveJourney].phase_type in [6] and (self.missionoptions.propagatorType == 1 or self.missionoptions.Journeys[self.missionoptions.ActiveJourney].propagatorType == 1)))):
            self.lbloverride_integration_step_size.Show(True)
            self.chkoverride_integration_step_size.Show(True)
            if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_integration_step_size:
                self.lblintegration_step_size.Show(True)
                self.txtintegration_step_size.Show(True)

                if self.missionoptions.mission_type == 8 or (self.missionoptions.mission_type == 9 and self.missionoptions.Journeys[self.missionoptions.ActiveJourney].phase_type == 8):
                    self.lblintegration_step_size.SetLabel('Integration step size (degrees)')
                else:
                    self.lblintegration_step_size.SetLabel('Integration step size (seconds)')
            else:
                self.lblintegration_step_size.Show(False)
                self.txtintegration_step_size.Show(False)
            self.lblCoastPhaseMatchPointFraction.Show(False)
            self.txtCoastPhaseMatchPointFraction.Show(False)
            self.lblCoastPhaseForwardIntegrationStepLength.Show(False)
            self.txtCoastPhaseForwardIntegrationStepLength.Show(False)
            self.lblCoastPhaseBackwardIntegrationStepLength.Show(False)
            self.txtCoastPhaseBackwardIntegrationStepLength.Show(False)

        else:
            self.lbloverride_integration_step_size.Show(False)
            self.chkoverride_integration_step_size.Show(False)
            self.lblintegration_step_size.Show(False)
            self.txtintegration_step_size.Show(False)
            self.lblCoastPhaseMatchPointFraction.Show(False)
            self.txtCoastPhaseMatchPointFraction.Show(False)
            self.lblCoastPhaseForwardIntegrationStepLength.Show(False)
            self.txtCoastPhaseForwardIntegrationStepLength.Show(False)
            self.lblCoastPhaseBackwardIntegrationStepLength.Show(False)
            self.txtCoastPhaseBackwardIntegrationStepLength.Show(False)

        #ellipsoid controls
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_class == 2: #ephemeris referenced
            self.lbldeparture_ellipsoid_axes.Show(True)
            self.txtdeparture_ellipsoid_axes.Show(True)
            self.btndeparture_ellipsoid_axes.Show(True)
        else:
            self.lbldeparture_ellipsoid_axes.Show(False)
            self.txtdeparture_ellipsoid_axes.Show(False)
            self.btndeparture_ellipsoid_axes.Show(False)

        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_class == 2: #ephemeris referenced
            self.lblarrival_ellipsoid_axes.Show(True)
            self.txtarrival_ellipsoid_axes.Show(True)
            self.btnarrival_ellipsoid_axes.Show(True)
        else:
            self.lblarrival_ellipsoid_axes.Show(False)
            self.txtarrival_ellipsoid_axes.Show(False)
            self.btnarrival_ellipsoid_axes.Show(False)

        #zero-turn flyby distance
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_class == 0 \
            and self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_type == 6: #ephemeris-pegged zero-turn flyby
            self.lblzero_turn_flyby_distance.Show(True)
            self.txtzero_turn_flyby_distance.Show(True)
        else:
            self.lblzero_turn_flyby_distance.Show(False)
            self.txtzero_turn_flyby_distance.Show(False)

        # show drag options?
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].perturb_drag == 1:
            self.lblSpacecraftDragArea.Show(True)
            self.txtSpacecraftDragArea.Show(True)
            self.lblSpacecraftDragCoefficient.Show(True)
            self.txtSpacecraftDragCoefficient.Show(True)
            self.lblAtmosphericDensityModel.Show(True)
            self.cmbAtmosphericDensityModel.Show(True)
            self.lblAtmosphericDensityModelDataFile.Show(True)
            self.atmosphericDensityModelDataFileBox.Show(True)
            self.txtAtmosphericDensityModelDataFile.Show(True)
            self.btnAtmosphericDensityModelDataFile.Show(True)
        else:
            self.lblSpacecraftDragArea.Show(False)
            self.txtSpacecraftDragArea.Show(False)
            self.lblSpacecraftDragCoefficient.Show(False)
            self.txtSpacecraftDragCoefficient.Show(False)
            self.lblAtmosphericDensityModel.Show(False)
            self.cmbAtmosphericDensityModel.Show(False)
            self.lblAtmosphericDensityModelDataFile.Show(False)
            self.atmosphericDensityModelDataFileBox.Show(False)
            self.txtAtmosphericDensityModelDataFile.Show(False)
            self.btnAtmosphericDensityModelDataFile.Show(False)
        
        self.Layout()
        if platform.system() == 'Windows':
           self.SetupScrolling(scrollToTop=False)

    #event handlers for journey options
    def ChangeJourneySelectBoxChoice(self, e):
        #determine which journey is "active"
        self.missionoptions.ActiveJourney = self.JourneySelectBox.GetSelection()

        self.parent.update()

    
    def ClickAddNewJourney(self, e):
        #add
        temp_JourneyOptions = JO.JourneyOptions()
        temp_JourneyOptions.sequence = []
        
        if self.missionoptions.mission_type == 9:
            temp_JourneyOptions.phase_type = 7
        else:
            temp_JourneyOptions.phase_type = self.missionoptions.mission_type
        temp_JourneyOptions.impulses_per_phase = 1
        temp_JourneyOptions.journey_enable_periapse_burns = 0
        self.missionoptions.Journeys.append(temp_JourneyOptions)
        self.missionoptions.number_of_journeys += 1
        self.JourneySelectBox.SetSelection(-1)
        self.parent.update()

        #determine which journey is "active"
        self.missionoptions.ActiveJourney = self.JourneySelectBox.GetSelection()

    def ClickDeleteJourney(self, e):
        #determine which journey is "active"
        self.missionoptions.ActiveJourney = self.JourneySelectBox.GetSelection()

        #delete
        self.missionoptions.Journeys.pop(self.missionoptions.ActiveJourney)
        self.missionoptions.number_of_journeys -= 1

        if self.missionoptions.ActiveJourney > self.missionoptions.number_of_journeys - 1:
            self.missionoptions.ActiveJourney -= 1

        self.JourneySelectBox.SetSelection(0)
        self.parent.update()

    def ClickMoveJourneyUp(self, e):
        #determine which journey is "active"
        self.missionoptions.ActiveJourney = self.JourneySelectBox.GetSelection()

        #move up
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney], self.missionoptions.Journeys[self.missionoptions.ActiveJourney-1] = self.missionoptions.Journeys[self.missionoptions.ActiveJourney-1], self.missionoptions.Journeys[self.missionoptions.ActiveJourney]
        self.missionoptions.ActiveJourney -= 1
        self.JourneySelectBox.SetSelection(self.missionoptions.ActiveJourney)
        self.parent.update()

    def ClickMoveJourneyDown(self, e):
        #determine which journey is "active"
        self.missionoptions.ActiveJourney = self.JourneySelectBox.GetSelection()

        #move down
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney], self.missionoptions.Journeys[self.missionoptions.ActiveJourney+1] = self.missionoptions.Journeys[self.missionoptions.ActiveJourney+1], self.missionoptions.Journeys[self.missionoptions.ActiveJourney]
        self.missionoptions.ActiveJourney += 1
        self.JourneySelectBox.SetSelection(self.missionoptions.ActiveJourney)
        self.parent.update()

    def Changejourney_name(self, e):
         e.Skip()
         namestring = self.txtjourney_name.GetValue().replace(' ', '_')
         self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_name = namestring
         self.update()        
        
    def ChangePhaseType(self, e):
        e.Skip()

        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].phase_type = self.cmbPhaseType.GetSelection()
                
        self.missionoptions.DisassembleMasterDecisionVector()
        self.missionoptions.ConvertDecisionVector()
        self.missionoptions.AssembleMasterDecisionVector()

        self.update()

    def Changeoverride_num_steps(self, e):
         e.Skip()
         self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_num_steps = int(self.chkoverride_num_steps.GetValue())
         self.update()

    def Changenumber_of_steps(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].number_of_steps = int(eval(self.txtnumber_of_steps.GetValue()))
        self.update()        

    def Changenum_interior_control_points(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].num_interior_control_points = int(self.txtnum_interior_control_points.GetValue())
        self.parent.update()
     
    def Changeimpulses_per_phase(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].impulses_per_phase = int(self.txtimpulses_per_phase.GetValue())
        self.parent.update()    

    def Changeforce_unit_magnitude_control(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].force_unit_magnitude_control = self.cmbforce_unit_magnitude_control.GetSelection()
        self.parent.update()  
        
    def Changeforce_fixed_inertial_control(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].force_fixed_inertial_control = int(self.chkforce_fixed_inertial_control.GetValue())
        self.parent.update()

    def Changejourney_central_body(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body = self.txtjourney_central_body.GetValue()

    def Changedestination_list(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list = eval(self.txtdestination_list.GetValue())
        self.parent.update()

    def Changefixed_ending_mass_increment(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].fixed_ending_mass_increment = eval(self.txtfixed_ending_mass_increment.GetValue())
        self.update()

    def Changefixed_starting_mass_increment(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].fixed_starting_mass_increment = eval(self.txtfixed_starting_mass_increment.GetValue())
        self.update()

    def Changeminimum_starting_mass_increment(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].minimum_starting_mass_increment = eval(self.txtminimum_starting_mass_increment.GetValue())
        self.update()
        
    def Changemaximum_starting_mass_increment(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].maximum_starting_mass_increment = eval(self.txtmaximum_starting_mass_increment.GetValue())
        self.update()

    def Changevariable_mass_increment(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].variable_mass_increment = int(self.chkvariable_mass_increment.GetValue())
        self.update()

    def Changeconstrain_initial_mass(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].constrain_initial_mass = int(self.chkconstrain_initial_mass.GetValue())
        self.update()
        
    def Changemaximum_initial_mass(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].maximum_initial_mass = eval(self.txtmaximum_initial_mass.GetValue())
        self.update()

    def Changewait_time_bounds_lower(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].wait_time_bounds[0] = eval(self.txtwait_time_bounds_lower.GetValue())
        self.update()

    def Changewait_time_bounds_upper(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].wait_time_bounds[1] = eval(self.txtwait_time_bounds_upper.GetValue())
        self.update()

    def Changetimebounded(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].timebounded = self.cmbtimebounded.GetSelection()
        self.parent.update()

    def Changeflight_time_bounds_lower(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].flight_time_bounds[0] = eval(self.txtflight_time_bounds_lower.GetValue())
        self.update()

    def Changeflight_time_bounds_upper(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].flight_time_bounds[1] = eval(self.txtflight_time_bounds_upper.GetValue())
        self.update()

    def Changearrival_date_bounds_lower(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_date_bounds[0] = eval(self.txtarrival_date_bounds_lower.GetValue())
                
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_date_bounds[0] > 2400000.5:
            self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_date_bounds[0] -= 2400000.5

        self.update()

    def ChangeArrivalDateLowerCalendar(self, e):
        date = self.ArrivalDateLowerCalendar.GetDate()
        date = date.FromTimezone(wx.DateTime.TimeZone(offset=0))
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_date_bounds[0] = date.GetMJD()
        self.update()

    def ChangeArrivalDateUpperCalendar(self, e):
        date = self.ArrivalDateUpperCalendar.GetDate()
        date = date.FromTimezone(wx.DateTime.TimeZone(offset=0))
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_date_bounds[1] = date.GetMJD()
        self.update()

    def Changearrival_date_bounds_upper(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_date_bounds[1] = eval(self.txtarrival_date_bounds_upper.GetValue())
        
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_date_bounds[1] > 2400000.5:
            self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_date_bounds[1] -= 2400000.5
        
        self.update()

    def Changebounded_departure_date(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].bounded_departure_date = self.chkbounded_departure_date.GetValue()
        self.update()

    def Changedeparture_date_bounds_lower(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_date_bounds[0] = eval(self.txtdeparture_date_bounds_lower.GetValue())

        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_date_bounds[0] > 2400000.5:
            self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_date_bounds[0] -= 2400000.5

        self.update()

    def ChangedepartureDateLowerCalendar(self, e):
        date = self.departureDateLowerCalendar.GetDate()
        date = date.FromTimezone(wx.DateTime.TimeZone(offset=0))
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_date_bounds[0] = date.GetMJD()
        self.update()

    def ChangedepartureDateUpperCalendar(self, e):
        date = self.departureDateUpperCalendar.GetDate()
        date = date.FromTimezone(wx.DateTime.TimeZone(offset=0))
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_date_bounds[1] = date.GetMJD()
        self.update()

    def Changedeparture_date_bounds_upper(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_date_bounds[1] = eval(self.txtdeparture_date_bounds_upper.GetValue())
        
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_date_bounds[1] > 2400000.5:
            self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_date_bounds[1] -= 2400000.5
        
        self.update()

    def Changeinitial_impulse_bounds_lower(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].initial_impulse_bounds[0] = eval(self.txtinitial_impulse_bounds_lower.GetValue())
        self.update()

    def Changeinitial_impulse_bounds_upper(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].initial_impulse_bounds[1] = eval(self.txtinitial_impulse_bounds_upper.GetValue())
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].initial_impulse_bounds[1] < 1.0e-8:
            self.missionoptions.Journeys[self.missionoptions.ActiveJourney].initial_impulse_bounds[1] = 1.0e-8

        self.update()
        
    def Changeforce_free_point_direct_insertion_along_velocity_vector(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].force_free_point_direct_insertion_along_velocity_vector = int(self.chkforce_free_point_direct_insertion_along_velocity_vector.GetValue())
        self.update()

    def Changedeparture_type(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_type = self.cmbdeparture_type.GetSelection()
        self.parent.update()

    def Changedeparture_class(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_class = self.cmbdeparture_class.GetSelection()
        self.parent.update()

    def Changedeparture_ellipsoid_axes(self, e):
        e.Skip()
        temp_axes = eval(self.txtdeparture_ellipsoid_axes.GetValue())
        
        while len(temp_axes) < 3:
            temp_axes.append(1.0e+6)
        while len(temp_axes) > 3:
            temp_axes.pop()

        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_ellipsoid_axes = temp_axes
        self.update()

    def ClickButtondeparture_ellipsoid_axes(self, e):
        e.Skip()
        self.universe = Universe.Universe(self.missionoptions.universe_folder + "\\" + self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body + ".emtg_universe")
        bodyIndex = self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list[0] - 1
        reference_sphere = 0.0
        if bodyIndex > 0: #body in the universe
            thisBody = self.universe.bodies[bodyIndex]
            reference_sphere = thisBody.SMA * (1.0 - thisBody.ECC) * (thisBody.mu / (3.0 * self.universe.mu)) ** (1.0 / 3.0)
        elif bodyIndex == 0: #central body
            reference_sphere = self.universe.central_body_radius
        else: #universe SOI
            reference_sphere = self.universe.r_SOI

        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_ellipsoid_axes = [reference_sphere] * 3
        self.update()

    def Changezero_turn_flyby_distance(self, e):
        e.Skip()

        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].zero_turn_flyby_distance = eval(self.txtzero_turn_flyby_distance.GetValue())
        self.update()

    def Changeoverride_flyby_altitude_bounds(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_flyby_altitude_bounds = self.chkoverride_flyby_altitude_bounds.GetValue()
        self.parent.update()
        
    def Changeflyby_altitude_boundsLower(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].flyby_altitude_bounds[0] = eval(self.txtflyby_altitude_boundsLower.GetValue())
        self.parent.update()

    def Changeflyby_altitude_boundsUpper(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].flyby_altitude_bounds[1] = eval(self.txtflyby_altitude_boundsUpper.GetValue())
        self.parent.update()

    def Changeescape_spiral_starting_radius(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].escape_spiral_starting_radius = eval(self.txtescape_spiral_starting_radius.GetValue())
        self.update()

    def Changeescape_spiral_final_radius(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].escape_spiral_final_radius = eval(self.txtescape_spiral_final_radius.GetValue())
        self.update()
        
    def Changearrival_type(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_type = self.cmbarrival_type.GetSelection()
        self.parent.update()

    def Changearrival_class(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_class = self.cmbarrival_class.GetSelection()
        self.parent.update()

    def Changearrival_ellipsoid_axes(self, e):
        e.Skip()
        temp_axes = eval(self.txtarrival_ellipsoid_axes.GetValue())
        
        while len(temp_axes) < 3:
            temp_axes.append(1.0e+6)
        while len(temp_axes) > 3:
            temp_axes.pop()

        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_ellipsoid_axes = temp_axes
        self.update()

    def ClickButtonarrival_ellipsoid_axes(self, e):
        e.Skip()
        self.universe = Universe.Universe(self.missionoptions.universe_folder + "\\" + self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body + ".emtg_universe")
        bodyIndex = self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list[1] - 1
        reference_sphere = 0.0
        if bodyIndex > 0: #body in the universe
            thisBody = self.universe.bodies[bodyIndex]
            reference_sphere = thisBody.SMA * (1.0 - thisBody.ECC) * (thisBody.mu / (3.0 * self.universe.mu)) ** (1.0 / 3.0)
        elif bodyIndex == 0: #central body
            reference_sphere = self.universe.central_body_radius
        else: #universe SOI
            reference_sphere = self.universe.r_SOI

        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_ellipsoid_axes = [reference_sphere] * 3
        self.update()

    def Changecapture_spiral_starting_radius(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].capture_spiral_starting_radius = eval(self.txtcapture_spiral_starting_radius.GetValue())
        self.update()

    def Changecapture_spiral_final_radius(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].capture_spiral_final_radius = eval(self.txtcapture_spiral_final_radius.GetValue())
        self.update()

    def Changeimpact_momentum_enhancement_factor(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].impact_momentum_enhancement_factor = eval(self.txtimpact_momentum_enhancement_factor.GetValue())

    def Changeephemeris_pegged_orbit_insertion_SMA(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].ephemeris_pegged_orbit_insertion_SMA = eval(self.txtephemeris_pegged_orbit_insertion_SMA.GetValue())
        self.update()

    def Changeephemeris_pegged_orbit_insertion_ECC(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].ephemeris_pegged_orbit_insertion_ECC = eval(self.txtephemeris_pegged_orbit_insertion_ECC.GetValue())
        self.update()

    def Changeforced_terminal_coast(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].forced_terminal_coast = eval(self.txtforced_terminal_coast.GetValue())
        self.update()
                
    def Changeforced_initial_coast(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].forced_initial_coast = eval(self.txtforced_initial_coast.GetValue())
        self.update()

    def Changefinal_velocity0(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].final_velocity[0] = eval(self.txtfinal_velocity0.GetValue())
        self.update()

    def Changefinal_velocity1(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].final_velocity[1] = eval(self.txtfinal_velocity1.GetValue())
        self.update()

    def Changefinal_velocity2(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].final_velocity[2] = eval(self.txtfinal_velocity2.GetValue())
        self.update()

    def ChangeFreePointArrival_print_target_spec(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].FreePointArrival_print_target_spec = int(self.chkFreePointArrival_print_target_spec.GetValue())
        self.update()

    def Changeoverride_duty_cycle(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_duty_cycle = int(self.chkoverride_duty_cycle.GetValue())
        self.update()
    
    def Changeduty_cycle(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].duty_cycle = eval(self.txtduty_cycle.GetValue())
        self.update()
        
    def Changeoverride_PropagatorType(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_PropagatorType = int(self.chkoverride_PropagatorType.GetValue())
        self.update()

    def ChangepropagatorType(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].propagatorType = self.cmbpropagatorType.GetSelection()
        self.parent.update()

    def Changeoverride_integration_step_size(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].override_integration_step_size = int(self.chkoverride_integration_step_size.GetValue())
        self.update()

    def Changeintegration_step_size(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].integration_step_size = eval(self.txtintegration_step_size.GetValue())
        self.update()

    def ChangeCoastPhaseMatchPointFraction(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].CoastPhaseMatchPointFraction = eval(self.txtCoastPhaseMatchPointFraction.GetValue())

    def ChangeCoastPhaseForwardIntegrationStepLength(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].CoastPhaseForwardIntegrationStepLength = eval(self.txtCoastPhaseForwardIntegrationStepLength.GetValue())

    def ChangeCoastPhaseBackwardIntegrationStepLength(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].CoastPhaseBackwardIntegrationStepLength = eval(self.txtCoastPhaseBackwardIntegrationStepLength.GetValue())

    def Changesequence(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].sequence = eval(self.txtsequence.GetValue())
        self.update()

    def Changeperiapse_burns(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].enable_periapse_burns = int(self.chkperiapse_burns.GetValue())
        self.update()
        
    def Changeperturbation_bodies(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].perturbation_bodies = eval(self.txtperturbation_bodies.GetValue())
        self.update()

    def Clickdestination_list(self, e):
        #call dialog to choose destination list
        self.universe = Universe.Universe(self.missionoptions.universe_folder + "\\" + self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body + ".emtg_universe")
        dlg = BodyPicker.DestinationPicker(self, -1,
                                           self.universe,
                                           self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list[0],
                                           self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list[1])

        dlg.ShowModal()

        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list[0] = dlg.destination1
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].destination_list[1] = dlg.destination2

        dlg.Destroy()
        self.update()

    def Clickjourney_central_body(self, e):
        #call dialog to choose destination list
        dlg = wx.FileDialog(self, "Choose an emtg_universe file", self.missionoptions.universe_folder, self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body+".emtg_universe", "*.emtg_universe", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
        else:
            filename = self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body+".emtg_universe"

        fileparts = filename.split(".")
        dlg.Destroy()

        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body = fileparts[0]

        self.update()

    def Clicksequence(self, e):
        #call dialog to choose sequence list
        self.universe = Universe.Universe(os.path.join(self.missionoptions.universe_folder, self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body + ".emtg_universe"))
        dlg = BodyPicker.SequencePicker(self, -1,
                                        self.universe,
                                        self.missionoptions,
                                        self.missionoptions.ActiveJourney)

        dlg.ShowModal()

        dlg.Destroy()

        self.parent.update()

    def Clickperturbation_bodies(self, e):
        #call dialog to choose perturbation list
        self.universe = Universe.Universe(os.path.join(self.missionoptions.universe_folder, self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_central_body + ".emtg_universe"))

        dlg = wx.MultiChoiceDialog(self, "", "Choose perturbation bodies", choices=self.universe.perturbation_menu)

        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].perturbation_bodies != [0]:
            selections = []
            for i in self.missionoptions.Journeys[self.missionoptions.ActiveJourney].perturbation_bodies:
                #get the index in the perturbation menu of the body of interest
                for j in range(0, len(self.universe.perturbation_indices)):
                    if self.universe.perturbation_indices[j] == i - 1:
                        selections.append(j)
            dlg.SetSelections(selections)

        if dlg.ShowModal() == wx.ID_OK:
            PertubationListIndices = dlg.GetSelections()
            if len(PertubationListIndices) > 0:
                self.missionoptions.Journeys[self.missionoptions.ActiveJourney].perturbation_bodies = []
                for i in PertubationListIndices:
                    self.missionoptions.Journeys[self.missionoptions.ActiveJourney].perturbation_bodies.append(self.universe.perturbation_indices[i]+1)
            else:
                self.missionoptions.Journeys[self.missionoptions.ActiveJourney].perturbation_bodies = []

        dlg.Destroy()
        self.update()

    def Changefreeze_decision_variables(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].freeze_decision_variables = self.chkfreeze_decision_variables.GetValue()
        self.update()

    def Changejourney_end_deltav(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_end_deltav = eval(self.txtjourney_end_deltav.GetValue())
        self.update()

    def Changejourney_end_TCM(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].journey_end_TCM = eval(self.txtjourney_end_TCM.GetValue())
        self.update()

    #####################
    # drag event handlers
    #####################
    
    def ChangeEnableDrag(self, e):
        # done
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].perturb_drag = int(self.chkEnableDrag.GetValue())
        self.update()
        return

    def ChangeSpacecraftDragArea(self, e):
        # done
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].spacecraft_drag_area = eval(self.txtSpacecraftDragArea.GetValue())
        self.update()
        return

    def ChangeSpacecraftDragCoefficient(self, e):
        # done
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].coefficient_of_drag = eval(self.txtSpacecraftDragCoefficient.GetValue())
        self.update()
        return

    def ChangeAtmosphericDensityModel(self, e):
        # done
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].AtmosphericDensityModelKey = self.cmbAtmosphericDensityModel.GetStringSelection()
        self.parent.update()
        return

    def ClickAtmosphericDensityModelDataFile(self, e):
        # based on MissionPanel.py ClickThrottleTableFileButton method
        #file load dialog to get name of 
        dlg = wx.FileDialog(self, "Select an atmospheric density model data file", "", "", '*', wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.missionoptions.Journeys[self.missionoptions.ActiveJourney].AtmosphericDensityModelDataFile = dlg.GetPath()
            self.txtAtmosphericDensityModelDataFile.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].AtmosphericDensityModelDataFile)
        dlg.Destroy()
        self.update()
        return

    def ChangeAtmosphericDensityModelDataFile(self, e):
        # done
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].AtmosphericDensityModelDataFile = self.txtAtmosphericDensityModelDataFile.GetValue()
        self.update()
        return