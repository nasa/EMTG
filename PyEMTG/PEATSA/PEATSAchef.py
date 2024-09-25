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
PEATSAchef.py
==================

File holds the PEATSAchef class, which creates the .emtgopt files for PEATSA.

"""


class PEATSAchef(object):
    """
    Creates the .emtgopt files for PEATSA
    
    Parameters
    ----------
    None.

    Returns
    -------
    None.
    
    """
    def __init__(self):
        
        # Nothing to do here
        pass  

    def parse_ephem_file(self,ephem_file):
        import State
        
        # Get a handle to the file
        file_handle = open(ephem_file)
    
        # Initialize the output data array
        states = []
    
        # Loop through all lines
        for line in file_handle:
        
            if line[0] == "#":
                continue
        
            states.append(State.state())
        
            states[-1].parse_ephem_line(line)
                    
        # Return the ephemeris data
        return states
        
    def parse_maneuver_spec(self,maneuver_spec_file):
        import Burn
        
        burns = []
    
        fileHandle = open(maneuver_spec_file,'r')
    
        ctr = -1
    
        all_lines = fileHandle.readlines()
            
        for line in all_lines:
            if line.startswith("<EVENTNAME>"):
                continue
            burns.append(Burn.burn())
            
            nBurns = burns[-1].parse_maneuver_spec_line(line,0)
    
            for burn_idx in range(1,nBurns):
                burns.append(Burn.burn())
                burns[-1].parse_maneuver_spec_line(line,burn_idx)
            
                if burns[-2].start_epoch + burns[-2].duration / 86400.0 > burns[-1].start_epoch:
                    burns[-1].start_epoch = burns[-2].start_epoch + burns[-2].duration / 86400.0
                    burns[-1].assumed_start = burns[-1].start_epoch
    
        return burns
        
    def process_reference(self,PEATSAorder):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import PEATSAbox

        # Grab the data from the reference mission
        PEATSAorder.reference = PEATSAbox.PEATSAbox(PEATSAorder,PEATSAorder.reference_files_location,PEATSAorder.reference_options_file_name,PEATSAorder.reference_mission_file_name)

        try:
            import PyHardware
            use_pyhardware = True
        except:
            print ("PyHardware not available. Might cause errors. Youve been warned")
            use_pyhardware = False

        # Move the emtg_spacecraftopt file over, if any
        if PEATSAorder.reference.PEATSAdough.SpacecraftModelInput == 1 and use_pyhardware:
            PEATSAorder.reference.PEATSAstone.write_output_file(PEATSAorder.hardware_dir + "/" + PEATSAorder.reference.PEATSAdough.SpacecraftOptionsFile)
            
            # Loop through all stages so that all throttle table files can be moved to the new folder
            for throttle_options in PEATSAorder.reference.PEATSAcooking_temperatures:
                throttle_options.print_throttle_file(PEATSAorder.hardware_dir + "/" + throttle_options.FileName.split("/")[-1])
        
        # Move the launch vehicle library file over, if any
        if PEATSAorder.reference.PEATSAdough.SpacecraftModelInput != 2 and use_pyhardware:
            PEATSAorder.reference.PEATSAovenmitt.write_output_file(PEATSAorder.hardware_dir + "/" + PEATSAorder.reference.PEATSAdough.LaunchVehicleLibraryFile,True)
        
        return PEATSAorder
    
    def create_PEATSAs(self,PEATSAorder):
        if PEATSAorder.PEATSA_type == 0:
            boxes = self.Custom(PEATSAorder)
        elif PEATSAorder.PEATSA_type == 1:
            boxes = self.MissedThrust(PEATSAorder)
        elif PEATSAorder.PEATSA_type == 2:
            boxes = self.TradeStudy(PEATSAorder)
        elif PEATSAorder.PEATSA_type == 3:
            boxes = self.TradeStudy(PEATSAorder)
        elif PEATSAorder.PEATSA_type == 4:
            boxes = self.BatchRun(PEATSAorder)
        elif PEATSAorder.PEATSA_type == 5:
            boxes = self.MonteCarlo(PEATSAorder)
        else:
            raise Exception("Unknown PEATSA type")
        return boxes
    
    def MonteCarlo(self,PEATSAorder,prefix = "",MObase = None,SObase = None):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        # sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG/SimpleMonteCarlo")
        # import ConOpsPeriod
        # import SimpleSEPMonteCarloStep
        import copy
        import PEATSAwaiter
        import PEATSAbox
        import PyHardware        
        
        # Create a waiter for splining things
        waiter = PEATSAwaiter.PEATSAwaiter()
                
        # Process the reference mission
        PEATSAorder = self.process_reference(PEATSAorder)
        
        # Iniitalize the stack of cases
        boxes2run = []
        all_boxes = []
        
        # Rename the EMTG objects for simplicity
        if PEATSAorder.thin_crust:
            RefM = Mission.Mission(PEATSAorder.reference.PEATSApath + "/" + PEATSAorder.reference.PEATSAcrust_path)
        else:
            RefM = PEATSAorder.reference.PEATSAcrust        
        RefMO = PEATSAorder.reference.PEATSAdough
        RefSO = PyHardware.SpacecraftOptions(RefMO.HardwarePath.rstrip("\r\n ") + RefMO.SpacecraftOptionsFile.rstrip("\r\n "))
        
        if MObase == None:
            MObase = RefMO
        if SObase == None:
            SObase = RefSO        
        
        burns  = self.parse_maneuver_spec(PEATSAorder.reference_files_location + "/" + PEATSAorder.maneuver_file)
        states = self.parse_ephem_file(PEATSAorder.reference_files_location + "/" + PEATSAorder.reference_spk_data_file)
        
        if PEATSAorder.stop_after_journey < 0:
            PEATSAorder.stop_after_journey += len(RefMO.Journeys)
            
        monte_carlo_depth_index = -1
        
        for jdx, RefJO in enumerate(RefMO.Journeys):
            
            if jdx > PEATSAorder.stop_after_journey:
                continue
                                        
            # Determine if this is a multi-phase Journey.
            # First, assume it only has one phase
            nPhases = 1
            # Loop through all of the sequence values in this Journey
            for seq_idx in range(0,len(RefMO.Journeys[jdx].sequence)):
                # Check if this sequence ID is non-zero, in which case there is a flyby and this is a multiphase journey
                if (RefMO.Journeys[jdx].sequence[seq_idx]):
                    # It is, update the counter
                    nPhases += 1
            
            # Check if the options override the default number of steps
            if RefJO.override_num_steps:
                # It did, grab the actual value
                reference_timesteps = RefJO.number_of_steps
            else:
                # Get the reference journeys number of timesteps
                reference_timesteps = RefMO.num_timesteps
            
            # If there are more than 9 phases it will break things, so throw an error
            if nPhases > 9:
                raise Exception("Phase number finding code in the decision vector re-creation cant handle more than 9 phases currently")
            
            # Set the start date for missed thrust cases
            journey_start_epoch = RefM.Journeys[jdx].missionevents[0].JulianDate
            
            # Loop through all of the phases
            for pdx in range(0,nPhases):
                    
                # Loop through all of the timesteps
                for nStep in range(reference_timesteps):                
                    
                    loc_prefix = "j" + str(jdx) + "p" + str(pdx) + "PSFB_Step" + str(nStep)
                                    
                    startBurn = None                
                    for burn in burns:
                        if loc_prefix in burn.name:
                            startBurn = burn
                            break              
                    if startBurn == None:
                        continue
                    elif startBurn.start_epoch < PEATSAorder.epoch_range[0]:
                        continue
                    elif startBurn.start_epoch > PEATSAorder.epoch_range[1]:
                        continue
                         
                    monte_carlo_depth_index += 1
                                    
                    # Loop through all of the samples
                    for nSample in range(PEATSAorder.nSamples):
            
                        # Create a copy of the reference options
                        MO = copy.deepcopy(MObase) 
                        # Create a copy of the reference stage options
                        SO = copy.deepcopy(SObase) 
                                               
                        for journey_removal_index in range(PEATSAorder.stop_after_journey+1,len(RefMO.Journeys)):
                            MO.number_of_journeys -= 1
                            MO.Journeys.pop(-1)                      
                        
                        # Create a counter for the number of stages to remove
                        stage_removal_counter = 0
                        
                        # Loop through previous journeys up to the current one
                        for journey_removal_index in range(0,jdx):
                        
                            # Check if staging occured after departure
                            if RefMO.Journeys[journey_removal_index].stage_after_departure:
                                stage_removal_counter += 1
                        
                            # Check if staging occured after arrival
                            if RefMO.Journeys[journey_removal_index].stage_after_arrival:
                                stage_removal_counter += 1
                            
                            # Remove the journey    
                            MO.Journeys.pop(0)
                            MO.number_of_journeys -= 1
                                                
                        # Check if staging occured after departure
                        if MO.Journeys[0].stage_after_departure:
                            stage_removal_counter += 1                            
                                                                        
                        # Loop through the number of stages to remove and remove them
                        for stage_idx in range(0,stage_removal_counter):                                
                            SO.remove_stage(0)
                                                
                        MO.Journeys[0].override_num_steps = 1
                        MO.Journeys[0].number_of_steps = reference_timesteps - nStep
                            
                        # Check if the options override the default number of steps
                        if RefJO.override_num_steps:
                            # It did, grab the actual value
                            reference_timesteps = RefJO.number_of_steps
                        else:
                            # Get the reference journeys number of timesteps
                            reference_timesteps = RefMO.num_timesteps
                                                
                        # Set the state as a free point in space
                        MO.Journeys[0].departure_type = 2
                        # Set the departure class to be a free point
                        MO.Journeys[0].departure_class = 1
                        # Set the state as a free point in space
                        MO.Journeys[0].destination_list[0] = -1
                        # Set the state for the start_date
                        MO.Journeys[0].departure_elements_reference_epoch = startBurn.start_epoch - 2400000.5
                        # Set the date of the event 
                        MO.launch_window_open_date = startBurn.start_epoch - 2400000.5
                    
                        # Turn on mbh in case its off
                        MO.run_inner_loop = 1
                        
                        if monte_carlo_depth_index < 2:
                            last_state = states[0]
                            for state in states:
                                if state.epoch == startBurn.start_epoch:
                                    initial_state = state
                                    break
                                elif state.epoch > startBurn.start_epoch:
                                    if abs(state.epoch-startBurn.start_epoch) < 1e-4:
                                        initial_state = state
                                        break
                                    elif abs(last_state.epoch - startBurn.start_epoch) < 1e-4:
                                        initial_state = state
                                        break
                                    else:
                                        raise Exception("Cant find first burn in ephemeris")
                                last_state = state
                            if monte_carlo_depth_index == 0:
                                MO.user_data["InitialStateWithAllErrors"] = initial_state.getList()
                                MO.user_data["InitialStateWithAllErrors"].append(initial_state.mass)
                            MO.user_data["InitialConditionsSet"] = 1 
                            MO.Journeys[0].departure_elements = initial_state.getList()
                            MO.maximum_mass = initial_state.mass
                            MO.run_inner_loop = 0
                            
                            MO.trialX = []
                            
                            add_all = False
                                                                      
                            for idx,entry in enumerate(RefM.DecisionVector):
                                descrip = RefM.Xdescriptions[idx]
                                
                                if "FBLT" in descrip or "MGALT" in descrip:
                                    print(descrip)
                                    raise Exception("Only PSFB is setup to work with this so far")
                                
                                jdx2 = int(descrip.split("p")[0].lstrip("j"))
                                if jdx2 > jdx and jdx2 <= PEATSAorder.stop_after_journey:
                                    MO.trialX.append(["j" + str(jdx2-jdx) + descrip.lstrip("j0123456789"),entry])
                                elif jdx2 == jdx:
                                    if "phase flight time" in descrip:
                                        add_all = True
                                        MO.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state x",MO.Journeys[0].departure_elements[0]])
                                        MO.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state y",MO.Journeys[0].departure_elements[1]])
                                        MO.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state z",MO.Journeys[0].departure_elements[2]])
                                        MO.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state xdot",MO.Journeys[0].departure_elements[3]])
                                        MO.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state ydot",MO.Journeys[0].departure_elements[4]])
                                        MO.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state zdot",MO.Journeys[0].departure_elements[5]])
                                        # Add the mass
                                        MO.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state mass",MO.maximum_mass])
                                        # Add the epoch
                                        MO.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state epoch",startBurn.start_epoch-2400000.5])
                                        # Add the flight time                                            
                                        MO.trialX.append(["j0p0PSFB: phase flight time",RefM.Journeys[jdx].missionevents[-1].JulianDate-startBurn.start_epoch])
                                                           
                                    elif "PSFB: virtual chemical" in descrip and "Step" not in descrip:
                                        add_all = False
                                    elif "PSFB: virtual electric" in descrip and "Step" not in descrip:
                                        # Add the chem propellant
                                        MO.trialX.append(["j0p0PSFB: virtual chemical fuel",0.0])
                                        # Add the ep propellant
                                        MO.trialX.append(["j0p0PSFB: virtual electric propellant",entry])
                                        virtual_ep_idx = len(MO.trialX) - 1
                                    elif add_all:
                                        MO.trialX.append(["j0" + descrip.lstrip("j0123456789"),entry])                                        
                                    elif "Step" in descrip:
                                        step2 = int(descrip.split("_Step")[1].split(":")[0])
                                        if step2 >= nStep:
                                            MO.trialX.append(["j0p" + str(pdx) + "PSFB_Step" + str(step2-nStep) + ":" + descrip.split(":")[1],entry])
                                            if "virtual electric" in descrip:
                                                if step2 == nStep:
                                                    virtual_ep_ctr = entry 
                                                    MO.user_data["PreviousPropUsage"] = entry
                                                    MO.trialX[-1][1] -= virtual_ep_ctr
                                                else:
                                                    MO.trialX[-1][1] -= virtual_ep_ctr
                            MO.trialX[virtual_ep_idx][1] -= virtual_ep_ctr   
                        else:
                            for idx,entry in enumerate(RefM.DecisionVector):
                                descrip = RefM.Xdescriptions[idx]
                                
                                if "FBLT" in descrip or "MGALT" in descrip:
                                    print(descrip)
                                    raise Exception("Only PSFB is setup to work with this so far")
                                
                                jdx2 = int(descrip.split("p")[0].lstrip("j"))
                                if jdx2 > jdx:
                                    raise Exception("We must have missed it")
                                elif jdx2 == jdx:                                   
                                    if "Step" in descrip:
                                        step2 = int(descrip.split("_Step")[1].split(":")[0])
                                        if step2 == nStep:
                                            if "virtual electric" in descrip:
                                                MO.user_data["PreviousPropUsage"] = entry
                                                break                                      
                                    
                        MO.allow_initial_mass_to_vary = 0
                        
                        # Allow the initial state to propagate (aka coast)
                        MO.Journeys[0].AllowJourneyFreePointDepartureToPropagate = 0    
                        # The state is fixed
                        MO.Journeys[0].departure_elements_vary_flag = [0]*6
                        # The state is cartesian
                        MO.Journeys[0].departure_elements_type = 0
                        # The state is equatorial
                        MO.Journeys[0].departure_elements_frame = 0
                        # Set the date requirement
                        MO.Journeys[0].wait_time_bounds = [0.0,1e-8]
                        # Turn off the initial coast
                        MO.Journeys[0].forced_initial_coast = 0
                        # Turn off journey departure date bounds
                        MO.Journeys[0].bounded_departure_date = 0
                            
                        # Use user data to assign the monte carlo trackers    
                        MO.user_data["MonteCarloSampleIndex"] = nSample   
                        MO.user_data["MonteCarloDepthIndex"] = monte_carlo_depth_index  
                        MO.user_data["MonteCarloOriginalJourney"] = jdx 
                        MO.user_data["InitialConditionsSet"] = 0                         
                        
                        # Create a list of the phase distance constraints that will be in the new problem
                        new_boundary_constraints = []
                    
                        # Loop through all of the existing constraints to remove any from previous phases
                        for constraint in MO.BoundaryConstraintDefinitions:
                        
                            # Split the string by underscores to get the prefix
                            constraint_prefix = constraint.split("_")[0]
                        
                            # Get the journey
                            constraint_journey = int(constraint_prefix.lstrip("j").split("p")[0])
                        
                            # Check if the constraint is from a previous journey
                            if constraint_journey < jdx:
                                # It is. Do not add it to the new list
                                continue                  
                            # We do need to add it
                            else:
                            
                                # We need to get the phase number
                                phase_number = int(constraint_prefix.split("p")[1])
                            
                                # Check if its still the same journey, but the phase is from a previous phase
                                if constraint_journey == jdx and phase_number < pdx:
                                    # We do not need to add it!
                                    continue
                                else:
                                    # We do need to add it
                                                            
                                    # Check if we need to decrement the phase
                                    if constraint_journey == jdx:
                                        # Decrement the phase
                                        phase_number -= pdx
                                                        
                                    # Decrement the journey number by how many were popped
                                    constraint_journey -= jdx
                            
                                    # Create new prefix and add it to the list
                                    new_boundary_constraints.append("j" + str(constraint_journey) + "p" + str(phase_number) + constraint.split(constraint_prefix)[1])
                                                                                    
                        # Overwrite the new boundary constraint list onto the old list
                        MO.BoundaryConstraintDefinitions = new_boundary_constraints
                    
                        # Create a list of the phase distance constraints that will be in the new problem
                        new_phase_distance_constraints = []
                    
                        # Loop through all of the existing constraints to remove any from previous phases
                        for constraint in MO.PhaseDistanceConstraintDefinitions:
                        
                            # Split the string by underscores to get the prefix
                            constraint_prefix = constraint.split("_")[0]
                        
                            # Get the journey
                            constraint_journey = int(constraint_prefix.lstrip("j").split("p")[0])
                        
                            # Check if the constraint is from a previous journey
                            if constraint_journey < jdx:
                                # It is. Do not add it to the new list
                                continue            
                            else:
                                # We might need to add it
                            
                                # We need to get the phase number
                                phase_number = int(constraint_prefix.split("p")[1])
                            
                                # Check if its still the same journey, but the phase is from a previous phase
                                if constraint_journey == jdx and phase_number < pdx:
                                    # We do not need to add it!
                                    continue
                                else:      
                                    # We do need to add it
                                                            
                                    # Check if we need to decrement the phase
                                    if constraint_journey == jdx:
                                        # Decrement the phase
                                        phase_number -= pdx
                                                        
                                    # Decrement the journey number by how many were popped
                                    constraint_journey -= jdx
                            
                                    # Create new prefix and add it to the list
                                    new_phase_distance_constraints.append("j" + str(constraint_journey) + "p" + str(phase_number) + constraint.split(constraint_prefix)[1])
                                                                                    
                        # Overwrite the new phase constraint list onto the old list
                        MO.PhaseDistanceConstraintDefinitions = new_phase_distance_constraints                        

                        # Create a list of the phase distance constraints that will be in the new problem
                        new_maneuver_constraints = []
                    
                        # Loop through all of the existing constraints to remove any from previous phases
                        for constraint in MO.ManeuverConstraintDefinitions:
                        
                            # Split the string by underscores to get the prefix
                            constraint_prefix = constraint.split("_")[0]
                        
                            # Get the journey
                            constraint_journey = int(constraint_prefix.lstrip("j").split("p")[0])
                        
                            # Check if the constraint is from a previous journey
                            if constraint_journey < jdx:
                                # It is. Do not add it to the new list
                                continue              
                            else:    
                                # We might need to add it
                            
                                # We need to get the phase number
                                phase_number = int(constraint_prefix.split("b")[0].split("p")[1])
                            
                                # Check if its still the same journey, but the phase is from a previous phase
                                if constraint_journey == jdx and phase_number < pdx:
                                    # We do not need to add it!
                                    continue
                                else:
                                
                                    # We need to get the maneuver number
                                    burn_number = int(constraint_prefix.split("b")[1])
                                
                                    # We do need to add it
                                    if constraint_journey == jdx and phase_number == pdx and burn_number < nStep:
                                        # We do not need to add it
                                        continue
                                    else:
                                        # We do need to add it
                                                      
                                        # Check if we need to decrement the burn:
                                        if phase_number == pdx and constraint_journey == jdx:
                                            # Decrement the burn
                                            burn_number -= nStep
                                                                
                                        # Check if we need to decrement the phase
                                        if constraint_journey == jdx:
                                            # Decrement the phase
                                            phase_number -= pdx                             
                                                        
                                        # Decrement the journey number by how many were popped
                                        constraint_journey -= jdx
                            
                                        # Create new prefix and add it to the list
                                        new_maneuver_constraints.append("j" + str(constraint_journey) + "p" + str(phase_number) + "b" + str(burn_number) + constraint.split(constraint_prefix)[1])
                                                                                    
                        # Overwrite the new phase constraint list onto the old list
                        MO.ManeuverConstraintDefinitions = new_maneuver_constraints
                            
                        # Set the mission name
                        if monte_carlo_depth_index >= 2:
                            stepStr = str(monte_carlo_depth_index-2)
                        elif monte_carlo_depth_index == 0:
                            stepStr = "n2"
                        elif monte_carlo_depth_index == 1:
                            stepStr = "n1"
                        MO.mission_name = prefix + 'J' + str(jdx) + "_Step" + stepStr + '_MC'+ str(nSample)
                    
                        # Set the spacecraft file
                        MO.SpacecraftOptionsFile = MO.mission_name + '.emtg_spacecraftopt'
                    
                        # Update the peatsa hardware path
                        MO.HardwarePath = PEATSAorder.hardware_dir
                        
                        # Turn on the working directory overrides
                        MO.override_working_directory = 1
                        MO.forced_working_directory = PEATSAorder.results_directory
                        MO.override_mission_subfolder = 1
                        MO.forced_mission_subfolder = "/"
                                
                        # Turn on seeding
                        MO.seed_MBH = 1                    
                    
                        # Quiet snopt and mbh
                        MO.quiet_basinhopping = 1
                        MO.quiet_NLP = 1
                    
                        # Set the max mbh run time
                        MO.MBH_max_run_time = PEATSAorder.MBH_max_run_time

                        # Apply the override options
                        for option in PEATSAorder.override_options:
                            try:
                                if "SMO" in option[0]:
                                    continue
                                # Check the condition
                                if eval(option[0]):
                                    if "SMO" in option[1]:
                                        continue
                    
                                    if option[1].startswith("MO"):
                                        varname = option[1].split("MO.")[1].split("=")[0].rstrip(" ")
                                        if "Journeys" in varname:
                                            varname = varname.split("].")[1]
                                            obj = MO.Journeys[0]
                                        else:
                                            obj = MO
                    
                                        if "[" in varname:
                                            varname = varname.split("[")[0]  
                                            
                                        if varname not in dir(obj):
                                            logging.info("WARNING: Override might not work, because " + varname + " doesnt seem to be in PyEMTG Mission Options or Journey Options")
                                            logging.info("Option: " + str(option))
                            
                                    # It evaluated true, so apply it
                                    exec(option[1])
                            except Exception as e:
                                print("WARNING: Unable to override an option because:" + str(e))
                                logging.info("Option: " + str(option))
                                if "text" in dir(e):
                                    print("*******************************************")
                                    print(e.text)
                                    print(' ' * e.offset + "^")
                                    print("*******************************************")
                                    
                                    
                        # Write the options to file!
                        MO.write_options_file(PEATSAorder.cases_directory + "/" + MO.mission_name + '.emtgopt')
                        SO.write_output_file(PEATSAorder.hardware_dir + "/" + MO.SpacecraftOptionsFile)

                        # Only the boxes to run will be put in the boxes_out folder. So to make sure that all of them end up in the study, write all of them directly to the output folder
                        MO.write_options_file(PEATSAorder.results_directory + "/" + MO.mission_name + '.emtgopt')
                                    
                        # Create a peatsa box object for this case
                        box = PEATSAbox.PEATSAbox(PEATSAorder,PEATSAorder.cases_directory,MO.mission_name + ".emtgopt")
                        
                        if monte_carlo_depth_index < 2:
                            # Add the peatsa box object to the stack
                            boxes2run.append(box)
        return boxes2run        
                    
    def BatchRun(self,PEATSAorder):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        import copy
        import PEATSAbox
        
        MO = MissionOptions.MissionOptions(PEATSAorder.reference_files_location + "/" + PEATSAorder.reference_options_file_name)
        
        # Even though the list will only have one file in it, create an empty list for the output boxes
        boxes_out = []
        
        # Turn on the working directory overrides
        MO.override_working_directory = 1
        MO.forced_working_directory = PEATSAorder.results_directory
        MO.override_mission_subfolder = 1
        MO.forced_mission_subfolder = "/"
        
        # Create a peatsa box object for this case
        box = PEATSAbox.PEATSAbox(PEATSAorder,PEATSAorder.cases_directory,MO.mission_name + ".emtgopt")

        # Add the peatsa box object to the stack
        boxes_out.append(box)

        # Write the options to file!
        MO.write_options_file(PEATSAorder.cases_directory + "/" + MO.mission_name + '.emtgopt')
        
        return boxes_out
    
    def ThreeD_Flybys_to_PatchedConic(self,RefMO):
        import copy
        import math
        
        MO = copy.deepcopy(RefMO)       
        
        
        newTrialX = []
        
        to_pop = []        
        
        removed = 0
            
        remove_ctr = []
        orig_jdx = []
        new_jdx = []
        for jdx,JO in enumerate(MO.Journeys):
            new_jdx.append(0)
            if JO.phase_type[0] == 7:
                to_pop.append(jdx)
                removed += 1
            else:
                remove_ctr.append(removed)
                new_jdx[-1] = len(orig_jdx)
                orig_jdx.append(jdx)
        to_pop.reverse()
        for jdx in to_pop:
            MO.Journeys.pop(jdx)
            MO.number_of_journeys -= 1   
                 
        # Loop through the decision vector to get the phase flight times and number of steps
        all_phase_flight_times = {}
        for local_decision_vector_index in range(0,len(RefMO.trialX)):
            
            # Get the transcription type
            if "FBLT" in RefMO.trialX[local_decision_vector_index][0]:
                transcription = "FBLT"
            elif "PSFB" in RefMO.trialX[local_decision_vector_index][0]:
                transcription = "PSFB"
            elif "MGALT" in RefMO.trialX[local_decision_vector_index][0]:
                transcription = "MGALT"
            elif "CoastPhase" in RefMO.trialX[local_decision_vector_index][0]:
                transcription = "CoastPhase"
            else:
                continue
                                      
            # Check if this is a phase flight time variable
            if "phase flight time" in RefMO.trialX[local_decision_vector_index][0]:
                # It is! 
                
                # Get the prefix
                phase_string = RefMO.trialX[local_decision_vector_index][0].split(transcription)[0]
                
                # Set the phase flight time in the dictionary
                all_phase_flight_times.update({phase_string : float(RefMO.trialX[local_decision_vector_index][1])})
        

        
        if RefMO.Journeys[0].phase_type[0] == 7:
            MO.Journeys[0].departure_type = 0
            MO.Journeys[0].departure_class = 0
            MO.Journeys[0].forced_initial_coast += all_phase_flight_times["j0p0"]
            MO.Journeys[0].destination_list[0] = 3
            for entry in RefMO.trialX:
                if entry[0] == "j0p0CoastPhasePeriapseLaunchOrImpulsiveDeparture: event left state epoch":
                    newTrialX.append(["j0p0PSFBEphemerisPeggedLaunchDirectInsertion: event left state epoch",entry[1]])
                    break
            for entry in RefMO.trialX:
                if entry[0] == "j0p0CoastPhaseEphemerisReferencedInterceptInterior: event interface state vMAG":
                    newTrialX.append(["j0p0PSFBEphemerisPeggedLaunchDirectInsertion: magnitude of outgoing velocity asymptote",entry[1]])
                    break
            for entry in RefMO.trialX:
                if entry[0] == "j0p0CoastPhaseEphemerisReferencedInterceptInterior: event interface state vRA":
                    newTrialX.append(["j0p0PSFBEphemerisPeggedLaunchDirectInsertion: RA of departure asymptote",entry[1]])
                    break
            for entry in RefMO.trialX:
                if entry[0] == "j0p0CoastPhaseEphemerisReferencedInterceptInterior: event interface state vDEC":
                    newTrialX.append(["j0p0PSFBEphemerisPeggedLaunchDirectInsertion: DEC of departure asymptote",entry[1]])
                    break
        
        for jdx,JO in enumerate(MO.Journeys):
            if orig_jdx[jdx] +1 != len(RefMO.Journeys):
                if RefMO.Journeys[orig_jdx[jdx] + 1].phase_type[0] == 7:
                    JO.arrival_class = 0
                    JO.forced_terminal_coast += all_phase_flight_times["j" + str(orig_jdx[jdx]+1) + "p0"]
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]) + "p0PSFBEphemerisReferencedInterceptExterior: event left state mass":
                            newTrialX.append(["j" + str(jdx) + "p0PSFBEphemerisPeggedIntercept: event left state mass",entry[1]])
                            break
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]) + "p0PSFBEphemerisReferencedInterceptExterior: event interface state vMAG":
                            vMag = float(entry[1])
                            break
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]) + "p0PSFBEphemerisReferencedInterceptExterior: event interface state vRA":
                            vRa = float(entry[1])
                            break
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]) + "p0PSFBEphemerisReferencedInterceptExterior: event interface state vDEC":
                            vDec = float(entry[1])
                            break
                    newTrialX.append(["j" + str(jdx) + "p0PSFBEphemerisPeggedIntercept: V_infinity_x",vMag * math.cos(vRa) * math.cos(vDec)])
                    newTrialX.append(["j" + str(jdx) + "p0PSFBEphemerisPeggedIntercept: V_infinity_y",vMag * math.sin(vRa) * math.cos(vDec)])
                    newTrialX.append(["j" + str(jdx) + "p0PSFBEphemerisPeggedIntercept: V_infinity_z",vMag * math.sin(vDec)])
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]) + "p0PSFBEphemerisReferencedInterceptExterior: virtual chemical fuel":
                            newTrialX.append(["j" + str(jdx) + "p0PSFBEphemerisPeggedIntercept: virtual chemical fuel",entry[1]])
                            break
            if jdx != 0:
                if RefMO.Journeys[orig_jdx[jdx] - 1 ].phase_type[0] == 7:
                    JO.departure_type = 3
                    JO.departure_class = 0
                    JO.forced_initial_coast += all_phase_flight_times["j" + str(orig_jdx[jdx]-1) + "p0"]
                    JO.destination_list[0] = MO.Journeys[jdx-1].destination_list[1]
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]-1) + "p0CoastPhaseEphemerisReferencedInterceptInterior: event interface state vMAG":
                            vMag = float(entry[1])
                            break
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]-1) + "p0CoastPhaseEphemerisReferencedInterceptInterior: event interface state vRA":
                            vRa = float(entry[1])
                            break
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]-1) + "p0CoastPhaseEphemerisReferencedInterceptInterior: event interface state vDEC":
                            vDec = float(entry[1])
                            break
                    newTrialX.append(["j" + str(jdx) + "p0PSFBEphemerisPeggedUnpoweredFlyby: V_infinity_x",vMag * math.cos(vRa) * math.cos(vDec)])
                    newTrialX.append(["j" + str(jdx) + "p0PSFBEphemerisPeggedUnpoweredFlyby: V_infinity_y",vMag * math.sin(vRa) * math.cos(vDec)])
                    newTrialX.append(["j" + str(jdx) + "p0PSFBEphemerisPeggedUnpoweredFlyby: V_infinity_z",vMag * math.sin(vDec)])
                    
            if JO.destination_list[0] in JO.perturbation_bodies:
                JO.perturbation_bodies.remove(JO.destination_list[0])  
                JO.number_of_perturbation_bodies -= 1 
            if JO.destination_list[1] in JO.perturbation_bodies:
                JO.perturbation_bodies.remove(JO.destination_list[1])   
                JO.number_of_perturbation_bodies -= 1 
        
        # Create a list of the phase distance constraints that will be in the new problem
        new_boundary_constraints = []
        
        # Loop through all of the existing constraints to remove any from previous phases
        for constraint in MO.BoundaryConstraintDefinitions:
            
            # Split the string by underscores to get the prefix
            constraint_prefix = constraint.split("_")[0]
            
            # Get the journey
            constraint_journey = int(constraint_prefix.lstrip("j").split("p")[0])
            
            # Check if the constraint is from a previous journey
            if constraint_journey in to_pop:
                # It is. Do not add it to the new list
                continue                  
            # We do need to add it
            else:
                
                # We need to get the phase number
                phase_number = int(constraint_prefix.split("p")[1])                                        
            
                # Create new prefix and add it to the list
                new_boundary_constraints.append("j" + str(new_jdx[constraint_journey]) + "p" + str(phase_number) + constraint.split(constraint_prefix)[1])
                                                                        
        # Overwrite the new boundary constraint list onto the old list
        MO.BoundaryConstraintDefinitions = new_boundary_constraints
        
        # Create a list of the phase distance constraints that will be in the new problem
        new_phase_distance_constraints = []
        
        # Loop through all of the existing constraints to remove any from previous phases
        for constraint in MO.PhaseDistanceConstraintDefinitions:
            
            # Split the string by underscores to get the prefix
            constraint_prefix = constraint.split("_")[0]
            
            # Get the journey
            constraint_journey = int(constraint_prefix.lstrip("j").split("p")[0])
            
            # Check if the constraint is from a previous journey
            if constraint_journey in to_pop:
                # It is. Do not add it to the new list
                continue            
            else:
                # We might need to add it
                
                # We need to get the phase number
                phase_number = int(constraint_prefix.split("p")[1])
                
                # Create new prefix and add it to the list
                new_phase_distance_constraints.append("j" + str(new_jdx[constraint_journey]) + "p" + str(phase_number) + constraint.split(constraint_prefix)[1])
                                                                        
        # Overwrite the new phase constraint list onto the old list
        MO.PhaseDistanceConstraintDefinitions = new_phase_distance_constraints                        

        # Create a list of the phase distance constraints that will be in the new problem
        new_maneuver_constraints = []
        
        # Loop through all of the existing constraints to remove any from previous phases
        for constraint in MO.ManeuverConstraintDefinitions:
            
            # Split the string by underscores to get the prefix
            constraint_prefix = constraint.split("_")[0]
            
            # Get the journey
            constraint_journey = int(constraint_prefix.lstrip("j").split("p")[0])
            
            # Check if the constraint is from a previous journey
            if constraint_journey in to_pop:
                # It is. Do not add it to the new list
                continue              
            else:    
                # We might need to add it
                
                # We need to get the phase number
                phase_number = int(constraint_prefix.split("b")[0].split("p")[1])
                    
                # We need to get the maneuver number
                burn_number = int(constraint_prefix.split("b")[1])
                
                    # Create new prefix and add it to the list
                new_maneuver_constraints.append("j" + str(new_jdx[constraint_journey])  + "p" + str(phase_number) + "b" + str(burn_number) + constraint.split(constraint_prefix)[1])
                                                                        
        # Overwrite the new phase constraint list onto the old list
        MO.ManeuverConstraintDefinitions = new_maneuver_constraints
            
        for entry in MO.trialX:
            # Get the transcription type
            if "FBLT" in entry[0]:
                transcription = "FBLT"
            elif "PSFB" in entry[0]:
                transcription = "PSFB"
            elif "MGALT" in entry[0]:
                transcription = "MGALT"
            elif "CoastPhase" in entry[0]:
                transcription = "CoastPhase"
                # It is! 
                
            # Get the prefix
            phase_string = entry[0].split(transcription)[0] 
            
            journey_no = int(phase_string.split("p")[0].lstrip("j"))
            
            if journey_no not in to_pop:
                if "FreePointFreeDirect" in entry[0] and MO.Journeys[new_jdx[journey_no]].departure_type != 2:
                    continue
                new_prefix = "j" + str(new_jdx[journey_no]) + "p0"
                newTrialX.append([new_prefix + transcription + entry[0].split(transcription)[1] ,entry[1]])
        MO.trialX = newTrialX   
                
        return MO      
    
    def MissedThrust(self,PEATSAorder,prefix = "",MObase = None,SObase = None):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import Mission
        import logging
        import copy
        import EMTGInterfaceReaders
        import PEATSAwaiter
        import PEATSAbox
        import PyHardware
        from numpy import floor
        
        logging.info("Creating Missed Thrust PEATSA")
        
        # Process the reference mission
        PEATSAorder = self.process_reference(PEATSAorder)
        
        # Create a waiter for splining things
        waiter = PEATSAwaiter.PEATSAwaiter()
        
        # Iniitalize the stack of cases
        boxes_out = []
                        
        # Initialize a counter for how many cases have been created
        missed_thrust_event_counter = 0
        
        # Rename the EMTG objects for simplicity
        if PEATSAorder.thin_crust:
            RefM = Mission.Mission(PEATSAorder.reference.PEATSApath + "/" + PEATSAorder.reference.PEATSAcrust_path)
        else:
            RefM = PEATSAorder.reference.PEATSAcrust      
        RefMO = PEATSAorder.reference.PEATSAdough
        RefSO = PyHardware.SpacecraftOptions(RefMO.HardwarePath.rstrip("\r\n ") + RefMO.SpacecraftOptionsFile.rstrip("\r\n "))
        
        if PEATSAorder.convert_3d_to_patched_conic:
            RefMO = self.ThreeD_Flybys_to_PatchedConic(RefMO)
        
        if MObase == None:
            MObase = RefMO
        if SObase == None:
            SObase = RefSO
        
        # Load the ephemeris data
        reference_ephem = EMTGInterfaceReaders.read_spk_data_to_list_of_lists(PEATSAorder.reference_files_location + "/" + PEATSAorder.reference_spk_data_file)
        
        # Get an index to keep track of where we are in the 
        # reference mission decision vector
        phase_start_decision_vector_index = 0
        
        # Get an index to keep track of where we are in the
        # reference mission ephemeris data
        reference_ephem_index = 0
        
        # Get an index to keep track of where the phase starts in the 
        # reference mission ephemeris data
        phase_start_ephem_index = 0
        
        # Start a counter for previous phase propellant usage
        net_ep_used = 0.0
        
        # Loop through the decision vector to get the phase flight times and number of steps
        all_phase_flight_times = {}
        for local_decision_vector_index in range(0,len(RefM.Xdescriptions)):
            
            # Get the transcription type
            if "FBLT" in RefM.Xdescriptions[local_decision_vector_index]:
                transcription = "FBLT"
            elif "PSFB" in RefM.Xdescriptions[local_decision_vector_index]:
                transcription = "PSFB"
            elif "MGALT" in RefM.Xdescriptions[local_decision_vector_index]:
                transcription = "MGALT"  
            elif "PSBI" in RefM.Xdescriptions[local_decision_vector_index]:
                transcription = "PSBI"
            elif "MGAnDSMs" in RefM.Xdescriptions[local_decision_vector_index]:
                transcription = "MGAnDSMs"
            elif "ProbeEntryPhase" in RefM.Xdescriptions[local_decision_vector_index]:
                transcription = "ProbeEntryPhase"
            elif "CoastPhase" in RefM.Xdescriptions[local_decision_vector_index]:
                transcription = "CoastPhase"
            else:
                continue
                                      
            # Check if this is a phase flight time variable
            if "phase flight time" in RefM.Xdescriptions[local_decision_vector_index]:
                # It is! 
                
                # Get the prefix
                phase_string = RefM.Xdescriptions[local_decision_vector_index].split(transcription)[0]
                
                # Set the phase flight time in the dictionary
                all_phase_flight_times.update({phase_string : RefM.DecisionVector[local_decision_vector_index]})
        
        # Loop through all of the journeys. 
        for jdx in range(len(RefM.Journeys)):
            
            logging.info("Entering Journey " + str(jdx))
        
            # Determine if this is a multi-phase Journey.
            # First, assume it only has one phase
            nPhases = 1
            # Loop through all of the sequence values in this Journey
            for seq_idx in range(0,len(RefMO.Journeys[jdx].sequence)):
                # Check if this sequence ID is non-zero, in which case there is a flyby and this is a multiphase journey
                if (RefMO.Journeys[jdx].sequence[seq_idx]):
                    # It is, update the counter
                    nPhases += 1
                    
            # If there are more than 9 phases it will break things, so throw an error
            if nPhases > 9:
                raise Exception("Phase number finding code in the decision vector re-creation cant handle more than 9 phases currently")
                                        
            # Grab the objects for this Journey and journey options
            J = RefM.Journeys[jdx]
            JO = RefMO.Journeys[jdx]
                                    
            # Set the start date for missed thrust cases
            current_epoch = J.missionevents[0].JulianDate
            
            # Loop through all of the phases
            for pdx in range(0,nPhases):      
            
                logging.info("Entering Phase " + str(pdx))  
            
                # Assume the wait time to start this phase is zero. 
                # If it isnt, we will find it and update it as we loop through the decision vector
                wait_time = 0
                
                # Assume the phase flight time is zero.
                # We will find it when searching through the decision vector, and if not, we can
                # throw an error
                phase_flight_time = 0
                
                # Assume the ep prop used is zero. 
                # If it isnt we will find it and update it as we loop through the decision vector
                reference_phase_prop_used = 0
                
                # Set the index of the control commands to -1 so we know it hasnt been found yet
                control_start_index = -1
            
                # Get the start of the phase decision variables:
                for local_decision_vector_index in range(phase_start_decision_vector_index,len(RefM.Xdescriptions)):
                    # Create a string for this phase
                    target_string = "j" + str(jdx) + "p" + str(pdx)
                    # Check if the phase prefix matches
                    if RefM.Xdescriptions[local_decision_vector_index].startswith(target_string):
                        # Cool, we now have the start of the phase's decision variables
                        phase_start_decision_vector_index = local_decision_vector_index
                        # kill the loop
                        break
                
                # Set a flag to indicate that we should be in the rhs boundary condition
                in_rhs_bc = False
                
                # Create an empty list to store the rhs boundary condition
                rhs_boundary_decision_variables = []
                
                # Initialize the transcription string
                transcription = ""
                
                # Create a flag to indicate that we have not found the rhs_bc_type
                found_rh_bc_type = False
                
                # Get the end of the phase decision variables, and indices to the variables that we will need
                for local_decision_vector_index in range(phase_start_decision_vector_index,len(RefM.Xdescriptions)):
                    
                    if "FBLT" in RefM.Xdescriptions[local_decision_vector_index]:
                        transcription = "FBLT"
                    elif "PSFB" in RefM.Xdescriptions[local_decision_vector_index]:
                        transcription = "PSFB"
                    elif "MGALT" in RefM.Xdescriptions[local_decision_vector_index]:
                        transcription = "MGALT"   
                    elif "PSBI" in RefM.Xdescriptions[local_decision_vector_index]:
                        transcription = "PSBI"  
                    elif "MGAnDSMs" in RefM.Xdescriptions[local_decision_vector_index]:
                        transcription = "MGAnDSMs"
                    elif "ProbeEntryPhase" in RefM.Xdescriptions[local_decision_vector_index]:
                        transcription = "ProbeEntryPhase"
                    elif "CoastPhase" in RefM.Xdescriptions[local_decision_vector_index]:
                        transcription = "CoastPhase"                         
                    
                    # Create a string for this phase
                    target_string = "j" + str(jdx) + "p" + str(pdx)
                    
                    # Check if we are in the rhs boundary
                    if in_rhs_bc:
                        # We are. 
                        
                        # Check if we know the bc type yet
                        if found_rh_bc_type == False:
                            # We have not. lets find it.
                            
                            # First split the decision variable by the transcription type
                            dv_split = RefM.Xdescriptions[local_decision_vector_index].split(transcription)
                            
                            # Check fo make sure the decision variable had the transcription type in it
                            if len(dv_split) == 1:
                                raise Exception("Cant find the rhs boundary condition type!")
                                
                            # Get the bc type
                            rhs_bc_type = dv_split[1].split(":")[0]
                            
                            # Update the flag
                            found_rh_bc_type = True
                                                    
                        # Check if the boundary variables are actually done
                        if rhs_bc_type not in RefM.Xdescriptions[local_decision_vector_index]:
                            in_rhs_bc = False
                        else:
                            # Add this decision variable to the list
                            rhs_boundary_decision_variables.append(["j0p0" + RefM.Xdescriptions[local_decision_vector_index].lstrip(target_string),RefM.DecisionVector[local_decision_vector_index]])
                    
                    # Check if this is the phase flight time
                    if "phase flight time" in RefM.Xdescriptions[local_decision_vector_index]:
                        phase_flight_time = RefM.DecisionVector[local_decision_vector_index]
                        in_rhs_bc = True
                        
                    # Check if this is the wait time
                    if "wait_time" in RefM.Xdescriptions[local_decision_vector_index]:
                        wait_time = RefM.DecisionVector[local_decision_vector_index]
                        
                    # Check if this is the strat of the control commands    
                    if "step 0 u_x" in RefM.Xdescriptions[local_decision_vector_index] or "Step0: substep0 u_x" in RefM.Xdescriptions[local_decision_vector_index]:
                        control_start_index = local_decision_vector_index
                        
                    # Check if this is the virtual electric propellant
                    if "virtual electric propellant" in RefM.Xdescriptions[local_decision_vector_index] and "Step" not in RefM.Xdescriptions[local_decision_vector_index] and "step" not in RefM.Xdescriptions[local_decision_vector_index]:
                        reference_phase_prop_used = RefM.DecisionVector[local_decision_vector_index]
                    
                    # Check if the phase prefix matches
                    if not RefM.Xdescriptions[local_decision_vector_index].startswith(target_string):
                        # Cool, we now have the end of the phase's decision variables
                        phase_end_decision_vector_index = local_decision_vector_index
                        # End the loop
                        break
                    # Make sure we arent at the end of the decision vector
                    elif local_decision_vector_index == len(RefM.Xdescriptions)-1:
                        # We are
                        phase_end_decision_vector_index = len(RefM.Xdescriptions)
                                                              
                # Make sure we found the phase flight time. If not, throw an error
                if phase_flight_time == 0:
                    raise Exception("Not able to find the phase flight time for " + target_string)    
                elif phase_flight_time != all_phase_flight_times[target_string]:
                    raise Exception("Two different phase flight times for phase " + target_string)                               
                
                # Get the length of the initial coast, if any
                if pdx == 0:
                    actual_departure_coast = JO.forced_initial_coast
                else:
                    actual_departure_coast = RefMO.forced_post_flyby_coast
                
                # We are going to use the departure coast for multiple things. Sometimes we want it to be zero
                departure_coast = 0
                            
                # Get the length of the arrival coast, if any
                if pdx == nPhases - 1:
                    final_coast = JO.forced_terminal_coast
                else:
                    final_coast = RefMO.forced_pre_flyby_coast 
                
                # We are going to use the departure coast for multiple things. Sometimes we want it to be zero
                actual_final_coast = final_coast
                if PEATSAorder.skip_forced_coast_periods == 0:
                    final_coast = 0
                
                # Check if the options override the default number of steps
                if RefMO.Journeys[jdx].override_num_steps:
                    # It did, grab the actual value
                    reference_timesteps = RefMO.Journeys[jdx].number_of_steps
                else:
                    # Get the reference journeys number of timesteps
                    reference_timesteps = RefMO.num_timesteps
                
                # Calculate the reference timestep
                reference_control_timestep = (phase_flight_time - actual_departure_coast - actual_final_coast) / reference_timesteps
                
                # Set the coast reduction 
                if PEATSAorder.missed_thrust_recovery_forced_terminal_coast_reduction > 0:
                    terminal_coast_reduction = PEATSAorder.missed_thrust_recovery_forced_terminal_coast_reduction * reference_control_timestep
                else:
                    terminal_coast_reduction = - PEATSAorder.missed_thrust_recovery_forced_terminal_coast_reduction
                                
                # Set the new final coast
                if actual_final_coast != 0.0:
                    new_final_coast = actual_final_coast - terminal_coast_reduction
                else:
                    new_final_coast = 0.0
                
                # Check to make sure we didnt reduce the coast too much
                if new_final_coast < 0:
                    raise Exception("Cannon reduce the forced terminal coast more than the reference terminal coast!")
                    
                # Get the arrival_epoch to the stop creating missed thrust cases
                arrival_epoch = current_epoch + wait_time + phase_flight_time - final_coast
                                    
                # Make sure we found the start of the control commands                                       
                if control_start_index == -1:
                    # This means there isnt any control. Check if we still want to create a missed thrust case
                    if PEATSAorder.skip_forced_coast_periods:
                        # no need for missed thrust in this phase. 
                        # Update the epoch and continue
                        current_epoch = arrival_epoch + final_coast
                        continue
                    else:
                        # Set a flag to indicate there is no control in this phase
                        no_control = True
                else:
                    # Set a flag to indicate there is no control in this phase
                    no_control = False
                    
                # Update the current_epoch to the start of this phase, and store that date for use with
                # initial guess creation
                phase_start_epoch = current_epoch + wait_time
                current_epoch = phase_start_epoch
                if PEATSAorder.skip_forced_coast_periods:
                    current_epoch += actual_departure_coast   
                control_start_epoch = current_epoch    
                
                # We need to find the mass at the beginning of the phase. This will be different if this is a subsequent phase or the initial
                if pdx == 0:
                    # Get the phase start mass
                    phase_start_mass = RefM.Journeys[jdx].missionevents[0].Mass
                else:
                                                         
                    # Create ephemeris splines aruond the phase start
                    phase_start_splines = waiter.CreateSplines(reference_ephem,phase_start_epoch)
  
                    # Get the phase start mass
                    phase_start_mass = phase_start_splines[6](phase_start_epoch)
                                
                # Loop through all dates of this phase
                while current_epoch < arrival_epoch:
                        
                    # Make sure we arent past an arrival date constraint if one exists            
                    if JO.timebounded == 2:
                        if current_epoch + .1 >= JO.arrival_date_bounds[1] + 2400000.5:
                            current_epoch = arrival_epoch + 1
                            continue
                                    
                    # Copy the mission object
                    MO = copy.deepcopy(MObase)
                    # SO = PyHardware.SpacecraftOptions(RefSO)
                    SO = copy.deepcopy(SObase)
                                 
                    # Create a counter for the number of stages to remove
                    stage_removal_counter = 0
                           
                    # Loop through previous journeys up to the current one
                    for journey_removal_index in range(0,jdx):
                        
                        # Check if staging occured after departure
                        if RefMO.Journeys[journey_removal_index].stage_after_departure:
                            stage_removal_counter += 1
                        
                        # Check if staging occured after arrival
                        if RefMO.Journeys[journey_removal_index].stage_after_arrival:
                            stage_removal_counter += 1
                            
                        # Remove the journey    
                        MO.Journeys.pop(0)
                        MO.number_of_journeys -= 1
                                                    
                    # Check if staging occured after departure
                    if MO.Journeys[0].stage_after_departure:
                        stage_removal_counter += 1                            
                                                                        
                    # Loop through the number of stages to remove and remove them
                    for stage_idx in range(0,stage_removal_counter):                                
                        SO.remove_stage(0)
                    
                    # Create the splines of the ephemeris around the current epoch
                    reference_splines = waiter.CreateSplines(reference_ephem,current_epoch)
                                           
                    # Turn off any forced coast
                    MO.Journeys[0].forced_post_launch_coast = 0
                    
                    # Handle mass
                                        
                    # Get the current mass at the missed thrust event date
                    current_mass = reference_splines[6](current_epoch)
                    
                    # Set the initial mass to be the current mass
                    MO.maximum_mass = current_mass
                    # And dont let that initial mass vary
                    MO.allow_initial_mass_to_vary = 0
                                        
                    # Determine how much electric propellant was used before the missed thrust event
                    # Initalize to the amount used in this phase so far
                    phase_electric_propellant_so_far = phase_start_mass-reference_splines[6](current_epoch)
                    
                    # Calculate the total ep usage up until now
                    used_electric_propellant = phase_electric_propellant_so_far + net_ep_used
                                                                                            
                    # Calculate the propellant margin in kg, up until the missed thrust event
                    if PEATSAorder.pre_missed_thrust_electric_propellant_margin >= 0:
                        # If the user options specify a margin greater than or equal to zero, use it. Otherwise, use the default
                        # value from the reference options
                        propellant_margin_used = used_electric_propellant * PEATSAorder.pre_missed_thrust_electric_propellant_margin
                    else:
                        propellant_margin_used = used_electric_propellant * RefMO.electric_propellant_margin
                                        
                    # Update the tank constraint
                    SO.setGlobalElectricPropellantTankCapacity(RefM.total_electric_propellant_including_margin - used_electric_propellant - propellant_margin_used)
                    
                    # Check if we are venting the propellant margin or not
                    if PEATSAorder.vent_pre_missed_thrust_electric_propellant_margin == 0:
                        # We are not, but it is still no longer usable, 
                        # so we need to increase the mass of any remaining stages
                        MO.final_mass_constraint_bounds[0] = RefM.spacecraft_dry_mass + propellant_margin_used
                        MO.final_mass_constraint_bounds[1] = MO.final_mass_constraint_bounds[0] + 1000.0
                    else:
                        MO.final_mass_constraint_bounds[0] = RefM.spacecraft_dry_mass
                        MO.final_mass_constraint_bounds[1] = MO.final_mass_constraint_bounds[0] + 1000.0
                        MO.maximum_mass = current_mass - propellant_margin_used                        
                    MO.constrain_dry_mass = 1
                    
                    # Determine the post-missed thrust prop margin. This variable can be a string or a float
                    if isinstance(PEATSAorder.post_missed_thrust_electric_propellant_margin,str):
                        
                        if PEATSAorder.post_missed_thrust_electric_propellant_margin == "same":
                            PEATSAorder.post_missed_thrust_electric_propellant_margin = MO.electric_propellant_margin
                        else:
                            # Remove any whitespace at the end of the string
                            PEATSAorder.post_missed_thrust_electric_propellant_margin = PEATSAorder.post_missed_thrust_electric_propellant_margin.rstrip(" ")
                            # Check if the units are percentage or kg
                            if PEATSAorder.post_missed_thrust_electric_propellant_margin.endswith("pct") or PEATSAorder.post_missed_thrust_electric_propellant_margin.endswith("%"):
                                # Units are percentage. Use the number as usual. 
                            
                                # Strip off the units and convert to a float
                                PEATSAorder.post_missed_thrust_electric_propellant_margin = float(PEATSAorder.post_missed_thrust_electric_propellant_margin.rstrip("pct% ")) / 100.0
                            elif PEATSAorder.post_missed_thrust_electric_propellant_margin.endswith("kg"):
                                # Units are in kg. You're welcome Bruno.
                            
                                # Get the value in kilograms, by stripping off the units and converting to a float
                                post_missed_thrust_electric_propellant_margin_kg = float(PEATSAorder.post_missed_thrust_electric_propellant_margin.rstrip("kg "))
                            
                                # Now calculate the margin as a percentage
                                PEATSAorder.post_missed_thrust_electric_propellant_margin = post_missed_thrust_electric_propellant_margin_kg / (MO.maximum_electric_propellant - post_missed_thrust_electric_propellant_margin_kg)
                            else:
                                try:
                                    PEATSAorder.post_missed_thrust_electric_propellant_margin = float(PEATSAorder.post_missed_thrust_electric_propellant_margin)
                                except:
                                    raise Exception("Unknown units on 'post_missed_thrust_electric_propellant_margin', use kg or pct")    
                    
                    # Set the post-missed thrust propellant margin
                    MO.electric_propellant_margin = PEATSAorder.post_missed_thrust_electric_propellant_margin
                           
                    # Create dictionaries ot store variables we will need later
                    all_remaining_phase_timesteps = {}
                    all_remaining_phase_num_timesteps = {}
                    all_remaining_nominal_terminal_coasts = {}
                    
                    # loop through all journeys and shorten the terminal coasts and shorten the mission level one also
                    for newJOjdx in range(len(MO.Journeys)):
                    
                        # Get the options for this journey
                        newJO = MO.Journeys[newJOjdx]
                        
                        # Make sure this journey has a terminal coast
                        if newJO.forced_terminal_coast != 0:
                            # Determine how we are going to reduce the terminal coasts
                            if PEATSAorder.missed_thrust_recovery_forced_terminal_coast_reduction > 0.0:        
                                                          
                                # Determine if this is a multi-phase Journey.
                                # First, assume it only has one phase
                                nPhases = 1
                                # Loop through all of the sequence values in this Journey
                                for seq_idx in range(0,len(newJO.sequence)):
                                    # Check if this sequence ID is non-zero, in which case there is a flyby and this is a multiphase journey
                                    if (newJO.sequence[seq_idx]):
                                        # It is, update the counter
                                        nPhases += 1
                        
                                # Loop through all of the phases
                                for newJOpdx in range(nPhases):
                            
                                    # Check if we can use the journey level terminal coast length
                                    if newJOpdx != nPhases - 1:
                                        # Nothing can be done about the lenght of the terminal coast. It is proscribed by
                                        # A mission level variable. 
                                        continue
                                    else:
                                
                                        # Get the length of the initial coast, if any
                                        if newJOpdx == 0:
                                            local_departure_coast = newJO.forced_initial_coast
                                        else:
                                            local_departure_coast = MO.forced_post_flyby_coast
                                    
                                        # We can fix the terminal coast exactly
                                        local_final_coast = newJO.forced_terminal_coast
                                                            
                                        # Get the prefix for this journey/phase   
                                        if newJOjdx != 0:                         
                                            target_string = "j" + str(newJOjdx+jdx) + "p" + str(newJOpdx)
                                        else:                         
                                            target_string = "j" + str(newJOjdx+jdx) + "p" + str(newJOpdx + pdx)
                                            
                                        # Get the new prefix
                                        new_prefix = "j" + str(newJOjdx) + "p" + str(newJOpdx)
                            
                                        all_remaining_nominal_terminal_coasts.update({new_prefix : local_final_coast})
                            
                                        # Check if the options override the default number of steps
                                        if newJO.override_num_steps:
                                            # It did, grab the actual value
                                            local_reference_timesteps = newJO.number_of_steps
                                        else:
                                            # Get the reference journeys number of timesteps
                                            local_reference_timesteps = MO.num_timesteps
                                        
                                        # Add the num timesteps to the dictionary
                                        all_remaining_phase_num_timesteps.update({new_prefix : local_reference_timesteps})
                                        
                                        # Calculate the reference timestep
                                        local_control_timestep = (all_phase_flight_times[target_string] - local_departure_coast - local_final_coast) / local_reference_timesteps
                                        
                                        # Add the timestep to the dictionary
                                        all_remaining_phase_timesteps.update({new_prefix : local_control_timestep})

                                        # Set the coast reduction 
                                        local_coast_reduction = PEATSAorder.missed_thrust_recovery_forced_terminal_coast_reduction * local_control_timestep
                                                                                                            
                                        # Appply the change
                                        newJO.forced_terminal_coast -= local_coast_reduction
                                        
                                        # Set this journey to override its number of steps
                                        newJO.override_num_steps = 1
                                        
                                        # Set the nubmer of timesteps
                                        import math
                                        newJO.number_of_steps = local_reference_timesteps + int(math.ceil(PEATSAorder.missed_thrust_recovery_forced_terminal_coast_reduction)) 
                            else:
                                # The exact terminal coast reduction is given, so just set it
                                local_coast_reduction = - PEATSAorder.missed_thrust_recovery_forced_terminal_coast_reduction
                    
                                # Appply the change
                                newJO.forced_terminal_coast -= local_coast_reduction
                    
                    # Reduce the mission level terminal coast
                    MO.forced_pre_flyby_coast -= terminal_coast_reduction       
                               
                    # Set the lhs boundary condition
                
                    # Set the state as a free point in space
                    MO.Journeys[0].departure_type = 2
                    # Set the departure class to be a free point
                    MO.Journeys[0].departure_class = 1
                    # Set the state as a free point in space
                    MO.Journeys[0].destination_list[0] = -1
                    # Get the state for this date
                    for state_idx in range(6):
                        MO.Journeys[0].departure_elements[state_idx] = reference_splines[state_idx](current_epoch)
                    # Get the state for this date
                    MO.Journeys[0].departure_elements_reference_epoch = current_epoch - 2400000.5
                    # Set the date of the event 
                    MO.launch_window_open_date = current_epoch - 2400000.5
                    # Allow the initial state to propagate (aka coast)
                    MO.Journeys[0].AllowJourneyFreePointDepartureToPropagate = 1    
                    # The state is fixed
                    MO.Journeys[0].departure_elements_vary_flag = [0]*6
                    # The state is cartesian
                    MO.Journeys[0].departure_elements_type = 0
                    # The state is equatorial
                    MO.Journeys[0].departure_elements_frame = 0
                    # Set the date requirement
                    MO.Journeys[0].wait_time_bounds = [0.0,PEATSAorder.peatsa_goal_objective + 1.0]
                    # Turn off the initial coast
                    MO.Journeys[0].forced_initial_coast = 0
                    # Turn off journey departure date bounds
                    MO.Journeys[0].bounded_departure_date = 0
                
                    # Now for the fun part, set up the new initial guess
                    new_decision_vector = []
                    
                    # First, handle the first journey decision vector LHS boundary condition
                    # Add the state
                    new_decision_vector.append(["j0p0" + transcription + "FreePointFreeDirectDeparture: event left state x",MO.Journeys[0].departure_elements[0]])
                    new_decision_vector.append(["j0p0" + transcription + "FreePointFreeDirectDeparture: event left state y",MO.Journeys[0].departure_elements[1]])
                    new_decision_vector.append(["j0p0" + transcription + "FreePointFreeDirectDeparture: event left state z",MO.Journeys[0].departure_elements[2]])
                    new_decision_vector.append(["j0p0" + transcription + "FreePointFreeDirectDeparture: event left state xdot",MO.Journeys[0].departure_elements[3]])
                    new_decision_vector.append(["j0p0" + transcription + "FreePointFreeDirectDeparture: event left state ydot",MO.Journeys[0].departure_elements[4]])
                    new_decision_vector.append(["j0p0" + transcription + "FreePointFreeDirectDeparture: event left state zdot",MO.Journeys[0].departure_elements[5]])
                    # Add the mass
                    new_decision_vector.append(["j0p0" + transcription + "FreePointFreeDirectDeparture: event left state mass",current_mass])
                    # Add the epoch
                    new_decision_vector.append(["j0p0" + transcription + "FreePointFreeDirectDeparture: event left state epoch",current_epoch - 2400000.5])
                    # Add the phase flight time
                    remaining_phase_flight_time = arrival_epoch - current_epoch + actual_final_coast
                                            
                    new_decision_vector.append(["j0p0" + transcription + ": phase flight time",remaining_phase_flight_time])
                    
                    # Add the rhs boundary variables
                    new_decision_vector += rhs_boundary_decision_variables
                                        
                    if not no_control:
                        # Add the chem propellant
                        new_decision_vector.append(["j0p0" + transcription + ": virtual chemical fuel",0.0])
                        # Add the ep propellant
                        new_decision_vector.append(["j0p0" + transcription + ": virtual electric propellant",reference_phase_prop_used - phase_electric_propellant_so_far])
                    
                    # Override the departure journey number of time steps
                    MO.Journeys[0].override_num_steps = 1
                                                                  
                    # Calculate the number of timesteps to use
                    MO.Journeys[0].number_of_steps = int(round((remaining_phase_flight_time - new_final_coast) / reference_control_timestep))
                    
                    # Make sure we have at least one step
                    if MO.Journeys[0].number_of_steps == 0:
                        MO.Journeys[0].number_of_steps = 1           
                                                             
                    # Calculate the control timestep
                    control_timestep = (remaining_phase_flight_time - new_final_coast) / MO.Journeys[0].number_of_steps 
                                                                                                             
                    # Loop through the control steps
                    if not no_control:
                        
                        # Set the number of variables per step in the control section
                        multiplier = 3
                        
                        # If the transcription is parallel shooting, there are actually 12, not 3 variables
                        if transcription == "PSFB":
                            multiplier = 12
                        
                        for control_index in range(0,MO.Journeys[0].number_of_steps):
                        
                            # Find the epoch of the midpoint of this control arc
                            if control_index == 0:
                                control_step_midpoint_epoch = current_epoch + .5 * control_timestep
                            else:
                                control_step_midpoint_epoch += control_timestep
                                        
                            # Figure out how long its been since the control began
                            time_since_control_start = control_step_midpoint_epoch - control_start_epoch
                    
                            # Determine how far we are fractionally through the phase
                            phase_fraction = time_since_control_start / (phase_flight_time - actual_departure_coast - actual_final_coast)
                                                
                            # Make sure we didnt overstep 
                            if phase_fraction > 1:
                                # We did. Check if thats okay
                                if PEATSAorder.missed_thrust_recovery_forced_terminal_coast_reduction and (time_since_control_start - (phase_flight_time - actual_departure_coast - actual_final_coast)) < terminal_coast_reduction:
                                    # We are past where we should be, but its okay because we are still inside of the reduced forced coast
                                    pass
                                elif control_index + 1 == MO.Journeys[0].number_of_steps and phase_fraction < 1 + reference_control_timestep / (2 * (phase_flight_time - actual_departure_coast - actual_final_coast)):
                                    # It should be okay
                                    pass
                                else:
                                    raise Exception("The closest control step is out of bounds and it shouldnt be")                                
                                                
                            # Find the closest control step
                            closest_control_step = int(floor(phase_fraction * reference_timesteps))

                            # Store the closest control step if this is the first step
                            if control_index == 0:
                                first_control_step = closest_control_step
                                                                
                            # If we are here, it should be okay if the closest control step > number of steps, we just need to not try to extract the decision vector for those
                            if closest_control_step >= reference_timesteps:
                                # We are past the previous mission's control phase, so we must be in a reduced forced terminal coast. Set control to zero and if needed
                                # get state from the splines
                                
                                # Adding control variables is different depending on transcription, so check what type is being used
                                if transcription == "FBLT":
                                    # Add the control command to the new decision vector
                                    new_decision_vector.append(["j0p0" + transcription + ": step " + str(control_index) + " u_x",0.0])
                                    new_decision_vector.append(["j0p0" + transcription + ": step " + str(control_index) + " u_y",0.0])
                                    new_decision_vector.append(["j0p0" + transcription + ": step " + str(control_index) + " u_z",0.0])
                                elif transcription == "PSFB":
                                    import math
                                    
                                    # Get the ephemeris
                                    reduced_coast_splines = waiter.CreateSplines(reference_ephem,control_step_midpoint_epoch - control_timestep / 2.0)
                                                                            
                                    # Get the cartesian state
                                    x = reduced_coast_splines[0](control_step_midpoint_epoch - control_timestep / 2.0)
                                    y = reduced_coast_splines[1](control_step_midpoint_epoch - control_timestep / 2.0)
                                    z = reduced_coast_splines[2](control_step_midpoint_epoch - control_timestep / 2.0)
                                    vx = reduced_coast_splines[3](control_step_midpoint_epoch - control_timestep / 2.0)
                                    vy = reduced_coast_splines[4](control_step_midpoint_epoch - control_timestep / 2.0)
                                    vz = reduced_coast_splines[5](control_step_midpoint_epoch - control_timestep / 2.0)
                                                                        
                                    # Convert to spherical
                                    # Page 53 http://gmat.sourceforge.net/doc/R2013a/GMATMathSpec.pdf
                                    r = math.sqrt(x*x + y*y + z*z)
                                    v = math.sqrt(vx*vx + vy*vy + vz*vz)
                                    RA = math.atan2(y,x)
                                    DEC = math.asin(z / r)
                                    FPA = math.acos((x*vx + y*vy + z*vz) / (r * v))
                                    
                                    xhat = [math.cos(DEC)*math.cos(RA),math.cos(DEC)*math.sin(RA),math.sin(DEC)]
                                    yhat = [math.cos(RA + math.pi/2),math.sin(RA + math.pi/2),0.0]
                                    zhat = [-math.sin(DEC)*math.cos(RA),-math.sin(DEC)*math.sin(RA),math.cos(DEC)]
                                    
                                    vy_prime = vx * yhat[0] + vy * yhat[1] + vz * yhat[2]
                                    vz_prime = vx * zhat[0] + vy * zhat[1] + vz * zhat[2]
                                    
                                    AZ = math.atan2(vy_prime,vz_prime)
                                                                        
                                    # Add the state
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state r",r])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state RA",RA])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state DEC",DEC])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state v",v])
                                    vRA = math.atan2(vy,vx)
                                    vDEC = math.asin(vz / v)
                                    if MO.ParallelShootingStateRepresentation == 1:
                                        new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state vRA",vRA])
                                        new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state vDEC",vDEC])
                                    elif MO.ParallelShootingStateRepresentation == 2:
                                        new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state AZ",AZ])
                                        new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state FPA",FPA])
                                    # Add the masses
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state mass",reduced_coast_splines[6](control_step_midpoint_epoch - control_timestep / 2.0)])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": virtual chemical fuel",0.0])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": virtual electric propellant",reference_phase_prop_used - phase_electric_propellant_so_far])
                                    # Add the control command to the new decision vector
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": substep0 u_x",0.0])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": substep0 u_y",0.0])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": substep0 u_z",0.0])
                                elif transcription == "MGALT":
                                    raise Exception("MGALT has not been added fully to Missed Thrust. It wont take long so just bug Jeremy about it")
                                    
                            else:
                                # This step is within the decision vector of hte initial guess, so just get those values and add them to the new decision vector
                                
                                # Find the index of this control step
                                closest_control_index = closest_control_step * multiplier + control_start_index
                            
                                # Adding control variables is different depending on transcription, so check what type is being used
                                if transcription == "FBLT":
                                    # Add the control command to the new decision vector
                                    new_decision_vector.append(["j0p0" + transcription + ": step " + str(control_index) + " u_x",RefM.DecisionVector[closest_control_index + 0]])
                                    new_decision_vector.append(["j0p0" + transcription + ": step " + str(control_index) + " u_y",RefM.DecisionVector[closest_control_index + 1]])
                                    new_decision_vector.append(["j0p0" + transcription + ": step " + str(control_index) + " u_z",RefM.DecisionVector[closest_control_index + 2]])
                                elif transcription == "PSFB":
                                    # Add the state
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state r",RefM.DecisionVector[closest_control_index - 9]])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state RA",RefM.DecisionVector[closest_control_index - 8]])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state DEC",RefM.DecisionVector[closest_control_index - 7]])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state v",RefM.DecisionVector[closest_control_index - 6]])
                                    if MO.ParallelShootingStateRepresentation == 1:
                                        new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state vRA",RefM.DecisionVector[closest_control_index - 5]])
                                        new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state vDEC",RefM.DecisionVector[closest_control_index - 4]])
                                    elif MO.ParallelShootingStateRepresentation == 2:
                                        new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state AZ",RefM.DecisionVector[closest_control_index - 5]])
                                        new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state FPA",RefM.DecisionVector[closest_control_index - 4]])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": left state mass",RefM.DecisionVector[closest_control_index - 3]])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": virtual chemical fuel",RefM.DecisionVector[closest_control_index - 2]])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": virtual electric propellant",RefM.DecisionVector[closest_control_index - 1] - phase_electric_propellant_so_far])
                                    # Add the control command to the new decision vector
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": substep0 u_x",RefM.DecisionVector[closest_control_index + 0]])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": substep0 u_y",RefM.DecisionVector[closest_control_index + 1]])
                                    new_decision_vector.append(["j0p0" + transcription + "_Step" + str(control_index) + ": substep0 u_z",RefM.DecisionVector[closest_control_index + 2]])
                                elif transcription == "MGALT":
                                    raise Exception("MGALT has not been added fully to Missed Thrust. It wont take long so just bug Jeremy about it")
                    
                    # Keep track of the epoch within the remaining decision vector
                    local_control_epoch = MO.launch_window_open_date + 2400000.5 + remaining_phase_flight_time
                    
                    # Loop through the rest of the decision vector
                    for local_decision_vector_index in range(phase_end_decision_vector_index,len(RefM.DecisionVector)):
                        # We need to update the journey and phase number prefixes
                        
                        # check if this is a time variable
                        if 'time' in RefM.Xdescriptions[local_decision_vector_index]:
                            local_control_epoch += RefM.DecisionVector[local_decision_vector_index]
                                                    
                        # Check if this is the virtual electric propellant
                        if "virtual electric propellant" in RefM.Xdescriptions[local_decision_vector_index] and "Step" not in RefM.Xdescriptions[local_decision_vector_index] and "step" not in RefM.Xdescriptions[local_decision_vector_index]:
                            local_reference_phase_prop_used = RefM.DecisionVector[local_decision_vector_index]
                        
                        # First split the nominal description by 'p' to isolate the journey number
                        nominal_description_phase_split = RefM.Xdescriptions[local_decision_vector_index].rstrip(" \r\n").split("p")
                                              
                        # Extract the journey number
                        journey_number = int(nominal_description_phase_split[0].lstrip("j"))
                        
                        # Check if the journey of this variable matches the one we are iterating through
                        if journey_number == jdx:
                            # It does match. Therefore this is a multiphase journey
                            
                            # We need to get the phase number
                            phase_number = int(nominal_description_phase_split[1][0])
                            
                            # Decrement the phase
                            phase_number -= pdx
                            
                            # Create new prefix 
                            new_prefix = "j0" + "p" + str(phase_number)
                            
                            # Add the prefix to the description
                            new_description = new_prefix + nominal_description_phase_split[1][1:]
                            
                            # Loop through the rest of the nominal description, re-adding the split p's
                            for split_string in nominal_description_phase_split[2:]:
                                new_description += "p" + split_string
                        else:
                            # It does not match, therefore we are in a subsequent journey
                            
                            # Decrement the journey number by how many were popped
                            journey_number -= jdx
                            
                            # Re-create the description
                            new_description = "j" + str(journey_number)
                            for split_string in nominal_description_phase_split[1:]:
                                new_description += "p" + split_string
                            
                            # We need to get the phase number
                            phase_number = int(nominal_description_phase_split[1][0])
                            
                            # Create the new prefix in case we need it
                            new_prefix = "j" + str(journey_number) + "p" + str(phase_number)
                                                                
                        # Add the decision variable with its updated prefix'd name
                        new_decision_vector.append([new_description,RefM.DecisionVector[local_decision_vector_index]])
                        
                        # See if we need to add a new set of control points here
                        if PEATSAorder.missed_thrust_recovery_forced_terminal_coast_reduction > 0:
                            # We might need to add poitns eventually
                                                        
                            # Check if this phase has a reduced terminal coast
                            if MO.Journeys[journey_number].forced_terminal_coast != 0:
                                # It does. still possible we need to add some control pts
                                                            
                                # Check if this is the last control point 
                                if "u_z" in RefM.Xdescriptions[local_decision_vector_index] and str(all_remaining_phase_num_timesteps[new_prefix]-1) in RefM.Xdescriptions[local_decision_vector_index].split("tep")[1]:
                                    # It is. Need to add points
                                
                                    # Get the transcription
                                    if "FBLT" in RefM.Xdescriptions[local_decision_vector_index]:
                                        transcription = "FBLT"
                                    elif "PSFB" in RefM.Xdescriptions[local_decision_vector_index]:
                                        transcription = "PSFB"
                                    elif "MGALT" in RefM.Xdescriptions[local_decision_vector_index]:
                                        transcription = "MGALT"
                                    
                                    import math                                    
                                    
                                    # Move the control epoch backwards
                                    local_control_epoch -= all_remaining_nominal_terminal_coasts[new_prefix]
                                                                                                    
                                    for control_index in range(int(math.ceil(PEATSAorder.missed_thrust_recovery_forced_terminal_coast_reduction))):
                                    
                                        # Adding control variables is different depending on transcription, so check what type is being used
                                        if transcription == "FBLT":
                                            # Add the control command to the new decision vector
                                            new_decision_vector.append([new_prefix + transcription + ": step " + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + " u_x",0.0])
                                            new_decision_vector.append([new_prefix + transcription + ": step " + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + " u_y",0.0])
                                            new_decision_vector.append([new_prefix + transcription + ": step " + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + " u_z",0.0])
                                        elif transcription == "PSFB":
                                                                    
                                            # Get the ephemeris
                                            reduced_coast_splines = waiter.CreateSplines(reference_ephem,local_control_epoch)
                                                                            
                                            # Get the cartesian state
                                            x = reduced_coast_splines[0](local_control_epoch)
                                            y = reduced_coast_splines[1](local_control_epoch)
                                            z = reduced_coast_splines[2](local_control_epoch)
                                            vx = reduced_coast_splines[3](local_control_epoch)
                                            vy = reduced_coast_splines[4](local_control_epoch)
                                            vz = reduced_coast_splines[5](local_control_epoch)
                                                                        
                                            # Convert to spherical
                                            # Page 53 http://gmat.sourceforge.net/doc/R2013a/GMATMathSpec.pdf
                                            r = math.sqrt(x*x + y*y + z*z)
                                            v = math.sqrt(vx*vx + vy*vy + vz*vz)
                                            RA = math.atan2(y,x)
                                            DEC = math.asin(z / r)
                                            FPA = math.acos((x*vx + y*vy + z*vz) / (r * v))
                                                                     
                                            vRA = math.atan2(vy,vx)
                                            vDEC = math.asin(vz/v)
                                                                                
                                            xhat = [math.cos(DEC)*math.cos(RA),math.cos(DEC)*math.sin(RA),math.sin(DEC)]
                                            yhat = [math.cos(RA + math.pi/2),math.sin(RA + math.pi/2),0.0]
                                            zhat = [-math.sin(DEC)*math.cos(RA),-math.sin(DEC)*math.sin(RA),math.cos(DEC)]
                                    
                                            vy_prime = vx * yhat[0] + vy * yhat[1] + vz * yhat[2]
                                            vz_prime = vx * zhat[0] + vy * zhat[1] + vz * zhat[2]
                                    
                                            AZ = math.atan2(vy_prime,vz_prime)
                                                                        
                                            # Add the state
                                            new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": left state r",r])
                                            new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": left state RA",RA])
                                            new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": left state DEC",DEC])
                                            new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": left state v",v])
                                            if MO.ParallelShootingStateRepresentation == 1:
                                                new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": left state vRA",vRA])
                                                new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": left state vDEC",vDEC])
                                            elif MO.ParallelShootingStateRepresentation == 2:
                                                new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": left state AZ",AZ])
                                                new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": left state FPA",FPA])
                                            # Add the masses
                                            new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": left state mass",reduced_coast_splines[6](local_control_epoch)])
                                            new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": virtual chemical fuel",0.0])
                                            new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": virtual electric propellant",local_reference_phase_prop_used])
                                            # Add the control command to the new decision vector
                                            new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": u_x",0.0])
                                            new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": u_y",0.0])
                                            new_decision_vector.append([new_prefix + transcription + "_Step" + str(control_index + all_remaining_phase_num_timesteps[new_prefix]) + ": u_z",0.0])
                                        elif transcription == "MGALT":
                                            raise Exception("MGALT has not been added fully to Missed Thrust. It wont take long so just bug Jeremy about it")
                                    
                                        # Move the control epoch forward
                                        local_control_epoch += all_remaining_phase_timesteps[new_prefix]
                                                                            
                                    # Move the control epoch forward
                                    local_control_epoch += (all_remaining_nominal_terminal_coasts[new_prefix] - PEATSAorder.missed_thrust_recovery_forced_terminal_coast_reduction * all_remaining_phase_timesteps[new_prefix])

                    MO.trialX = new_decision_vector

                    # Create a list of the phase distance constraints that will be in the new problem
                    new_boundary_constraints = []
                    
                    # Loop through all of the existing constraints to remove any from previous phases
                    for constraint in MO.BoundaryConstraintDefinitions:
                        
                        # Split the string by underscores to get the prefix
                        constraint_prefix = constraint.split("_")[0]
                        
                        # Get the journey
                        constraint_journey = int(constraint_prefix.lstrip("j").split("p")[0])
                        
                        # Check if the constraint is from a previous journey
                        if constraint_journey < jdx:
                            # It is. Do not add it to the new list
                            continue                  
                        # We do need to add it
                        else:
                            
                            # We need to get the phase number
                            phase_number = int(constraint_prefix.split("p")[1])
                            
                            # Check if its still the same journey, but the phase is from a previous phase
                            if constraint_journey == jdx and phase_number < pdx:
                                # We do not need to add it!
                                continue
                            else:
                                # We do need to add it
                                                            
                                # Check if we need to decrement the phase
                                if constraint_journey == jdx:
                                    # Decrement the phase
                                    phase_number -= pdx
                                                        
                                # Decrement the journey number by how many were popped
                                constraint_journey -= jdx
                            
                                # Create new prefix and add it to the list
                                new_boundary_constraints.append("j" + str(constraint_journey) + "p" + str(phase_number) + constraint.split(constraint_prefix)[1])
                                                                                    
                    # Overwrite the new boundary constraint list onto the old list
                    MO.BoundaryConstraintDefinitions = new_boundary_constraints
                    
                    # Create a list of the phase distance constraints that will be in the new problem
                    new_phase_distance_constraints = []
                    
                    # Loop through all of the existing constraints to remove any from previous phases
                    for constraint in MO.PhaseDistanceConstraintDefinitions:
                        
                        # Split the string by underscores to get the prefix
                        constraint_prefix = constraint.split("_")[0]
                        
                        # Get the journey
                        constraint_journey = int(constraint_prefix.lstrip("j").split("p")[0])
                        
                        # Check if the constraint is from a previous journey
                        if constraint_journey < jdx:
                            # It is. Do not add it to the new list
                            continue            
                        else:
                            # We might need to add it
                            
                            # We need to get the phase number
                            phase_number = int(constraint_prefix.split("p")[1])
                            
                            # Check if its still the same journey, but the phase is from a previous phase
                            if constraint_journey == jdx and phase_number < pdx:
                                # We do not need to add it!
                                continue
                            else:      
                                # We do need to add it
                                                            
                                # Check if we need to decrement the phase
                                if constraint_journey == jdx:
                                    # Decrement the phase
                                    phase_number -= pdx
                                                        
                                # Decrement the journey number by how many were popped
                                constraint_journey -= jdx
                            
                                # Create new prefix and add it to the list
                                new_phase_distance_constraints.append("j" + str(constraint_journey) + "p" + str(phase_number) + constraint.split(constraint_prefix)[1])
                                                                                    
                    # Overwrite the new phase constraint list onto the old list
                    MO.PhaseDistanceConstraintDefinitions = new_phase_distance_constraints                        

                    # Create a list of the phase distance constraints that will be in the new problem
                    new_maneuver_constraints = []
                    
                    # Loop through all of the existing constraints to remove any from previous phases
                    for constraint in MO.ManeuverConstraintDefinitions:
                        
                        # Split the string by underscores to get the prefix
                        constraint_prefix = constraint.split("_")[0]
                        
                        # Get the journey
                        constraint_journey = int(constraint_prefix.lstrip("j").split("p")[0])
                        
                        # Check if the constraint is from a previous journey
                        if constraint_journey < jdx:
                            # It is. Do not add it to the new list
                            continue              
                        else:    
                            # We might need to add it
                            
                            # We need to get the phase number
                            phase_number = int(constraint_prefix.split("b")[0].split("p")[1])
                            
                            # Check if its still the same journey, but the phase is from a previous phase
                            if constraint_journey == jdx and phase_number < pdx:
                                # We do not need to add it!
                                continue
                            else:
                                
                                # We need to get the maneuver number
                                burn_number = int(constraint_prefix.split("b")[1])
                                
                                # We do need to add it
                                if constraint_journey == jdx and phase_number == pdx and burn_number < first_control_step:
                                    # We do not need to add it
                                    continue
                                else:
                                    # We do need to add it
                                                      
                                    # Check if we need to decrement the burn:
                                    if phase_number == pdx:
                                        # Decrement the burn
                                        burn_number -= first_control_step
                                                                
                                    # Check if we need to decrement the phase
                                    if constraint_journey == jdx:
                                        # Decrement the phase
                                        phase_number -= pdx                             
                                                        
                                    # Decrement the journey number by how many were popped
                                    constraint_journey -= jdx
                            
                                    # Create new prefix and add it to the list
                                    new_maneuver_constraints.append("j" + str(constraint_journey) + "p" + str(phase_number) + "b" + str(burn_number) + constraint.split(constraint_prefix)[1])
                                                                                    
                    # Overwrite the new phase constraint list onto the old list
                    MO.ManeuverConstraintDefinitions = new_maneuver_constraints

                    # Set the missed thrust recovery mission name
                    MO.mission_name = prefix + 'J' + str(jdx) + '_MT'+ str(missed_thrust_event_counter)
                    
                    # Set the spacecraft file
                    MO.SpacecraftOptionsFile = MO.mission_name + '.emtg_spacecraftopt'
                    
                    # Update the peatsa hardware path
                    MO.HardwarePath = PEATSAorder.hardware_dir
                        
                    # Turn on the working directory overrides
                    MO.override_working_directory = 1
                    MO.forced_working_directory = PEATSAorder.results_directory
                    MO.override_mission_subfolder = 1
                    MO.forced_mission_subfolder = "/"
                
                    # Update the objective function type
                    MO.objective_type = 4
                    MO.objective_journey = 0
                
                    # Turn on seeding
                    MO.seed_MBH = 1
                    
                    # Turn on mbh in case its off
                    MO.run_inner_loop = 1
                    
                    # Quiet snopt and mbh
                    MO.quiet_basinhopping = 1
                    MO.quiet_NLP = 1
                    
                    # Set the max mbh run time
                    MO.MBH_max_run_time = PEATSAorder.MBH_max_run_time

                    # Apply the override options
                    for option in PEATSAorder.override_options:
                        try:
                            if "SMO" in option[0]:
                                continue
                            # Check the condition
                            if eval(option[0]):
                                if "SMO" in option[1]:
                                    continue
                
                                if option[1].startswith("MO"):
                                    varname = option[1].split("MO.")[1].split("=")[0].rstrip(" ")
                                    if "Journeys" in varname:
                                        varname = varname.split("].")[1]
                                        obj = MO.Journeys[0]
                                    else:
                                        obj = MO
                
                                    if "[" in varname:
                                        varname = varname.split("[")[0]  
                                        
                                    if varname not in dir(obj):
                                        logging.info("WARNING: Override might not work, because " + varname + " doesnt seem to be in PyEMTG Mission Options or Journey Options")
                                        logging.info("Option: " + str(option))
                        
                                # It evaluated true, so apply it
                                exec(option[1])
                        except Exception as e:
                            print("WARNING: Unable to override an option because:" + str(e))
                            logging.info("Option: " + str(option))
                            if "text" in dir(e):
                                print("*******************************************")
                                print(e.text)
                                print(' ' * e.offset + "^")
                                print("*******************************************")
                                    
                    # Write the options to file!
                    MO.write_options_file(PEATSAorder.cases_directory + "/" + MO.mission_name + '.emtgopt')
                    SO.write_output_file(PEATSAorder.hardware_dir + "/" + MO.SpacecraftOptionsFile)
                                    
                    # Create a peatsa box object for this case
                    box = PEATSAbox.PEATSAbox(PEATSAorder,PEATSAorder.cases_directory,MO.mission_name + ".emtgopt")
                
                    # Add the peatsa box object to the stack
                    boxes_out.append(box)  
                                
                    # Setup the next iteration
                                    
                    # Increment the date counter
                    if PEATSAorder.missed_thrust_timestep > 1:
                        current_epoch += PEATSAorder.missed_thrust_timestep
                    else:
                        current_epoch += PEATSAorder.missed_thrust_timestep * reference_control_timestep             
                                
                    # Cycle the counter
                    missed_thrust_event_counter += 1
                    
                # Update the current epoch to the exact end date of this phase
                current_epoch = arrival_epoch + final_coast
                
                # Increase the net propellant usage
                net_ep_used = reference_phase_prop_used
                
        return boxes_out
    
    def TradeStudy(self,PEATSAorder):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        import copy
        import PEATSAbox
        
        logging.info("Creating Trade Study PEATSA")
        
        # Initialize the output array
        boxes_out = []
        
        # Initialize a counter
        options_counter = PEATSAorder.starting_case_index
                
        # Loop through all of the options files
        for options_file_tuple in PEATSAorder.trade_study_options_files:
            
            # Check how this options file is formatted
            if options_file_tuple[1] == 0:
                # Call the type zero trade study peatsa chef
                options_list = self.TradeStudyTypeZero(PEATSAorder,options_file_tuple[0])
            elif options_file_tuple[1] == 1:
                # Call the type one trade study peatsa chef
                options_list = self.TradeStudyTypeOne(PEATSAorder,options_file_tuple[0])
            elif options_file_tuple[1] == 2:
                # Call the type one trade study peatsa chef
                options_list = self.TradeStudyTypeTwo(PEATSAorder,options_file_tuple[0])                    

            # Loop through all of the options and check if any use spacecraft options or stage options
            # Assume we dont need them
            need_SO = False
            for option in options_list:
                # Loop through all of the variables
                for variable_tuple in option:
                    if variable_tuple[0].startswith("SO."):
                        need_SO = True
                        break
                    if variable_tuple[0].startswith("StgO"):
                        need_SO = True
                        break
                        
            # Loop through all of the trade study cases
            for option in options_list:
                
                # Copy the mission options object
                MO = copy.deepcopy(PEATSAorder.reference.PEATSAdough)
        
                # If we need stage or spaceraft options, then load the spacecraft options 
                if need_SO:
                    # Import the pyhardware module
                    import PyHardware
                    SO = PyHardware.SpacecraftOptions(MO.HardwarePath + "/" + MO.SpacecraftOptionsFile)
        
                # Loop through the variables:
                for variable_tuple in option:
                    
                    # Check if this is misison data being traded
                    if variable_tuple[0].startswith("MO.") or variable_tuple[0].startswith("PEATSAorder."):
                        # It is. Easy
                        
                        # Apply the value to the variable
                        exec(variable_tuple[0] + " = " + variable_tuple[1])          
                    # Check if the variable being traded is spacecraft option data
                    elif variable_tuple[0].startswith("SO."):
                        # It is. We need to use get/set methods
                                                
                        # Apply the variable using the set function
                        exec(variable_tuple[0] + "(" + variable_tuple[1] + ")")
                    # Check if the variable being traded is a stage 
                    elif variable_tuple[0].startswith("StgO"):                 
                        # It is. We need to use get/set methods
                        
                        # Get the stage index
                        stage_idx = variable_tuple[0].split("[")[1].split("]")[0] 
                        
                        # Get the stage
                        StgO = SO.getStageOptions(int(stage_idx))
                        
                        # Get the set function name
                        set_function = variable_tuple[0].lstrip("StgO[").lstrip(stage_idx).lstrip("].")
                        
                        # Execute the options
                        exec("StgO." + set_function + "(" + variable_tuple[1] + ")")
                        
                        # Put the stage options back on the spacecraft
                        SO.setStageOptions(int(stage_idx),StgO)
                        
                    else:
                        raise Exception("Unknown option to be traded: " + variable_tuple[0])
             
                # Check if this is a missed thrust trade or regular trade
                if PEATSAorder.PEATSA_type == 2:
                 
                    # Set the trade study mission name
                    MO.mission_name = "TradeStudy_Case" + str(options_counter)
                    
                    # Turn on the working directory overrides
                    MO.override_working_directory = 1
                    MO.forced_working_directory = PEATSAorder.results_directory
                    MO.override_mission_subfolder = 1
                    MO.forced_mission_subfolder = "/"
                    
                    # Update the max run time set from the peatsa options
                    MO.MBH_max_run_time = PEATSAorder.MBH_max_run_time
                            
                    # Apply the override options
                    for option in PEATSAorder.override_options:
                        try:
                            if "SMO" in option[0]:
                                continue
                            # Check the condition
                            if eval(option[0]):
                                if "SMO" in option[1]:
                                    continue
                
                                if option[1].startswith("MO"):
                                    varname = option[1].split("MO.")[1].split("=")[0].rstrip(" ")
                                    if "Journeys" in varname:
                                        varname = varname.split("].")[1]
                                        obj = MO.Journeys[0]
                                    else:
                                        obj = MO
                
                                    if "[" in varname:
                                        varname = varname.split("[")[0]  
                                        
                                    if varname not in dir(obj):
                                        logging.info("WARNING: Override might not work, because " + varname + " doesnt seem to be in PyEMTG Mission Options or Journey Options")
                                        logging.info("Option: " + str(option))
                        
                                # It evaluated true, so apply it
                                exec(option[1])
                        except Exception as e:
                            logging.info("WARNING: Unable to override an option because:" + str(e))
                            logging.info("Option: " + str(option))
                            if "text" in dir(e):
                                logging.info("*******************************************")
                                logging.info(e.text)
                                logging.info(' ' * e.offset + "^")
                                logging.info("*******************************************")
                    
                    # Check if we need a new spacecraft option file
                    if need_SO:
                        # set the new name
                        SO.setName("TradeStudy_Case" + str(options_counter))
                        
                        # Write the options to file
                        SO.write_output_file(PEATSAorder.hardware_dir + "/" + "TradeStudy_Case" + str(options_counter) + ".emtg_spacecraftopt")
                    
                        # Update the spacecraft options file in the mission options
                        MO.HardwarePath = PEATSAorder.hardware_dir
                        MO.SpacecraftOptionsFile = "TradeStudy_Case" + str(options_counter) + ".emtg_spacecraftopt"                        
                    
                    # assemble scripted constraints
                    MO.AssembleMasterConstraintVectors()
                    
                    # Write the options to file!
                    MO.write_options_file(PEATSAorder.cases_directory + "/" + MO.mission_name + '.emtgopt')
                                
                    # Create a peatsa box object for this case
                    box = PEATSAbox.PEATSAbox(PEATSAorder,PEATSAorder.cases_directory,MO.mission_name + ".emtgopt")
            
                    # Add the peatsa box object to the stack
                    boxes_out.append(box)
            
                elif PEATSAorder.PEATSA_type == 3:
                                            
                    # Check if we need a new spacecraft option file
                    if need_SO:
                        # Write the options to file
                        SO.write_output_file(PEATSAorder.hardware_dir + "/MissedThrust_Trade" + str(options_counter) + ".emtg_spacecraftopt")

                        # Update the spacecraft options file in the mission options
                        MO.HardwarePath = PEATSAorder.hardware_dir
                        MO.SpacecraftOptionsFile = "MissedThrust_Trade" + str(options_counter) + ".emtg_spacecraftopt"
                        
                    # Write the options to a fiel that the Missed Thrust module can use
                    # MO.write_options_file(PEATSAorder.reference_files_location + "/TempReference.emtgopt")
                    
                    # Update the PEATSAorder's peatsa box to use the temporary options
                    # PEATSAorder.reference = PEATSAbox.PEATSAbox(PEATSAorder,PEATSAorder.reference_files_location,"TempReference.emtgopt",PEATSAorder.reference_mission_file_path)
                    
                    # Call the missed thrust method and add these boxes to the list
                    if need_SO:
                        boxes_out += self.MissedThrust(PEATSAorder,"TS" + str(options_counter) + "_",MO,SO)
                    else:
                        boxes_out += self.MissedThrust(PEATSAorder,"TS" + str(options_counter) + "_",MO)
                  
                # Increment the options counter
                options_counter += 1                    
        
        # If this is a missed thrust trade, put the default mission options object back where it belongs
        if PEATSAorder.PEATSA_type == 3:
            # Recreate the box
            PEATSAorder.reference = PEATSAbox.PEATSAbox(PEATSAorder,PEATSAorder.reference_files_location,PEATSAorder.reference_options_file_name,PEATSAorder.reference_mission_file_name)
                                   
        # Return the list of boxes so the oven can run them all
        return boxes_out
    
    def cleanup_line_endings(self,filename):
        # Determine the line endings of the file and strip them leaving a list of csv lines
        
        # Open file
        file_handle = open(filename)
        
        # Grab all the lines of the file
        all_lines = file_handle.readlines()
                
        # First, python might have barfed and put the whole file in one string
        if len(all_lines) == 1:        
            # It did. 
            
            # Find the line endings from within that long string
            if "\r\n" in all_lines[0]:
                newline = "\r\n"
            elif "\n" in all_lines[0]:
                newline = "\n"
            elif "\r" in all_lines[0]:
                newline = "\r"
            else:
                raise Exception("Unable to determine line endings in options file: " + filename)
                
            # Now split the file based on the new line character into a list of csv strings
            all_lines = all_lines[0].split(newline)
        elif len(all_lines) == 0:
            # The file has no lines
            raise Exception("Options file empty: " + filename)
        else:
            # Python split the lines into a list already, but those lines might have line endings. Just to be safe
            # Loop through them all and strip them
            for line_idx in range(len(all_lines)):
                all_lines[line_idx] = all_lines[line_idx].rstrip("\r\n")
                
        return all_lines
    
    def TradeStudyTypeZero(self,PEATSAorder,filename,returntype = 0):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        
        logging.info("Opening Trade Study file: " + filename)
        
        # Clean up the line endings
        all_lines = self.cleanup_line_endings(filename)
            
        # # Initialize the dictionary
        all_options = {}
 
        # Set the reference
        PEATSAorder.reference_files_location = all_lines[0].split(",")[0]
        PEATSAorder.reference_options_file_name = all_lines[1].split(",")[0]
        
        # Process the reference mission
        PEATSAorder = self.process_reference(PEATSAorder)
        
        # Split the first line:
        variables = all_lines[2].split(",")

        # Find default values
        defaults = all_lines[3].split(",")
        
        # Loop through the variable strings to create an empty list in the dictionary for each variable
        for var in variables:
            all_options.update({var : []})
                    
        # Loop through the rest of the options file lines    
        for line in all_lines[4:]:
                        
            # Split the line by comma
            linesplit = line.split(",")
            
            newlinesplit = []
            for entry in linesplit:
                entry = entry.replace("|",",")
                newlinesplit.append(entry)
            linesplit = newlinesplit
            
            # Loop through the entries in the split line, keeping in mind
            # That there might not be as many entries as there are variables
            # Because empty columns to the right will drop out eventually
            for idx in range(len(linesplit)):
                if linesplit[idx] != "":
                    # Add this option to the dictionary
                    all_options[variables[idx]].append(linesplit[idx])
                    
        # Initialize the final trade list            
        list_of_trades = []
        
        # Initialize a list of counters for each variable
        nOptions = {}
        
        # Loop through all the variables
        for var_idx in range(len(variables)):
            
            # Add a new counter
            nOptions.update({variables[var_idx]:[defaults[var_idx],0]})
            
            # Loop through all of the options for that variable
            for val in all_options[variables[var_idx]]:
                
                # update the counter
                nOptions[variables[var_idx]][1] += 1
                
                # Initialize an empty list for this set of options
                list_of_trades.append([])
                
                # Loop through the variables again, this time to add their values to the trade study list
                for var_jdx in range(len(variables)):
                    
                    # If this variable is not the one being traded, then it is going to get a default value
                    if var_idx != var_jdx:
                        
                        # It is not, give it the default value
                        list_of_trades[-1].append((variables[var_jdx],defaults[var_jdx]))
                    
                    else:
                        
                        # It is the trade variable, give it the new option
                        list_of_trades[-1].append((variables[var_jdx],val.rstrip("\r\n")))
                        
        if returntype == 0:
            # Return the full list
            return list_of_trades  
        elif returntype == 1:
            # Return the number of options for each variable
            return nOptions        
    
    def TradeStudyTypeOne(self,PEATSAorder,filename,returntype = 0):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        import itertools as it
        
        logging.info("Opening Trade Study file: " + filename)
        
        # Clean up the line endings
        all_lines = self.cleanup_line_endings(filename)
            
        # # Initialize the dictionary
        all_options = {}
 
        # Set the reference
        PEATSAorder.reference_files_location = all_lines[0].split(",")[0]
        PEATSAorder.reference_options_file_name = all_lines[1].split(",")[0]
        
        # Process the reference mission
        PEATSAorder = self.process_reference(PEATSAorder)
        
        # Split the first line:
        variables = all_lines[2].split(",")
        
        # Initialize a variable counter
        var_ctr = 0
        
        # Loop through the variable strings to create an empty list in the dictionary for each variable
        for var in variables:
            if var != "":
                var_ctr += 1
                all_options.update({var : []})
            else:
                break   
                    
        # Loop through the rest of the options file lines    
        for line in all_lines[3:]:
            # Split the line by comma
            linesplit = line.split(",")
            
            newlinesplit = []
            for entry in linesplit:
                entry = entry.replace("|",",")
                newlinesplit.append(entry)
            linesplit = newlinesplit
                
            # Loop through the entries in the split line, keeping in mind
            # That there might not be as many entries as there are variables
            # Because empty columns to the right will drop out eventually
            for idx in range(0,var_ctr):
                if linesplit[idx] != "":
                    # Add this option to the dictionary
                    all_options[variables[idx]].append(linesplit[idx])
                    
        # Initialize the final trade list            
        list_of_trades = []
        
        # Initialize a list of counters for each variable
        nOptions = {}
        
        # Loop through all the variables
        for var in variables:
            if var != "":
                    
                # Initialize an empty list for this set of options
                list_of_trades.append([])
            
                # Loop through all of the options for that variable
                for val in all_options[var]:
                                                            
                    # It is the trade variable, give it the new option
                    list_of_trades[-1].append((var,val.rstrip("\r\n")))
                
            else:
                break    
                
        # Use the python iterator class to quickly parse out all combinations of these options        
        list_of_trades = list(it.product(*list_of_trades))
             
        # Initialize a new list                   
        new_list = []           
                       
        # Loop through the list and convert the tuples to lists
        for lot in list_of_trades:    
            new_list.append(list(lot))    
        
        if returntype == 0:
            # Return the full list
            return new_list 
        elif returntype == 1:
            # Return the number of options for each variable
            return nOptions       
    
    def TradeStudyTypeTwo(self,PEATSAorder,filename):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        import itertools as it
        
        logging.info("Opening Trade Study file: " + filename)
        
        # Clean up the line endings
        all_lines = self.cleanup_line_endings(filename)
            
        # Set the reference
        PEATSAorder.reference_files_location = all_lines[0].split(",")[0]
        PEATSAorder.reference_options_file_name = all_lines[1].split(",")[0]
        
        # Process the reference mission
        PEATSAorder = self.process_reference(PEATSAorder)
        
        # Split the first line:
        variables = all_lines[2].split(",")
                            
        # Initialize the final trade list            
        list_of_trades = []
                    
        # Loop through the rest of the options file lines    
        for line in all_lines[3:]:
            # Split the line by comma
            linesplit = line.split(",")
            
            newlinesplit = []
            for entry in linesplit:
                entry = entry.replace("|",",")
                newlinesplit.append(entry)
            linesplit = newlinesplit
                
            # Initialize a new trade study option
            list_of_trades.append([])
            
            # Loop through the entries in the split line
            for idx in range(len(linesplit)):
                # Get the value
                val = linesplit[idx]
                        
                # Add this variable/value to the trade option      
                list_of_trades[-1].append((variables[idx],val.rstrip("\r\n")))
                                    
        # Return the full list
        return list_of_trades                
            
    def Custom(self):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        import importlib
                
        logging.info("Creating Custom PEATSA")
        
        # Import the custom script
        custom_module = importlib.import_module(PEATSAorder.custom_PEATSA_chef_path)
        
        # Calling the custom PEATSA chef
        boxes_out = custom_module.CustomPEATSAchef(PEATSAorder)
                
        logging.info("Custom PEATSAs successfully generated: " + str(len(boxes_out)))
        
        # Return the generated boxes
        return boxes_out
        
