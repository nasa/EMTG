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

#PEATSAwaiter class
#script to analyze the PEATSAboxes (determine if they need to go back in the oven) and 
#if so, finds them a new initial guess
#Used to be step 4
#Jeremy Knittel 1-24-2018

class PEATSAwaiter(object):
    def propulate(self,PEATSAorder,all_boxes):
        import os
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG/SimpleMonteCarlo/")
        import MissionOptions
        import maneuverExecutionError
        import logging
        
        if PEATSAorder.thin_crust:
            raise Exception("Propulator code not set to work with thin crust peatsas yet. Sorry.")
        
        # Create a list to hold the cases 2 run
        boxes_out = []
        propulated_boxes = []
        box1ups = []       
        total_time_to_propagate = [] 
        duty_cycles = []
        propagated_times = []
        og_journeys = []
        unit_thrust_dir = []

        prop_in_file = open(PEATSAorder.cases_directory + "/PropulationIn0.csv", 'w')
        prop_in_file.write('#JD, journey, duty cycle, step size [sec], x [km/s], y [km/s], z [km/s], xdot [km/s], ydot [km/s], zdot [km/s], mass [kg], [thrust_x, thrust_y, thrust_z,] propagation time [sec]\n')
                
        # Loop through all of the samples
        for nSample in range(PEATSAorder.nSamples):
            
            # Get all of the cases that have this index
            boxes = [box for box in all_boxes if box.PEATSAdough.user_data["MonteCarloSampleIndex"] == nSample]
            
            # Sort the cases by their depth into the mission
            boxes = sorted(boxes,key = lambda box: box.PEATSAdough.user_data["MonteCarloDepthIndex"])
            
            # Initialize an index to double check we dont have duplicates
            last_index = -1
            
            # Loop through each box with the relevant index
            for box in boxes:
                
                # Make sure that this box doesnt have the same index as the previous
                if box.PEATSAdough.user_data["MonteCarloDepthIndex"] == last_index:
                    raise Exception("Multiple versions of sample index: " + str(nSample))
                else:
                    last_index = box.PEATSAdough.user_data["MonteCarloDepthIndex"]
                
                # Check what the peatsa objective type is
                if PEATSAorder.objective_type == 0 and PEATSAorder.max_or_min == "max":
                    # The peatsa objective type is simply to have an objective
                    # function value greater than the target objective value
                
                    # Check if the objective value is greater or not
                    if box.objective_value < PEATSAorder.peatsa_goal_objective:
                        # It is less, so we need to add it to the outgoing boxes
                        logging.info("Needs to be re-run: " + box.mission_name)
                    else:
                        logging.info("Completed: " + box.mission_name)
                        continue
                elif PEATSAorder.objective_type == 0 and PEATSAorder.max_or_min == "min":
                    # The peatsa objective type is simply to have an objective
                    # function value less than the target objective value
                
                    # Check if the objective value is greater or not
                    if box.objective_value > PEATSAorder.peatsa_goal_objective:
                        # It is greater, so we need to add it to the outgoing boxes
                        logging.info("Needs to be re-run: " + box.mission_name)
                    else:
                        logging.info("Completed: " + box.mission_name)      
                        continue  
                
                # If we are here, then this case needs to be re-run
                
                # Check if we need to propulate it's initial conditions
                if box.PEATSAdough.user_data["InitialConditionsSet"]:
                    # We dont
                    boxes_out.append(box)
                else:
                
                    # First we need to find the case that has the burn that error will be applied to
                    box2up = None
                    for box2 in boxes:
                        if box.PEATSAdough.user_data["MonteCarloDepthIndex"] == box2.PEATSAdough.user_data["MonteCarloDepthIndex"]:
                            raise Exception("We've gone too far.")
                        elif box.PEATSAdough.user_data["MonteCarloDepthIndex"] - 2 == box2.PEATSAdough.user_data["MonteCarloDepthIndex"]:
                            box2up = box2
                            break
                    if box2up == None:
                        raise Exception("Box2Up not found")
                    
                    
                    # Second we need to find the case that has the control that will be simply propagated
                    box1up = None
                    for box1 in boxes:
                        if box.PEATSAdough.user_data["MonteCarloDepthIndex"] == box1.PEATSAdough.user_data["MonteCarloDepthIndex"]:
                            raise Exception("We've gone too far.")
                        elif box.PEATSAdough.user_data["MonteCarloDepthIndex"] - 1 == box1.PEATSAdough.user_data["MonteCarloDepthIndex"]:
                            box1up = box1
                            break
                    if box1up == None:
                        raise Exception("Box1Up not found")

                    
                    nominal_thrust = box2up.PEATSAcrust.Journeys[0].missionevents[1].Thrust
                    summation = 0
                    for entry in nominal_thrust:
                        summation += entry * entry
                    import math
                    norm_thrust = math.sqrt(summation)
                    unit_thrust = []
                    for entry in nominal_thrust:
                         unit_thrust.append(entry / norm_thrust)
                    
                    control_with_error_obj = maneuverExecutionError.maneuverExecutionError(unit_thrust,PEATSAorder.thrusting_standard_deviations[0],0.0,PEATSAorder.thrusting_standard_deviations[1],0.0,PEATSAorder.thrusting_standard_deviations[2])

                    # Write the input propulator file    
                    prop_in_file.write(str(box2up.PEATSAdough.launch_window_open_date + 2400000.5) + ",")
                    prop_in_file.write(str(box2up.PEATSAdough.user_data["MonteCarloOriginalJourney"]) + ",")
                    if PEATSAorder.override_propulated_duty_cycle != 0.0:
                        duty_cycles.append(PEATSAorder.override_propulated_duty_cycle)
                    else:
                        if box2up.PEATSAdough.Journeys[0].override_duty_cycle:
                            duty_cycles.append(box2up.PEATSAdough.Journeys[0].duty_cycle)
                        else:
                            duty_cycles.append(box2up.PEATSAdough.engine_duty_cycle)
                    prop_in_file.write(str(1.0) + ",")       
                    prop_in_file.write("86400,")
                    for state in box2up.PEATSAdough.user_data["InitialStateWithAllErrors"]:
                        prop_in_file.write(str(state) + ",")
                    for control in control_with_error_obj.uRandomInertial:
                        prop_in_file.write(str(control) + ",")
                    prop_in_file.write(str(86400.0) + "\n")
                
                    propulated_boxes.append(box)
                    box1ups.append(box1up)
                    total_time_to_propagate.append(box2up.PEATSAcrust.Journeys[0].missionevents[1].TimestepLength*86400.0/duty_cycles[-1])
                    propagated_times.append(86400.0)
                    og_journeys.append(box2up.PEATSAdough.user_data["MonteCarloOriginalJourney"])
                    unit_thrust_dir.append(unit_thrust)
                
                # Now we need to break, because we dont want to go any deeper into the monte carlo until the current case is complete
                break
        
        if len(box1ups): 
            
            max_time_to_propagate = 86400.0

            ctr = 0

            notdone = True
            
            lastbox1ups = box1ups

            while notdone:

                prop_in_file.close()

                # Run the propulator
                os.system(PEATSAorder.emtg_root_directory + "/PropulatorDriver" + ' ' + PEATSAorder.reference_files_location + "/" + PEATSAorder.reference_options_file_name + " " + PEATSAorder.cases_directory + "/PropulationIn" + str(ctr) + ".csv " + ' ' + PEATSAorder.cases_directory + "/PropulationOut" + str(ctr) + ".csv > /dev/null\n")
                                
                import time
                time.sleep(.5)
                  
                prop_out_file_lines = open(PEATSAorder.cases_directory + "/PropulationOut" + str(ctr) + ".csv", 'r').readlines()
                
                ctr += 1

                line_start = 0            
                for line in prop_out_file_lines:
                    line_start += 1
                    if line.startswith("JD"):
                        break
                
                notdone = False
                
                newbox1ups = []
                
                prop_in_file = open(PEATSAorder.cases_directory + "/PropulationIn" + str(ctr) +".csv", 'w')
                prop_in_file.write('#JD, journey, duty cycle, step size [sec], x [km/s], y [km/s], z [km/s], xdot [km/s], ydot [km/s], zdot [km/s], mass [kg], [thrust_x, thrust_y, thrust_z,] propagation time [sec]\n')
                
                for idx,box1up in enumerate(lastbox1ups):
                    
                    linesplit = prop_out_file_lines[idx+line_start].split(",")
                    out_state = []
                    for state_idx in range(8):
                        out_state.append(float(linesplit[state_idx]))
            
                    if total_time_to_propagate[idx] - propagated_times[idx] > max_time_to_propagate:
                        time_to_propagate = max_time_to_propagate
                    elif total_time_to_propagate[idx] > propagated_times[idx] :
                        time_to_propagate = total_time_to_propagate[idx] - propagated_times[idx]
                    else:
                        box1up.PEATSAdough.user_data["InitialStateWithAllErrors"] = out_state[1:]
                        box1up.PEATSAdough.write_options_file(box1up.PEATSApath + "/" + box1up.PEATSAdough_path)
                        continue
                    
                    notdone = True
                    
                    coast_time = total_time_to_propagate[idx] * (1-duty_cycles[idx])
                    coast_start = total_time_to_propagate[idx] - coast_time

                    if propagated_times[idx] + time_to_propagate < coast_start:
                        duty_cycle = 1.0
                    elif propagated_times[idx] >= coast_start:
                        duty_cycle = 0.0
                    else:
                        duty_cycle = 1.0
                        time_to_propagate = coast_start - propagated_times[idx]       
                                 
                    control_with_error_obj = maneuverExecutionError.maneuverExecutionError(unit_thrust_dir[idx],PEATSAorder.thrusting_standard_deviations[0],0.0,PEATSAorder.thrusting_standard_deviations[1],0.0,PEATSAorder.thrusting_standard_deviations[2])
                    
                    prop_in_file.write(str(out_state[0]) + ",")
                    prop_in_file.write(str(og_journeys[idx]) + ",")
                    prop_in_file.write(str(duty_cycle) + ",")       
                    prop_in_file.write("86400,")
                    for state in out_state[1:]:
                        prop_in_file.write(str(state) + ",")
                    for control in control_with_error_obj.uRandomInertial:
                        prop_in_file.write(str(control) + ",")
                    prop_in_file.write(str(time_to_propagate) + "\n")
                    
                    propagated_times[idx] += time_to_propagate
                    newbox1ups.append(box1up)
                lastbox1ups = newbox1ups
                
            prop_in_file_final = open(PEATSAorder.cases_directory + "/PropulationInFinal.csv", 'w')
            prop_in_file_final.write('#JD, journey, duty cycle, step size [sec], x [km/s], y [km/s], z [km/s], xdot [km/s], ydot [km/s], zdot [km/s], mass [kg], [thrust_x, thrust_y, thrust_z,] propagation time [sec]\n')
            
            for idx,box1up in enumerate(box1ups):      
                
                nominal_thrust = box1up.PEATSAcrust.Journeys[0].missionevents[1].Thrust
                summation = 0
                for entry in nominal_thrust:
                    summation += entry * entry
                import math
                norm_thrust = math.sqrt(summation)
                unit_thrust = []
                for entry in nominal_thrust:
                     unit_thrust.append(entry / norm_thrust)
        
                # Write the input propulator file    
                prop_in_file_final.write(str(box1up.PEATSAdough.launch_window_open_date + 2400000.5) + ",")
                prop_in_file_final.write(str(box1up.PEATSAdough.user_data["MonteCarloOriginalJourney"]) + ",")

                if PEATSAorder.override_propulated_duty_cycle != 0.0:
                    dc = PEATSAorder.override_propulated_duty_cycle
                else:
                    if box1up.PEATSAdough.Journeys[0].override_duty_cycle:
                        dc = box1up.PEATSAdough.Journeys[0].duty_cycle
                    else:
                        dc = box1up.PEATSAdough.engine_duty_cycle 
                prop_in_file_final.write("1.0,")
                prop_in_file_final.write("86400,")
                for state in box1up.PEATSAdough.user_data["InitialStateWithAllErrors"]:
                    prop_in_file_final.write(str(state) + ",")
                for control in unit_thrust:
                    prop_in_file_final.write(str(control) + ",")
                prop_in_file_final.write(str(box1up.PEATSAcrust.Journeys[0].missionevents[1].TimestepLength*86400.0) + "\n")
            
            prop_in_file_final.close()
         
            # Run the propulator
            os.system(PEATSAorder.emtg_root_directory + "/PropulatorDriver" + ' ' + PEATSAorder.reference_files_location + "/" + PEATSAorder.reference_options_file_name + " " + PEATSAorder.cases_directory + "/PropulationInFinal.csv " + ' ' + PEATSAorder.cases_directory + "/PropulationOutFinal.csv > /dev/null\n")
        
            import time
            time.sleep(5)
            
            prop_in_file_final2 = open(PEATSAorder.cases_directory + "/PropulationInFinal2.csv", 'w')
            prop_in_file_final2.write('#JD, journey, duty cycle, step size [sec], x [km/s], y [km/s], z [km/s], xdot [km/s], ydot [km/s], zdot [km/s], mass [kg], [thrust_x, thrust_y, thrust_z,] propagation time [sec]\n')
            
            prop_out_file_lines = open(PEATSAorder.cases_directory + "/PropulationOutFinal.csv", 'r').readlines()
                                           
            line_start = 0            
            for line in prop_out_file_lines:
                line_start += 1
                if line.startswith("JD"):
                    break
                    
            for idx,box1up in enumerate(box1ups):      
                
                linesplit = prop_out_file_lines[idx+line_start].split(",")
             
                out_state = []
                for state_idx in range(7):
                    out_state.append(float(linesplit[1+state_idx]))
                        
                # Write the input propulator file    
                prop_in_file_final2.write(linesplit[0] + ",")
                prop_in_file_final2.write(str(box1up.PEATSAdough.user_data["MonteCarloOriginalJourney"]) + ",")
                if PEATSAorder.override_propulated_duty_cycle != 0.0:
                    dc = PEATSAorder.override_propulated_duty_cycle
                else:
                    if box1up.PEATSAdough.Journeys[0].override_duty_cycle:
                        dc = box1up.PEATSAdough.Journeys[0].duty_cycle
                    else:
                        dc = box1up.PEATSAdough.engine_duty_cycle 
                prop_in_file_final2.write("0.0,")
                prop_in_file_final2.write("86400,")
                for state in out_state:
                    prop_in_file_final2.write(str(state) + ",")
                prop_in_file_final2.write("0.0,0.0,0.0,")
                prop_in_file_final2.write(str(box1up.PEATSAcrust.Journeys[0].missionevents[1].TimestepLength*86400.0*(1/dc-1)) + "\n")
            
            prop_in_file_final2.close()
         
            # Run the propulator
            os.system(PEATSAorder.emtg_root_directory + "/PropulatorDriver" + ' ' + PEATSAorder.reference_files_location + "/" + PEATSAorder.reference_options_file_name + " " + PEATSAorder.cases_directory + "/PropulationInFinal2.csv " + ' ' + PEATSAorder.cases_directory + "/PropulationOutFinal2.csv > /dev/null\n")
        
            import time
            time.sleep(5)
    
            prop_out_file_lines2 = open(PEATSAorder.cases_directory + "/PropulationOutFinal2.csv", 'r').readlines()
                                           
            line_start = 0            
            for line in prop_out_file_lines2:
                line_start += 1
                if line.startswith("JD"):
                    break
                            
            for idx,box in enumerate(propulated_boxes):
                box1up = box1ups[idx]
            
                linesplit = prop_out_file_lines2[idx+line_start].split(",")
             
                out_state = []
                for state_idx in range(7):
                    out_state.append(float(linesplit[1+state_idx]))
                box.PEATSAdough.Journeys[0].departure_elements = out_state[:6]
                box.PEATSAdough.maximum_mass = out_state[6]
                box.PEATSAdough.user_data["InitialConditionsSet"] = 1
            
                box.PEATSAdough.trialX = []
            
                add_all = False
            
                virtual_ep_ctr = 0         
                              
                for idx,entry in enumerate(box1up.PEATSAcrust.DecisionVector):
                    descrip = box1up.PEATSAcrust.Xdescriptions[idx]
                
                    jdx2 = int(descrip.split("p")[0].lstrip("j"))
                    if jdx2 > 0:
                        if box.PEATSAdough.user_data["MonteCarloOriginalJourney"] == box1up.PEATSAdough.user_data["MonteCarloOriginalJourney"]:
                            box.PEATSAdough.trialX.append(["j" + str(jdx2) + descrip.lstrip("j0123456789"),entry])
                        else:
                            raise Exception("Havent figured out how to do this yet")
                    elif jdx2 == 0:
                        if box.PEATSAdough.user_data["MonteCarloOriginalJourney"] != box1up.PEATSAdough.user_data["MonteCarloOriginalJourney"]:
                            raise Exception("Havent figured out how to do this yet")
                        if "phase flight time" in descrip:
                            add_all = True
                            box.PEATSAdough.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state x",box.PEATSAdough.Journeys[0].departure_elements[0]])
                            box.PEATSAdough.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state y",box.PEATSAdough.Journeys[0].departure_elements[1]])
                            box.PEATSAdough.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state z",box.PEATSAdough.Journeys[0].departure_elements[2]])
                            box.PEATSAdough.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state xdot",box.PEATSAdough.Journeys[0].departure_elements[3]])
                            box.PEATSAdough.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state ydot",box.PEATSAdough.Journeys[0].departure_elements[4]])
                            box.PEATSAdough.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state zdot",box.PEATSAdough.Journeys[0].departure_elements[5]])
                            # Add the mass
                            box.PEATSAdough.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state mass",box.PEATSAdough.maximum_mass])
                            # Add the epoch
                            box.PEATSAdough.trialX.append(["j0p0PSFBFreePointFreeDirectDeparture: event left state epoch",box.PEATSAdough.launch_window_open_date])
                            # Add the flight time                                            
                            box.PEATSAdough.trialX.append(["j0p0PSFB: phase flight time",box1up.PEATSAcrust.Journeys[0].missionevents[-1].JulianDate-2400000.5-box.PEATSAdough.launch_window_open_date])
                                           
                        elif "PSFB: virtual chemical" in descrip and 'Step' not in descrip:
                            add_all = False
                            # Add the chem propellant
                            box.PEATSAdough.trialX.append(["j0p0PSFB: virtual chemical fuel",0.0])
                        elif "PSFB: virtual electric" in descrip and 'Step' not in descrip:
                            # Add the ep propellant
                            box.PEATSAdough.trialX.append(["j0p0PSFB: virtual electric propellant",entry])
                            virtual_ep_idx = len(box.PEATSAdough.trialX)-1
                        elif add_all:
                            box.PEATSAdough.trialX.append(["j0" + descrip.lstrip("j0123456789"),entry])                                       
                        elif "Step" in descrip:
                            step2 = int(descrip.split("_Step")[1].split(":")[0])
                            if step2 >= 1:
                                box.PEATSAdough.trialX.append(["j0p0PSFB_Step" + str(step2-1) + ":" + descrip.split(":")[1],entry])
                                if "virtual electric" in descrip:
                                    if step2 == 1:
                                        virtual_ep_ctr = entry 
                                        box.PEATSAdough.trialX[-1][1] -= virtual_ep_ctr
                                    else:
                                        box.PEATSAdough.trialX[-1][1] -= virtual_ep_ctr    
                    
                box.PEATSAdough.trialX[virtual_ep_idx][1] -= virtual_ep_ctr
                
                if PEATSAorder.override_propulated_duty_cycle != 0.0:
                    MO.Journeys[0].ManeuverConstraintDefinitions.append('p0b0_dutycycle_' + str(PEATSAorder.override_propulated_duty_cycle))
                    MO.AssembleMasterConstraintVectors()
                
                boxes_out.append(box)
        
        return boxes_out, all_boxes

    def check_only_run_if(self,PEATSAorder,PEATSAboxes_in,override = False):                
        # Check if non-optimality criteria is being applied to determine if cases
        # should be re-run
        if len(PEATSAorder.only_rerun_if):
        
            # Initialize the list of peatsa boxes (cases) that need to be re-run
            boxes_out = []
            
            # Loop through all of the cases
            for box in PEATSAboxes_in:
            
                # Assume we will not re-run it
                re_run = False

                # Grab the mission object so that it can be used in evaluation
                if hasattr(box,"PEATSAcrust"):
                    M = box.PEATSAcrust
            
                MO = box.PEATSAdough
            
                # Loop through all criteria.
                for condition in PEATSAorder.only_rerun_if:
                                                    
                    # Check condition
                    try:
                        if eval(condition):
                            # It evaluated true, change the flag
                            re_run = True
                    except:
                        # Condition cant be evaluated, so just re-run
                        re_run = True
            
                # If nothing has overriden the re-run = false condition, then go to the next box
                if re_run == False:
                    continue
                else: 
                    boxes_out.append(box)
        
            if len(PEATSAboxes_in) > 1:
                return boxes_out
            else:
                if re_run == False:
                    return False
                else:
                    return True      
        else:
            if len(PEATSAboxes_in) > 1 or override:
                return PEATSAboxes_in
            else:
                return True     
    
    def check_optimality(self,PEATSAorder,PEATSAboxes_in):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        logging.info("Checking which cases are complete")
        
        # Initialize the list of peatsa boxes (cases) that need to be re-run
        boxes_out = []
                
        # Loop through all of the cases
        for box in PEATSAboxes_in:
        
            re_run = self.check_only_run_if(PEATSAorder,[box])
            
            if re_run == False:
                continue
            
            # Check what the peatsa objective type is
            if PEATSAorder.objective_type == 0 and PEATSAorder.max_or_min == "max":
                # The peatsa objective type is simply to have an objective
                # function value greater than the target objective value
                
                # Check if the objective value is greater or not
                if box.objective_value < PEATSAorder.peatsa_goal_objective:
                    # It is less, so we need to add it to the outgoing boxes
                    boxes_out.append(box)
                    logging.info("Needs to be re-run: " + box.mission_name)
                else:
                    logging.info("Completed: " + box.mission_name)
            elif PEATSAorder.objective_type == 0 and PEATSAorder.max_or_min == "min":
                # The peatsa objective type is simply to have an objective
                # function value less than the target objective value
                
                # Check if the objective value is greater or not
                if box.objective_value > PEATSAorder.peatsa_goal_objective:
                    # It is greater, so we need to add it to the outgoing boxes
                    boxes_out.append(box)
                    logging.info("Needs to be re-run: " + box.mission_name)
                else:
                    logging.info("Completed: " + box.mission_name)
            
            elif PEATSAorder.objective_type == 1:
                # We need to make a polynomial fit of the seed cases and then
                # check how well the current case does by comparison
                
                # First, find the seed cases
                box = self.find_seed_cases(PEATSAorder,box,PEATSAboxes_in)
                
                # Assume it needs to be re-run
                needs_to_be_rerun = False
                
                # Loop through the different seed criteria
                for seed_criteria in PEATSAorder.seed_criteria:
                    
                    # Check if there are enough cases to satisfy the 
                    # needs of the polynomial order requested
                    if len(box.seed_boxes[seed_criteria[0]]) > PEATSAorder.polyfit_order:
                    
                        # Polyfit is a member of numpy, so import it
                        import numpy
                    
                        # Get the arrays of data that make the polyfit
                        xData = [case[1] for case in box.seed_boxes[seed_criteria[0]]]
                        yData = [case[2] for case in box.seed_boxes[seed_criteria[0]]]
                        minVal = min(yData)
                        yData = [y - minVal for y in yData]
                    
                        # Get the coefficients of the polynomial fit
                        coeffs = numpy.polyfit(xData,yData,PEATSAorder.polyfit_order)
                    
                        # Create the polynomial function object
                        poly = numpy.poly1d(coeffs)
                    
                        # Check if this is a minimization or maximization problem
                        if PEATSAorder.max_or_min == "max":
                            # It is maximization
                        
                            # Check if the case exceeds the polynomial fit's expected value
                            if box.objective_value - minVal < poly(box.seed_criteria[seed_criteria[0]]) * (1.0 - PEATSAorder.polyfit_margin):
                                # It is not, so we can stop looping
                                needs_to_be_rerun = True
                                
                                box.re_run_using_this_seed[seed_criteria[0]] = True
                            else:
                                box.re_run_using_this_seed[seed_criteria[0]] = False
                                # It is, so check th enext criteria
                        elif PEATSAorder.max_or_min == "min":
                            # It is minimization
                        
                            # Check if the case is less than tshe polynomial fit's expected value
                            if box.objective_value - minVal > poly(box.seed_criteria[seed_criteria[0]]) * (1.0 + PEATSAorder.polyfit_margin):
                                # It is not, so we can stop looping
                                needs_to_be_rerun = True
                                box.re_run_using_this_seed[seed_criteria[0]] = True
                            else:
                                box.re_run_using_this_seed[seed_criteria[0]] = False
                                # It is, so check the next criteria
                    else:
                        box.re_run_using_this_seed[seed_criteria[0]] = True
                        needs_to_be_rerun = True
                        # Not enough seed_boxes, assume it must be re-run
                        logging.info("Needs to be re-run: " + box.mission_name)
                
                if needs_to_be_rerun:
                    boxes_out.append(box)
                    logging.info("Needs to be re-run: " + box.mission_name)
                else:
                    logging.info("Completed: " + box.mission_name)
                    
                
        # All done, return the boxes
        return boxes_out
        
    def find_all_seed_cases(self,PEATSAorder,PEATSAboxes2run,PEATSAboxes,external_seed_boxes):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        logging.info("Finding all seed cases for all cases")
        
        # Initialize the output list
        boxes_out = []
        
        # Loop through all cases
        for box in PEATSAboxes2run:
            
            # Find the seed cases for this case
            out_box = self.find_seed_cases(PEATSAorder,box,PEATSAboxes + external_seed_boxes)
            
            # append the returned box. The returned box is the same, except that the seed cases have been set
            boxes_out.append(out_box)
            
        # Return the full list
        return boxes_out
        
    def find_seed_cases(self,PEATSAorder,PEATSAbox,PEATSAboxes,include_infeasible = False,solo_criteria = None):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        logging.info("Finding all seeed cases for: " + PEATSAbox.mission_name)
        
                
        if solo_criteria != None:
            loop_criteria = [solo_criteria]
        else:
            # Re-set the list of seeds
            PEATSAbox.seed_boxes = {}
        
            loop_criteria = PEATSAorder.seed_criteria       
                
        # Loop through all the seed criteria
        for seed_criteria in loop_criteria:
            
            if seed_criteria[0] not in PEATSAbox.seed_boxes.keys():
                # Add a new empty list for each seed criteria
                PEATSAbox.seed_boxes.update({seed_criteria[0]:[]})
        
            # Get the target fingerprint
            target_fingerprint = PEATSAbox.fingerprint[seed_criteria[0]][:-1]

            # The boxes should be sorted, so once we find the section of boxes that match, we dont have to keep looking.
            # So, create a flag to indicate when that section of boxes is found
            found_matching_section = 0
        
            # Loop through all the other cases
            for possible_seed in PEATSAboxes:
                    
                # Get the possible seeds fingerprint
                possible_seed_fingerprint = possible_seed.fingerprint[seed_criteria[0]][:-1]

                # Check if the cases have the same fingerprint, excluding the seed criteria of note
                if target_fingerprint == possible_seed_fingerprint:
                
                    # Update the flag
                    found_matching_section = 1
                
                    # Check if the seed is to far away from the current case or if we care:
                    if seed_criteria[1] <= 0 and possible_seed.seed_criteria[seed_criteria[0]] - PEATSAbox.seed_criteria[seed_criteria[0]] < seed_criteria[1]:

                         # Too far, try the next
                         continue
                         
                    if seed_criteria[2] >= 0 and possible_seed.seed_criteria[seed_criteria[0]] - PEATSAbox.seed_criteria[seed_criteria[0]] > seed_criteria[2]:

                         # Too far, try the next
                         continue
                
                    # Check if we are allowe to seed ourselves
                    if abs(PEATSAbox.seed_criteria[seed_criteria[0]] - possible_seed.seed_criteria[seed_criteria[0]]) == 0 and possible_seed.in_study == True and PEATSAorder.allow_cases_to_seed_themselves == 0:
                        # we are not and this the same case
                        continue
                
                    # Check if the potential seed converged, and can therefore be used as a seed, or if we dont care
                    if possible_seed.converged == 1 or include_infeasible:
                        # It did, or we dont care
                    
                        # Make sure that if this case hasnt converged, that it it actually has an infeasible .emtg file written out
                        if include_infeasible and "fake" in possible_seed.PEATSAcrust_path:
                            # It doesnt.
                            continue
                    
                        # Check if we need to determine if the seed case met the objective
                        if PEATSAorder.seed_from_cases_that_havent_met_target == 0:
                            # We do
                        
                            # Check if this is a minimization or maximizatino problem
                            if PEATSAorder.max_or_min == "max":
                                # It is a maximization
                            
                                # Check if the case met the goal criteria
                                if possible_seed.objective_value < PEATSAorder.peatsa_goal_objective:
                                    
                                    # It did not. Continue to the next potential seed
                                    continue
                            elif PEATSAorder.max_or_min == "min":
                                # It is a maximization
                            
                                # Check if the case met the goal criteria
                                if possible_seed.objective_value > PEATSAorder.peatsa_goal_objective:
                                    
                                    # It did not. Continue to the next potential seed
                                    continue
                                    
                        logging.info("Found a seeder. " + PEATSAbox.mission_name + " seeded by " + possible_seed.mission_name)
                        
                        if PEATSAorder.thin_crust:
                            import Mission
                            SM = Mission.Mission(possible_seed.PEATSApath + possible_seed.PEATSAcrust_path)
                        else:
                            SM = possible_seed.PEATSAcrust
                        seed_code = (possible_seed.mission_name,possible_seed.seed_criteria[seed_criteria[0]],possible_seed.objective_value,SM.Xdescriptions,SM.DecisionVector,possible_seed.PEATSApath + possible_seed.PEATSAcrust_path)                

                        PEATSAbox.seed_boxes[seed_criteria[0]].append(seed_code)
                
                if PEATSAorder.seed_from_infeasible_cases_if_no_feasible_seeds_found and len(PEATSAbox.seed_boxes[seed_criteria[0]]) == 0 and solo_criteria == None:
                    PEATSAbox = self.find_seed_cases(PEATSAorder,PEATSAbox,PEATSAboxes,True,seed_criteria)
                    
                # # Check if we can stop searching because this is the first criterion in the criteria list (true also if its the only one)
                # elif (seed_criteria[0] == PEATSAorder.seed_criteria[0][0]) and found_matching_section == 1:
                #
                #     # We are no longer in the matching section, so we can stop searching
                #     break
             
        # return the PEATSAbox with the seed list set
        return PEATSAbox
    
    def reseed(self,PEATSAorder,PEATSAboxes,history):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        logging.info("Selecting a new seed for all cases")
        
        # Initialize the output array
        boxes_out = []
                
        # Loop through all the cases
        for box in PEATSAboxes:  
        
            # Create a tracker to note if an unseeded case has been written yet
            unseeded_written = False
            
            # Check if we are overriding the local seed_criteria 
            if PEATSAorder.only_one_seed_in_all_seed_directions > 0:
                # We must accumulate a list of all seed criteria
                net_seed_cases = []
                                
                # Loop through the different seed criterion
                for seed_criteria in PEATSAorder.seed_criteria:
                    net_seed_cases += box.seed_boxes[seed_criteria[0]]
                
                # Make sure there are some cases
                if len(net_seed_cases) > 0:        
                    # Call the seed from all routine with the net seed case list
                    if PEATSAorder.only_one_seed_in_all_seed_directions == 1:
                        # seed from best
                        out_box = self.seed_from_best(PEATSAorder,box,net_seed_cases,history)
                    elif PEATSAorder.only_one_seed_in_all_seed_directions == 2:
                        # seed from random
                        out_box = self.seed_from_random(PEATSAorder,box,net_seed_cases,history)
                    else:
                        # invalid input. tell the user and write an unseeded case
                        logging.info("Invalid input value for PEATSAorder.only_one_seed_in_all_seed_directions.")
                        out_box = self.write_unseeded_case(PEATSAorder,box,history)
                    # Add the returned case to the stack
                    boxes_out.append(out_box)
                elif PEATSAorder.wait_until_seeded == 0:
                    # No seed were found, but write an unseeded case 
                    out_box = self.write_unseeded_case(PEATSAorder,box,history)
            
                    # Add the returned case to the stack
                    boxes_out.append(out_box)                
            else:
                # Loop through the different seed criterion
                for seed_criteria in PEATSAorder.seed_criteria:
                
                    # Check if there are seeds available
                    if len(box.seed_boxes[seed_criteria[0]]) == 0 and unseeded_written == False:
                        # There are not
                
                        # Update the flag so we dont write anymore unseededs
                        unseeded_written = True
                
                        logging.info("No seed for " + box.mission_name)
                        
                        # Check if we want to wait until this case has a seed?
                        if PEATSAorder.wait_until_seeded:
                            # We do. dont write an unseeded case
                            continue
                            
                        # Write an unseeded case
                        out_box = self.write_unseeded_case(PEATSAorder,box,history)
                
                        # Add the returned case to the stack
                        boxes_out.append(out_box)
                    elif len(box.seed_boxes[seed_criteria[0]]):
                        # There are seed cases available
                
                        if PEATSAorder.objective_type == 1:
                            if box.re_run_using_this_seed[seed_criteria[0]] == False:
                                continue
                                
                        # Determine the seed criteria
                        if seed_criteria[3] == 0:
                            out_box = self.seed_from_closest(PEATSAorder,box,box.seed_boxes[seed_criteria[0]],seed_criteria[0],history)
                
                            # Add the returned case to the stack
                            boxes_out.append(out_box)
                        elif seed_criteria[3] == 1:
                            out_boxes = self.seed_from_all(PEATSAorder,box,box.seed_boxes[seed_criteria[0]],history)
                
                            # Add the returned cases to the stack
                            for out_box in out_boxes:
                                boxes_out.append(out_box)
                        elif seed_criteria[3] == 2:
                            out_box = self.seed_from_best(PEATSAorder,box,box.seed_boxes[seed_criteria[0]],history)
                
                            # Add the returned case to the stack
                            boxes_out.append(out_box)
                        elif seed_criteria[3] == 3:
                            out_box = self.seed_from_random(PEATSAorder,box,box.seed_boxes[seed_criteria[0]],history)
                            
                            # Add the returned case to the stack
                            boxes_out.append(out_box)
                        elif seed_criteria[3] == 4:
                            out_box = self.seed_from_curve_fit(PEATSAorder,box,box.seed_boxes[seed_criteria[0]],seed_criteria[0],history)
                            
                            boxes_out.append(out_box)
        
        return boxes_out
    
    def seed_from_closest(self,PEATSAorder,box,seed_cases,key,history):
        
        # Sort the seeds in the relevant order
        seeds = sorted(seed_cases, key = lambda case: abs(case[1] - box.seed_criteria[key]))
        
        # Write the new case, seeded from the closest 
        new_box = self.write_seeded_case(PEATSAorder,box,seeds[0],history)
        
        # Return the box with the new options
        return new_box 
        
    def seed_from_all(self,PEATSAorder,box,seed_cases,history):
        
        # Initialize the output array
        out_boxes = []
        
        # Loop through all seed cases
        for seed in seed_cases:
        
            # Write the new case, seeded from the closest 
            new_box = self.write_seeded_case(PEATSAorder,box,seed,history)
        
            # Add the new box to the stack
            out_boxes.append(new_box)
        
        # Return the box with the new options
        return out_boxes
    
    def seed_from_curve_fit(self,PEATSAorder,box,seed_cases,key,history):
        import logging
        import copy
        import Mission
        
        # Make sure we have enough cases
        if len(seed_cases) > PEATSAorder.polyfit_order:
            polyfit_order = PEATSAorder.polyfit_order
        elif len(seed_cases) > 1:
            logging.info("Lowering polyfit order for case: " + box.mission_name)
            polyfit_order = len(seed_cases)-1
        else:
            logging.info("Not enough cases for a curve fit for case: " + box.mission_name)
            return self.seed_from_closest(PEATSAorder,box,seed_cases,key,history)
            
        # Save the nomianl initial guess type
        nominal_initial_guess_type = PEATSAorder.initial_guess_type
        
        # Override the initial guess type
        PEATSAorder.initial_guess_type = "copy"
                
        # Overwrite the mission name
        dummy_seed_mission_name = "Average_of"
            
        # Update the dummy seeds objective_value
        dummy_seed_objective_value = 0.0
        
        # Initialize a list of usable seeds
        usable_seeds = []
        
        # Loop through all the seed cases
        for seed in seed_cases:
            
            # Update the name
            dummy_seed_mission_name += "_" + seed[0].split("_seeded_by")[0].split("_v")[0]
            
            # Update the objective value
            dummy_seed_objective_value += seed[2]
            
            # Create a flag if the dummy needs to be interpolated
            needs_interp = False
            
            # Make sure the seed is the same length
            if len(seed[4]) != len(box.PEATSAdough.trialX):
                # It isnt, so we need to interpolate the timesteps
                needs_interp = True
            else:
                # Loop through all decision vector variables
                for idx,descrip in enumerate(seed[3]):
                    # Make sure all of the decision variable descriptions are the same
                    if descrip.lstrip(" ").rstrip(" \r\n") != box.PEATSAdough.trialX[idx][0].lstrip(" ").rstrip(" \r\n"):
                        needs_interp = True
                        break
                    
            if needs_interp:        
                    
                logging.info("Creating an interpolated initial guess: " + seed[0] + " for " + box.mission_name)
                
                # Create a holder for the number of timesteps per journey
                nSteps = []
                
                # Loop through each journey
                for JO in box.PEATSAdough.Journeys:
                    # Check if this journey overrides the global number of steps or not and append the number of steps
                    if JO.override_num_steps:
                        nSteps.append(JO.number_of_steps)
                    else:
                        nSteps.append(box.PEATSAdough.num_timesteps)
                
                # Call the interpolator to get a new mission optiosn structure with the improved initial guess
                SM = Mission.Mission(seed[5])
                newMO = self.modified_timestep_initial_guess_generator(SM,copy.deepcopy(box.PEATSAdough),"n_steps",nSteps)
                
                # Copy the seed since we are going to overwrite its decision vector
                temp_seed = (seed[0],seed[1],seed[2],[dv[0] for dv in newMO.trialX],[dv[1] for dv in newMO.trialX])
                
                # Loop through all decision vector variables
                for idx,descrip in enumerate(temp_seed[3]):
                    # Make sure all of the decision variable descriptions are the same
                    if descrip.lstrip(" ").rstrip(" \r\n") != box.PEATSAdough.trialX[idx][0].lstrip(" ").rstrip(" \r\n"):
                        raise Exception("This case has the wrong decisionv ector structure: " + descrip + " vs " + box.PEATSAdough.trialX[idx][0])
                        
            else:
                # Dont need to do anything, just copy the seed into a new variable
                temp_seed = seed              
                    
            # Add this seed to the usable seed list
            usable_seeds.append(temp_seed)
        
        # Initialize the new decision vector
        new_dv_values = []       
            
        # Loop through the decision vector entries
        for idx in range(len(box.PEATSAdough.trialX)):
            
            # Polyfit is a member of numpy, so import it
            import numpy
        
            # Get the arrays of data that make the polyfit
            xData = [seed[1] for seed in usable_seeds]
            yData = [seed[4][idx] for seed in usable_seeds]
        
            # Get the coefficients of the polynomial fit
            coeffs = numpy.polyfit(xData,yData,polyfit_order)
        
            # Create the polynomial function object
            poly = numpy.poly1d(coeffs)            
            
            # Assign the new decision vector value
            new_dv_values.append( poly(box.seed_criteria[key]) )
                
        # Create the new dummy seed        
        dummy_seed = (dummy_seed_mission_name,0,dummy_seed_objective_value,[dv[0] for dv in box.PEATSAdough.trialX],new_dv_values)        
                    
        # Write the new case, seeded fromt he dummy 
        new_box = self.write_seeded_case(PEATSAorder,box,dummy_seed,history)
                
        # Put the initial guess type back
        PEATSAorder.initial_guess_type = nominal_initial_guess_type
        
        # Return the box with the new options
        return new_box
        
    def seed_from_best(self,PEATSAorder,box,seed_cases,history):
        
        # Sort the seeds in the relevant order
        seeds = sorted(seed_cases, key = lambda case: case[2])
        
        # Write the new case, seeded from the closest 
        if PEATSAorder.max_or_min == "max":
            new_box = self.write_seeded_case(PEATSAorder,box,seeds[-1],history)
        elif PEATSAorder.max_or_min == "min":
            new_box = self.write_seeded_case(PEATSAorder,box,seeds[0],history)
        
        # Return the box with the new options
        return new_box
        
    def seed_from_random(self,PEATSAorder,box,seed_cases,history):
        import random
        
        # Write the new case, seeded from the closest 
        new_box = self.write_seeded_case(PEATSAorder,box,seed_cases[random.randint(0,len(seed_cases)-1)],history)
        
        # Return the box with the new options
        return new_box
           
    def write_seeded_case(self,PEATSAorder,box,seed,history):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import Mission
        import logging
        import copy
        import PEATSAbox
        
        # NH added a try/except block on 1/4/2022 because I was tired of getting
        # errors like 
        # MO.trialX.append([seed[3][idx],seed[4][idx]])
        # IndexError: list index out of range
        #
        # I wasn't sure what the cause actually was, so I just put the whole
        # thing in a try/except block and told it to use create un unseeded case
        # if creating the seeded case failed
        try:
            logging.info("Writing a new seeded options file for case " + box.mission_name.split("_seeded_by")[0].split("_v")[0] + ", using seed " + seed[0].split("_seeded_by")[0].split("_v")[0])
            
            # Copy the existing options object
            MO = copy.deepcopy(box.PEATSAdough)
            
            # Get the new mission name
            MO.mission_name = box.mission_name.split("_seeded_by")[0].split("_v")[0] + "_seeded_by_" + seed[0].split("_seeded_by")[0].split("_v")[0]
            
            # Update the working directory
            MO.forced_working_directory = PEATSAorder.results_directory
            
            # Turn seeding on, if it wasnt already
            MO.seed_MBH = 1
            
            # Move the decision vector over
            if PEATSAorder.initial_guess_type == 'copy':
                MO.trialX = []
                for idx in range(len(seed[4])):
                    MO.trialX.append([seed[3][idx],seed[4][idx]])
            else:
                SM = Mission.Mission(seed[5])
                MO = self.modified_timestep_initial_guess_generator(SM,MO,PEATSAorder.initial_guess_type,PEATSAorder.new_timestep)       
            
            # Get a tuple of the seed name and its objective value
            seed_code = (seed[0].split("_seeded_by")[0].split("_v")[0],seed[2])
    
            # Assume that we should not skip the first nlp run
            MO.skip_first_nlp_run = 0
            
            # Check the history to see if these seed has been used before
            key = box.mission_name.split("_seeded_by")[0].split("_v")[0]
            # Loop through the iterations
            for iter_idx in range(2,PEATSAorder.iteration+1):
                # Keep a flag if we need to break out of this loop
                break_out = False
                # Check if this case was run in the previous iteration
                if key in history["seed"][-iter_idx].keys():
                    # It was. Check if the seed was used in the previous iteration
                    if seed_code in history["seed"][-iter_idx][key]:
                        # Make sure the case actually ran, because the iteration might have been killed
                        if key in history["run"][-iter_idx].keys():
                            # It was. Loop through the run history from that iteration
                            for prev_entry in history["run"][-iter_idx][key]:
                                # Check if this run is the one with the relevant seed
                                if prev_entry[3] == seed_code:
                                    # It is. 
                                
                                    # Check if the run from the previous generation usign this seed improved
                                    if prev_entry[0] == 0:
                                        # It didnt improve but it also didnt finish, so dont stop searching backwards
                                        pass
                                    elif prev_entry[2] == 0:
                                        # IT did not improve. We should skip the first nlp run
                                        MO.skip_first_nlp_run = 1
                                        # Because we found the most recent time this case has been run with this seed we can stop looping now
                                        break_out = True
                                    else:
                                        # It did. And it did improve so we should keep skip first nlp run = 0 
                                        # Because we found the most recent time this case has been run with this seed we can stop looping now
                                        break_out = True
                                    break
                if break_out:
                    break
                    
            # Apply the override options
            for option in PEATSAorder.override_options:
                try:
                    if "SMO" in option[0]:
                        SMO = MissionOptions.MissionOptions(seed[6])
                    # Check the condition
                    if eval(option[0]):
                        if "SMO" in option[1] and "SMO" not in option[0]:
                            SMO = MissionOptions.MissionOptions(seed[6])
                        
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
                             
            if len(PEATSAorder.constraint_walking):    
                # Keep a counter of how many times this case has been run since last improvement
                nRuns = 0
                # Check the history to see if these seed has been used before
                key = box.mission_name.split("_seeded_by")[0].split("_v")[0]
                # Loop through the iterations
                for iter_idx in range(2,PEATSAorder.iteration+1):
                    # Keep a flag if we need to break out of this loop
                    break_out = False
                    # Check if this case was run in the previous iteration
                    if key in history["objective"][-iter_idx].keys():
                        # It was. Check if the seed was used in the previous iteration
                        if history["objective"][-iter_idx][key][3]:
                            break
                        else:
                            nRuns += 1
                logging.info("Last improvement: " + str(nRuns) + " iterations ago")      
            
            for idx,walk in enumerate(PEATSAorder.constraint_walking):
                if nRuns > 0 and round(float(nRuns)/walk[1]) == float(nRuns)/walk[1]:
                    logging.info("Walking constraint: " + walk[0])
                    exec(walk[0] + " += " + str(walk[2]))
            
            # Write the options file out
            MO.write_options_file(PEATSAorder.cases_directory + MO.mission_name + ".emtgopt")
            
            # Create a new peatsa box to hold the new case
            new_box = PEATSAbox.PEATSAbox(PEATSAorder,PEATSAorder.cases_directory,MO.mission_name + ".emtgopt")
        
            # Make a note of the mission name used to seed and the objective value of the seed
            new_box.used_seed = seed_code
        
            # Return the new box
            return new_box
        except:
            # creating a seeded case failed
            # create an unseeded case instead
            logging.info("WARNING: Created a seeded options file failed.")
            logging.info("Creating an unseeded options file instead.")
            new_box = self.write_unseeded_case(PEATSAorder, box, history)
            
            # Return the new box
            return new_box
            
    def write_unseeded_case(self,PEATSAorder,box,history=None):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import Mission
        import logging
        import copy
        import PEATSAbox
        logging.info("Writing a new unseeded options file for case " + box.mission_name)
        
        # Copy the existing options object
        MO = copy.deepcopy(box.PEATSAdough)
        
        # Get the new mission name
        MO.mission_name = box.mission_name.split("_seeded_by")[0].split("_v")[0]
        
        # Update the working directory
        MO.forced_working_directory = PEATSAorder.results_directory
        
        if PEATSAorder.seed_from_infeasible_self_when_unseeded:
            if hasattr(box,'PEATSAcrust_path') and box.PEATSAcrust_path != None and 'fake' not in box.PEATSAcrust_path:
                MO.trialX = []
                if PEATSAorder.thin_crust:
                    M = Mission.Mission(box.PEATSApath + box.PEATSAcrust_path)
                else:
                    M = box.PEATSAcrust
                for idx in range(len(M.DecisionVector)):
                    MO.trialX.append([M.Xdescriptions[idx],M.DecisionVector[idx]])
                MO.seed_MBH = 1
            elif PEATSAorder.turn_off_seed_MBH_when_unseeded:
                # Turn seeding off, if it wasnt already
                MO.seed_MBH = 0
            else:
                MO.seed_MBH = 1
        elif PEATSAorder.turn_off_seed_MBH_when_unseeded:
            # Turn seeding off, if it wasnt already
            MO.seed_MBH = 0
        else:
            MO.seed_MBH = 1
        
        # Apply the override options
        for option in PEATSAorder.override_options:
            try:
                if "SMO" in option[0]:
                    SMO = MissionOptions.MissionOptions(seed[6])
                # Check the condition
                if eval(option[0]):
                    if "SMO" in option[1] and "SMO" not in option[0]:
                        SMO = MissionOptions.MissionOptions(seed[6])
                    
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
        
        if len(PEATSAorder.constraint_walking):   
            # Keep a counter of how many times this case has been run since last improvement
            nRuns = 0
            if history != None: 
                # Check the history to see if these seed has been used before
                key = box.mission_name.split("_seeded_by")[0].split("_v")[0]
                # Loop through the iterations
                for iter_idx in range(2,PEATSAorder.iteration+1):
                    # Keep a flag if we need to break out of this loop
                    break_out = False
                    # Check if this case was run in the previous iteration
                    if key in history["objective"][-iter_idx].keys():
                        # It was. Check if the seed was used in the previous iteration
                        if history["objective"][-iter_idx][key][3]:
                            break
                        else:
                            nRuns += 1
                logging.info("Last improvement: " + str(nRuns) + " iterations ago")      
        
        for idx,walk in enumerate(PEATSAorder.constraint_walking):
            if nRuns > 0 and round(float(nRuns)/walk[1]) == float(nRuns)/walk[1]:
                logging.info("Walking constraint: " + walk[0])
                exec(walk[0] + " += " + str(walk[2]))
            
        # Write the options file out
        MO.write_options_file(PEATSAorder.cases_directory + MO.mission_name + ".emtgopt")
        
        # Create a new peatsa box to hold the new case
        new_box = PEATSAbox.PEATSAbox(PEATSAorder,PEATSAorder.cases_directory,MO.mission_name + ".emtgopt")
    
        # Return the new box
        return new_box
    
    def CreateSplines(self,reference_ephem,current_epoch,turn_off_point_density_check = 0):
        from scipy.interpolate import interp1d
        
        # Star the index at zero
        reference_ephem_index = 0
        
        # Find the most recent ephemeris
        # Get the date of the ephemeris at the current index:
        ephem_epoch = reference_ephem[0][reference_ephem_index]
        # Loop through the ephemeris data until we are past the current_epoch
        while ephem_epoch < current_epoch:
            # Increment the index
            reference_ephem_index += 1
            # Update the date
            ephem_epoch = reference_ephem[0][reference_ephem_index]
        
        # Check if the lhs point is closer to the current epoch than the rhs
        if reference_ephem_index > 0 and ephem_epoch - current_epoch > current_epoch - reference_ephem[0][reference_ephem_index-1]:
            # It is, so lets go back to it
            reference_ephem_index -= 1
            # Update the date
            ephem_epoch = reference_ephem[0][reference_ephem_index]
                            
        # We should now have the index into the reference ephemeris data
        # But lets make sure that we arent too far off, timewise
        if turn_off_point_density_check == 0:
            if abs(ephem_epoch-current_epoch) > 2.0:
                # We should always have ephemeris closer to the current epoch than half a day, 
                # because we should always have ephemeris at at least 1 day intervals. Because we dont
                # I dont trust this ephemeris, so throw an error
                raise Exception("Not able to find an ephemeris data point within .5 days of the desired missed thrust date (JD" + str(current_epoch) + " vs JD" + str(ephem_epoch) + "). Either the ephemeris is out of order, or it is not dense enough. Current index into the ephemeris = " + str(reference_ephem_index))
    
        # Get the first index to use for the spline. Note we are making a spline with only +- 5 points around the missed thrust date because it is too costly to load the full journey let alone mission
        if reference_ephem_index < 5:
            spline_start_index = 0
        else:
            spline_start_index = reference_ephem_index - 5
                        
        # Check that the phase start isnt within 10 entries of the end of the data
        if spline_start_index + 11 > len(reference_ephem[7]):
            spline_start_index = len(reference_ephem[7]) - 10                
                            
        # Initialize the spline list
        reference_splines = []
                               
        # Create the splines and add them to the list
        for spline_idx in range(7):
            try:
                reference_splines.append(interp1d(reference_ephem[0][spline_start_index:spline_start_index+10],reference_ephem[1+spline_idx][spline_start_index:spline_start_index+10],kind='cubic'))
            except:
                reference_splines.append(interp1d(reference_ephem[0][spline_start_index:spline_start_index+10],reference_ephem[1+spline_idx][spline_start_index:spline_start_index+10],kind='linear'))
                    
                        
        return reference_splines
     
    def get_ephem_lists_from_emtg_mission(self,M):
        
        # Initialize the lists:
        ephem_lists = [[] for i in range(8)]
        
        # loop through the journeys
        for J in M.Journeys:
            
            # Loop through the events
            for e in J.missionevents:
                
                if e.JulianDate not in ephem_lists[0]:
                    # Get the epoch
                    ephem_lists[0].append(e.JulianDate)
                
                    # Loop through the states
                    for sdx in range(6):
                    
                        # Append the state
                        ephem_lists[1+sdx].append(e.SpacecraftState[sdx])
                
                    # APpend the mass
                    ephem_lists[-1].append(e.Mass)
        
        # Send back the output
        return ephem_lists
            
    def modified_timestep_initial_guess_generator(self,SM,MO,timestep_type,timestep_input):
        # SM = seed mission. the initial guess source. it's journeys have n timesteps per journey (or n1, n2, n3...)
        # MO = destination mission options. where to put the mission options. It will be given an initial guess and 
        #      m timesteps of ~x days each
        from math import floor
                
        # Do some error checking
        if len(SM.Journeys) != len(MO.Journeys):
            raise Exception("Sorry, this initial guess generator only works with the same number of journeys")
        
        # Initialize the decision vector
        new_decision_vector = []
        
        # Get an index to keep track of where we are in the 
        # seed mission decision vector
        decision_vector_index_tracker = 0
        
        # Get the reference ephem in case we need it
        reference_ephem = self.get_ephem_lists_from_emtg_mission(SM)
        
        # Loop through the journeys
        for jdx in range(len(MO.Journeys)):
                        
            # Loop through the phases
            # Determine if this is a multi-phase Journey.
            # First, assume it only has one phase
            nPhases = 1
            # Loop through all of the sequence values in this Journey
            for seq_idx in range(0,len(MO.Journeys[jdx].sequence)):
                # Check if this sequence ID is non-zero, in which case there is a flyby and this is a multiphase journey
                if (MO.Journeys[jdx].sequence[seq_idx]):
                    # It is, update the counter
                    nPhases += 1
                    
            # If there are more than 9 phases it will break things, so throw an error
            if nPhases > 9:
                raise Exception("Phase number finding code in the decision vector re-creation cant handle more than 9 phases currently")
                                        
            # # Grab the objects for this Journey and journey options
            J = SM.Journeys[jdx]
            JO = MO.Journeys[jdx]
            
            # Loop through all of the phases
            for pdx in range(0,nPhases):
                
                # We need to figure out if there is a departure and/or arrival coast
                # Start by assuming both are zero
                sm_departure_coast = 0
                sm_arrival_coast = 0
                duty_cycle_coast = 0
                sm_control_timestep = -1
                # Then determine whether the journey is single phase or multiphase and whetehr
                # we are starting in the relevant phase or not
                if pdx == 0:
                    in_phase = True
                    in_thrusting = False
                else:
                    in_phase = False
                # Loop through all the events in the journey
                for event in J.missionevents[1:-1]:
                    
                    # Check if this is a flyby
                    if 'flyby' in event.EventType:
                        # It is. A phase is starting or ending
                        
                        # If we were in the phase, we cna stop looping
                        if in_phase:
                            break
                        else:
                            in_phase = True 
                            in_thrusting = False
                            continue
                    
                    # Check if we are in the relevant phase or not
                    if in_phase == 0:
                        # We arent, so move along to the next event
                        continue
                        
                    # Check if this is a trhust or coast
                    if 'force-coast' in event.EventType:
                        # It is the forced coast. So just need to determine if it is a departure or arrival
                        if in_thrusting == False:
                            sm_departure_coast = event.TimestepLength
                        elif in_thrusting == True:
                            sm_arrival_coast = event.TimestepLength
                    elif in_thrusting == False and 'thrust' in event.EventType:
                        # It isnt a coast, but it is a thrust, so make sure the thrusting flag is set
                        in_thrusting = True    
                        
                        # Get the timestep
                        sm_control_timestep = event.TimestepLength            
                    elif 'coast' == event.EventType or "nav-coast" == event.EventType:
                        # Check if this is a realistic duty cycle coast
                        if event.TimestepLength < sm_control_timestep :
                            # It must be
                            duty_cycle_coast = event.TimestepLength
                            
                # If this is a full coast phase (even though it COULD thrust, EMTG chose not to)
                if sm_control_timestep == -1:
                    for event in J.missionevents[1:-1]:
                        if event.EventType == "coast":
                            sm_control_timestep = event.TimestepLength
                        elif "nav-coast" == event.EventType:
                            # It must be
                            duty_cycle_coast = event.TimestepLength
                
                # Add the duty cycle coast ot the control timestep if any
                sm_control_timestep += duty_cycle_coast    
                
                sm_num_steps = ( J.missionevents[-1].JulianDate - J.missionevents[0].JulianDate - sm_departure_coast - sm_arrival_coast ) / sm_control_timestep               
                                
                # Get the length of the initial coast, if any
                if pdx == 0:
                    departure_coast = JO.forced_initial_coast
                    phase_start_epoch = J.missionevents[0].JulianDate
                else:
                    departure_coast = MO.forced_post_flyby_coast
                    phase_start_epoch = J.missionevents[int(len(J.missionevents)/2)].JulianDate
                                            
                # Get the length of the arrival coast, if any
                if pdx == nPhases - 1:
                    arrival_coast = JO.forced_terminal_coast
                else:
                    arrival_coast = MO.forced_pre_flyby_coast
                
                # Create the decision vector prefix
                prefix = "j" + str(jdx) + "p" + str(pdx)
                                
                # Initialize the transcription string
                transcription = ""
                
                # Initialize the phase flight time to make sure its been found
                phase_flight_time = 0
                
                # Note that we havent found the control start index
                control_start_index = -1
                
                # Loop through the decision vector
                for dv_idx in range(decision_vector_index_tracker,len(SM.Xdescriptions)):
                    
                    currDVvar = SM.Xdescriptions[dv_idx]
                                        
                    # Make sure we are in the correct phase
                    if prefix not in currDVvar:
                        if transcription != "":
                            # Update the tracker
                            decision_vector_index_tracker = dv_idx
                        
                            # Break the loop
                            break
                    
                        else:
                            # Not to the right phase yet
                            continue
                            
                    
                    # Get the transcription
                    if transcription == "":
                        if "FBLT" in currDVvar:
                            transcription = "FBLT"
                        elif "PSFB" in currDVvar:
                            transcription = "PSFB"
                        elif "MGALT" in currDVvar:
                            transcription = "MGALT"
                        elif "CoastPhase" in currDVvar:
                            transcription = "Coast" 
                        else:
                            raise Exception("This transcription is not in the initial guess utility yet. Sorry")
                          
                    # Check if we are in the steps yet
                    if ((transcription == "MGALT" or transcription == "FBLT") and "step" in currDVvar) or (transcription == "PSFB" and "Step" in currDVvar):
                        # We are in the steps.
                        
                        # Check if we have already added the control for this phase
                        if control_start_index == -1:
                            # We have not
                            
                            # Note that we have found the control start
                            control_start_index = dv_idx
                        
                            # Make sure we have the phase flight time
                            if phase_flight_time == 0:
                                raise Exception("Phase flight time wasnt found for " + prefix)
                            
                            # Set the number of variables per step in the control section
                            multiplier = 3
                        
                            # If the transcription is parallel shooting, there are actually 12, not 3 variables
                            if transcription == "PSFB":
                                multiplier = 9 + 3 * JO.num_interior_control_points
                                if JO.num_interior_control_points > 1:
                                    raise Exception("Time interpolation not setup to handle multiple interior contorl points yet. Use 'copy' instead of 'x_days' or 'n_steps'")
                            
                            # Override the departure journey number of time steps
                            JO.override_num_steps = 1
                                                                  
                            # Calculate the number of timesteps to use from inputs
                            if timestep_type == "x_days":
                                if isinstance(timestep_input,list):
                                    if len(timestep_input) < len(MO.Journeys):
                                        raise Exception("Not enough entries in new_timestep")
                                    else:
                                        if isinstance(timestep_input[jdx-len(MO.Journeys)],list):
                                            local_list = timestep_input[jdx-len(MO.Journeys)]
                                            if len(local_list) != 3:
                                                raise Exception("Wrong entry size")
                                            if sm_control_timestep > local_list[0] and sm_control_timestep < local_list[2]:
                                                nDays = sm_control_timestep
                                            else:
                                                nDays = local_list[1]
                                                print("Changing to " + str(nDays))
                                        else:
                                            nDays = timestep_input[jdx-len(MO.Journeys)]
                                        
                                elif timestep_input > 0:
                                    nDays = timestep_input
                                elif timestep_input == -2:
                                    nDays = sm_control_timestep
                                    
                                JO.number_of_steps = int(round((phase_flight_time - departure_coast - arrival_coast) / nDays))
                                                                
                                # Make sure we have at least one timestep
                                if JO.number_of_steps <= 0:
                                    JO.number_of_steps = 1
                            elif timestep_type == "n_steps":
                                if isinstance(timestep_input,list):
                                    if len(timestep_input) < len(MO.Journeys):
                                        raise Exception("Not enough entries in new_timestep")
                                    else:
                                        JO.number_of_steps = timestep_input[jdx-len(MO.Journeys)]
                                elif timestep_input > 0:
                                    JO.number_of_steps = timestep_input
                                elif timestep_input == -1:
                                    JO.number_of_steps = MO.num_timesteps     
                                elif timestep_input == -2:
                                    JO.number_of_steps = sm_num_steps                               
                        
                            # Timestep wont be exactly x. Calculate what it is exactly
                            control_timestep = (phase_flight_time - arrival_coast - departure_coast) / JO.number_of_steps 
                        
                            # if abs(control_timestep - sm_control_timestep)
                        
                            # Set the current epoch
                            control_epoch = phase_start_epoch + departure_coast
                        
                            # Set the epoch when sm control started and ended
                            sm_control_start_epoch = phase_start_epoch + sm_departure_coast
                            sm_control_end_epoch = phase_start_epoch + phase_flight_time - sm_arrival_coast
                        
                            # Loop through the new control steps
                            for control_index in range(0,JO.number_of_steps):
                        
                                # Find the epoch of the midpoint of this control arc
                                if control_index == 0:
                                    control_epoch += .5 * control_timestep
                                else:
                                    control_epoch += control_timestep
                                        
                                # Figure out how long its been since the control began
                                time_since_sm_control_start = control_epoch - sm_control_start_epoch
                    
                                # Determine how far we are fractionally through the phase
                                phase_fraction = time_since_sm_control_start / (phase_flight_time - sm_departure_coast - sm_arrival_coast)
                                                
                                # Make sure we didnt overstep 
                                if phase_fraction > 1:
                                    # We did. Check if thats okay
                                    if arrival_coast < sm_arrival_coast and control_epoch - sm_control_end_epoch < sm_arrival_coast - arrival_coast:
                                        # We are past where we should be, but its okay because we are still inside of the reduced forced coast
                                        pass
                                    else:
                                        raise Exception("The closest control step is out of bounds and it shouldnt be")      
                                elif phase_fraction < 0:
                                    # we are before control starts
                                    if departure_coast < sm_departure_coast and sm_control_start_epoch - control_epoch < sm_departure_coast - departure_coast:
                                        # we are before control starts, but its okay because we are just in the reduced forced coast
                                        pass
                                    else:
                                        raise Exception("The closest control step is out of bounds and it shouldnt be")                                         
                                                
                                # Find the closest control step
                                closest_control_step = int(round(phase_fraction * sm_control_steps))

                                # Store the closest control step if this is the first step
                                if control_index == 0:
                                    first_control_step = closest_control_step
                                                                
                                # If we are here, it should be okay if the closest control step > number of steps, we just need to not try to extract the decision vector for those
                                if closest_control_step > sm_control_steps or closest_control_step < 0:
                                    # We are past the previous mission's control phase, so we must have reduced the length of the arrival coast
                                    from IPython import embed
                                    embed()
                                    stop
                                    # Adding control variables is different depending on transcription, so check what type is being used
                                    if transcription == "FBLT" or transcription == "MGALT":
                                        # Add the control command to the new decision vector
                                        new_decision_vector.append([prefix + transcription + ": step " + str(control_index) + " u_x",0.0])
                                        new_decision_vector.append([prefix + transcription + ": step " + str(control_index) + " u_y",0.0])
                                        new_decision_vector.append([prefix + transcription + ": step " + str(control_index) + " u_z",0.0])
                                    elif transcription == "PSFB":
                                        import math
                                                                        
                                        # Get the ephemeris
                                        reduced_coast_splines = self.CreateSplines(reference_ephem,control_epoch - control_timestep / 2.0,1)
                                                                            
                                        # Get the cartesian state
                                        x = reduced_coast_splines[0](control_epoch - control_timestep / 2.0)
                                        y = reduced_coast_splines[1](control_epoch - control_timestep / 2.0)
                                        z = reduced_coast_splines[2](control_epoch - control_timestep / 2.0)
                                        vx = reduced_coast_splines[3](control_epoch - control_timestep / 2.0)
                                        vy = reduced_coast_splines[4](control_epoch - control_timestep / 2.0)
                                        vz = reduced_coast_splines[5](control_epoch - control_timestep / 2.0)
                                                                        
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
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state r",r])
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state RA",RA])
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state DEC",DEC])
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state v",v])
                                        if MO.PSFBstateRepresentation == 1:
                                            new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state vRA",vRA])
                                            new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state vDEC",vDEC])                                            
                                        elif MO.PSFBstateRepresentation == 2:
                                            new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state AZ",AZ])
                                            new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state FPA",FPA])
                                        # Add the masses
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state mass",reduced_coast_splines[6](control_epoch - control_timestep / 2.0)])
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": virtual chemical fuel",0.0])
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": virtual electric propellant",net_ep_used])
                                        for ss_idx in range(JO.num_interior_control_points):
                                        # Add the control command to the new decision vector
                                            new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": substep" + str(ss_idx) +" u_x",0.0])
                                            new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": substep" + str(ss_idx) +" u_y",0.0])
                                            new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": substep" + str(ss_idx) +" u_z",0.0])
                                    
                                else:
                                    # This step is within the decision vector of hte initial guess, so just get those values and add them to the new decision vector
                                
                                    # Find the index of this control step
                                    closest_control_index = closest_control_step * multiplier + control_start_index
                            
                                    print(closest_control_step)
                                    print(control_index)
                                    if control_index != closest_control_step:
                                        from IPython import embed
                                        embed()
                                        stop
                                    print("-----------------")
                            
                                    # Adding control variables is different depending on transcription, so check what type is being used
                                    if transcription == "FBLT" or transcription == "MGALT":
                                        # Add the control command to the new decision vector
                                        new_decision_vector.append([prefix + transcription + ": step " + str(control_index) + " u_x",SM.DecisionVector[closest_control_index + 0]])
                                        new_decision_vector.append([prefix + transcription + ": step " + str(control_index) + " u_y",SM.DecisionVector[closest_control_index + 1]])
                                        new_decision_vector.append([prefix + transcription + ": step " + str(control_index) + " u_z",SM.DecisionVector[closest_control_index + 2]])
                                    elif transcription == "PSFB":
                                        # Add the state
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state r",SM.DecisionVector[closest_control_index]])
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state RA",SM.DecisionVector[closest_control_index + 1]])
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state DEC",SM.DecisionVector[closest_control_index + 2]])
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state v",SM.DecisionVector[closest_control_index + 3]])
                                        if "AZ" in SM.Xdescriptions[closest_control_index + 4] and MO.PSFBstateRepresentation == 2:
                                            new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state AZ",SM.DecisionVector[closest_control_index + 4]])
                                            new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state FPA",SM.DecisionVector[closest_control_index + 5]])
                                        elif "vRA" in SM.Xdescriptions[closest_control_index + 4] and MO.PSFBstateRepresentation == 1:
                                            new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state vRA",SM.DecisionVector[closest_control_index + 4]])
                                            new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state vDEC",SM.DecisionVector[closest_control_index + 5]])
                                        else:
                                            raise Exception("I dont do psfb state conversions yet. sorry")
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": left state mass",SM.DecisionVector[closest_control_index + 6]])
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": virtual chemical fuel",SM.DecisionVector[closest_control_index + 7]])
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": virtual electric propellant",SM.DecisionVector[closest_control_index + 8]])
                                        # Add the control command to the new decision vector
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": substep0 u_x",SM.DecisionVector[closest_control_index + 9]])
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": substep0 u_y",SM.DecisionVector[closest_control_index + 10]])
                                        new_decision_vector.append([prefix + transcription + "_Step" + str(control_index) + ": substep0 u_z",SM.DecisionVector[closest_control_index + 11]])
                    else:
                        
                        # Check if this variable is the phase flight time
                        if "phase flight time" in currDVvar:
                            # It is
                            
                            # Set the phase flight time
                            phase_flight_time = SM.DecisionVector[dv_idx]
                        
                            # Get the number of timesteps
                            sm_control_steps = (phase_flight_time - sm_departure_coast - sm_arrival_coast) / sm_control_timestep
                        
                            # Make sure we got an integer
                            if abs(round(sm_control_steps) - sm_control_steps) > 1e-6:
                                 if phase_flight_time < 1:
                                     # This is probably okay
                                     sm_control_steps = 1
                                 else:
                                     raise Exception("sm_control steps is not an integer. Somethign went wrong finding the control timestep lengths")
                            else:
                                sm_control_steps = int(sm_control_steps)
                            
                        # Check if this variable is the ep    
                        elif "virtual electric propellant" in currDVvar:
                            # It is
                            
                            # Store the ep used in this phase
                            net_ep_used = SM.DecisionVector[dv_idx]
                        
                        # Add this variable to the decision vector
                        new_decision_vector.append([currDVvar,SM.DecisionVector[dv_idx]])
        
        MO.trialX = new_decision_vector
        
        # return the estination options structure with the new initial guess
        return MO
    
    
    
    
    
    
    
    
    
    
    
        
        
        
