



def MGAnDSMs_3D_Flybys_to_PatchedConic(RefMO):
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

        # remove any white space from the TrialX descriptions
        for entry in RefMO.trialX:
            entry[0] = entry[0].strip()
                 
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
            elif "MGAnDSMs" in RefMO.trialX[local_decision_vector_index][0]:
                transcription = "MGAnDSMs"
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
                    newTrialX.append(["j0p0MGAnDSMsEphemerisPeggedLaunchDirectInsertion: event left state epoch",entry[1]])
                    break
            for entry in RefMO.trialX:
                if entry[0] == "j0p0CoastPhaseEphemerisReferencedInterceptInterior: event interface state vMAG":
                    newTrialX.append(["j0p0MGAnDSMsEphemerisPeggedLaunchDirectInsertion: magnitude of outgoing velocity asymptote",entry[1]])
                    break
            for entry in RefMO.trialX:
                if entry[0] == "j0p0CoastPhaseEphemerisReferencedInterceptInterior: event interface state vRA":
                    newTrialX.append(["j0p0MGAnDSMsEphemerisPeggedLaunchDirectInsertion: RA of departure asymptote",entry[1]])
                    break
            for entry in RefMO.trialX:
                if entry[0] == "j0p0CoastPhaseEphemerisReferencedInterceptInterior: event interface state vDEC":
                    newTrialX.append(["j0p0MGAnDSMsEphemerisPeggedLaunchDirectInsertion: DEC of departure asymptote",entry[1]])
                    break
        
        for jdx,JO in enumerate(MO.Journeys):
            if orig_jdx[jdx] +1 != len(RefMO.Journeys):
                if RefMO.Journeys[orig_jdx[jdx] + 1].phase_type[0] == 7:
                    JO.arrival_class = 0
                    JO.forced_terminal_coast += all_phase_flight_times["j" + str(orig_jdx[jdx]+1) + "p0"]
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]) + "p0MGAnDSMsEphemerisReferencedInterceptExterior: event left state mass":
                            newTrialX.append(["j" + str(jdx) + "p0MGAnDSMsEphemerisPeggedIntercept: event left state mass",entry[1]])
                            break
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]) + "p0MGAnDSMsEphemerisReferencedInterceptExterior: event interface state vMAG":
                            vMag = float(entry[1])
                            break
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]) + "p0MGAnDSMsEphemerisReferencedInterceptExterior: event interface state vRA":
                            vRa = float(entry[1])
                            break
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]) + "p0MGAnDSMsEphemerisReferencedInterceptExterior: event interface state vDEC":
                            vDec = float(entry[1])
                            break
                    newTrialX.append(["j" + str(jdx) + "p0MGAnDSMsEphemerisPeggedIntercept: V_infinity_x",vMag * math.cos(vRa) * math.cos(vDec)])
                    newTrialX.append(["j" + str(jdx) + "p0MGAnDSMsEphemerisPeggedIntercept: V_infinity_y",vMag * math.sin(vRa) * math.cos(vDec)])
                    newTrialX.append(["j" + str(jdx) + "p0MGAnDSMsEphemerisPeggedIntercept: V_infinity_z",vMag * math.sin(vDec)])
                    for entry in RefMO.trialX:
                        if entry[0] == "j" + str(orig_jdx[jdx]) + "p0MGAnDSMsEphemerisReferencedInterceptExterior: virtual chemical fuel":
                            newTrialX.append(["j" + str(jdx) + "p0MGAnDSMsEphemerisPeggedIntercept: virtual chemical fuel",entry[1]])
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
                    newTrialX.append(["j" + str(jdx) + "p0MGAnDSMsEphemerisPeggedUnpoweredFlyby: V_infinity_x",vMag * math.cos(vRa) * math.cos(vDec)])
                    newTrialX.append(["j" + str(jdx) + "p0MGAnDSMsEphemerisPeggedUnpoweredFlyby: V_infinity_y",vMag * math.sin(vRa) * math.cos(vDec)])
                    newTrialX.append(["j" + str(jdx) + "p0MGAnDSMsEphemerisPeggedUnpoweredFlyby: V_infinity_z",vMag * math.sin(vDec)])
                    
            if JO.destination_list[0] in JO.perturbation_bodies:
                JO.perturbation_bodies.remove(JO.destination_list[0])
            if JO.destination_list[1] in JO.perturbation_bodies:
                JO.perturbation_bodies.remove(JO.destination_list[1])
        
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
    
            if "EphemerisReferenced" in entry[0]:
                continue    

            # Get the transcription type
            if "FBLT" in entry[0]:
                transcription = "FBLT"
            elif "PSFB" in entry[0]:
                transcription = "PSFB"
            elif "MGALT" in entry[0]:
                transcription = "MGALT"
            elif "CoastPhase" in entry[0]:
                transcription = "CoastPhase"
            elif "MGAnDSMs" in entry[0]:
                transcription = "MGAnDSMs"
                
            # Get the prefix
            phase_string = entry[0].split(transcription)[0] 

            nums = phase_string.split("p")
            journeyIndex = int(nums[0].lstrip('j'))
            phaseIndex = int(nums[1])

            if journeyIndex not in to_pop:
                if "FreePointFreeDirect" in entry[0] and MO.Journeys[new_jdx[journeyIndex]].departure_type != 2:
                    continue
                new_prefix = "j" + str(new_jdx[journeyIndex]) + "p0"
                newTrialX.append([new_prefix + transcription + entry[0].split(transcription)[1] ,entry[1]])

        # now sort the TrialX so everything is grouped together
        trialX_list = []
        for entry in newTrialX:
        
             # Get the transcription type
            if "FBLT" in entry[0]:
                transcription = "FBLT"
            elif "PSFB" in entry[0]:
                transcription = "PSFB"
            elif "MGALT" in entry[0]:
                transcription = "MGALT"
            elif "CoastPhase" in entry[0]:
                transcription = "CoastPhase"
            elif "MGAnDSMs" in entry[0]:
                transcription = "MGAnDSMs"
                
            # Get the prefix
            phase_string = entry[0].split(transcription)[0]

            nums = phase_string.split("p")
            journeyIndex = int(nums[0].lstrip('j'))
            phaseIndex = int(nums[1])

            trialX_list.append([journeyIndex, phaseIndex, entry])

        trialX_list = sorted(trialX_list, key=lambda x: (x[0], x[1]))

        newTrialX = []
        for entry in trialX_list:
            newTrialX.append(entry[2])

        MO.trialX = newTrialX   
                
        return MO


if __name__ == '__main__':
    import MissionOptions

    # First read in an existing options file that has high-fidelity flybys
    high_fidelity_options_file = "C:/Discovery/MAGIC/high_fidelity/Multiple_Maneuvers/Launch_to_Callisto_SOI_HighFidelity_July_2025_AtlasV401_10172018_161254/Launch_to_Callisto_SOI_HighFidelity_July_2025_AtlasV401.emtgopt"
    MO_high_fidelity = MissionOptions.MissionOptions(high_fidelity_options_file)

    # convert the high-fidelity script to zero-radius SOI flybys
    MO_low_fidelity = MGAnDSMs_3D_Flybys_to_PatchedConic(MO_high_fidelity)

    low_fidelity_options_file = "C:/Discovery/MAGIC/high_fidelity/Multiple_Maneuvers/SOI_to_GCCCCCCCCC_zero_radius.emtgopt"
    MO_low_fidelity.write_options_file(low_fidelity_options_file)