#EMTG: Evolutionary Mission Trajectory Generator
#An open-source global optimization tool for preliminary mission design
#Provided by NASA Goddard Space Flight Center
#
#Copyright (c) 2014 - 2020 United States Government as represented by the
#Administrator of the National Aeronautics and Space Administration.
#All Other Rights Reserved.
#
#Licensed under the NASA Open Source License (the "License"); 
#You may not use this file except in compliance with the License. 
#You may obtain a copy of the License at:
#https://opensource.org/licenses/NASA-1.3
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
#express or implied.   See the License for the specific language
#governing permissions and limitations under the License.

#script to convert multi-phase journeys to single-phase journeys

import sys
sys.path.append("c:/emtg/PyEMTG/")
import Mission
import MissionOptions
import Universe
import copy


def convert_to_single_phase_journeys(originalOptions):

    #copy the old options file
    newOptions = copy.deepcopy(originalOptions)

    #clear the Journeys list
    newOptions.Journeys = []

    #walk through the journeys in the old options file and split them if needed
    for oldJourney in originalOptions.Journeys:
        #Step 1: does this journey have only one phase? If so just copy it and move on
        if oldJourney.sequence == []:
            newOptions.Journeys.append(copy.deepcopy(oldJourney))
            continue
            
        #Step 2: if we made it this far then there are multiple phases. We'll need to make them into journeys and add them one at a time.
        for phaseIndex in range(0, len(oldJourney.sequence) + 1):
            #Step 2.2: copy most settings from the old journey
            newJourney = copy.deepcopy(oldJourney)
            
            #Step 2.3: set the new destination list
            #if this is the first phase then we don't need to change the departure body, otherwise...
            if phaseIndex > 0: #not the first phase
                newJourney.destination_list[0] = oldJourney.sequence[phaseIndex - 1]
                
            #if this is the last phase then we don't need to change the arrival body otherwise...
            if phaseIndex < len(oldJourney.sequence):
                newJourney.destination_list[1] = oldJourney.sequence[phaseIndex]
                
            #Step 2.4: clear the new journey's flyby sequence
            newJourney.sequence = []
            
            #Step 2.5: set the new departure and arrival event
            #if this is the first phase then we don't need to change the departure event, otherwise...
            if phaseIndex > 0: #not the first phase
                newJourney.departure_class = 0 #ephemeris-pegged
                newJourney.departure_type = 3 #flyby
                
            #if this is the last phase then we don't need to change the arrival event, otherwise...
            if phaseIndex < len(oldJourney.sequence):
                newJourney.arrival_class = 0 #ephemeris-pegged
                newJourney.arrival_type = 2 #intercept with bounded v-infinity
                #we are intentionally not adjusting the v-infinity magnitude constraint. Users can adjust this if they want to.
            
            #Step 2.6: set the new initial guess
            prefix = 'p' + str(phaseIndex)
            newJourney.trialX = []
            for Xentry in oldJourney.trialX:
                if prefix in Xentry[0]:
                    #copy the entry
                    newEntry = Xentry
                    
                    #rename to phase 0 for a single-phase journey
                    newEntry[0] = newEntry[0].replace(prefix, 'p0')
                    
                    #replace EphemerisPeggedFlybyIn (only exists in multi-phase journeys) with EphemerisPeggedIntercept
                    newEntry[0] = newEntry[0].replace('EphemerisPeggedFlybyIn', 'EphemerisPeggedIntercept')
                    
                    #add the new entry to the decision vectors
                    newJourney.trialX.append(newEntry)
            
            #Step 2.7: set the new scripted constraints
            newJourney.ManeuverConstraintDefinitions = []
            newJourney.BoundaryConstraintDefinitions = []
            newJourney.PhaseDistanceConstraintDefinitions = []
            
            #maneuver constraints
            for entry in oldJourney.ManeuverConstraintDefinitions:
                if prefix in entry:
                    #copy the entry
                    newEntry = entry
                    
                    #rename to phase 0 for a single-phase journey
                    newEntry = newEntry.replace(prefix, 'p0')
                    
                    #add the new entry to the decision vectors
                    newJourney.ManeuverConstraintDefinitions.append(newEntry)
                    
            #boundary constraints
            for entry in oldJourney.BoundaryConstraintDefinitions:
                if prefix in entry:
                    #copy the entry
                    newEntry = entry
                    
                    #rename to phase 0 for a single-phase journey
                    newEntry = newEntry.replace(prefix, 'p0')
                    
                    #add the new entry to the decision vectors
                    newJourney.BoundaryConstraintDefinitions.append(newEntry)
            
            #phase distance constraints
            for entry in oldJourney.PhaseDistanceConstraintDefinitions:
                if prefix in entry:
                    #copy the entry
                    newEntry = entry
                    
                    #rename to phase 0 for a single-phase journey
                    newEntry = newEntry.replace(prefix, 'p0')
                    
                    #add the new entry to the decision vectors
                    newJourney.PhaseDistanceConstraintDefinitions.append(newEntry)

            #step 2.8: journey time bounds, forced coasts, and post-journey delta-v
            #forced coasts
            if phaseIndex == 0:
                if newJourney.departure_type == 0 and newJourney.departure_class == 0:
                    newJourney.forced_initial_coast = originalOptions.forced_post_launch_coast
            else:
                newJourney.forced_initial_coast = originalOptions.forced_post_flyby_coast

            if phaseIndex < len(oldJourney.sequence):                
                    newJourney.forced_terminal_coast = originalOptions.forced_pre_flyby_coast

            #time bounds and post-journey delta-v
            if phaseIndex == len(oldJourney.sequence):
                if newJourney.timebounded == 1:# if the journey has time bounds, as opposed to arrival date bounds or aggregate time bounds
                    newJourney.timebounded = 0
            else: #if this isn't the last phase of the journey, turn off all time bounds and post-journey delta-v
                newJourney.timebounded = 0
                newJourney.journey_end_TCM = 0.0
                newJourney.journey_end_deltav = 0.0
            
            #Step 2.9: rename the journey
            newJourney.journey_name = oldJourney.journey_name + '_phase' + str(phaseIndex)

            #Step 2.10: staging
            if phaseIndex == len(oldJourney.sequence):
                newJourney.stage_before_arrival = 0
                newJourney.stage_after_arrival = 0
            else:
                newJourney.stage_after_departure = 0
            
            #Step 2.11: attach the new journey to the mission
            newOptions.Journeys.append(copy.deepcopy(newJourney))
                
    #get rid of the mission-level forced coasts since they are now captured at the journey level in the new options object
    newOptions.forced_post_flyby_coast = 0.0
    newOptions.forced_pre_flyby_coast = 0.0
    newOptions.forced_post_launch_coast = 0.0

    #set the number of journeys
    newOptions.number_of_journeys = len(newOptions.Journeys)

    #assemble the master decision and constraint vectors
    newOptions.AssembleMasterDecisionVector()
    newOptions.AssembleMasterConstraintVectors()

    #rename the mission
    newOptions.mission_name = originalOptions.mission_name + '_singlePhase'
    
    return newOptions
    
    
if __name__ == '__main__':
    # Ensure the correct number of command line options were provided
    if len(sys.argv) < 2:
        raise Exception("Unknown number of command line options!\n\
                         Syntax: python convert_to_single_phase_journeys.py source_emtgopt [source_emtg]\n")
                         
    source_emtgopt = sys.argv[1]
    
    originalOptions = []

    #try to ingest the old mission
    try:
        originalOptions = MissionOptions.MissionOptions(source_emtgopt)
    except:
        print(source_emtgopt, ' not found. Exiting')
        raise SystemExit

    if len(sys.argv) == 3:
        source_emtg = sys.argv[2]
        
        try:
            originalMission = Mission.Mission(source_emtg)
        except:
            print(source_emtg, ' not found. Exiting')
            raise SystemExit

        #seed the old mission from the old results file otherwise the seeder won't work
        originalOptions.trialX = []
        for Xindex in range(0, len(originalMission.DecisionVector)):
            originalOptions.trialX.append([originalMission.Xdescriptions[Xindex].replace('\n','').replace('\r',''), originalMission.DecisionVector[Xindex]])

        originalOptions.DisassembleMasterDecisionVector()
        originalOptions.ConvertDecisionVector()
        originalOptions.AssembleMasterDecisionVector()

    newOptions = convert_to_single_phase_journeys(originalOptions)

    #write the mission
    newOptions.write_options_file(source_emtgopt.replace(originalOptions.mission_name + '.emtgopt', newOptions.mission_name + '.emtgopt'))