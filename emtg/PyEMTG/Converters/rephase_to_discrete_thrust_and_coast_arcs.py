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

#script to convert from regular MGALT/FBLT with the opportunity for partial thrust
#to individual journeys that are either thrust full on or thrust full off
#uses CoastPhase whenever the spacecraft was not thrusting more than 10%
#
#could later be updated to work with PSFB/PSBI. But this is more work so I am not doing it right now.
#
#assumes single-phase journeys. Run convert_to_single_phase_journeys() before running this script

import sys
sys.path.append("c:/emtg/PyEMTG/")
import Mission
import MissionOptions
import Universe
import copy
import numpy as np
from math import fabs


class phaseBlock(object):
    def __init__(self, oldJourney):
        self.phaseType = 'CoastPhase'
        self.beginState = [0.0]*6
        self.beginEpoch = 0.0
        self.beginMass = 0.0
        self.endState = [0.0]*6
        self.endMass = 0.0
        self.endEpoch = 0.0
        self.phaseFlightTime = 0.0
        self.controlVector = [] #becomes a vector of 3-vectors
        self.state_frame = oldJourney.state_frame
        self.alpha0 = float(oldJourney.alpha0)
        self.delta0 = float(oldJourney.delta0)
        
    def Rx(self, angle):
        return np.array([[1.0, 0.0          ,  0.0           ],
                         [0.0, np.cos(angle), -np.sin(angle)],
                         [0.0, np.sin(angle),  np.cos(angle) ]])
                        
    def Ry(self, angle):
        return np.array([[ np.cos(angle), 0.0, np.sin(angle)],
                         [ 0.0          , 1.0, 0.0          ],
                         [-np.sin(angle), 0.0, np.cos(angle)]])
                        
    def Rz(self, angle):
        return np.array([[np.cos(angle), -np.sin(angle), 0.0],
                         [np.sin(angle),  np.cos(angle), 0.0],
                         [0.0          ,  0.0          , 1.0]])
    
    
    def rotate_from_local_to_ICRF(self, three_vector):
        if self.state_frame == 'ICRF':
            return three_vector
        else:
            #construct rotation matrix
            Rx = self.Rx(np.pi / 2 - self.delta0)
            Rz = self.Rz(np.pi / 2 + self.alpha0)
            R_from_local_to_ICRF = Rz.dot(Rx)
                        
            return R_from_local_to_ICRF.dot(three_vector.transpose())
        
    def makeJourney(self, baseJourney, phaseIndex = 0, isFirst = False, isLast = False, mission_time_step_size = 86400.0):
        #Step 1: copy the base journey
        newJourney = copy.deepcopy(baseJourney)
        newJourney.trialX = []
        phaseIndex = phaseIndex + 1
        
        basePhaseType = []
        if baseJourney.phase_type == 2:
            basePhaseType = 'MGALT'
        elif baseJourney.phase_type == 3:
            basePhaseType = 'FBLT'
        elif baseJourney.phase_type == 4:
            basePhaseType = 'PSBI'
        elif baseJourney.phase_type == 5:
            basePhaseType = 'PSFB'
       
        
        #Step 2: what is the phase type?
        if self.phaseType == 'CoastPhase':
            newJourney.phase_type = 7
        elif self.phaseType == 'MGALT':
            newJourney.phase_type = 2
            newJourney.force_unit_magnitude_control = 1
        elif self.phaseType == 'FBLT':
            newJourney.phase_type = 3
            newJourney.force_unit_magnitude_control = 1
        elif self.phaseType == 'PSBI':
            newJourney.phase_type = 4
            newJourney.force_unit_magnitude_control = 1
        elif self.phaseType == 'PSFB':
            newJourney.phase_type = 5
            newJourney.force_unit_magnitude_control = 1
            
        newJourney.journey_name = newJourney.journey_name + '_' + self.phaseType + '_' + str(phaseIndex)
        
        #Step 3: boundary conditions and their initial guesses
        
        #Step 3.1: departure
        #if this is the first thrust or coast phase in the new journey, then it keeps the original left boundary condition
        #otherwise it's a free point free direct departure
        if isFirst:
            #no changes to the departure
            #copy the initial guess from the existing journey
            for entry in baseJourney.trialX:
                if 'EphemerisPeggedLaunchDirectInsertion' in entry[0] \
                    or 'EphemerisPeggedFreeDirectDeparture' in entry[0] \
                    or 'EphemerisPeggedZeroTurnFlyby' in entry[0] \
                    or 'EphemerisPeggedUnpoweredFlyby' in entry[0] \
                    or 'EphemerisPeggedPoweredFlyby' in entry[0] \
                    or 'EphemerisPeggedSpiralDeparture' in entry[0] \
                    or 'FreePointDirectInsertion' in entry[0] \
                    or 'FreePointFreeDirectDeparture' in entry[0] \
                    or 'EphemerisReferencedFreeDirectDepartureExterior' in entry[0] \
                    or 'EphemerisReferencedFreeDirectDepartureInterior' in entry[0] \
                    or 'PeriapseLaunchOrImpulsiveDeparture' in entry[0]:
                    
                    #this decision variable is attached to the journey departure event. Add it to the new journey
                    newJourney.trialX.append([entry[0].replace(basePhaseType, self.phaseType), entry[1]])
        else:
            newJourney.departure_class = 1 #free point
            newJourney.departure_type = 2 #free direct departure
            newJourney.destination_list[0] = -1
            newJourney.AllowJourneyFreePointDepartureToPropagate = 0
            newJourney.wait_time_bounds = [0.0, 0.0]
            
            #initial guess for mass
            newJourney.trialX.append(['p0' + self.phaseType + 'FreePointFreeDirectDeparture: event left state mass', self.beginMass])
            
        #Step 3.2: arrival
        #if this is the last thrust or coast phase in the new journey, then it keeps the original right boundary condition
        #otherwise it's a free point rendezvous with no maneuver
        if isLast:
            #no changes to the arrival
            #copy the initial guess from the existing journey
            for entry in baseJourney.trialX:
                if 'EphemerisPeggedLTRendezvous' in entry[0] \
                    or 'EphemerisPeggedFlybyIn' in entry[0] \
                    or 'EphemerisPeggedIntercept' in entry[0] \
                    or 'EphemerisPeggedChemRendezvous' in entry[0] \
                    or 'EphemerisPeggedOrbitInsertion' in entry[0] \
                    or 'EphemerisPeggedMomentumTransfer' in entry[0] \
                    or 'EphemerisPeggedSpiralArrival' in entry[0] \
                    or 'FreePointLTRendezvous' in entry[0] \
                    or 'FreePointChemRendezvous' in entry[0] \
                    or 'FreePointIntercept' in entry[0] \
                    or 'EphemerisReferencedLTRendezvousExterior' in entry[0] \
                    or 'EphemerisReferencedInterceptExterior' in entry[0] \
                    or 'EphemerisReferencedLTRendezvousInterior' in entry[0] \
                    or 'EphemerisReferencedInterceptInterior' in entry[0] \
                    or 'PeriapseFlybyIn' in entry[0]:
                    
                    #this decision variable is attached to the journey arrival event. Add it to the new journey
                    newJourney.trialX.append([entry[0].replace(basePhaseType, self.phaseType), entry[1]])
        else:
            newJourney.arrival_class = 1 #free point
            newJourney.arrival_type = 3 #rendezvous with no maneuver
            newJourney.destination_list[1] = -1
            newJourney.AllowJourneyFreePointArrivalToPropagate = 0
            newJourney.arrival_elements_state_representation = 0 #Cartesian
            newJourney.arrival_elements_vary_flag = [1, 1, 1, 1, 1, 1]
            
            rotR = self.rotate_from_local_to_ICRF(np.array(self.endState[0:3]))
            rotV = self.rotate_from_local_to_ICRF(np.array(self.endState[3:6]))
            
            newJourney.arrival_elements_bounds = [-2.0 * fabs(max(rotR[0:3], key=abs)), 2.0 * fabs(max(rotR[0:3], key=abs)),\
                                                  -2.0 * fabs(max(rotR[0:3], key=abs)), 2.0 * fabs(max(rotR[0:3], key=abs)),\
                                                  -2.0 * fabs(max(rotR[0:3], key=abs)), 2.0 * fabs(max(rotR[0:3], key=abs)),\
                                                  -2.0 * fabs(max(rotV[0:3], key=abs)), 2.0 * fabs(max(rotV[0:3], key=abs)),\
                                                  -2.0 * fabs(max(rotV[0:3], key=abs)), 2.0 * fabs(max(rotV[0:3], key=abs)),\
                                                  -2.0 * fabs(max(rotV[0:3], key=abs)), 2.0 * fabs(max(rotV[0:3], key=abs))]
                                                  
            #inital guess for the arrival condition
            newJourney.trialX.append(['p0' + self.phaseType + 'FreePointLTRendezvous: event left state x', rotR[0]])
            newJourney.trialX.append(['p0' + self.phaseType + 'FreePointLTRendezvous: event left state y', rotR[1]])
            newJourney.trialX.append(['p0' + self.phaseType + 'FreePointLTRendezvous: event left state z', rotR[2]])
            newJourney.trialX.append(['p0' + self.phaseType + 'FreePointLTRendezvous: event left state vx', rotV[0]])
            newJourney.trialX.append(['p0' + self.phaseType + 'FreePointLTRendezvous: event left state vy', rotV[1]])
            newJourney.trialX.append(['p0' + self.phaseType + 'FreePointLTRendezvous: event left state vz', rotV[2]])
            newJourney.trialX.append(['p0' + self.phaseType + 'FreePointLTRendezvous: event left state mass', self.endMass])
        
    
        #Step 4: flight time initial guess
        newJourney.trialX.append(['p0' + self.phaseType + ': phase flight time', self.phaseFlightTime])
        
        #Step 5: virtual propellant tank initial guess
        newJourney.trialX.append(['p0' + self.phaseType + ': virtual chemical fuel', 0.0])
        if self.phaseType in ['MGALT','FBLT','PSBI','PSFB']:
            newJourney.trialX.append(['p0' + self.phaseType + ': virtual electric propellant', self.beginMass - self.endMass])
        
        #only do these steps if this is a thrust phase of some kind
        if self.phaseType in ['MGALT','FBLT','PSBI','PSFB']:
            #Step 6: number of control steps
            newJourney.override_num_steps = 1
            newJourney.number_of_steps = len(self.controlVector)
            
            #Step 7: control initial guess
            #there are separate versions for MGALT/FBLT and PSBI/PSFB because of interior control points and states
            if self.phaseType in ['MGALT','FBLT']:
                stepIndex = 0
                for controlVectorRow in self.controlVector:
                    rotControl = self.rotate_from_local_to_ICRF(np.array(controlVectorRow))
                    
                    newJourney.trialX.append(['p0' + self.phaseType + ': step ' + str(stepIndex) + ' u_x', rotControl[0]])
                    newJourney.trialX.append(['p0' + self.phaseType + ': step ' + str(stepIndex) + ' u_y', rotControl[1]])
                    newJourney.trialX.append(['p0' + self.phaseType + ': step ' + str(stepIndex) + ' u_z', rotControl[2]])
                    stepIndex += 1
                    
                #if we have an odd number of control entries, copy the last one
                if len(self.controlVector) % 2:
                    rotControl = self.rotate_from_local_to_ICRF(np.array(self.controlVector[-1]))
                    newJourney.number_of_steps += 1
                    newJourney.trialX.append(['p0' + self.phaseType + ': step ' + str(stepIndex) + ' u_x', rotControl[0]])
                    newJourney.trialX.append(['p0' + self.phaseType + ': step ' + str(stepIndex) + ' u_y', rotControl[1]])
                    newJourney.trialX.append(['p0' + self.phaseType + ': step ' + str(stepIndex) + ' u_z', rotControl[2]])
                    stepIndex += 1
                    
            
            if self.phaseType in ['PSBI','PSFB']:
                print('this thing doesn\'t make initial guesses for parallel shooting phase types yet')
                
                
        #step 8: remove any forced coasts
        newJourney.forced_initial_coast = 0.0
        newJourney.forced_terminal_coast = 0.0
        
        #step 9: journey time bounds and post-journey delta-v
        if isLast:
            if newJourney.timebounded == 1:# if the journey has time bounds, as opposed to arrival date bounds or aggregate time bounds
                newJourney.timebounded = 0
        else: #if this isn't the last phase of the journey, turn off all time bounds and post-journey delta-v
            newJourney.timebounded = 0
            newJourney.journey_end_TCM = 0.0
            newJourney.journey_end_deltav = 0.0
            
        #step 10: integration time step size for coast phases - has its own control fields
        if self.phaseType == 'CoastPhase':
            if newJourney.override_integration_step_size:
                newJourney.CoastPhaseForwardIntegrationStepLength = newJourney.integration_step_size
                newJourney.CoastPhaseBackwardIntegrationStepLength = newJourney.integration_step_size
            else:
                newJourney.CoastPhaseForwardIntegrationStepLength = mission_time_step_size
                newJourney.CoastPhaseBackwardIntegrationStepLength = mission_time_step_size
                
        #Step 11: staging
        if not isLast:
            newJourney.stage_before_arrival = 0
            newJourney.stage_after_arrival = 0
        if not isFirst:
            newJourney.stage_after_departure = 0

        #step 12: please don't print so many free point boundaries in the target spec file...
        #if the old journey did not end with a free point, turn off the boundary target printer
        if baseJourney.arrival_class != 1:
            newJourney.FreePointArrival_print_target_spec = 0
            
        return newJourney, phaseIndex
        

def rephase_to_discrete_thrust_and_coast_arcs(originalOptions, originalMission):

    #copy the old options file
    newOptions = copy.deepcopy(originalOptions)

    #clear the Journeys list
    newOptions.Journeys = []

    #walk through the journeys in the old options file and split them if needed
    for oldJourneyIndex in range(0, len(originalOptions.Journeys)):
        oldJourneyOptions = originalOptions.Journeys[oldJourneyIndex]
        oldJourney = originalMission.Journeys[oldJourneyIndex]
        
        #Step 1: is this journey *not* an MGALT or FBLT journey?
        if oldJourneyOptions.phase_type not in [2, 3]:
            newOptions.Journeys.append(copy.deepcopy(oldJourneyOptions))
            continue
            
        #Step 2: if we made it this far then there is thrusting. We'll need to break it up into discrete thrust arcs
        #let's scan through the journey and look for thrust events.
        #Make a basket of minimal thrust event information
        phaseStack = []
        
        #there is always a first phase
        currentPhase = phaseBlock(oldJourney)
        currentPhase.beginState = oldJourney.missionevents[0].SpacecraftState
        currentPhase.beginEpoch = oldJourney.missionevents[0].JulianDate
        currentPhase.beginMass = oldJourney.missionevents[0].Mass
        
        #does the first phase start with a coast or a thrust? set the phase type appropriately
        firstEvent = oldJourney.missionevents[1]
        if 'SFthrust' in firstEvent.EventType and np.linalg.norm(np.array(firstEvent.Thrust) / firstEvent.AvailableThrust) > 0.5:
            currentPhase.phaseType = 'MGALT'
        elif 'FBLTthrust' in firstEvent.EventType and np.linalg.norm(np.array(firstEvent.Thrust) / firstEvent.AvailableThrust) > 0.5:
            currentPhase.phaseType = 'FBLT'
        elif 'PSBIthrust' in firstEvent.EventType and np.linalg.norm(np.array(firstEvent.Thrust) / firstEvent.AvailableThrust) > 0.5:
            currentPhase.phaseType = 'PSBI'
        elif 'PSFBthrust' in firstEvent.EventType and np.linalg.norm(np.array(firstEvent.Thrust) / firstEvent.AvailableThrust) > 0.5:
            currentPhase.phaseType = 'PSFB'
        #the default case is already CoastPhase
        
        #scan through the journey. We do this by index instead of by list comprehension because we need to be able to look backward in time
        for eventIndex in range(1, len(oldJourney.missionevents) - 1): #we know that the last event is an endpoint
            #make an alias to the current event
            event = oldJourney.missionevents[eventIndex]
            
            #4 cases:
            #case 1: this is a thrust event and we are in a thrust phase - add the control vector and keep going
            #case 2: this is a coast event and we are in a thrust phase - stop the thrust phase and save the end point. Start a coast phase
            #case 3: this is a thrust event and we are in a coast phase - stop the coast phase and save the end point. Start a thrust phase and save the control vector
            #case 4: this is a coast event and we are in a coast phase. Continue without doing anything
            #
            #<now begin the cases>
            
            #case 1: this is a thrust event and we are in a thrust phase - add the control vector and keep going
            if 'thrust' in event.EventType and currentPhase.phaseType in ['MGALT','FBLT','PSBI','PSFB'] \
                and np.linalg.norm(np.array(event.Thrust) / event.AvailableThrust) > 0.5:
                
                currentPhase.controlVector.append(np.array(event.Thrust) / event.AvailableThrust) #note that dutyCycle is already built into AvailableThrust
            
            #case 2: this is a coast event and we are in a thrust phase - stop the thrust phase and save the end point. Start a coast phase
            #set the final mass and epoch based on the end of the thrust phase
            if ('coast' in event.EventType or ('thrust' in event.EventType and np.linalg.norm(np.array(event.Thrust) / event.AvailableThrust) < 0.5)) \
                and currentPhase.phaseType in ['MGALT','FBLT','PSBI','PSFB']:
                
                #alias to the previous event - this is safe because you can't get here in the first event of the mission
                previousEvent = oldJourney.missionevents[eventIndex - 1]
                
                #finish off the current phase
                currentPhase.endEpoch = previousEvent.JulianDate + previousEvent.TimestepLength / 2
                currentPhase.endMass = previousEvent.Mass
                currentPhase.phaseFlightTime = currentPhase.endEpoch - currentPhase.beginEpoch
                currentPhase.endState = previousEvent.SpacecraftState
                
                #add the current phase to the stack
                phaseStack.append(copy.deepcopy(currentPhase))
                
                #reset the current phase to a coast phase (which is the default setting anyway)
                currentPhase = phaseBlock(oldJourney)
                currentPhase.beginEpoch = event.JulianDate - event.TimestepLength / 2
                currentPhase.beginMass = event.Mass
                currentPhase.beginState = event.SpacecraftState
            
            #case 3: this is a thrust event and we are in a coast phase - stop the coast phase and save the end point. Start a thrust phase and save the control vector
            if 'thrust' in event.EventType and currentPhase.phaseType == 'CoastPhase':
                if np.linalg.norm(np.array(event.Thrust) / event.AvailableThrust) > 0.5:                      
                    #alias to the previous event - this is safe because you can't get here in the first event of the mission
                    previousEvent = oldJourney.missionevents[eventIndex - 1]
                    
                    #finish off the current phase
                    currentPhase.endEpoch = previousEvent.JulianDate + previousEvent.TimestepLength / 2
                    currentPhase.endMass = previousEvent.Mass
                    currentPhase.phaseFlightTime = currentPhase.endEpoch - currentPhase.beginEpoch
                    currentPhase.endState = previousEvent.SpacecraftState
                    
                    #add the current phase to the stack
                    phaseStack.append(copy.deepcopy(currentPhase))
                    
                    #reset the current phase to a coast phase - but which type of coast phase?
                    currentPhase = phaseBlock(oldJourney)
                    currentPhase.beginEpoch = event.JulianDate - event.TimestepLength / 2
                    currentPhase.beginMass = event.Mass
                    currentPhase.beginState = event.SpacecraftState
                    
                    if 'SFthrust' in event.EventType:
                        currentPhase.phaseType = 'MGALT'
                    elif 'FBLTthrust' in event.EventType:
                        currentPhase.phaseType = 'FBLT'
                    elif 'PSBIthrust' in event.EventType:
                        currentPhase.phaseType = 'PSBI'
                    elif 'PSFBthrust' in event.EventType:
                        currentPhase.phaseType = 'PSFB'
                        
                    #add the current control vector
                    currentPhase.controlVector.append(np.array(event.Thrust) / event.AvailableThrust) #note that dutyCycle is already built into AvailableThrust
                
            #case 4: this is a coast event and we are in a coast phase. Continue without doing anything
            #no need for any code here!
        
        #Step 3: finish of the last phase and add it to the stack
        lastEvent = oldJourney.missionevents[-1]
        currentPhase.endEpoch = lastEvent.JulianDate
        currentPhase.endMass = lastEvent.Mass
        currentPhase.phaseFlightTime = currentPhase.endEpoch - currentPhase.beginEpoch
        currentPhase.endState = lastEvent.SpacecraftState
        
        phaseStack.append(copy.deepcopy(currentPhase))
        
        #Step 4: we now have a stack of phase objects. Let's make them into journeys.
        phaseIndex = 0
        #Step 4.1: the first coast/thrust phase of the journey (which may also be the last
        firstJourney, phaseIndex = phaseStack[0].makeJourney(baseJourney = oldJourneyOptions, phaseIndex = phaseIndex, isFirst = True, mission_time_step_size = originalOptions.integration_time_step_size)
        newOptions.Journeys.append(firstJourney)
        
        #Step 4.2: intermediate coast/thrust phases (if the journey has more than two)
        if len(phaseStack) > 2:
            for phase in phaseStack[1:-1]:
                thisJourney, phaseIndex = phase.makeJourney(baseJourney = oldJourneyOptions, phaseIndex = phaseIndex, mission_time_step_size = originalOptions.integration_time_step_size)
                newOptions.Journeys.append(thisJourney)
        
        #Step 4.3: the last coast/thrust phase of the journey (if there was more than one)
        if len(phaseStack) > 1:
            lastJourney, phaseIndex = phaseStack[-1].makeJourney(baseJourney = oldJourneyOptions, phaseIndex = phaseIndex, isLast = True, mission_time_step_size = originalOptions.integration_time_step_size)
            newOptions.Journeys.append(lastJourney)
        

    #set the number of journeys
    newOptions.number_of_journeys = len(newOptions.Journeys)

    #assemble the master decision and constraint vectors
    newOptions.AssembleMasterDecisionVector()
    newOptions.AssembleMasterConstraintVectors()

    #rename the mission
    newOptions.mission_name = originalOptions.mission_name + '_discreteThrust'
    
    #set to variable phase type
    newOptions.mission_type = 9
    
    #remove forced coasts
    newOptions.forced_post_launch_coast = 0.0
    newOptions.forced_pre_flyby_coast = 0.0
    newOptions.forced_post_flyby_coast = 0.0

    #return the new options file
    return newOptions
    
    
if __name__ == '__main__':
    # Ensure the correct number of command line options were provided
    if len(sys.argv) < 2:
        raise Exception("Unknown number of command line options!\n\
                         Syntax: python rephase_to_discrete_thrust_and_coast_arcs.py source_emtgopt source_emtg\n")
                         
    source_emtgopt = sys.argv[1]
    source_emtg = sys.argv[2]
    
    originalOptions = []
    originalMission = []

    #try to ingest the old mission
    try:
        originalOptions = MissionOptions.MissionOptions(source_emtgopt)
    except:
        print(source_emtgopt, ' not found. Exiting')
        raise SystemExit
        
    
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

    newOptions = rephase_to_discrete_thrust_and_coast_arcs(originalOptions, originalMission)

    #write the mission
    newOptions.write_options_file(source_emtgopt.replace(originalOptions.mission_name + '.emtgopt', newOptions.mission_name + '.emtgopt'))