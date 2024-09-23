#EMTG: Evolutionary Mission Trajectory Generator
#An open-source global optimization tool for preliminary mission design
#Provided by NASA Goddard Space Flight Center
#
#Copyright (c) 2014 - 2018 United States Government as represented by the
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

#script to convert a two-point shooting low-thrust (MGALT or FBLT) mission to PSFB
#this code assumes that all journeys have a single phase - otherwise bad things happen!
#idealy this is for use with a "high-fidelity" .emtgopt and corresponding .emtg

import sys
sys.path.append("c:/emtg/PyEMTG/")
import Mission
import MissionOptions
import Universe
import copy

from math import atan2, acos, asin, cos, sin, pi
import numpy

class MissionConverter_TPSLT_to_PSFB(object):
    
    def __init__(self, originalMissionPath, originalOptionsPath, OriginalTranscription):
        self.originalMission = Mission.Mission(originalMissionPath)
        self.originalOptions = MissionOptions.MissionOptions(originalOptionsPath)

        self.newOptions = []
        self.HighFidelityJourneys = []
        
        Transcriptions = {'MGALTS': 0,
                          'FBLTS': 1,
                          'MGALT': 2,
                          'FBLT': 3,
                          'PSBI': 4,
                          'PSFB': 5,
                          'MGAnDSMs': 6,
                          'CoastPhase': 7}
        self.OriginalTranscription = OriginalTranscription#[OriginalTranscription,Transcriptions[OriginalTranscription[0]]]

        
        self.NewTranscription = ['PSFB', 5]

    def getMissionOptions(self):
        return self.newOptions    

    def CreateMission(self):
        self.newOptions = copy.deepcopy(self.originalOptions)

        for journeyIndex in range(0, len(self.newOptions.Journeys)):
            if self.newOptions.Journeys[journeyIndex].phase_type == self.OriginalTranscription[1]:
                self.newOptions.Journeys[journeyIndex].phase_type = self.NewTranscription[1]

    def CreateInitialGuess(self):
        #Step 1: construct the old initial guess
        DecisionVector = []
        for Xindex in range(0, len(self.originalMission.DecisionVector)):
            DecisionVector.append([self.originalMission.Xdescriptions[Xindex].replace('\n','').replace('\r',''), self.originalMission.DecisionVector[Xindex]])

        controlDV = []

        #Step 2: first convert all mentions of the original transcription to PSFB
        newDecisionVector = []
        for entry in DecisionVector:
            line = copy.deepcopy(entry)

            if self.OriginalTranscription[0] not in line[0]:
                newDecisionVector.append(line)
                continue

            line = [line[0].replace(self.OriginalTranscription[0], 'PSFB'), line[1]]
            
            journeyIndex = int(line[0][1:line[0].find('p')])
            stepIndex = []

            if 'u_' in line[0]:
                stepIndex = int(line[0][line[0].find('step'):line[0].find('u')].strip('step '))
                line[0] = line[0].replace('PSFB: step ' + str(stepIndex),'PSFB_Step' + str(stepIndex) + ': substep0')
            
                controlDV.append(line)
            else:
                newDecisionVector.append(line)
            
            #find the missionevent offset
            originalJourney = self.originalMission.Journeys[journeyIndex]
            originalJourneyOptions = self.originalOptions.Journeys[journeyIndex]
            eventStartIndex = []
            eventMatchPointIndex = []
            for event in originalJourney.missionevents:
                if event.EventType in ['coast','Scoast','SFthrust','FBLTthrust','SSFthrust','FBLTSthrust']:
                    eventStartIndex = event.EventNumber
                    break
            for event in originalJourney.missionevents:
                if event.EventType in ['match_point']:
                    eventMatchPointIndex = event.EventNumber
                    break

            originalUniverse = Universe.Universe(self.originalOptions.universe_folder + '/' + originalJourneyOptions.journey_central_body + '.emtg_universe')

            
            delta0 = float(originalUniverse.reference_angles[2]) * pi / 180.0
            alpha0 = float(originalUniverse.reference_angles[0]) * pi / 180.0

            Rx = numpy.matrix([[1.0, 0.0, 0.0],
                  [0.0, cos(pi/2.0 - delta0), -sin(pi/2 -delta0)],
                  [0.0, sin(pi/2.0 - delta0), cos(pi/2 - delta0)]])
            
            Rz = numpy.matrix([[cos(pi/2 + alpha0), -sin(pi/2 + alpha0), 0.0],
                               [sin(pi/2 + alpha0), cos(pi/2 + alpha0), 0.0],
                               [0.0, 0.0, 1.0]])

            R = Rz * Rx

            #insert step states after u_z
            if 'Step' in line[0] and 'u_z' in line[0]:
                eventIndex = eventStartIndex + stepIndex
                if eventIndex >= eventMatchPointIndex:
                    eventIndex += 1

                event = []
                for oldEvent in originalJourney.missionevents:
                    if oldEvent.EventNumber == eventIndex:
                        event = copy.deepcopy(oldEvent)              
                        break

                #create the initial guess (note: rotation may be needed later...)
                position = numpy.matrix([event.SpacecraftState[0],
                                         event.SpacecraftState[1],
                                         event.SpacecraftState[2]]).T
                velocity = numpy.matrix([event.SpacecraftState[3],
                                         event.SpacecraftState[4],
                                         event.SpacecraftState[5]]).T

                Rposition = R * position
                Rvelocity = R * velocity

                x = float(Rposition[0])
                y = float(Rposition[1])
                z = float(Rposition[2])
                xdot = float(Rvelocity[0])
                ydot = float(Rvelocity[1])
                zdot = float(Rvelocity[2])

                prefix = 'j' + str(journeyIndex) + 'p0PSFB_Step' + str(stepIndex) + ': '

                if self.originalOptions.ParallelShootingStateRepresentation in [1, 2]:
                    r = (x**2 + y**2 + z**2)**0.5
                    v = (xdot**2 + ydot**2 + zdot**2)**0.5
                    RA = atan2(y, x)
                    DEC = asin(z / r)
                    FPA = acos( (x*xdot + y*ydot + z*zdot) / r / v )
        
                    #azimuth is complicated
                    xhat = numpy.matrix([cos(RA)*cos(DEC), sin(RA)*cos(DEC), sin(DEC)]).T
                    yhat = numpy.matrix([cos(RA + pi / 2.0), sin(RA + pi / 2.0), 0.0]).T
                    zhat = numpy.matrix([-cos(RA)*sin(DEC), -sin(RA)*sin(DEC), cos(DEC)]).T
                    R = numpy.hstack([xhat, yhat, zhat]).T
                    V = numpy.matrix([xdot, ydot, zdot]).T
                    Vprime = R * V
                    AZ = atan2(Vprime[1], Vprime[2])
                  
                    #vRA and vDEC
                    vRA = atan2(ydot, xdot)                                                                                 
                    vDEC = asin(zdot / v)   

                    newDecisionVector.append([prefix + 'left state r', r])
                    newDecisionVector.append([prefix + 'left state RA', RA])
                    newDecisionVector.append([prefix + 'left state DEC', DEC])
                    newDecisionVector.append([prefix + 'left state v', v])
                    if self.originalOptions.ParallelShootingStateRepresentation == 2: #spherical AZFPA
                        newDecisionVector.append([prefix + 'left state AZ', AZ])
                        newDecisionVector.append([prefix + 'left state FPA', FPA])
                    else: #spherical RADEC
                        newDecisionVector.append([prefix + 'left state vRA', vRA])
                        newDecisionVector.append([prefix + 'left state vDEC', vDEC])
                else: #cartesian! this is easier...
                    newDecisionVector.append([prefix + 'left state x', x])
                    newDecisionVector.append([prefix + 'left state y', y])
                    newDecisionVector.append([prefix + 'left state z', z])
                    newDecisionVector.append([prefix + 'left state vx', xdot])
                    newDecisionVector.append([prefix + 'left state vy', ydot])
                    newDecisionVector.append([prefix + 'left state vz', zdot])
                newDecisionVector.append([prefix + 'left state mass', event.Mass])
                newDecisionVector.append([prefix + 'virtual chemical fuel', 0.0])
                newDecisionVector.append([prefix + 'virtual electric propellant', 0.0])
                
                newDecisionVector += controlDV
                
                controlDV = []
            
        return newDecisionVector
