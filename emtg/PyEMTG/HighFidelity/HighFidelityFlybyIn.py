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

#base class for high-fidelity flyby IN
#Jacob Englander 1-23-2018

import HighFidelityHalfFlyby
import Universe

import numpy
from math import asin, atan2, acos, cos, sin, pi, tan

class HighFidelityFlybyIn(HighFidelityHalfFlyby.HighFidelityHalfFlyby):
    def __init__(self, parent, UniversePath, originalJourney, originalMissionOptions, originalJourneyOptions, nextJourney, nextJourneyOptions, Model='kepler'):
        HighFidelityHalfFlyby.HighFidelityHalfFlyby.__init__(self, parent, UniversePath, originalJourney, originalMissionOptions, originalJourneyOptions, Model)

        self.nextJourney = nextJourney
        self.nextJourneyOptions = nextJourneyOptions

    def createEvent(self):
        HighFidelityHalfFlyby.HighFidelityHalfFlyby.createEvent(self)

        #periapse state is the first periapse of the NEXT journey
        self.periapseState = self.nextJourney.flyby_periapse_states[0]

        self.newJourney.journey_name = self.nextJourney.missionevents[0].Location + 'GAi'
        self.newJourney.journey_central_body = self.nextJourney.missionevents[0].Location
        self.myUniverse = Universe.Universe(self.UniversePath + '/' + self.newJourney.journey_central_body + '.emtg_universe')
        
        #set the left boundary - free point feeding from the previous journey's ephemeris-referenced arrival
        self.newJourney.departure_class = 1 #free point
        self.newJourney.destination_list[0] = -1 #doesn't matter, -1 by convention
        self.newJourney.departure_type = 2 #"free direct departure"
        self.newJourney.AllowJourneyFreePointDepartureToPropagate = 0
        self.newJourney.wait_time_bounds = [0.0, 0.0]
        
        #set the right boundary
        self.newJourney.arrival_class = 3 #periapse
        self.newJourney.destination_list[1] = -1 #doesn't matter, -1 by convention
        self.newJourney.arrival_type = 2 #intercept
        self.newJourney.arrival_elements_type = 0 
        self.newJourney.arrival_elements_frame = 1#J2000 BCI
        self.newJourney.arrival_elements = self.periapseState
        self.newJourney.arrival_elements_reference_epoch = self.nextJourney.missionevents[0].JulianDate

        #set the v-infinity
        self.vMAG = self.nextJourney.missionevents[0].C3 ** 0.5

        #set the altitude constraint
        self.newJourney.PeriapseArrival_override_altitude = 1

        if self.nextJourneyOptions.override_flyby_altitude_bounds:
            self.newJourney.PeriapseArrival_altitude_bounds = self.nextJourneyOptions.flyby_altitude_bounds
        else:
            self.newJourney.PeriapseArrival_altitude_bounds = [self.myUniverse.minimum_safe_distance - self.myUniverse.central_body_radius, 1.0e+6]

        self.computeTOF()

        #set up the coast phase match point
        if self.Model == 'sundman':
            self.newJourney.CoastPhaseMatchPointFraction = 0.5
        else:
            self.newJourney.CoastPhaseMatchPointFraction = 0.9
        self.newJourney.CoastPhaseForwardIntegrationStepLength = 600.0
        self.newJourney.CoastPhaseBackwardIntegrationStepLength = 60.0

    def CreateInitialGuess(self, originalInitialGuess):
        incomingVinfinity = []
        
        for Xindex in range(0, len(originalInitialGuess)):
            if 'Intercept: V_infinity_x' in originalInitialGuess[Xindex][0] \
                and self.originalPrefix in originalInitialGuess[Xindex][0]:

                incomingVinfinity = [float(originalInitialGuess[Xindex][1]), \
                                     float(originalInitialGuess[Xindex + 1][1]), \
                                     float(originalInitialGuess[Xindex + 2][1])]
                break
        #get the cartesian state, then convert it
        x = self.periapseState[0]
        y = self.periapseState[1]
        z = self.periapseState[2]
        vx = self.periapseState[3]
        vy = self.periapseState[4]
        vz = self.periapseState[5]
        
        prefix = self.prefix
        if self.Model in ['kepler', 'integrated']:
            prefix = prefix + 'CoastPhase'
        elif self.Model == 'sundman':
            prefix = prefix + 'SundmanCoastPhase'

        self.InitialGuess.append([prefix + 'FreePointFreeDirectDeparture: event left state mass', self.nextJourney.missionevents[0].Mass])
        self.InitialGuess.append([prefix + ': phase flight time', self.TOF / 86400.0])
        if self.Model == 'sundman':
            self.InitialGuess.append([prefix + ': phase Sundman independent variable', self.DeltaTA])
            
        self.InitialGuess.append([prefix + ': virtual chemical fuel', 0.0])
        self.InitialGuess.append([prefix + 'PeriapseFlybyIn: event left state x', x])
        self.InitialGuess.append([prefix + 'PeriapseFlybyIn: event left state y', y])
        self.InitialGuess.append([prefix + 'PeriapseFlybyIn: event left state z', z])
        self.InitialGuess.append([prefix + 'PeriapseFlybyIn: event left state vx', vx])
        self.InitialGuess.append([prefix + 'PeriapseFlybyIn: event left state vy', vy])
        self.InitialGuess.append([prefix + 'PeriapseFlybyIn: event left state vz', vz])
        self.InitialGuess.append([prefix + 'PeriapseFlybyIn: event left state mass', self.nextJourney.missionevents[0].Mass])

        #convert the initial guess from Cartesian to the state representation that the user wants
        from StateConverter import StateConverter
        myStateConverter = StateConverter()

        stateRepresentationNames = ["Cartesian", "SphericalRADEC", "SphericalAZFPA", "COE", "MEE", "IncomingBplane", "OutgoingBplane"]

        self.InitialGuess = myStateConverter.convertDecisionVector(self.InitialGuess,                                                                  
                                                                   stateRepresentationNames[self.originalMissionOptions.PeriapseBoundaryStateRepresentation],                 
                                                                   ["PeriapseFlybyIn"],                                                                  
                                                                   self.myUniverse.mu)    

        return self.InitialGuess