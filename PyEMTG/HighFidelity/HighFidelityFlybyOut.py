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

#base class for high-fidelity flyby OUT
#Jacob Englander 1-23-2018

import HighFidelityHalfFlyby
import Universe
import kepler

from math import asin, atan2, acos, cos, sin, pi, tan

class HighFidelityFlybyOut(HighFidelityHalfFlyby.HighFidelityHalfFlyby):
    def __init__(self, parent, UniversePath, originalJourney, originalMissionOptions, originalJourneyOptions, Model='kepler'):
        HighFidelityHalfFlyby.HighFidelityHalfFlyby.__init__(self, parent, UniversePath, originalJourney, originalMissionOptions, originalJourneyOptions, Model)

    def createEvent(self):
        HighFidelityHalfFlyby.HighFidelityHalfFlyby.createEvent(self)

        #periapse state is the first periapse of the current journey
        self.periapseState = self.originalJourney.flyby_periapse_states[0]

        self.newJourney.journey_name = self.originalJourney.missionevents[0].Location + 'GAo'
        self.newJourney.journey_central_body = self.originalJourney.missionevents[0].Location
        self.myUniverse = Universe.Universe(self.UniversePath + '/' + self.newJourney.journey_central_body + '.emtg_universe')
        
        #set the left boundary - free point feeding from the previous journey's ephemeris-referenced arrival
        self.newJourney.departure_class = 1 #free point
        self.newJourney.destination_list[0] = -1 #really doesn't matter, choose -1 by convention
        self.newJourney.departure_type = 2 #"free direct departure"
        self.newJourney.AllowJourneyFreePointDepartureToPropagate = 0
        self.newJourney.wait_time_bounds = [0.0, 0.0]
        
        #set the right boundary
        self.newJourney.arrival_class = 2 #ephemeris-referenced
        self.newJourney.destination_list[1] = -1 #SOI
        self.newJourney.arrival_type = 2 #intercept
        self.newJourney.final_velocity = [1.0e-8, 50.0, 0.0]#lower, upper, null
        self.newJourney.arrival_ellipsoid_axes = [self.myUniverse.r_SOI]*3

        #set the v-infinity
        self.vMAG = self.originalJourney.missionevents[0].C3 ** 0.5

        self.computeTOF()

        #set up the coast phase match point
        if self.Model == 'sundman':
            self.newJourney.CoastPhaseMatchPointFraction = 0.5
        else:
            self.newJourney.CoastPhaseMatchPointFraction = 0.1

        self.newJourney.CoastPhaseForwardIntegrationStepLength = 60.0
        self.newJourney.CoastPhaseBackwardIntegrationStepLength = 600.0

    def ComputeSphericalState(self):
        pass
        
    def CreateInitialGuess(self, originalInitialGuess):
        outgoingVinfinity = []
        
        for Xindex in range(0, len(originalInitialGuess)):
            if 'UnpoweredFlyby: V_infinity_x' in originalInitialGuess[Xindex][0] \
                and self.originalPrefix in originalInitialGuess[Xindex][0]:

                outgoingVinfinity = [float(originalInitialGuess[Xindex][1]), \
                                     float(originalInitialGuess[Xindex + 1][1]), \
                                     float(originalInitialGuess[Xindex + 2][1])]
                break

        
        
        #before we do anything else, we need to look through the OLD trialX and find the outgoing v-infinity vector
        vMAG = (outgoingVinfinity[0]**2 + outgoingVinfinity[1]**2 + outgoingVinfinity[2]**2)**0.5
        vRA = atan2(outgoingVinfinity[1], outgoingVinfinity[0])
        vDEC = asin(outgoingVinfinity[2] / vMAG)

        
        
        #generate the point at infinity
        elements = kepler.cart2kep(r=self.periapseState[0:3], v=self.periapseState[3:6], mu=self.myUniverse.mu)
        elements[5] = self.DeltaTA
        Rinfinity_star, Vinfinity_star = kepler.coe2rv(elements, self.myUniverse.mu)

        RA_SOI = atan2(Rinfinity_star[1], Rinfinity_star[0])
        DEC_SOI = asin(Rinfinity_star[2] / self.myUniverse.r_SOI)

        prefix = self.prefix
        if self.Model in ['kepler', 'integrated']:
            prefix = prefix + 'CoastPhase'
        elif self.Model == 'sundman':
            prefix = prefix + 'SundmanCoastPhase'

        self.InitialGuess.append([prefix + 'FreePointFreeDirectDeparture: event left state mass', self.originalJourney.missionevents[0].Mass])
        self.InitialGuess.append([prefix + ': phase flight time', self.TOF / 86400.0])
        if self.Model == 'sundman':
            self.InitialGuess.append([prefix + ': phase Sundman independent variable', self.DeltaTA])
            
        self.InitialGuess.append([prefix + ': virtual chemical fuel', 0.0])
        self.InitialGuess.append([prefix + 'EphemerisReferencedInterceptInterior: event interface state vMAG', vMAG])
        self.InitialGuess.append([prefix + 'EphemerisReferencedInterceptInterior: event interface state vRA', vRA])
        self.InitialGuess.append([prefix + 'EphemerisReferencedInterceptInterior: event interface state vDEC', vDEC])
        self.InitialGuess.append([prefix + 'EphemerisReferencedInterceptInterior: event interface state RA', RA_SOI])
        self.InitialGuess.append([prefix + 'EphemerisReferencedInterceptInterior: event interface state DEC', DEC_SOI])
        self.InitialGuess.append([prefix + 'EphemerisReferencedInterceptInterior: event left state mass', self.originalJourney.missionevents[0].Mass])

        return self.InitialGuess