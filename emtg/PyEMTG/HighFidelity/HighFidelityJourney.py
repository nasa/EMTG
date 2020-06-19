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

#base class for high-fidelity journey
#Jacob Englander 1-23-2018

import sys
sys.path.append("c:/EMTG/PyEMTG")
import JourneyOptions
import Journey

import HighFidelityLaunch
import HighFidelityFlybyOut
import HighFidelityFlybyIn
import kepler

import copy
from math import atan2, asin, pi, acos

class HighFidelityJourney(object):
    def __init__(self, parent, UniversePath, originalMissionOptions, originalJourney, originalJourneyOptions, nextJourney=None, nextJourneyOptions=None, Transcription=['MGALT', 2], useSundman=False):
        self.parent = parent
        self.UniversePath = UniversePath
        self.originalMissionOptions = originalMissionOptions
        self.originalJourney = originalJourney
        self.originalJourneyOptions = originalJourneyOptions
        self.nextJourney = nextJourney
        self.nextJourneyOptions = nextJourneyOptions

        self.newJourneyOptions = []
        self.newJourneyOptionsSet = []
        self.InitialGuessSet = [] #vector of guess vectors
        self.journeyInitialGuess = []

        self.originalJourneyIndex = self.originalJourney.journey_number
        self.originalPrefix = 'j' + str(self.originalJourneyIndex) + 'p0'

        self.FlybyModel = []

        # Adopt the main mission transcription type, unless we are an existing coast phase
        if self.originalJourneyOptions.phase_type == 7:
            self.Transcription = ['CoastPhase', 7]
        else:
            self.Transcription = Transcription

        if Transcription[0] in ['MGALT', 'MGALTS'] or (Transcription[0] in ['MGAnDSMs'] and self.originalMissionOptions.propagatorType == 0):
            self.FlybyModel = 'kepler'
        elif useSundman == True:
            self.FlybyModel = 'sundman'
        else:
            self.FlybyModel = 'integrated'

        self.DepartureEvent = None
        
        if self.originalJourney.missionevents[0].EventType == 'launch' and self.originalJourney.missionevents[0].Location != 'periapse':
            self.DepartureEvent = HighFidelityLaunch.HighFidelityLaunch(self,
                                                                        self.UniversePath,
                                                                        self.originalJourney,
                                                                        self.originalMissionOptions,
                                                                        self.originalJourneyOptions,
                                                                        self.FlybyModel)
        elif self.originalJourney.missionevents[0].EventType == 'upwr_flyby':
            self.DepartureEvent = HighFidelityFlybyOut.HighFidelityFlybyOut(self,
                                                                            self.UniversePath,
                                                                            self.originalJourney, 
                                                                            self.originalMissionOptions,
                                                                            self.originalJourneyOptions,
                                                                            self.FlybyModel)


        self.ArrivalEvent = None
        
        if self.originalJourney.missionevents[-1].EventType == 'intercept':
            if self.nextJourney != None: #i.e., there exists a next journey
                if self.nextJourney.missionevents[0].EventType in ['upwr_flyby']:
                #if self.nextJourney.missionevents[0].EventType in ['upwr_flyby', 'zeroflyby']:
                    self.ArrivalEvent = HighFidelityFlybyIn.HighFidelityFlybyIn(self,
                                                                                self.UniversePath,
                                                                                self.originalJourney,
                                                                                self.originalMissionOptions,
                                                                                self.originalJourneyOptions,
                                                                                self.nextJourney,
                                                                                self.nextJourneyOptions,
                                                                                self.FlybyModel)


    def CreateJourney(self):
        self.newJourneyOptionsSet = []

        #create a new journey, has all of the properties of the original
        self.newJourneyOptions = copy.deepcopy(self.originalJourneyOptions)
        self.newJourneyOptions.phase_type = self.Transcription[1]

        #does the journey begin with a 3D event?
        if not self.DepartureEvent == None:
            self.DepartureEvent.createEvent()

            #refactor the current journey to start with a free point free direct departure
            self.newJourneyOptions.departure_class = 1 #free point
            self.newJourneyOptions.destination_list[0] = -1
            self.newJourneyOptions.departure_type = 2 #"free direct departure"
            self.newJourneyOptions.AllowJourneyFreePointDepartureToPropagate = 0
            self.newJourneyOptions.wait_time_bounds = [0.0, 0.0]

            #adjust the forced initial coast
            if self.newJourneyOptions.forced_initial_coast > 0.0:
                self.newJourneyOptions.forced_initial_coast -= self.DepartureEvent.getTOF() / 86400.0

                if self.newJourneyOptions.forced_initial_coast < 0.0:
                    self.newJourneyOptions.forced_initial_coast = 0.0

        #does the journey end with a 3D flyby?
        if not self.ArrivalEvent == None:

            #create the half-flyby
            self.ArrivalEvent.createEvent()

            #refactor the current journey to start with an ephemeris-referenced intercept
            self.newJourneyOptions.arrival_class = 2 #free point
            #don't change the destination list - it already tells us whose SOI we are going to
            self.newJourneyOptions.arrival_type = 2 #"intercept"
            #set the SOI boundary
            self.newJourneyOptions.arrival_ellipsoid_axes = [self.ArrivalEvent.myUniverse.r_SOI]*3

            #adjust the forced terminal coast
            if self.newJourneyOptions.forced_terminal_coast > 0.0:
                self.newJourneyOptions.forced_terminal_coast -= self.ArrivalEvent.getTOF() / 86400.0

                if self.newJourneyOptions.forced_terminal_coast < 0.0:
                    self.newJourneyOptions.forced_terminal_coast = 0.0


        #create the options set
        if not self.DepartureEvent == None:
            self.newJourneyOptionsSet.append(self.DepartureEvent.getJourneyOptions())

        self.newJourneyOptionsSet.append(self.newJourneyOptions)

        if not self.ArrivalEvent == None:
            self.newJourneyOptionsSet.append(self.ArrivalEvent.getJourneyOptions())

        return self.newJourneyOptionsSet

    def setJourneyIndex(self, journeyIndex):
        if not self.DepartureEvent == None:
            journeyIndex = self.DepartureEvent.setJourneyIndex(journeyIndex)

        self.journeyIndex = journeyIndex
        self.prefix = 'j' + str(self.journeyIndex) + 'p0'
        journeyIndex += 1

        if not self.ArrivalEvent == None:
            journeyIndex = self.ArrivalEvent.setJourneyIndex(journeyIndex)

        return journeyIndex

    def CreateInitialGuess(self, originalInitialGuess):
        self.InitialGuessSet = [] #vector of guess vectors
        self.journeyInitialGuess = []

        #departure flyby?
        if not self.DepartureEvent == None:
            self.InitialGuessSet.append(self.DepartureEvent.CreateInitialGuess(originalInitialGuess))

        #main journey
        #********************************************************************************
        #TODO we have to rename the boundary guesses appropriately!
        #********************************************************************************
        for entry in originalInitialGuess:
            if self.originalPrefix in entry[0]:
                #TODO departure stuff

                #main phase stuff
                newEntry = copy.deepcopy(entry)
                newEntry[0] = newEntry[0].replace(self.originalPrefix, self.prefix)

                #adjust the flight time guess after removing the flight time for each flyby
                if 'phase flight time' in entry[0]:
                    if not self.DepartureEvent == None:
                        newEntry[1] -= self.DepartureEvent.getTOF() / 86400.0
                    if not self.ArrivalEvent == None:
                        newEntry[1] -= self.ArrivalEvent.getTOF() / 86400.0

                self.journeyInitialGuess.append(newEntry)

                #TODO arrival stuff

        #if the journey starts with a 3D flyby, we need to encode left mass
        if not self.DepartureEvent == None:
            self.journeyInitialGuess.append([self.prefix + self.Transcription[0] + 'FreePointFreeDirectDeparture: event left state mass', self.originalJourney.missionevents[0].Mass])

        #if the journey ends in a 3D flyby, we need to encode the EphemerisReferencedInterceptExterior
        if not self.ArrivalEvent == None:
            #first we need to construct the v-infinity vector at the end of the journey
            
            for Xindex in range(0, len(originalInitialGuess)):
                if 'Intercept: V_infinity_x' in originalInitialGuess[Xindex][0] \
                    and self.originalPrefix in originalInitialGuess[Xindex][0]:

                    incomingVinfinity = [float(originalInitialGuess[Xindex][1]), \
                                         float(originalInitialGuess[Xindex + 1][1]), \
                                         float(originalInitialGuess[Xindex + 2][1])]
                    break

            vMAG = (incomingVinfinity[0]**2 + incomingVinfinity[1]**2 + incomingVinfinity[2]**2)**0.5
            vRA = atan2(incomingVinfinity[1], incomingVinfinity[0])
            vDEC = asin(incomingVinfinity[2] / vMAG)

            #generate the point at infinity
            elements = kepler.cart2kep(r=self.ArrivalEvent.periapseState[0:3], v=self.ArrivalEvent.periapseState[3:6], mu=self.ArrivalEvent.myUniverse.mu)
            elements[5] = -self.ArrivalEvent.getDeltaTA()
            
            Rinfinity_star, Vinfinity_star = kepler.coe2rv(elements, mu=self.ArrivalEvent.myUniverse.mu)
            RA_SOI = atan2(Rinfinity_star[1], Rinfinity_star[0])
            DEC_SOI = asin(Rinfinity_star[2] / self.ArrivalEvent.myUniverse.r_SOI)

            self.journeyInitialGuess.append([self.prefix + self.Transcription[0] + 'EphemerisReferencedInterceptExterior: event interface state vMAG', vMAG])
            self.journeyInitialGuess.append([self.prefix + self.Transcription[0] + 'EphemerisReferencedInterceptExterior: event interface state vRA', vRA])
            self.journeyInitialGuess.append([self.prefix + self.Transcription[0] + 'EphemerisReferencedInterceptExterior: event interface state vDEC', vDEC])
            self.journeyInitialGuess.append([self.prefix + self.Transcription[0] + 'EphemerisReferencedInterceptExterior: event interface state RA', RA_SOI])
            self.journeyInitialGuess.append([self.prefix + self.Transcription[0] + 'EphemerisReferencedInterceptExterior: event interface state DEC', DEC_SOI])
            self.journeyInitialGuess.append([self.prefix + self.Transcription[0] + 'EphemerisReferencedInterceptExterior: event left state mass', self.originalJourney.missionevents[-1].Mass])
            self.journeyInitialGuess.append([self.prefix + self.Transcription[0] + ': virtual chemical fuel', 0.0])


        #add the guess to the stack
        self.InitialGuessSet.append(self.journeyInitialGuess)

        #arrival flyby?
        if not self.ArrivalEvent == None:
            self.InitialGuessSet.append(self.ArrivalEvent.CreateInitialGuess(originalInitialGuess))

        return self.InitialGuessSet

    def createConstraints(self, constraintList):
        newConstraintList = []
        for constraint in constraintList:
            if self.originalPrefix in constraint:
                newConstraint = constraint.replace(self.originalPrefix, self.prefix)
                newConstraintList.append(newConstraint)

        return newConstraintList



def computeEphemerisReferencedGuess(periapse_state, SOI_state, rSOI, mu_flyby, incoming_flyby):
    import numpy as np
    vMAG = (SOI_state[3]**2 + SOI_state[4]**2 + SOI_state[5]**2)**0.5
    vRA = atan2(SOI_state[4], SOI_state[3])
    vDEC = asin(SOI_state[5] / vMAG)
    
    r_periapse = (periapse_state[0] ** 2 + periapse_state[1] ** 2 + periapse_state[2] ** 2) ** 0.5
    
    SMA = -mu_flyby / vMAG / vMAG
    ECC = 1.0 - r_periapse / SMA
    H = np.arccosh(1.0 / ECC * (rSOI / -SMA + 1.0))
    N = ECC * np.sinh(H) - H
    TOF = N / np.sqrt(mu_flyby / -(SMA*SMA*SMA))
        
    DeltaTA = 2.0 * np.arctan(((ECC + 1.0) / (ECC - 1.0))**0.5 * np.tanh(0.5 * H))

    #generate the point at infinity
    elements = kepler.cart2kep(r=periapse_state[0:3], v=periapse_state[3:6], mu=mu_flyby)

    if incoming_flyby:
        elements[5] = -DeltaTA
    else:
        elements[5] = DeltaTA
            
    Rinfinity_star, Vinfinity_star = kepler.coe2rv(elements, mu=mu_flyby)
    
    # instead of using Rinfinity_star, let's just take the actual position data directly from the
    # acceleration model output
    #RA_shell = atan2(Rinfinity_star[1], Rinfinity_star[0])
    #DEC_shell = asin(Rinfinity_star[2] / rSOI)
    RA_shell = atan2(SOI_state[1], SOI_state[0])
    DEC_shell = asin(SOI_state[2] / rSOI)

    return vMAG, vRA, vDEC, RA_shell, DEC_shell

def computePeriapseGuess(periapse_state):
    import numpy as np
    r = np.sqrt(periapse_state[0]**2 + periapse_state[1]**2 + periapse_state[2]**2)
    v = np.sqrt(periapse_state[3]**2 + periapse_state[4]**2 + periapse_state[5]**2)
    RA = np.arctan2(periapse_state[1], periapse_state[0])
    DEC = np.arcsin(periapse_state[2] / r)
    vRA = np.arctan2(periapse_state[4], periapse_state[3])
    vDEC = np.arcsin(periapse_state[5] / v)

    return r, v, RA, DEC, vRA, vDEC


if __name__ == '__main__':    
    
    # spacecraft state in the frame of the flyby body
    SOI_in_state   = [ -361828.0502141565084457,
                       -162347.8509645648300648,
                       471368.8083165474236012,
                       1.2216846701685835,
                       0.4520179826917605,
                       -1.4335821899866534,
                        1247.6492,
                        2461868.7890625]
    periapse_state = [-10885.98899246,
                       1163.66576399,
                       7083.78495716,
                       -2.07946164,
                       5.59034390,
                       -4.11394311,
                       1258.4936,
                       2461652.06074018]
    SOI_out_state  = [355275.7003868818283081,
                      192525.2492393213324249,
                      -464943.7380269346758723,
                      1.2212572322499553,
                      0.4905238333513218,
                      -1.4157222865481440,
                      1258.4936,
                      2461655.140763889]

    rSOI = 616000 # Callisto SOI in km
    mu_flyby = 324858.592079 # Callisto mu
    
    vMAG, vRA, vDEC, RA_shell, DEC_shell = computeEphemerisReferencedGuess(periapse_state, SOI_in_state, rSOI, mu_flyby, True)
    prefix = 'p0'
    transcription = 'MGAnDSMs'
    incoming_SOI_guess = []
    incoming_SOI_guess.append([prefix + transcription + 'EphemerisReferencedInterceptExterior: event interface state vMAG', vMAG])
    incoming_SOI_guess.append([prefix + transcription + 'EphemerisReferencedInterceptExterior: event interface state vRA', vRA])
    incoming_SOI_guess.append([prefix + transcription + 'EphemerisReferencedInterceptExterior: event interface state vDEC', vDEC])
    incoming_SOI_guess.append([prefix + transcription + 'EphemerisReferencedInterceptExterior: event interface state RA', RA_shell])
    incoming_SOI_guess.append([prefix + transcription + 'EphemerisReferencedInterceptExterior: event interface state DEC', DEC_shell])
    incoming_SOI_guess.append([prefix + transcription + 'EphemerisReferencedInterceptExterior: event left state mass', SOI_in_state[6]])
    
    r, v, RA, DEC, vRA, vDEC = computePeriapseGuess(periapse_state)
    transcription = 'CoastPhase'
    periapse_state_guess = []
    periapse_state_guess.append([prefix + transcription + ': phase flight time', periapse_state[7] - SOI_in_state[7]])
    periapse_state_guess.append([prefix + transcription + 'PeriapseFlybyIn: event left state r', r])
    periapse_state_guess.append([prefix + transcription + 'PeriapseFlybyIn: event left state RA', RA])
    periapse_state_guess.append([prefix + transcription + 'PeriapseFlybyIn: event left state DEC', DEC])
    periapse_state_guess.append([prefix + transcription + 'PeriapseFlybyIn: event left state v', v])
    periapse_state_guess.append([prefix + transcription + 'PeriapseFlybyIn: event left state vRA', vRA])
    periapse_state_guess.append([prefix + transcription + 'PeriapseFlybyIn: event left state vDEC', vDEC])
    periapse_state_guess.append([prefix + transcription + 'PeriapseFlybyIn: event left state mass', periapse_state[6]])
    periapse_state_guess.append([prefix + transcription + ': virtual chemical fuel', 0.0])

    vMAG, vRA, vDEC, RA_shell, DEC_shell = computeEphemerisReferencedGuess(periapse_state, SOI_out_state, rSOI, mu_flyby, False)
    outgoing_SOI_guess = []
    outgoing_SOI_guess.append([prefix + transcription + ': phase flight time', SOI_out_state[7] - periapse_state[7]])
    outgoing_SOI_guess.append([prefix + transcription + 'EphemerisReferencedInterceptInterior: event interface state vMAG', vMAG])
    outgoing_SOI_guess.append([prefix + transcription + 'EphemerisReferencedInterceptInterior: event interface state vRA', vRA])
    outgoing_SOI_guess.append([prefix + transcription + 'EphemerisReferencedInterceptInterior: event interface state vDEC', vDEC])
    outgoing_SOI_guess.append([prefix + transcription + 'EphemerisReferencedInterceptInterior: event interface state RA', RA_shell])
    outgoing_SOI_guess.append([prefix + transcription + 'EphemerisReferencedInterceptInterior: event interface state DEC', DEC_shell])
    outgoing_SOI_guess.append([prefix + transcription + 'EphemerisReferencedInterceptInterior: event left state mass', periapse_state[6]])
    outgoing_SOI_guess.append([prefix + transcription + ': virtual chemical fuel', 0.0])

    for entry in incoming_SOI_guess:
        print(entry[0]+','+str(entry[1]))

    print('\n')
    for entry in periapse_state_guess:
        print(entry[0]+','+str(entry[1]))

    print('\n')
    for entry in outgoing_SOI_guess:
        print(entry[0]+','+str(entry[1]))

    print("yay")