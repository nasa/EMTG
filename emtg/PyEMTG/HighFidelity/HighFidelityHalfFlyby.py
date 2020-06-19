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

#base class for high-fidelity flyby
#Jacob Englander 1-23-2018

import JourneyOptions
import Journey
import Universe

from math import sqrt, acosh, sinh, atan, tanh, acos

class HighFidelityHalfFlyby(object):
    def __init__(self, parent, UniversePath, originalJourney, originalMissionOptions, originalJourneyOptions, Model='kepler'):
        self.parent = parent
        self.UniversePath = UniversePath
        self.originalJourney = originalJourney
        self.originalMissionOptions = originalMissionOptions
        self.originalJourneyOptions = originalJourneyOptions
        self.Model = Model.lower()

        self.InitialGuess = []
        self.newJourney = []
        
        self.periapseState = []

        self.TOF = 0.0
        self.vMAG = 0.0

        self.originalJourneyIndex = self.originalJourney.journey_number
        self.originalPrefix = 'j' + str(self.originalJourneyIndex) + 'p0'
        
    def setJourneyIndex(self, journeyIndex):
        self.journeyIndex = journeyIndex
        self.prefix = 'j' + str(self.journeyIndex) + 'p0'
        journeyIndex += 1

        return journeyIndex

    def computeTOF(self):        
        mu = self.myUniverse.mu
        r = (self.periapseState[0] ** 2 + self.periapseState[1] ** 2 + self.periapseState[2] ** 2) ** 0.5

        rSOI = self.myUniverse.r_SOI
        SMA = -mu / self.vMAG / self.vMAG
        ECC = 1 - r / SMA
        H = acosh(1.0 / ECC * (rSOI / -SMA + 1.0))
        N = ECC * sinh(H) - H
        self.TOF = N / sqrt(mu / -(SMA*SMA*SMA))
        
        self.DeltaTA = 2.0 * atan(((ECC + 1.0) / (ECC - 1.0))**0.5 * tanh(0.5 * H))

        
        DeltaTA2 = acos(1.0 / ECC * (SMA * (1 - ECC**2) / rSOI - 1.0))
        #print(DeltaTA2 - self.DeltaTA)

    def getTOF(self):
        return self.TOF

    def getDeltaTA(self):
        return self.DeltaTA

    def getJourneyOptions(self):
        return self.newJourney

    def createEvent(self):
        if self.Model in ['kepler', 'integrated']:
            self.newJourney = JourneyOptions.JourneyOptions()
            self.newJourney.phase_type = 7 # CoastPhase
        elif self.Model == 'sundman':
            self.newJourney = JourneyOptions.JourneyOptions()
            self.newJourney.phase_type = 8 # SundmanCoastPhase
    def CreateInitialGuess(originalInitialGuess):
        pass

    def ComputeSphericalState(self):
        pass