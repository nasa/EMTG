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

#base class for high-fidelity launch
#Jacob Englander 2-8-2018

import JourneyOptions
import Journey
import Universe

from math import asin, atan, atan2, acos, cos, sin, pi, tan, tanh, acosh, sinh, sqrt, fabs
from scipy.optimize import minimize_scalar
import numpy as np

class HighFidelityLaunch(object):
    def __init__(self, parent, UniversePath, originalJourney, originalMissionOptions, originalJourneyOptions, Model='kepler'):
        self.parent = parent
        self.UniversePath = UniversePath
        self.originalMissionOptions = originalMissionOptions
        self.originalJourney = originalJourney
        self.originalJourneyOptions = originalJourneyOptions
        self.Model = Model.lower()

        self.InitialGuess = []
        self.newJourney = []
        
        self.periapseState = []

        self.C3 = 0.0
        self.TOF = 0.0
        self.PeriapseManeuverMagnitude = 0.0

        self.originalJourneyIndex = self.originalJourney.journey_number
        self.originalPrefix = 'j' + str(self.originalJourneyIndex) + 'p0'



    def createEvent(self):
        if self.Model in ['kepler', 'integrated']:
            self.newJourney = JourneyOptions.JourneyOptions()
            self.newJourney.phase_type = 7 # CoastPhase
        elif self.Model == 'sundman':
            self.newJourney = JourneyOptions.JourneyOptions()
            self.newJourney.phase_type = 8 # SundmanCoastPhase

        self.newJourney.journey_name = self.originalJourney.missionevents[0].Location + 'Launch'
        self.newJourney.journey_central_body = self.originalJourney.missionevents[0].Location
        self.myUniverse = Universe.Universe(self.UniversePath + '/' + self.newJourney.journey_central_body + '.emtg_universe')
        
        #set the left boundary - free point feeding from the previous journey's ephemeris-referenced arrival
        self.newJourney.departure_class = 3 #periapse
        self.newJourney.destination_list[0] = -1 #really doesn't matter, choose -1 by convention
        self.newJourney.departure_type = 0 #"launch or impulsive maneuver"
        self.newJourney.AllowJourneyFreePointDepartureToPropagate = 0
        self.newJourney.wait_time_bounds = self.originalJourneyOptions.wait_time_bounds
        
        #set the right boundary
        self.newJourney.arrival_class = 2 #ephemeris-referenced
        self.newJourney.destination_list[1] = -1 #SOI
        self.newJourney.arrival_type = 2 #intercept
        self.newJourney.final_velocity = [1.0e-8, 50.0, 0.0]#lower, upper, null
        self.newJourney.arrival_ellipsoid_axes = [self.myUniverse.r_SOI]*3

        # set the C3
        self.C3 = self.originalJourney.missionevents[0].C3

        if self.C3 < 1.0e-07:
            self.C3 = 1.0e-07

        #set the v-infinity
        self.vMAG = self.C3 ** 0.5

        self.computeTOF()

        #set up the coast phase match point
        if self.Model == 'sundman':
            self.newJourney.CoastPhaseMatchPointFraction = 0.5
        else:
            self.newJourney.CoastPhaseMatchPointFraction = 0.1

        self.newJourney.CoastPhaseForwardIntegrationStepLength = 60.0
        self.newJourney.CoastPhaseBackwardIntegrationStepLength = 600.0

    def setJourneyIndex(self, journeyIndex):
        self.journeyIndex = journeyIndex
        self.prefix = 'j' + str(self.journeyIndex) + 'p0'
        journeyIndex += 1

        return journeyIndex

    def computeTOFparabolic(self):
        # compute the TOF along a parabolic trajectory
        mu = self.myUniverse.mu
        r_periapse = self.myUniverse.central_body_radius + self.originalJourneyOptions.PeriapseDeparture_altitude_bounds[0]
        rSOI = self.myUniverse.r_SOI

        self.DeltaTA = np.arccos(2.0 * r_periapse / rSOI - 1.0)

        TA_initial = 0.0 # injection from periapse
        TA_final = TA_initial + self.DeltaTA
        D0 = np.tan(TA_initial / 2.0)
        Df = np.tan(TA_final / 2.0)
    
        self.TOF = np.sqrt(2.0 * r_periapse**3 / mu) * (Df + (1.0/3.0) * Df**3 - D0 - (1.0/3.0) * D0**3)

    def computeTOF(self):        
        mu = self.myUniverse.mu
        r = self.myUniverse.central_body_radius + self.originalJourneyOptions.PeriapseDeparture_altitude_bounds[0]

        rSOI = self.myUniverse.r_SOI
        SMA = -mu / self.vMAG / self.vMAG
        ECC = 1 - r / SMA
        H = acosh(1.0 / ECC * (rSOI / -SMA + 1.0))
        N = ECC * sinh(H) - H
        self.TOF = N / sqrt(mu / -(SMA*SMA*SMA))
        
        #TA that we passed through since periapse
        self.DeltaTA = 2.0 * atan(((ECC + 1.0) / (ECC - 1.0))**0.5 * tanh(0.5 * H))

        #what is the time since periapse passage at halfway through the propagation?
        self.HalfwayTime = (ECC - sinh(H / 2.0) - H/2.0) / sqrt(mu / -(SMA*SMA*SMA))
        
        

    def getTOF(self):
        return self.TOF

    def getJourneyOptions(self):
        return self.newJourney
        
    def CreateInitialGuess(self, originalInitialGuess):
        import kepler
        #hard code inclination of parking orbit for now - but in the long run it needs to come out of the launch vehicle file
        #this should include a decision on prograde vs retrograde launch
        direction = 'prograde'
        incPark = 28.5 * pi / 180.0
        mu = self.myUniverse.mu

        vRA = self.originalJourney.missionevents[0].RightAscension * pi / 180.0
        vDEC = self.originalJourney.missionevents[0].Declination * pi / 180.0
        
        rPark = self.myUniverse.central_body_radius + self.originalJourneyOptions.PeriapseDeparture_altitude_bounds[0]
        vPark = sqrt(self.myUniverse.mu / rPark)

        v_periapse = sqrt(self.originalJourney.missionevents[0].C3 + 2.0 * self.myUniverse.mu / rPark)
        v_escape = sqrt(2.0 * self.myUniverse.mu / rPark)

        deltav = v_periapse - vPark

        if (self.originalMissionOptions.LaunchVehicleKey == 'Fixed_Initial_mass'):
            self.newJourney.initial_impulse_bounds = [1e-8, sqrt(self.originalJourney.missionevents[0].C3)]      
        else:
            self.newJourney.initial_impulse_bounds = [1e-8, 20.0]      

        #compute both state-at-infinity and parking orbit state from 
        #periapsis velocity vector [km/s], inclination
        INC = incPark
        SMA = -mu / self.vMAG / self.vMAG
        ECC = 1 - rPark / SMA

        #so now we have SMA, ECC, INC, TA
        #we need to find the RAAN, AOP that combine with these and give the right terminal velocity vector
        
        vRA = self.originalJourney.missionevents[0].RightAscension * pi / 180.0
        vDEC = self.originalJourney.missionevents[0].Declination * pi / 180.0

        Vinfinity = np.array([self.C3 ** 0.5 * cos(vRA) * cos(vDEC),
                              self.C3 ** 0.5 * sin(vRA) * cos(vDEC),
                              self.C3 ** 0.5 * sin(vDEC)])

        def computeVinfinityError(inputs):
            import kepler
            RAAN = inputs[0]
            AOP = inputs[1]
            elements = [SMA, ECC, INC, RAAN, AOP, self.DeltaTA]
            Rinfinity_star, Vinfinity_star = kepler.coe2rv(elements, mu)

            VinfinityErrorVector = Vinfinity - Vinfinity_star

            return np.linalg.norm(VinfinityErrorVector)

        from scipy.optimize import minimize
        result = minimize(fun=computeVinfinityError, x0=[0.0, 0.0], tol=1.0e-8)#, options={'disp': True, 'gtol': 1.0e-6}
        
        RAAN = result.x[0]
        AOP = result.x[1]


        #compute the states at infinity
        import StateConverter

        myStateConverter = StateConverter.StateConverter()

        elements = [SMA, ECC, INC, RAAN, AOP, self.DeltaTA]
        [r_SOI, RA_SOI, DEC_SOI, v_SOI, vRA_SOI, vDEC_SOI] = myStateConverter.COEtoSphericalRADEC(elements, mu)
        

        [VINF, RHA, DHA, BRADIUS, BTHETA, TA] = myStateConverter.COEtoOutgoingBplane([SMA, ECC, INC, RAAN, AOP, 0.0], mu)

        prefix = self.prefix
        if self.Model in ['kepler', 'integrated']:
            prefix = prefix + 'CoastPhase'
        elif self.Model == 'sundman':
            prefix = prefix + 'SundmanCoastPhase'

        self.InitialGuess.append([prefix + 'PeriapseLaunch: event left state epoch', self.originalJourney.missionevents[0].JulianDate - 2400000.5])
        self.InitialGuess.append([prefix + 'PeriapseLaunch: event left state mass', self.originalJourney.missionevents[0].Mass])
        self.InitialGuess.append([prefix + 'PeriapseLaunch: event left state VINFout', VINF])
        self.InitialGuess.append([prefix + 'PeriapseLaunch: event left state RHAout', RHA])
        self.InitialGuess.append([prefix + 'PeriapseLaunch: event left state DHAout', DHA])
        self.InitialGuess.append([prefix + 'PeriapseLaunch: event left state BRADIUSout', BRADIUS])
        self.InitialGuess.append([prefix + 'PeriapseLaunch: event left state BTHETAout', BTHETA])
        self.InitialGuess.append([prefix + 'PeriapseLaunch: event left state TAout', TA])

        self.InitialGuess.append([prefix + ': phase flight time', self.TOF / 86400.0])
        if self.Model == 'sundman':
            self.InitialGuess.append([prefix + ': phase Sundman independent variable', self.DeltaTA])
        
        self.InitialGuess.append([prefix + ': virtual chemical fuel', 0.0])
        self.InitialGuess.append([prefix + 'EphemerisReferencedInterceptInterior: event interface state vMAG', self.vMAG])
        self.InitialGuess.append([prefix + 'EphemerisReferencedInterceptInterior: event interface state vRA', vRA])
        self.InitialGuess.append([prefix + 'EphemerisReferencedInterceptInterior: event interface state vDEC', vDEC])
        self.InitialGuess.append([prefix + 'EphemerisReferencedInterceptInterior: event interface state RA', RA_SOI])
        self.InitialGuess.append([prefix + 'EphemerisReferencedInterceptInterior: event interface state DEC', DEC_SOI])
        self.InitialGuess.append([prefix + 'EphemerisReferencedInterceptInterior: event left state mass', self.originalJourney.missionevents[0].Mass])

        return self.InitialGuess