import MissionEvent
import EOM
import AstroFunctions
import ThrottleTable

import math
import copy
import datetime
import numpy as np
import os
import warnings

class Journey(object):
    def __init__(self, parent):
        self.parent = parent

        self.missionevents = []
        self.journey_name = "AJourney"
        self.journey_number = 0
        self.central_body = "Sun"
        self.central_body_radius = 4.379e+6
        self.mu = 132712440017.99
        self.LU = 1.49597870691e+8
        self.TU = 5022642.890912973
        self.boundary_states = []
        self.flyby_periapse_states = [[0.0]*6]
        self.thruster_duty_cycle = 1.0
        self.journey_post_arrival_deltav = 0.0
        self.journey_post_arrival_propellant_used = 0.0
        self.journey_initial_mass_increment = 0.0
        self.journey_final_mass_increment = 0.0

        #default values for entry things
        self.journey_entryinterface_velocity_with_respect_to_rotating_atmosphere = 0.0
        self.journey_entry_landing_interface_point = [0.0]*6
        self.entry_latitude = 0.0
        self.entry_longitude = 0.0
        self.entry_sun_angle = 0.0
        self.journey_arrival_spacecraft_sun_Earth_angle = 180.0
        self.state_frame = ''
        self.alpha0 = 0.0
        self.delta0 = 0.0

    def PlotJourney(self, JourneyAxes, PlotOptions, CustomLabels = None):
        #plot each mission event
        BeforeMatchPoint = True
        if PlotOptions.ShowMissionEvents:
            for event in self.missionevents:
                #do we have a custom label?
                if CustomLabels != None:
                    for label in CustomLabels: #this is a list of lists, [eventNumber, labelstring]
                        if label[0] == event.EventNumber:
                            BeforeMatchPoint = event.PlotEvent(JourneyAxes, self.LU, self.TU, self.mu, PlotOptions, BeforeMatchPoint, label[1])
                        else:
                            BeforeMatchPoint = event.PlotEvent(JourneyAxes, self.LU, self.TU, self.mu, PlotOptions, BeforeMatchPoint)
                else:
                    BeforeMatchPoint = event.PlotEvent(JourneyAxes, self.LU, self.TU, self.mu, PlotOptions, BeforeMatchPoint)

        

    def PlotJourneyBoundaryOrbits(self, JourneyAxes, PlotOptions):
        from scipy.integrate import ode

        for boundaryStateIndex in range(0, len(self.boundary_states)):
            if not ((boundaryStateIndex == 0 and self.missionevents[0].Location == 'free point') \
                or (boundaryStateIndex == (len(self.boundary_states) - 1) and self.missionevents[-1].Location == 'free point')) \
                or PlotOptions.ShowFreePointBoundaryOrbits:
                boundarystate = self.boundary_states[boundaryStateIndex]
                BoundaryStateScaled = np.array(boundarystate) / self.LU
                BoundaryStateScaled[3:6] *= self.TU
                r = np.linalg.norm(BoundaryStateScaled[0:3])
                v = np.linalg.norm(BoundaryStateScaled[3:6])

                if r == 0.0:
                    continue
                
                a = r / (2.0 - r*v*v)
                if a > 0.0:
                    T = 2*math.pi*math.sqrt(a**3)
                else:
                    T = 2*math.pi*self.LU / self.TU

                StateIntegrateObject = ode(EOM.EOM_inertial_2body, jac=EOM.EOM_jacobian_intertial_2body).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                StateIntegrateObject.set_initial_value(BoundaryStateScaled).set_f_params(1.0).set_jac_params(1.0)

                dt = T / 10000
                StateHistory = []

                epoch = (self.missionevents[0].JulianDate - 2451544.5) * 86400.0
                while StateIntegrateObject.successful() and StateIntegrateObject.t < T * 1.01:
                    epoch += dt * self.TU
                    stateArray = StateIntegrateObject.integrate(StateIntegrateObject.t + dt)
                    state = [stateArray[0] * self.LU, stateArray[1] * self.LU, stateArray[2] * self.LU, epoch]
                    StateHistory.append(state)

                X = []
                Y = []
                Z = []
                for StateLine in StateHistory:
                    if PlotOptions.PlotCentralBody != 'Journey central body':
                        import spiceypy
                        journey_central_body = self.central_body
                            
                        if journey_central_body in ['Mars','Jupiter','Uranus','Neptune','Pluto']:
                            journey_central_body = journey_central_body + ' Barycenter'

                        if journey_central_body != PlotOptions.PlotCentralBody:
                            body_state, LT = spiceypy.spkezr(journey_central_body, StateLine[-1], 'J2000', "None", PlotOptions.PlotCentralBody)
                            
                            X.append(StateLine[0] + body_state[0])
                            Y.append(StateLine[1] + body_state[1])
                            Z.append(StateLine[2] + body_state[2])
                        else:
                            X.append(StateLine[0])
                            Y.append(StateLine[1])
                            Z.append(StateLine[2])
                    else:
                        X.append(StateLine[0])
                        Y.append(StateLine[1])
                        Z.append(StateLine[2])

                JourneyAxes.plot(X, Y, Z, lw=1, c='0.75')

    def PlotJourneyCentralBodyOrbits(self, JourneyAxes, PlotOptions):
        if PlotOptions.PlotCentralBody != 'Journey central body':
            import spiceypy
            journey_central_body = self.central_body
                            
            if journey_central_body in ['Mars','Jupiter','Uranus','Neptune','Pluto']:
                journey_central_body = journey_central_body + ' Barycenter'

            if journey_central_body != PlotOptions.PlotCentralBody:
                body_state, LT = spiceypy.spkezr(journey_central_body, (self.missionevents[0].JulianDate - 2451545.0) * 86400.0, 'J2000', "None", PlotOptions.PlotCentralBody)
                
                
                r = np.linalg.norm(body_state[0:3])
                v = np.linalg.norm(body_state[3:6])

                mu = {'Mercury':                   22031.780000,
                      'Venus':                    324858.592000,
                      'Earth':                    398600.435436,
                      'Mars':                      42828.375214,
                      'Jupiter Barycenter':    126712764.800000,
                      'Saturn Barycenter':      37940585.200000,
                      'Uranus Barycenter':       5794548.600000,
                      'Neptune Barycenter':      6836527.100580,
                      'Pluto Barycenter':            977.000000,
                      'Sun':                132712440041.939400,
                      'Moon':                       4902.800066,
                      'Earth-Moon Barycenter':    403503.235502}

                elements = spiceypy.oscltx(body_state, (self.missionevents[0].JulianDate - 2451545.0) * 86400.0, mu[PlotOptions.PlotCentralBody])
                a = elements[-2]
                #a = r / (2.0 - r*v*v)
                if a > 0.0:
                    T = 2*math.pi  *math.sqrt(a**3 / mu[PlotOptions.PlotCentralBody]) 
                else:
                    T = 2*math.pi * self.TU

                t0 = (self.missionevents[0].JulianDate - 2451545.0) * 86400.0
                tf = t0 + T
                X = []
                Y = []
                Z = []
                for epoch in np.linspace(t0, tf, num=100):
                    body_state, LT = spiceypy.spkezr(journey_central_body, epoch, 'J2000', "None", PlotOptions.PlotCentralBody)
                            
                    X.append(body_state[0])
                    Y.append(body_state[1])
                    Z.append(body_state[2])

                JourneyAxes.plot(X, Y, Z, lw=1, c='0.75')

    def UpdateLabelPosition(self, Figure, Axes):
        for event in self.missionevents:
            event.UpdateLabelPosition(Figure, Axes)

    def PlotPhaseBoundariesOnDataPlot(self, DataAxes, PlotOptions, firstpass, Ybounds):

        #create vertical lines at important events
        date_string_vector = []
        boundarylegendflag = True
        burnlegendflag = True
        for event in self.missionevents:
            if event.EventType in ['upwr_flyby', 'pwr_flyby', 'LT_rndzvs', 'rendezvous', 'intercept', 'insertion', 'match-vinf', 'launch', 'departure', "begin_spiral", "end_spiral"]:
                event_epoch = datetime.datetime.strptime(event.GregorianDate,'%m/%d/%Y').date()
                if firstpass and boundarylegendflag:
                    DataAxes.plot([event_epoch]*2, Ybounds, c='k', marker='+', ls = '-.', lw=3)#, label='Phase boundary')
                    boundarylegendflag = False
                else:
                    DataAxes.plot([event_epoch]*2, Ybounds, c='k', marker='+', ls = '-.', lw=3)

            if event.EventType == 'chem_burn':
                event_epoch = datetime.datetime.strptime(event.GregorianDate,'%m/%d/%Y').date()
                if firstpass and burnlegendflag:
                    DataAxes.plot([event_epoch]*2, Ybounds, c='r', marker='+', ls = '-.', lw=3, label='Deep-Space Maneuver')
                    burnlegendflag = False
                else:
                    DataAxes.plot([event_epoch]*2, Ybounds, c='r', marker='+', ls = '-.', lw=3)


    def GenerateJourneyDataPlot(self, DataAxesLeft, DataAxesRight, PlotOptions, firstpass):
        import math
        import astropy
        #generate a vector of dates
        date_string_vector = []
        for event in self.missionevents:
            if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby' and event.EventType != 'LT_rndzvs':
                date_string_vector.append(event.GregorianDate)

        date_vector = [datetime.datetime.strptime(d,'%m/%d/%Y').date() for d in date_string_vector]

        #dummy line across the bottom so that neither axes object crashes
        DataAxesLeft.plot(date_vector, np.zeros_like(date_vector), c='w', lw = 0.1)
        DataAxesRight.plot(date_vector, np.zeros_like(date_vector), c='w', lw = 0.1)

        #plot distance from central body
        if PlotOptions.PlotR:
            Rvector = []
            for event in self.missionevents:
                if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby' and event.EventType != 'LT_rndzvs':
                    Rvector.append(math.sqrt(event.SpacecraftState[0]**2 + event.SpacecraftState[1]**2 + event.SpacecraftState[2]**2) / self.LU)
            if firstpass:
                unitstring = 'LU'
                if 'sun' in self.central_body.lower() and math.fabs(self.LU - 149597870.69100001) < 1.0:
                    unitstring = 'AU'
                labelstring = 'Distance from ' + self.central_body + ' (' + unitstring + ')'
                DataAxesLeft.plot(date_vector, Rvector, c='k', lw=2, label=labelstring)
            else:
                DataAxesLeft.plot(date_vector, Rvector, c='k', lw=2)

        #plot velocity
        if PlotOptions.PlotV:
            Vvector = []
            for event in self.missionevents:
                if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby' and event.EventType != 'LT_rndzvs':
                    Vvector.append(math.sqrt(event.SpacecraftState[3]**2 + event.SpacecraftState[4]**2 + event.SpacecraftState[5]**2) / self.LU * self.TU)
            if firstpass:
                DataAxesLeft.plot(date_vector, Vvector, c='k', lw=2, ls='-.', label='Velocity magnitude (LU/TU)')
            else:
                DataAxesLeft.plot(date_vector, Vvector, c='k', lw=2, ls='-.')

        #plot Thrust
        if PlotOptions.PlotThrust:
            Thrustvector = []
            shortDateVector = []

            for event in self.missionevents:
                epoch = astropy.time.Time(event.JulianDate, format='jd', scale='tdb')
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'PSFBthrust']:
                    Thrustvector.append(math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2) * 10.0)
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Thrustvector.append(math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2) * 10.0)
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)
                elif event.EventType == 'LT_spiral':
                    Thrustvector.append(event.AvailableThrust * self.thruster_duty_cycle*10.0)
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Thrustvector.append(event.AvailableThrust * self.thruster_duty_cycle*10.0)
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)
                elif event.EventType != 'match_point':
                    Thrustvector.append(0.0)
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Thrustvector.append(0.0)
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)


            if firstpass:
                DataAxesLeft.plot(shortDateVector, Thrustvector, c='r', lw=2, ls='--', label='Applied thrust (0.1 N)')
            else:
                DataAxesLeft.plot(shortDateVector, Thrustvector, c='r', lw=2, ls='--')

        #plot Isp
        if PlotOptions.PlotIsp:
            Ispvector = []
            shortDateVector = []

            for event in self.missionevents:
                epoch = astropy.time.Time(event.JulianDate, format='jd', scale='tdb')
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'PSFBthrust', 'LT_spiral']:
                    Ispvector.append(event.Isp / 1000.0)
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Ispvector.append(event.Isp / 1000.0)
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)
                elif event.EventType != 'match_point':
                    Ispvector.append(0.0)
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Ispvector.append(0.0)
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)

            if firstpass:
                DataAxesLeft.plot(shortDateVector, Ispvector, 'dodgerblue', lw=4, ls='-.', label='Isp (1000 s)')
            else:
                DataAxesLeft.plot(shortDateVector, Ispvector, 'dodgerblue', lw=4, ls='-.')

        #plot mass flow rate
        if PlotOptions.PlotMdot:
            Mdotvector = []
            shortDateVector = []

            for event in self.missionevents:
                epoch = astropy.time.Time(event.JulianDate, format='jd', scale='tdb')
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'PSFBthrust', 'LT_spiral']:
                    Mdotvector.append(event.MassFlowRate * 1.0e6)
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Mdotvector.append(event.MassFlowRate * 1.0e6)
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)
                elif event.EventType != 'match_point':
                    Mdotvector.append(0.0)
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Mdotvector.append(0.0)
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)

            if firstpass:
                DataAxesLeft.plot(shortDateVector, Mdotvector, c='brown', lw=2, ls='-', label='Mass flow rate (mg/s)')
            else:
                DataAxesLeft.plot(shortDateVector, Mdotvector, c='brown', lw=2, ls='-')

        #plot Efficiency
        if PlotOptions.PlotEfficiency:
            Efficiencyvector = []
            shortDateVector = []

            for event in self.missionevents:
                epoch = astropy.time.Time(event.JulianDate, format='jd', scale='tdb')
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'PSFBthrust', 'LT_spiral']:
                    Efficiencyvector.append(event.AvailableThrust * event.Isp * 9.80665 / (2000 * event.ActivePower))
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Efficiencyvector.append(event.AvailableThrust * event.Isp * 9.80665 / (2000 * event.ActivePower))
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)
                elif event.EventType != 'match_point':
                    Efficiencyvector.append(0.0)
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Efficiencyvector.append(0.0)
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)


            if firstpass:
                DataAxesLeft.plot(shortDateVector, Efficiencyvector, c='DarkGreen', lw=2, ls='-', label='Propulsion system efficiency')
            else:
                DataAxesLeft.plot(shortDateVector, Efficiencyvector, c='DarkGreen', lw=2, ls='-')

        #plot control magnitude
        if PlotOptions.PlotThrottle:
            Throttlevector = []
            shortDateVector = []

            for event in self.missionevents:
                epoch = astropy.time.Time(event.JulianDate, format='jd', scale='tdb')
                if event.EventType in ['SFthrust', 'SSFthrust', 'PSBIthrust']:
                    Throttlevector.append(math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2) / (event.AvailableThrust))
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Throttlevector.append(math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2) / (event.AvailableThrust))
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)
                elif event.EventType in ['FBLTthrust','FBLTSthrust','PSFBthrust']:
                    Throttlevector.append(event.DVmagorThrottle)
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Throttlevector.append(event.DVmagorThrottle)
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)
                elif event.EventType == 'LT_spiral':
                    Throttlevector.append(1.0)
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Throttlevector.append(1.0)
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)
                elif event.EventType != 'match_point':
                    Throttlevector.append(0.0)
                    shortDateVector.append((epoch - event.TimestepLength / 2).datetime)
                    Throttlevector.append(0.0)
                    shortDateVector.append((epoch + event.TimestepLength / 2).datetime)

            if firstpass:
                DataAxesLeft.plot(shortDateVector, Throttlevector, c='r', lw=2, ls='--', label='Control magnitude')
            else:
                DataAxesLeft.plot(shortDateVector, Throttlevector, c='r', lw=2, ls='--')

        #plot power
        if PlotOptions.PlotPower:
            Powervector = []
            for event in self.missionevents:
                if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby' and event.EventType != 'LT_rndzvs':
                    if event.AvailablePower > 0.0:
                        Powervector.append(event.AvailablePower)
                    else:
                        Powervector.append(0.0)
            if firstpass:
                DataAxesLeft.plot(date_vector, Powervector, c='Navy', lw=2, ls='-.', label='Power available for propulsion (kW)')
            else:
                DataAxesLeft.plot(date_vector, Powervector, c='Navy', lw=2, ls='-.')

        #plot gamma
        if PlotOptions.PlotGamma:
            gammavector = []
            shortDateVector = []

            for event in self.missionevents:
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'PSFBthrust', 'LT_spiral']:
                    gammavector.append(math.atan2(event.Thrust[1], event.Thrust[0]) * 180.0 / math.pi)
                    shortDateVector.append(datetime.datetime.strptime(event.GregorianDate,'%m/%d/%Y').date())

            if firstpass:
                DataAxesRight.plot(shortDateVector, gammavector, c='DarkGreen', lw=2, ls='--', label=r'$\gamma$ (degrees)')
            else:
                DataAxesRight.plot(shortDateVector, gammavector, c='DarkGreen', lw=2, ls='--')

        #plot delta
        if PlotOptions.PlotDelta:
            deltavector = []
            shortDateVector = []

            for event in self.missionevents:
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'PSFBthrust', 'LT_spiral']:
                    AppliedThrust = math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2)
                    deltavector.append(math.asin(event.Thrust[2] / AppliedThrust) * 180.0 / math.pi)
                    shortDateVector.append(datetime.datetime.strptime(event.GregorianDate,'%m/%d/%Y').date())

            if firstpass:
                DataAxesRight.plot(shortDateVector, deltavector, c='LightGreen', lw=2, ls='--', label=r'$\delta$ (degrees)')
            else:
                DataAxesRight.plot(shortDateVector, deltavector, c='LightGreen', lw=2, ls='--')


        #plot central body to thrust vector angle
        if PlotOptions.PlotArray_Thrust_Angle:
            CBthrustvector = []
            shortDateVector = []

            for event in self.missionevents:
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'PSFBthrust', 'LT_spiral']:
                    r = math.sqrt(event.SpacecraftState[0]**2 + event.SpacecraftState[1]**2 + event.SpacecraftState[2]**2)
                    AppliedThrust = math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2)
                    rdotT = event.SpacecraftState[0]*event.Thrust[0] + event.SpacecraftState[1]*event.Thrust[1] + event.SpacecraftState[2]*event.Thrust[2]
                    CBthrustvector.append( 90.0 - math.acos( rdotT / (r * AppliedThrust) ) * 180.0 / math.pi)
                    shortDateVector.append(datetime.datetime.strptime(event.GregorianDate,'%m/%d/%Y').date())

            if firstpass:
                DataAxesRight.plot(shortDateVector, CBthrustvector, c='Salmon', marker='o', ls='None', label='Array to thrust angle (degrees)')
            else:
                DataAxesRight.plot(shortDateVector, CBthrustvector, c='Salmon', marker='o', ls='None')

        #plot mass
        if PlotOptions.PlotMass:
            mass = []
            for event in self.missionevents:
                if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby' and event.EventType != 'LT_rndzvs':
                    mass.append(event.Mass * 1.0e-3)
            if firstpass:
                DataAxesLeft.plot(date_vector, mass, c='DarkGrey', lw=2, ls='-', label='Mass (1000 kg)')
            else:
                DataAxesLeft.plot(date_vector, mass, c='DarkGrey', lw=2, ls='-')

        #plot number of active thrusters
        if PlotOptions.PlotNumberOfEngines:
            numberofengines = []
            shortDateVector = []

            for event in self.missionevents:
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'PSFBthrust', 'LT_spiral']:
                    numberofengines.append(event.ActivePower)
                    shortDateVector.append(datetime.datetime.strptime(event.GregorianDate,'%m/%d/%Y').date())

            if firstpass:
                DataAxesLeft.plot(shortDateVector, numberofengines, 'k*', mew=2, markersize=7, label='Number of active thrusters')
            else:
                DataAxesLeft.plot(shortDateVector, numberofengines, 'k*', mew=2, markersize=7)

        #plot power actively used by the thrusters
        if PlotOptions.PlotActivePower:
            activepowervector = []
            shortDateVector = []

            for event in self.missionevents:
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'PSFBthrust', 'LT_spiral']:
                    activepowervector.append(event.ActivePower)
                    shortDateVector.append(datetime.datetime.strptime(event.GregorianDate,'%m/%d/%Y').date())

            if firstpass:
                DataAxesLeft.plot(shortDateVector, activepowervector, c='Navy', lw=2, ls='-', label='Power used by the propulsion system (kW)')
            else:
                DataAxesLeft.plot(shortDateVector, activepowervector, c='Navy', lw=2, ls='-')

        #plot waste heat from the propulsion system
        if PlotOptions.PlotWasteHeat:
            WasteHeatvector = []
            shortDateVector = []

            for event in self.missionevents:
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'PSFBthrust', 'LT_spiral']:
                    WasteHeatvector.append( (1 - event.AvailableThrust * event.Isp * 9.80665 / (2000 * event.ActivePower)) * event.ActivePower )
                    shortDateVector.append(datetime.datetime.strptime(event.GregorianDate,'%m/%d/%Y').date())

            if firstpass:
                DataAxesLeft.plot(shortDateVector, WasteHeatvector, c='Crimson', lw=2, ls='--', label='Waste heat from propulsion system (kW)')
            else:
                DataAxesLeft.plot(shortDateVector, WasteHeatvector, c='Crimson', lw=2, ls='--')

        #plot Earth distance in LU and SPE angle
        if PlotOptions.PlotEarthDistance or PlotOptions.PlotSunSpacecraftEarthAngle or PlotOptions.PlotSpacecraftViewingAngle:
            if not 'sun' in self.central_body.lower():
                print('getting Earth distance and Sun-Spacecraft-Earth angle is only supported if the central body is the sun')
            else:
                import de423
                from jplephem import Ephemeris
                eph = Ephemeris(de423)

                EarthDistanceVector = []
                SunSpacecraftEarthAngle = []
                SpacecraftViewingAngle = []
                for event in self.missionevents:
                    if event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby' and event.EventType != 'LT_rndzvs':
                    
                        #get the position of the central body relative to the solar system barycenter
                        cb_relative_to_solar_system_barycenter_in_ICRF = eph.position('sun', event.JulianDate)

                        #get the position of the Earth relative to the central body
                        embarycenter_in_ICRF = eph.position('earthmoon', event.JulianDate)
                        moonvector_in_ICRF = eph.position('moon', event.JulianDate)
                        Earth_relative_to_solar_system_barycenter_in_ICRF = embarycenter_in_ICRF - moonvector_in_ICRF * eph.earth_share
                        Earth_relative_to_cb_in_ICRF = Earth_relative_to_solar_system_barycenter_in_ICRF + cb_relative_to_solar_system_barycenter_in_ICRF

                        #rotate the spacecraft state relative to the central body into the ICRF frame
                        spacecraft_state_relative_to_central_body_in_ICRF = AstroFunctions.rotate_from_ecliptic_to_equatorial3(np.array(event.SpacecraftState[0:3]))

                        #get the vector from the Earth to the Spacecraft
                        Earth_Spacecraft_Vector = spacecraft_state_relative_to_central_body_in_ICRF - np.transpose(Earth_relative_to_cb_in_ICRF)[0]

                        #return the magnitude of the Earth-centered position vector
                        EarthDistanceVector.append(np.linalg.norm(Earth_Spacecraft_Vector) / self.LU)

                        #if applicable, compute sun-spacecraft-Earth angle
                        if PlotOptions.PlotSunSpacecraftEarthAngle:
                            cosAngle = np.dot(-Earth_Spacecraft_Vector, -spacecraft_state_relative_to_central_body_in_ICRF)/np.linalg.norm(spacecraft_state_relative_to_central_body_in_ICRF)/np.linalg.norm(Earth_Spacecraft_Vector)
                            SunSpacecraftEarthAngle.append(np.arccos(cosAngle) * 180.0/math.pi)
                            
                        if PlotOptions.PlotSpacecraftViewingAngle:
                            tanAngle = Earth_Spacecraft_Vector[2]/math.sqrt(Earth_Spacecraft_Vector[0]*Earth_Spacecraft_Vector[0]+Earth_Spacecraft_Vector[1]*Earth_Spacecraft_Vector[1])
                            SpacecraftViewingAngle.append(np.arctan(tanAngle) * 180.0/math.pi)  
                
                if PlotOptions.PlotEarthDistance:
                    if firstpass:
                        unitstring = 'LU'
                        if 'sun' in self.central_body.lower() and math.fabs(self.LU - 149597870.69100001) < 1.0:
                            unitstring = 'AU'
                        labelstring = 'Distance from Earth (' + unitstring + ')'
                        DataAxesLeft.plot(date_vector, EarthDistanceVector, c='g', lw=2, label=labelstring)
                    else:
                        DataAxesLeft.plot(date_vector, EarthDistanceVector, c='g', lw=2)
                if PlotOptions.PlotSunSpacecraftEarthAngle:
                    if firstpass:
                        DataAxesRight.plot(date_vector, SunSpacecraftEarthAngle, c='orangered', lw=2, ls = '-.', label='Sun-Spacecraft-Earth Angle')
                    else:
                        DataAxesRight.plot(date_vector, SunSpacecraftEarthAngle, c='orangered', lw=2, ls = '-.')
                if PlotOptions.PlotSpacecraftViewingAngle:
                    if firstpass:
                        DataAxesRight.plot(date_vector, SpacecraftViewingAngle, c='orangered', lw=2, ls = '-.', label='Latitude of Earth-Spacecraft Vector')
                    else:
                        DataAxesRight.plot(date_vector, SpacecraftViewingAngle, c='orangered', lw=2, ls = '-.')

        #plot throttle level
        if PlotOptions.PlotThrottleLevel:            
            ThrottleLevelVector = []
            shortDateVector = []
        
            for event in self.missionevents:
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'PSFBthrust', 'LT_spiral']:
                    ThrottleLevelVector.append(event.ThrottleLevel)
                    shortDateVector.append(datetime.datetime.strptime(event.GregorianDate,'%m/%d/%Y').date())
                    
            if firstpass:
                DataAxesLeft.scatter(shortDateVector, ThrottleLevelVector, c='midnightblue', marker='h' , label='Throttle level')
            else:                                                                                   
                DataAxesLeft.scatter(shortDateVector, ThrottleLevelVector, c='midnightblue', marker='h' )
                
        #plot sun-boresight angle
        if PlotOptions.PlotSunBoresightAngle:
            SunBoresightAngle = []
            for event in self.missionevents:
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSFBthrust', "PSBIthrust"]:
                    r = math.sqrt(event.SpacecraftState[0]**2 + event.SpacecraftState[1]**2 + event.SpacecraftState[2]**2)
                    AppliedThrust = math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2)
                    rdotT = -event.SpacecraftState[0]*event.Thrust[0] - event.SpacecraftState[1]*event.Thrust[1] - event.SpacecraftState[2]*event.Thrust[2]
                    SunBoresightAngle.append(math.acos( rdotT / (r * AppliedThrust) ) * 180.0 / math.pi)
                elif event.EventType != 'match_point' and event.EventType != 'upwr_flyby' and event.EventType != 'pwr_flyby' and event.EventType != 'LT_rndzvs':
                    SunBoresightAngle.append(0.0)
            if firstpass:
                DataAxesRight.plot(date_vector, SunBoresightAngle, c='purple', marker='o', ls='None', label='Sun-to-Boresight angle (degrees)')
            else:
                DataAxesRight.plot(date_vector, SunBoresightAngle, c='purple', marker='o', ls='None')
        
        #DataAxesLeft.set_ylim(bottom=-0.7)
        #DataAxesRight.set_ylim(bottom=-0.7)
        #DataAxesLeft.set_ylim([-0.2,2.5])
        #DataAxesLeft.set_ylim([-0.2,2.5])

    def WriteDateReport(self, PlotOptions, reportfile):
        #write header
        reportfile.write('JD, MM-DD-YYYY, Event duration (days), Event type, Event location')
        if PlotOptions.IncludeStateVectorInReport:
            reportfile.write(',x (km), y (km), z (km), xdot (km/s), ydot (km/s), zdot (km/s)')
        if PlotOptions.PlotR:
            unitstring = 'LU'
            if 'sun' in self.central_body.lower() and math.fabs(self.LU - 149597870.69100001) < 1.0:
                unitstring = 'AU'
            reportfile.write(', Distance from ' + self.central_body + ' (' + unitstring + ')')
        if PlotOptions.PlotV:
            reportfile.write(', Velocity magnitude (LU/TU)')
        if PlotOptions.PlotThrust:
            reportfile.write(', Applied thrust (N)')
        if PlotOptions.PlotIsp:
            reportfile.write(', Isp (s)')
        if PlotOptions.PlotMdot:
            reportfile.write(', Mass flow rate (kg/s)')
        if PlotOptions.PlotEfficiency:
            reportfile.write(', Propulsion system efficiency')
        if PlotOptions.PlotThrottle:
            reportfile.write(', Control magnitude')
        if PlotOptions.PlotPower:
            reportfile.write(', Power produced by spacecraft (kW)')
        if PlotOptions.PlotGamma:
            reportfile.write(', in-plane control angle (degrees)')
        if PlotOptions.PlotDelta:
            reportfile.write(', out-of-plane control angle (degrees)')
        if PlotOptions.PlotArray_Thrust_Angle:
            reportfile.write(', Array to thrust angle (degrees)')
        if PlotOptions.PlotMass:
            reportfile.write(', Mass (kg)')
        if PlotOptions.PlotNumberOfEngines:
            reportfile.write(', Number of active thrusters')
        if PlotOptions.PlotActivePower:
            reportfile.write(', Power used by the propulsion system (kW)')
        if PlotOptions.PlotWasteHeat:
            reportfile.write(', Waste heat from propulsion system (kW)')  
        if PlotOptions.PlotEarthDistance:
            unitstring = 'LU'
            if 'sun' in self.central_body.lower() and math.fabs(self.LU - 149597870.69100001) < 1.0:
                unitstring = 'AU'
            reportfile.write(', Distance from Earth (' + unitstring + ')')
        if PlotOptions.PlotSunSpacecraftEarthAngle:
            reportfile.write(', Sun-Spacecraft-Earth Angle')
        if PlotOptions.PlotSpacecraftViewingAngle:
            reportfile.write(', Latitude of Earth-Spacecraft Vector')
        if PlotOptions.PlotThrottleLevel:
            reportfile.write(', Throttle level')
        if PlotOptions.PlotSunBoresightAngle:
            reportfile.write(', Sun-Boresight Angle')
        reportfile.write('\n')

        #load stuff if necessary
        if PlotOptions.PlotEarthDistance or PlotOptions.PlotSunSpacecraftEarthAngle or PlotOptions.PlotSpacecraftViewingAngle:
            if not 'sun' in self.central_body.lower():
                print('getting Earth distance and Sun-Spacecraft-Earth angle is only supported if the central body is the sun')
            else:
                import de423
                from jplephem import Ephemeris
                self.eph = Ephemeris(de423)
        if PlotOptions.PlotThrottleLevel:
            #instantiate a throttle table object
            thruster_throttle_table = ThrottleTable.ThrottleTable(PlotOptions.throttletablefile)
            self.high_Mdot_set, self.low_Mdot_set = thruster_throttle_table.create_non_dominated_sets()

        #now print the report
        for event in self.missionevents:
            reportfile.write(str(event.JulianDate))
            reportfile.write(', ' + event.GregorianDate)
            reportfile.write(', ' + str(event.TimestepLength))
            reportfile.write(', ' + event.EventType)
            reportfile.write(', ' + event.Location)
            if PlotOptions.IncludeStateVectorInReport:
                for state in event.SpacecraftState:
                    reportfile.write(', ' + str(state))
            if PlotOptions.PlotR:
                reportfile.write(', ' + str(math.sqrt(event.SpacecraftState[0]**2 + event.SpacecraftState[1]**2 + event.SpacecraftState[2]**2) / self.LU))
            if PlotOptions.PlotV:
                reportfile.write(', ' + str(math.sqrt(event.SpacecraftState[3]**2 + event.SpacecraftState[4]**2 + event.SpacecraftState[5]**2) / self.LU * self.TU))
            if PlotOptions.PlotThrust:
                if event.EventType not in ['match_point','upwr_flyby','pwr_flyby','LT_rndzvs','LT_spiral']:
                    reportfile.write(', ' + str(math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2)))
                elif event.EventType == 'LT_spiral':
                    reportfile.write(', ' + str(event.AvailableThrust * self.thruster_duty_cycle))
                else:
                    reportfile.write(', 0.0')
            if PlotOptions.PlotIsp:
                reportfile.write(', ' + str(event.Isp))
            if PlotOptions.PlotMdot:
                reportfile.write(', ' + str(event.MassFlowRate))
            if PlotOptions.PlotEfficiency:
                if event.EventType in ['SFthrust', 'SSFthrust','FBLTSthrust','FBLTthrust','PSBIthrust','PSFBthrust','LT_spiral']:
                    reportfile.write(', ' + str(event.AvailableThrust * event.Isp * 9.80665 / (2000 * event.ActivePower)))
                else:
                    reportfile.write(', 0.0')
            if PlotOptions.PlotThrottle:
                if event.EventType in ['SFthrust', 'SSFthrust','FBLTSthrust','FBLTthrust','PSBIthrust','PSFBthrust']:
                    reportfile.write(', ' + str(math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2) / (event.AvailableThrust * self.thruster_duty_cycle)))
                elif event.EventType == 'LT_spiral':
                    reportfile.write(', ' + str(1.0))
                elif event.EventType not in ['match_point','upwr_flyby','pwr_flyby','LT_rndzvs']:
                    reportfile.write(', ' + str(0.0))
            if PlotOptions.PlotPower:
                if event.EventType not in ['match_point','upwr_flyby','pwr_flyby','LT_rndzvs']:
                    if (event.AvailablePower > 0.0):
                        reportfile.write(', ' + str(event.AvailablePower))
                    else:
                        reportfile.write(', ' + str(0.0))
                else:
                    reportfile.write(', ' + str(0.0))
            if PlotOptions.PlotGamma:
                if event.EventType in ['SFthrust', 'SSFthrust','FBLTSthrust','FBLTthrust','PSBIthrust','PSFBthrust','LT_spiral']:
                    reportfile.write(', ' + str(math.atan2(event.Thrust[1], event.Thrust[0]) * 180.0 / math.pi))
                else:
                    reportfile.write(', ' + str(0.0))
            if PlotOptions.PlotDelta:
                if event.EventType in ['SFthrust', 'SSFthrust','FBLTSthrust','FBLTthrust','PSBIthrust','PSFBthrust','LT_spiral']:
                    AppliedThrust = math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2)
                    reportfile.write(', ' + str(math.asin(event.Thrust[2] / AppliedThrust) * 180.0 / math.pi))
                else:
                    reportfile.write(', ' + str(0.0))
            if PlotOptions.PlotArray_Thrust_Angle:
                if event.EventType in ['SFthrust', 'SSFthrust','FBLTSthrust','FBLTthrust','PSBIthrust','PSFBthrust','LT_spiral']:
                    r = math.sqrt(event.SpacecraftState[0]**2 + event.SpacecraftState[1]**2 + event.SpacecraftState[2]**2)
                    AppliedThrust = math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2)
                    rdotT = event.SpacecraftState[0]*event.Thrust[0] + event.SpacecraftState[1]*event.Thrust[1] + event.SpacecraftState[2]*event.Thrust[2]
                    reportfile.write(', ' + str(math.acos( rdotT / (r * AppliedThrust) ) * 180.0 / math.pi))
                else:
                    reportfile.write(', ' + str(0.0))
            if PlotOptions.PlotMass:
                reportfile.write(', ' + str(event.Mass))
            if PlotOptions.PlotNumberOfEngines:
                if event.EventType in ['SFthrust', 'SSFthrust','FBLTSthrust','FBLTthrust','PSBIthrust','PSFBthrust','LT_spiral']:
                    reportfile.write(', ' + str(event.Number_of_Active_Engines))
                else:
                    reportfile.write(', ' + str(0))
            if PlotOptions.PlotActivePower:
                if event.EventType in ['SFthrust', 'SSFthrust','FBLTSthrust','FBLTthrust','PSBIthrust','PSFBthrust','LT_spiral']:
                    reportfile.write(', ' + str(event.ActivePower))
                else:
                    reportfile.write(', ' + str(0.0))
            if PlotOptions.PlotWasteHeat:
                if event.EventType in ['SFthrust', 'SSFthrust','FBLTSthrust','FBLTthrust','PSBIthrust','PSFBthrust','LT_spiral']:
                    reportfile.write(', ' + str((1.0 - event.AvailableThrust * event.Isp * 9.80665 / (2000.0 * event.ActivePower)) * event.ActivePower))
                else:
                    reportfile.write(', ' + str(0.0))
            if PlotOptions.PlotEarthDistance or PlotOptions.PlotSunSpacecraftEarthAngle or PlotOptions.PlotSpacecraftViewingAngle:
                if not 'sun' in self.central_body.lower():
                    print('getting Earth distance and Sun-Spacecraft-Earth angle is only supported if the central body is the sun')
                else:
                    #get the position of the central body relative to the solar system barycenter
                    cb_relative_to_solar_system_barycenter_in_ICRF = self.eph.position('sun', event.JulianDate)

                    #get the position of the Earth relative to the central body
                    embarycenter_in_ICRF = self.eph.position('earthmoon', event.JulianDate)
                    moonvector_in_ICRF = self.eph.position('moon', event.JulianDate)
                    Earth_relative_to_solar_system_barycenter_in_ICRF = embarycenter_in_ICRF - moonvector_in_ICRF * self.eph.earth_share
                    Earth_relative_to_cb_in_ICRF = Earth_relative_to_solar_system_barycenter_in_ICRF + cb_relative_to_solar_system_barycenter_in_ICRF

                    #rotate the spacecraft state relative to the central body into the ICRF frame
                    spacecraft_state_relative_to_central_body_in_ICRF = AstroFunctions.rotate_from_ecliptic_to_equatorial3(np.array(event.SpacecraftState[0:3]))

                    #get the vector from the Earth to the Spacecraft
                    Earth_Spacecraft_Vector = spacecraft_state_relative_to_central_body_in_ICRF - np.transpose(Earth_relative_to_cb_in_ICRF)[0]

                    #return the magnitude of the Earth-centered position vector
                    EarthDistance = np.linalg.norm(Earth_Spacecraft_Vector) / self.LU

                    if PlotOptions.PlotEarthDistance:
                        reportfile.write(', ' + str(EarthDistance))
                    if PlotOptions.PlotSunSpacecraftEarthAngle:
                        cosAngle = np.dot(-Earth_Spacecraft_Vector, -spacecraft_state_relative_to_central_body_in_ICRF)/np.linalg.norm(spacecraft_state_relative_to_central_body_in_ICRF)/np.linalg.norm(Earth_Spacecraft_Vector)
                        reportfile.write(', ' + str(np.arccos(cosAngle) * 180.0/math.pi))
                    if PlotOptions.PlotSpacecraftViewingAngle:
                        tanAngle = Earth_Spacecraft_Vector[2]/math.sqrt(Earth_Spacecraft_Vector[0]*Earth_Spacecraft_Vector[0]+Earth_Spacecraft_Vector[1]*Earth_Spacecraft_Vector[1])
                        reportfile.write(', ' + str(np.arctan(tanAngle) * 180.0/math.pi))   
                        
            if PlotOptions.PlotThrottleLevel:
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust','PSFBthrust', 'end_spiral']:
                    if PlotOptions.throttlesetmode == 0:
                        reportfile.write(', ' + str(int(self.high_Mdot_set.find_nearest_throttle_setting_1D(event.ActivePower / event.Number_of_Active_Engines).TL.strip('TL'))))
                    elif PlotOptions.throttlesetmode == 1:
                        reportfile.write(', ' + str(int(self.low_Mdot_set.find_nearest_throttle_setting_1D(event.ActivePower / event.Number_of_Active_Engines).TL.strip('TL'))))
                    elif PlotOptions.throttlesetmode == 2:
                        print('2D throttle matching not currently implemented')
                else:
                    reportfile.write(', ' + str(0.0))
            #plot sun-boresight angle
            if PlotOptions.PlotSunBoresightAngle:
                if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', "PSBIthrust",'PSFBthrust']:
                    r = math.sqrt(event.SpacecraftState[0]**2 + event.SpacecraftState[1]**2 + event.SpacecraftState[2]**2)
                    AppliedThrust = math.sqrt(event.Thrust[0]**2 + event.Thrust[1]**2 + event.Thrust[2]**2)
                    rdotT = event.SpacecraftState[0]*event.Thrust[0] + event.SpacecraftState[1]*event.Thrust[1] + event.SpacecraftState[2]*event.Thrust[2]
                    reportfile.write(', ' + str(math.acos( rdotT / (r * AppliedThrust) ) * 180.0 / math.pi))
                else:
                    reportfile.write(', ' + str(0.0))
            reportfile.write('\n')

    def AutoTableJourney(self, TableFile, PlotOptions, skipNext):
        for event in self.missionevents:
            skipNext = event.AutoTableLine(TableFile, PlotOptions, skipNext)
        return skipNext

    def getManeuver(self, phaseIndex = 0, maneuverIndex = 0):
        #function to get the nth DSM/TCM in a phase
        #if that maneuver does not exist, return None
        
        #Step 1: search through the journey to find the right phase. We do this by counting the number of boundary events until we get to the end of the journey or the phase that we want
        #since phase boundaries are always upwr_flyby or pwr_flyby, this isn't too hard!
        phaseCount = 0
        maneuverCount = 0
        for event in self.missionevents:
            if event.EventType in ['pwr_flyby','upwr_flyby']:
                phaseCount += 1
                maneuverCount = 0

            if event.EventType == 'chem_burn':
                if maneuverCount == maneuverIndex:
                    #this is the maneuver that we want, so return it
                    return event
                
                #if we didn't just return, increment the maneuver count
                maneuverCount += 1

        #if we got this far then the event that the user asked for does not exist, so return "none"
        return None

    def getPeriapseDistance(self, scale = None):
        r_min = 1.0e+10

        if scale == None:
            scale = self.LU

        for event in self.missionevents:
            r = (event.SpacecraftState[0]**2 + event.SpacecraftState[1]**2 + event.SpacecraftState[2]**2 ) ** 0.5

            if r < r_min:
                r_min = r

        return r_min / scale