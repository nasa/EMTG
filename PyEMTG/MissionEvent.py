import EOM

import math
import numpy as np
import copy


class MissionEvent(object):
    lock = None

    def __init__(self, parent, inputcell):
        self.parent = parent

        self.EventNumber = -1 #index of event
        self.JulianDate = -1 #Julian Ephemeris Date of event
        self.GregorianDate = 'banana' #gregorian date (MM/DD/YYYY) of event
        self.EventType = 'shampoo' #what class of event? Launch? Thrust? Coast? Flyby? etc
        self.Location = 'bottle' #where did the event take place
        self.TimestepLength = 1.0e-4 #how long is the time step in days?
        self.Altitude = -1 #altitude for flybys
        self.BdotR = -1 #for flybys
        self.BdotT = -1 #for flybys
        self.RightAscension = 0.0 #RA for maneuvers
        self.Declination = 0.0 #DEC for maneuvers
        self.C3 = 0.0 #C3 for arrivals and departures
        self.SpacecraftState = [0.0]*6 #6-vector
        self.DeltaVorThrustVectorControl = [0.0]*3 #3-vector
        self.Thrust = [0.0]*3 #3-vector
        self.DVmagorThrottle = 0.0 #for impulses or thrust arcs
        self.AvailableThrust = 0.0 #for thrust arcs
        self.Isp = 1.0 #for all propulsive maneuvers
        self.AvailablePower = 0.0 #for thrust arcs
        self.MassFlowRate = 0.0 #kg/s
        self.Mass = 1.0 #Mass after event occurs
        self.Number_of_Active_Engines = 0 #number of thrusters firing at the center of this arc
        self.ActivePower = 0.0 #how much power is currently being used by the thrust system
        self.ThrottleLevel = "none"
        self.eventlabel = None

        self.parse_input_line(inputcell)

    def parse_input_line(self, inputcell):
        for i in range(0, len(inputcell)):
            inputcell[i] = inputcell[i].strip(' ')
            if inputcell[i] == '-' or inputcell[i] == '-\n':
                inputcell[i] = 0.0

        

        self.EventNumber = int(inputcell[0])
        self.JulianDate = float(inputcell[1])
        self.GregorianDate = inputcell[2]
        self.EventType = inputcell[3]
        self.Location = inputcell[4]
        self.TimestepLength = float(inputcell[5])
        self.Altitude = float(inputcell[6])
        self.BdotR = float(inputcell[7])
        self.BdotT = float(inputcell[8])
        self.RightAscension = float(inputcell[9])
        self.Declination = float(inputcell[10])
        self.C3 = float(inputcell[11])
        self.SpacecraftState = [float(inputcell[12]), float(inputcell[13]), float(inputcell[14]), float(inputcell[15]), float(inputcell[16]), float(inputcell[17])]
        self.DeltaVorThrustVectorControl = [float(inputcell[18]), float(inputcell[19]), float(inputcell[20])]
        self.Thrust = [float(inputcell[21]), float(inputcell[22]), float(inputcell[23])]
        self.DVmagorThrottle = float(inputcell[24])
        if inputcell[25] == "impulse":
            self.AvailableThrust = 0.0
        else:
            self.AvailableThrust = float(inputcell[25])
        if inputcell[26] == "LV-supplied":
            self.Isp = 0.0
        elif inputcell[26] == "UNHANDLED EVENT TYPE":
            self.Isp = 0.0
        else:
            self.Isp = float(inputcell[26])
            
        self.AvailablePower = float(inputcell[27])
        self.MassFlowRate = float(inputcell[28])
        self.Mass = float(inputcell[29])

        if len(inputcell) >= 31:
            self.Number_of_Active_Engines = int(inputcell[30])
        
        if len(inputcell) >= 32:
            self.ActivePower = float(inputcell[31])

        if len(inputcell) > 33:
            self.ThrottleLevel = inputcell[32]


    def PlotEvent(self, GraphicsObject, LU, TU, mu, PlotOptions, BeforeMatchPoint, CustomLabel=None):
        from scipy.integrate import ode
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d import proj3d
        from matplotlib.patches import FancyArrowPatch
        #switch between different event types
        pX = self.SpacecraftState[0]
        pY = self.SpacecraftState[1]
        pZ = self.SpacecraftState[2]

        if PlotOptions.PlotCentralBody != 'Journey central body':
            import spiceypy
            journey_central_body = self.parent.central_body
            if journey_central_body in ['Mars','Jupiter','Uranus','Neptune','Pluto']:
                journey_central_body = journey_central_body + ' Barycenter'

            if journey_central_body != PlotOptions.PlotCentralBody:
                body_state, LT = spiceypy.spkezr(journey_central_body, (self.JulianDate - 2451545.0) * 86400.0, 'J2000', "None", PlotOptions.PlotCentralBody)

                pX += body_state[0]
                pY += body_state[1]
                pZ += body_state[2]

        if self.EventType == "launch" or self.EventType == "departure":
            GraphicsObject.scatter(pX, pY, pZ, s=4, c='g', marker='^')
            BeforeMatchPoint = True
        elif self.EventType == "begin_spiral" or self.EventType == "end_spiral":
            GraphicsObject.scatter(pX, pY, pZ, s=4, c='orange', marker='^')
            BeforeMatchPoint = True
        elif self.EventType == "upwr_flyby":
            GraphicsObject.scatter(pX, pY, pZ, s=4, c='b', marker=r'$\circlearrowleft$')
            BeforeMatchPoint = True
        elif self.EventType == "pwr_flyby":
            GraphicsObject.scatter(pX, pY, pZ, s=4, c='r', marker=r'$\circlearrowleft$')
            BeforeMatchPoint = True
        elif self.EventType == "chem_burn" and self.DVmagorThrottle > 0.001:
            GraphicsObject.scatter(pX, pY, pZ, s=4, c='r', marker='o')
        
        if self.EventType in ["SFthrust", 'SSFthrust', "PSBIthrust", "coast", 'nav-coast', 'Scoast', "force-coast", "LT_spiral", 'end_spiral', 'begin_spiral']:
            color = ''
            if self.EventType in ["coast", 'Scoast', 'nav-coast']:
                linestyle = '-'
                color = 'k'
                weight = 0.5
            elif self.EventType in ['LT_spiral', 'end_spiral', 'begin_spiral']:
                linestyle = '-'
                color = 'g'
                weight = 2
            elif self.EventType in ["SFthrust", 'SSFthrust', "PSBIthrust"]:
                linestyle = '-'
                color = 'k'
                weight = 1
                if PlotOptions.ShowThrustVectors:
                    ControlVector = np.array(self.Thrust) / self.AvailableThrust * LU * 0.1
                    T_length = np.linalg.norm(self.Thrust) / self.AvailableThrust * LU * 0.1
                    X = self.SpacecraftState[0] + np.array([0.0, ControlVector[0]])
                    Y = self.SpacecraftState[1] + np.array([0.0, ControlVector[1]])
                    Z = self.SpacecraftState[2] + np.array([0.0, ControlVector[2]])

                    if PlotOptions.PlotCentralBody != 'Journey central body':
                        import spiceypy
                        journey_central_body = self.parent.central_body
                        if journey_central_body in ['Mars','Jupiter','Uranus','Neptune','Pluto']:
                            journey_central_body = journey_central_body + ' Barycenter'

                        if journey_central_body != PlotOptions.PlotCentralBody:
                            body_state, LT = spiceypy.spkezr(journey_central_body, (self.JulianDate - 2451545.0) * 86400.0, 'J2000', "None", PlotOptions.PlotCentralBody)

                            X += body_state[0]
                            Y += body_state[1]
                            Z += body_state[2]

                    GraphicsObject.plot(X, Y, Z, c='violet', lw=1)
            else:
                linestyle = '--'
                color = 'b'
                weight = 1

            # comment out to turn off event dots
            #GraphicsObject.scatter(self.SpacecraftState[0], self.SpacecraftState[1], self.SpacecraftState[2], s=1, c=color, marker='o')

            if PlotOptions.ShowPropagatedTrajectory:
                if self.EventType in ['begin_spiral', 'LT_spiral']:
                    #plot ONLY a forward propagation 
                    ManeuverPointState = np.zeros(6)
                    ManeuverPointState[0:6] = np.array(self.SpacecraftState) / LU
                    ManeuverPointState[3:6] *= TU
                    ForwardIntegrateObject = ode(EOM.EOM_inertial_2body).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                    ForwardIntegrateObject.set_initial_value(ManeuverPointState).set_f_params(1.0).set_jac_params(1.0)

                    dt = self.TimestepLength * 86400.0 / TU / 100
                    StateHistoryForward = []
                    while ForwardIntegrateObject.successful() and ForwardIntegrateObject.t < self.TimestepLength * 86400.0 / TU:
                        ForwardIntegrateObject.integrate(ForwardIntegrateObject.t + dt)
                        StateHistoryForward.append(ForwardIntegrateObject.y * LU)
                    StateHistoryForward.reverse()

                    X = []
                    Y = []
                    Z = []
                    for StateLine in StateHistoryForward:
                        X.append(StateLine[0])
                        Y.append(StateLine[1])
                        Z.append(StateLine[2])
                else: #thrusting!
                    CenterPointState = np.zeros(6)
                    CenterPointState[0:6] = np.array(self.SpacecraftState) / LU
                    CenterPointState[3:6] *= TU
                    CenterPointStateAfterManeuver = copy.deepcopy(CenterPointState)
                    CenterPointStateAfterManeuver[3:6] += np.array(self.DeltaVorThrustVectorControl) * TU /LU

                    ForwardIntegrateObject = ode(EOM.EOM_inertial_2body).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                    ForwardIntegrateObject.set_initial_value(CenterPointStateAfterManeuver).set_f_params(1.0).set_jac_params(1.0)

                    dt = self.TimestepLength * 86400.0 / TU / 100
                    StateHistoryForward = []
                    epoch = (self.JulianDate - 2451545.0) * 86400.0
                    while ForwardIntegrateObject.successful() and ForwardIntegrateObject.t < self.TimestepLength * 86400.0 / TU / 2.0:
                        epoch += dt * TU
                        stateArray = ForwardIntegrateObject.integrate(ForwardIntegrateObject.t + dt)
                        state = [stateArray[0] * LU, stateArray[1] * LU, stateArray[2] * LU, epoch]
                        StateHistoryForward.append(state)

                    BackwardIntegrateObject = ode(EOM.EOM_inertial_2body).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                    BackwardIntegrateObject.set_initial_value(CenterPointState).set_f_params(1.0).set_jac_params(1.0)

                    dt = self.TimestepLength * 86400.0 / TU / 100
                    StateHistoryBackward = []
                    epoch = (self.JulianDate - 2451545.0) * 86400.0
                    while BackwardIntegrateObject.successful() and BackwardIntegrateObject.t > -self.TimestepLength * 86400.0 / TU / 2.0:   
                        epoch -= dt * TU
                        stateArray = BackwardIntegrateObject.integrate(BackwardIntegrateObject.t - dt)
                        state = [stateArray[0] * LU, stateArray[1] * LU, stateArray[2] * LU, epoch]
                        StateHistoryBackward.append(state)

                    StateHistoryBackward.reverse()

                    X = []
                    Y = []
                    Z = []
                    for StateLine in StateHistoryBackward:
                        if PlotOptions.PlotCentralBody != 'Journey central body':
                            import spiceypy
                            journey_central_body = self.parent.central_body
                            
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
                    for StateLine in StateHistoryForward:
                        if PlotOptions.PlotCentralBody != 'Journey central body':
                            import spiceypy
                            journey_central_body = self.parent.central_body
                            if journey_central_body in ['Mars','Jupiter','Uranus','Neptune','Pluto']:
                                journey_central_body = journey_central_body + ' Barycenter'
                                
                            if journey_central_body != PlotOptions.PlotCentralBody:
                                body_state, LT = spiceypy.spkezr(journey_central_body, StateLine[-1], 'J2000', "NONE", PlotOptions.PlotCentralBody)
                            
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

                GraphicsObject.plot(X, Y, Z, lw=weight, c=color, ls=linestyle)

        elif self.EventType in ["FBLTthrust", "FBLTSthrust", 'PSFBthrust']:
            GraphicsObject.scatter(pX, pY, pZ, s=2, c='k', marker='o')
            
            if PlotOptions.ShowThrustVectors:
                ControlVector = np.array(self.Thrust) / self.AvailableThrust * LU * 0.1
                X = self.SpacecraftState[0] + np.array([0.0, ControlVector[0]])
                Y = self.SpacecraftState[1] + np.array([0.0, ControlVector[1]])
                Z = self.SpacecraftState[2] + np.array([0.0, ControlVector[2]])

                if PlotOptions.PlotCentralBody != 'Journey central body':
                    import spiceypy
                    journey_central_body = self.parent.central_body
                    if journey_central_body in ['Mars','Jupiter','Uranus','Neptune','Pluto']:
                        journey_central_body = journey_central_body + ' Barycenter'

                    body_state, LT = spiceypy.spkezr(journey_central_body, (self.JulianDate - 2451545.0) * 86400.0, 'J2000', "None", PlotOptions.PlotCentralBody)

                    X += body_state[0]
                    Y += body_state[1]
                    Z += body_state[2]
                GraphicsObject.plot(X, Y, Z, c='violet', lw=2)


            if PlotOptions.ShowPropagatedTrajectory:
                CenterPointState = np.zeros(7)
                CenterPointState[0:6] = np.array(self.SpacecraftState) / LU
                CenterPointState[3:6] *= TU
                CenterPointState[6] = 1.0
                ScaledThrust = np.array(self.Thrust) * TU * TU / 1000.0 / self.Mass / LU
                ScaledMdot = self.MassFlowRate / self.Mass * TU

                ForwardIntegrateObject = ode(EOM.EOM_inertial_2bodyconstant_thrust).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                ForwardIntegrateObject.set_initial_value(CenterPointState).set_f_params(ScaledThrust, ScaledMdot, 1.0).set_jac_params(ScaledThrust, ScaledMdot, 1.0)

                dt = self.TimestepLength * 86400.0 / TU / 10
                StateHistoryForward = []
                epoch = (self.JulianDate - 2451545.0) * 86400.0
                while ForwardIntegrateObject.successful() and ForwardIntegrateObject.t < self.TimestepLength * 86400.0 / TU / 2.0:
                    epoch += dt * TU
                    stateArray = ForwardIntegrateObject.integrate(ForwardIntegrateObject.t + dt)
                    state = [stateArray[0] * LU, stateArray[1] * LU, stateArray[2] * LU, epoch]
                    StateHistoryForward.append(state)

                BackwardIntegrateObject = ode(EOM.EOM_inertial_2bodyconstant_thrust).set_integrator('dop853', atol=1.0e-8, rtol=1.0e-8)
                BackwardIntegrateObject.set_initial_value(CenterPointState).set_f_params(ScaledThrust, ScaledMdot, 1.0).set_jac_params(ScaledThrust, ScaledMdot, 1.0)

                dt = self.TimestepLength * 86400.0 / TU / 10
                StateHistoryBackward = []
                epoch = (self.JulianDate - 2451545.0) * 86400.0
                while BackwardIntegrateObject.successful() and BackwardIntegrateObject.t > -self.TimestepLength * 86400.0 / TU / 2.0:
                    epoch -= dt * TU
                    stateArray = BackwardIntegrateObject.integrate(BackwardIntegrateObject.t - dt)
                    state = [stateArray[0] * LU, stateArray[1] * LU, stateArray[2] * LU, epoch]
                    StateHistoryBackward.append(state)

                StateHistoryBackward.reverse()

                X = []
                Y = []
                Z = []
                for StateLine in StateHistoryBackward:
                    if PlotOptions.PlotCentralBody != 'Journey central body':
                        import spiceypy
                        journey_central_body = self.parent.central_body
                        if journey_central_body in ['Mars','Jupiter','Uranus','Neptune','Pluto']:
                            journey_central_body = journey_central_body + ' Barycenter'
                            
                        if journey_central_body != PlotOptions.PlotCentralBody:
                            body_state, LT = spiceypy.spkezr(journey_central_body, (self.JulianDate - 2451545.0) * 86400.0 + StateLine[-1], 'J2000', "None", PlotOptions.PlotCentralBody)
                            
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
                for StateLine in StateHistoryForward:
                    if PlotOptions.PlotCentralBody != 'Journey central body':
                        import spiceypy
                        journey_central_body = self.parent.central_body
                        if journey_central_body in ['Mars','Jupiter','Uranus','Neptune','Pluto']:
                            journey_central_body = journey_central_body + ' Barycenter'

                        if journey_central_body != PlotOptions.PlotCentralBody:
                            body_state, LT = spiceypy.spkezr(journey_central_body, (self.JulianDate - 2451545.0) * 86400.0 + StateLine[-1], 'J2000', "None", PlotOptions.PlotCentralBody)
                            
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

                GraphicsObject.plot(X, Y, Z, lw=2, c='k')

        elif self.EventType == "insertion":
            GraphicsObject.scatter(pX, pY, pZ, s=4, c='r', marker='v')
        elif self.EventType == "LT_rndzvs":
            GraphicsObject.scatter(pX, pY, pZ, s=4, c='m', marker='o')
        elif self.EventType == "intercept" or self.EventType == 'momtransfer':
            GraphicsObject.scatter(pX, pY, pZ, s=4, c='m', marker='d')
        elif self.EventType == "rendezvous":
            GraphicsObject.scatter(pX, pY, pZ, s=4, c='r', marker='d')
        elif self.EventType == "match-vinf":
            GraphicsObject.scatter(pX, pY, pZ, s=4, c='c', marker='s')
        elif self.EventType == "match_point":
            GraphicsObject.scatter(pX, pY, pZ, s=0, c='k', marker='s')
            BeforeMatchPoint = False


        nonPlottableEvents = ['zeroflyby',
                              'coast', 
                              'nav-coast', 
                              'Scoast', 
                              'force-coast',
                              'SFthrust',
                              'SSFthrust',
                              'FBLTthrust', 
                              'PSBIthrust',
                              'PSFBthrust',
                              'match_point',
                              'waiting', 
                              'LT_spiral',
                              'TCM']

        if (PlotOptions.LabelEvents and (not self.EventType in nonPlottableEvents)) or CustomLabel != None:
            print(self.EventNumber, CustomLabel)
            self.LabelEvent(GraphicsObject, PlotOptions, CustomLabel)
        
        return BeforeMatchPoint

    def LabelEvent(self, GraphicsObject, PlotOptions, labelString = None):        
        import matplotlib
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d import proj3d

        description = ''
        if labelString == None:
            if PlotOptions.ShowTextDescriptions:
                EventTypeFormatted = self.EventType
                if EventTypeFormatted == 'upwr_flyby':
                    EventTypeFormatted = 'unpowered flyby'
                elif EventTypeFormatted == 'pwr_flyby':
                    EventTypeFormatted = 'powered flyby'
                elif EventTypeFormatted == 'chem_burn':
                    EventTypeFormatted = 'chemical burn'
                elif EventTypeFormatted == 'LT_rndzvs':
                    EventTypeFormatted = 'rendezvous'
                elif EventTypeFormatted == 'begin_spiral':
                    EventTypeFormatted = 'begin spiral'
                elif EventTypeFormatted == 'end_spiral':
                    EventTypeFormatted = 'end spiral'
                elif EventTypeFormatted == 'mission_end':
                    EventTypeFormatted = 'end of mission'
                elif EventTypeFormatted == 'momtransfer':
                    EventTypeFormatted = 'momentum transfer'

                if PlotOptions.NumberEventLabels:
                    EventNumber = 'Event # ' + str(PlotOptions.EventCounter) + ':\n'
                else:
                    EventNumber = ''

                description = EventNumber + EventTypeFormatted.capitalize() + '\n' + self.Location
        
                if PlotOptions.DisplayEventDates:
                    description += '\n' + self.GregorianDate
        
                if PlotOptions.DisplayEventSpecs:
                    #for launches a C3 and DLA are needed
                    if self.EventType == 'launch':
                        description += '\n$C_3$ = ' + "{0:.3f}".format(self.C3) + ' $km^2/s^2$'
                        #add the LV to the description?
                        description += '\nDLA = ' + "{0:.1f}".format(self.Declination) + '$^{\circ}$'

                    #for non-launch departures only the C3 is needed
                    if self.EventType == 'departure' and self.C3 > 0.0:
                        description += '\n$C_3$ = ' + "{0:.3f}".format(self.C3) + ' $km^2/s^2$'

                    #for other events, output v-infinity and DLA
                    if self.EventType in ['upwr_flyby', 'pwr_flyby', 'intercept', 'interface', 'insertion', 'momtransfer']:
                        description += '\n$v_\infty$ = ' + "{0:.3f}".format(math.sqrt(self.C3)) + ' $km/s$'
                        description += '\nDEC = ' + "{0:.1f}".format(self.Declination) + '$^{\circ}$'

                    #for flybys, altitude should be outputed
                    if self.EventType in ['upwr_flyby', 'pwr_flyby', 'periapse']:
                        description += '\naltitude = ' + "{0:.0f}".format(self.Altitude) + ' $km$'



                    #for propulsive events, a deltaV is needed
                    if (self.EventType == 'departure' or self.EventType == 'pwr_flyby' or self.EventType == 'insertion' or self.EventType == 'chem_burn' or self.EventType == 'rendezvous') and self.DVmagorThrottle > 0.0:
                            description += '\n$\Delta v$ = ' + "{0:.3f}".format(self.DVmagorThrottle) + ' $km/s$'

                    #always append the spacecraft Mass
                    if PlotOptions.DisplayEventMass:
                        if self.Mass < 1.0e+5:
                            description += '\nm = ' + "{0:.0f}".format(self.Mass) + ' $kg$'
                        else:
                            description += '\nm = ' + "{0:.4e}".format(self.Mass) + ' $kg$'

                    if PlotOptions.DisplayArrivalPhaseAngle and self.EventType == 'intercept':
                        R = self.SpacecraftState[0:3]
                        V = self.DeltaVorThrustVectorControl
                        r = math.sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2])
                        v = math.sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2])
                        try:
                            PhaseAngle = math.acos(-(R[0]*V[0] + R[1]*V[1] + R[2]*V[2]) / (r * v)) * 180.0 / math.pi
                            description += '\n' + r'$\beta $' + ' = ' + "{0:.1f}".format(PhaseAngle) + ' degrees'
                        except:
                            print("Bad velocity vector while calculating phase angle, event" + str(self.EventNumber))

            else:
                description = str(PlotOptions.EventCounter)
        else:
            description = labelString

        #draw the text
        #note, do not draw anything for chemical burns below 10 m/s
        if not (self.EventType == "chem_burn" and self.DVmagorThrottle < 0.0001):
            x2D, y2D, _ = proj3d.proj_transform(self.SpacecraftState[0],self.SpacecraftState[1],self.SpacecraftState[2], GraphicsObject.get_proj())
            
            if PlotOptions.ShowTextDescriptions:
                self.eventlabel = plt.annotate(description, xycoords = 'data', xy = (x2D, y2D), xytext = (20, (PlotOptions.EventCounter + 20)), textcoords = 'offset points', ha = 'left', va = 'bottom',
                #self.eventlabel = plt.annotate(description, xycoords = 'data', xy = (self.SpacecraftState[0],self.SpacecraftState[1]), xytext = (20, (PlotOptions.EventCounter + 20)), textcoords = 'offset points', ha = 'left', va = 'bottom',
                                               bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 1.0), arrowprops = dict(arrowstyle = '->',
                                               connectionstyle = 'arc3,rad=0'), size=PlotOptions.FontSize, family='serif') # changed to usurp drag-fu
            else:
                self.eventlabel = plt.annotate(description, xycoords = 'data', xy = (x2D, y2D), ha = 'left', va = 'bottom',
                                               bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 1.0), arrowprops = dict(arrowstyle = '->',
                                               connectionstyle = 'arc3,rad=0'), size=PlotOptions.FontSize, family='serif') # changed to usurp drag-fu

            self.AnnotationHelper = self.eventlabel.draggable(use_blit=True)
            self.pcid = GraphicsObject.figure.canvas.mpl_connect('button_press_event', self.ClickAnnotation)
            self.rcid = GraphicsObject.figure.canvas.mpl_connect('button_release_event', self.ReleaseAnnotation)

            #update the event counter only if the event is printable
            PlotOptions.EventCounter += 1

    def UpdateLabelPosition(self, Figure, Axes):        
        from mpl_toolkits.mplot3d import proj3d
        if self.eventlabel != None:
            print('Moving label #' + str(self.EventNumber))
            x2, y2, _ = proj3d.proj_transform(self.SpacecraftState[0],self.SpacecraftState[1],self.SpacecraftState[2], Axes.get_proj())
            self.eventlabel.xy = x2,y2
            try:
                self.eventlabel.update_positions(Figure.canvas.renderer)
            except:
                print('Figure.canvas.renderer is not initialized. Don\'t worry, everything will probably work fine. Don\'t you love bubbles?')

    def ClickAnnotation(self, event):
        if event.inaxes != self.eventlabel.axes: return
        if MissionEvent.lock is not None: return
        contains, attrd = self.eventlabel.contains(event)
        if not contains: return
        self.eventlabel.axes.disable_mouse_rotation()
        MissionEvent.lock = self

    def ReleaseAnnotation(self, event):
        if MissionEvent.lock is not self:
            return
        MissionEvent.lock = None

        if event.inaxes != self.eventlabel.axes: return

        contains, attrd = self.eventlabel.contains(event)
        if not contains: return
        self.eventlabel.axes.mouse_init()

    def AutoTableLine(self, TableFile, PlotOptions, skipNext):
        EventTypeFormatted = self.EventType
        LocationFormatted = self.Location.replace('_',' ')
        if EventTypeFormatted == 'upwr_flyby':
            EventTypeFormatted = 'unpowered flyby'
        elif EventTypeFormatted == 'pwr_flyby':
            EventTypeFormatted = 'powered flyby'
        elif EventTypeFormatted == 'chem_burn':
            EventTypeFormatted = 'chemical burn'
        elif EventTypeFormatted == 'LT_rndzvs':
            EventTypeFormatted = 'rendezvous'
        elif EventTypeFormatted == 'begin_spiral':
            EventTypeFormatted = 'begin spiral'
        elif EventTypeFormatted == 'end_spiral':
            EventTypeFormatted = 'end spiral'
        elif EventTypeFormatted == 'mission_end':
            EventTypeFormatted = 'end of mission'

        plotMe = True

        if  skipNext\
            or (self.EventType in ['chem_burn', 'TCM'] and self.DVmagorThrottle < 0.00001)\
            or self.EventType in ['zeroflyby','coast','nav-coast','SFthrust','SSFthrust','FBLTthrust','FBLTSthrust','PSBIthrust', 'PSFBthrust','match_point','force-coast']\
            or (self.EventType == 'TCM' and PlotOptions.AutoTableTCMcolumn == False):
            plotMe = False

        if " BE" in LocationFormatted:
            plotMe = False
            skipNext = True
        elif EventTypeFormatted == 'periapse':
            EventTypeFormatted = 'unpowered flyby'
            LocationFormatted = self.parent.central_body
            skipNext = True
        else:
            skipNext = False

        #columns are date, event type, event location, vinfinity, deltav, altitude, mass (if applicable)
        if plotMe:
            TableFile.write('       ' + str(PlotOptions.EventCounter) + ' & ' + self.GregorianDate + ' & ' + EventTypeFormatted + ' & ' + LocationFormatted + '&')
            if EventTypeFormatted in ['launch','unpowered flyby','powered flyby','rendezvous','intercept','interface','insertion'] and self.C3 > 0.0:
                TableFile.write("{0:.3f}".format(self.C3 ** 0.5) + ' & ')
            else:
                TableFile.write(' & ')

            if PlotOptions.AutoTableDeltavColumn:
                if self.EventType in ['chem_burn','rendezvous','pwr_flyby','insertion','TCM','departure']:
                    TableFile.write(" {0:.5f}".format(self.DVmagorThrottle) + ' & ')
                    TableFile.write(" {0:.0f}".format(self.Isp) + ' & ')
                else:
                    TableFile.write(' & ')
                    TableFile.write(' & ')
                            
            if PlotOptions.AutoTableAltitudeColumn:
                if EventTypeFormatted in ['unpowered flyby','powered flyby']:
                    TableFile.write(" {0:.0f}".format(self.Altitude) + ' & ')
                else:
                    TableFile.write(' & ')

            if PlotOptions.DisplayEventMass:
                TableFile.write( "{0:.1f}".format(self.Mass) + ' & ')

            if PlotOptions.DisplayArrivalPhaseAngle:
                if self.EventType == 'intercept':
                    R = self.SpacecraftState[0:3]
                    V = self.DeltaVorThrustVectorControl
                    r = math.sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2])
                    v = math.sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2])
                    PhaseAngle = math.acos(-(R[0]*V[0] + R[1]*V[1] + R[2]*V[2]) / (r * v)) * 180.0 / math.pi
                    TableFile.write("{0:.1f}".format(PhaseAngle))

            TableFile.write('\\\\\n')

            PlotOptions.EventCounter += 1

        return skipNext

    def getLaunchHyperbolicCoordinates(self):
        from math import atan2, asin, pi
        rVec = np.array(self.SpacecraftState[0:3])
        vVec = np.array(self.SpacecraftState[3:6])
        hVec = np.cross(rVec, vVec)
        mu = self.parent.mu
        
        r = np.linalg.norm(rVec)
        v = np.linalg.norm(vVec)

        h = np.linalg.norm(hVec)

        #hyperbolic C3 may be different from what is recorded in the event's C3 entry, or you may be computing C3 from a line that doesn't have a C3 entry
        #we're going to return it as an arguement rather than change the event's object data`
        C3 = v**2 - 2 * mu / r

        #now we need to find the v-infinity vector
        #first we need the eccentricity vector
        eVec = ((v**2 - mu / r) * rVec - np.dot(rVec, vVec) * vVec) / mu

        #outgoing asymptote unit vector
        sHat = 1 / (1.0 + C3 * (h / mu)**2) * (C3**0.5 / mu * np.cross(hVec, eVec) - eVec)

        #RLA and DLA
        RLA = atan2(sHat[1], sHat[0]) * 180.0 / pi

        DLA = asin(sHat[2]) * 180.0 / pi

        return C3, RLA, DLA