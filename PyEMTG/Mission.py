import Journey
import MissionEvent
import ThrottleTable
import os
import math
import kepler

import numpy as np

class Mission(object):
    #class to hold mission information

    def __init__(self, input_file_name):
        #default values, set as garbage so that if a .emtg file is not fully populated for some reason then bad cases will be thrown out by any filter
        self.mission_name = 'bob'
        self.total_deterministic_deltav = 1.0e+20
        self.total_statistical_deltav = 1.0e+20
        self.total_flight_time_years = 1.0e+20
        self.first_thrust_event = 1      
        self.spacecraft_dry_mass = 1.0e-20
        self.final_mass_including_propellant_margin = 1.0e-20
        self.total_electric_propellant_including_margin = 1.0e+20
        self.chemical_fuel_with_margin = 1.0e+20
        self.chemical_fuel_used = 1.0e+20
        self.chemical_oxidizer_with_margin = 1.0e+20
        self.chemical_oxidizer_used = 1.0e+20
        self.total_electric_propellant_used = 1.0e+20
        self.spacecraft_dry_mass_margin = -1.0e+20
        self.DecisionVector = []
        self.Xdescriptions = []
        self.Xlowerbounds = []
        self.Xupperbounds = []
        self.user_data = []
        self.worst_constraint = ""
        self.worst_violation = 0.0
        self.objective_value = 1.0e+20

        #now read the file
        self.parse_mission_file(input_file_name)

    def parse_mission_file(self, input_file_name):
        #Step 1: open the file
        if os.path.isfile(input_file_name):
            inputfile = open(input_file_name, "r")
            self.success = 1
        else:
            print("Unable to open " + input_file_name)
            self.success = 0
            return

        #Step 2: parse the file
        self.ActiveJourney = -1
        LastJourney = -1
        self.ChosenJourneys = []
        self.Journeys = []
        self.mu = 0
        linenumber = 0
        for line in inputfile:
            linenumber += 1
            #split the line by colons
            linecell = line.rstrip(" \r\n").split(':')

            if linecell[0] == "Mission":
                self.mission_name = linecell[1].strip('\n')

            elif linecell[0] == "Decision Vector":
                DecisionCell = []
                if ',' in linecell[1]:
                    DecisionCell = linecell[1].split(',')
                else:
                    DecisionCell = linecell[1].split(' ')
                if "" in DecisionCell:
                    DecisionCell.remove('')
                self.DecisionVector = []
                for entry in DecisionCell:
                    self.DecisionVector.append(float(entry))

            elif "Xdescriptions" in linecell[0]:
                linecellcomma = line.split(',')
                self.Xdescriptions = []
                for entry in linecellcomma[1:]:
                    self.Xdescriptions.append(entry.rstrip("\r\n"))

            elif line[0:12] == 'Xupperbounds':
                if ',' in line:
                    boundcell = line.split(',')[1:]
                else:
                    boundcell = line.split(' ')[1:]
                self.Xupperbounds = []
                for entry in boundcell:
                    self.Xupperbounds.append(float(entry))

            elif line[0:12] == 'Xlowerbounds':
                if ',' in line:
                    boundcell = line.split(',')[1:]
                else:
                    boundcell = line.split(' ')[1:]
                self.Xlowerbounds = []
                for entry in boundcell:
                    self.Xlowerbounds.append(float(entry))   
                    
            elif "Constraint_Vector" in linecell[0]:
                DecisionCell = []
                if ',' in linecell[0]:
                    DecisionCell = line.rstrip("/r/n ").split(',')[1:]
                else:
                    DecisionCell = line.rstrip("/r/n ").split(' ')[1:]
                if "" in DecisionCell:
                    DecisionCell.remove('')
                self.ConstraintVector = []
                for entry in DecisionCell:
                    self.ConstraintVector.append(float(entry))
                    
            elif "Fdescriptions" in linecell[0]:
                linecellcomma = line.split(',')
                self.Fdescriptions = []
                for entry in linecellcomma[1:]:
                    self.Fdescriptions.append(entry.rstrip("\r\n"))                 

            elif line[0:12] == 'Fupperbounds':
                if ',' in line:
                    boundcell = line.split(',')[1:]
                else:
                    boundcell = line.split(' ')[1:]
                self.Fupperbounds = []
                for entry in boundcell:
                    self.Fupperbounds.append(float(entry))                

            elif line[0:12] == 'Flowerbounds':
                if ',' in line:
                    boundcell = line.split(',')[1:]
                else:
                    boundcell = line.split(' ')[1:]
                self.Flowerbounds = []
                for entry in boundcell:
                    self.Flowerbounds.append(float(entry))      

            elif linecell[0] == "Journey":
                LastJourney += 1
                self.ActiveJourney = LastJourney
                self.Journeys.append(Journey.Journey(self))
            
            elif linecell[0] == "user_data":
                full_notes = line.lstrip("user_data:").lstrip(" ").rstrip("\r\n ")
                                        
                self.user_data = dict()
                
                if full_notes != '':
                    if ":" in full_notes:
                        full_notes = full_notes.split(":")
                        for note in full_notes:
                            var = note.lstrip('("').split(",")[0].rstrip('"')
                            val = eval(note.lstrip("(").lstrip(var + '"').lstrip(", ").rstrip(") "))
                            self.user_data.update({var:val})
                    else:
                        var = full_notes.lstrip('("').split(",")[0].rstrip('"')
                        val = eval(full_notes.lstrip("(").lstrip(var + '"').lstrip(", ").rstrip(") "))
                        self.user_data.update({var:val})

            if self.ActiveJourney >= 0:
                self.Journeys[self.ActiveJourney].journey_number = self.ActiveJourney
                if linecell[0] == "Journey name":
                    self.Journeys[self.ActiveJourney].journey_name = linecell[1].strip()
                elif linecell[0] == "Central Body":
                    self.Journeys[self.ActiveJourney].central_body = linecell[1].strip()
                elif linecell[0] == "Thruster duty cycle":
                    self.Journeys[self.ActiveJourney].thruster_duty_cycle = eval(linecell[1])
                elif linecell[0] == "Radius (km)":
                    self.Journeys[self.ActiveJourney].central_body_radius = eval(linecell[1])
                elif linecell[0] == "mu (km^3/s^2)":
                    self.Journeys[self.ActiveJourney].mu = eval(linecell[1])
                elif linecell[0] == "Characteristic length unit (km)":
                    self.Journeys[self.ActiveJourney].LU = eval(linecell[1])
                    self.Journeys[self.ActiveJourney].TU = math.sqrt(self.Journeys[self.ActiveJourney].LU**3 / self.Journeys[self.ActiveJourney].mu)
                elif linecell[0] == "Frame":
                    self.Journeys[self.ActiveJourney].state_frame = linecell[1].strip()
                elif linecell[0] == 'alpha0':
                    self.Journeys[self.ActiveJourney].alpha0 = eval(linecell[1])
                elif linecell[0] == 'delta0':
                    self.Journeys[self.ActiveJourney].delta0 = eval(linecell[1])
                elif linecell[0] == "Boundary":
                    boundary_elements = linecell[1].split(' ')
                    boundary_state = [float(boundary_elements[2]), float(boundary_elements[3]), float(boundary_elements[4]), float(boundary_elements[5]), float(boundary_elements[6]), float(boundary_elements[7])]
                    self.Journeys[self.ActiveJourney].boundary_states.append(boundary_state)
                elif linecell[0] == "Flyby":
                    flyby_elements = linecell[1].split(' ')
                    flyby_state = [float(flyby_elements[3]), float(flyby_elements[4]), float(flyby_elements[5]), float(flyby_elements[6]), float(flyby_elements[7]), float(flyby_elements[8])]
                    self.Journeys[self.ActiveJourney].flyby_periapse_states = []
                    self.Journeys[self.ActiveJourney].flyby_periapse_states.append(flyby_state)
                elif linecell[0] == "End journey":
                    self.ActiveJourney = -1
                elif linecell[0] == "Journey entry/landing velocity with respect to atmosphere (km/s)":
                    self.journey_entryinterface_velocity_with_respect_to_rotating_atmosphere = float(linecell[1])
                elif linecell[0] == "Journey initial mass increment":
                    self.Journeys[self.ActiveJourney].journey_mass_increment = float(linecell[1].split(' ')[1])
                elif linecell[0] == "Journey final mass increment":
                    self.Journeys[self.ActiveJourney].journey_mass_increment = float(linecell[1].split(' ')[1])

                elif "Journey post-arrival delta-v" in linecell[0]:
                    self.Journeys[self.ActiveJourney].journey_post_arrival_deltav = float(linecell[1].split(' ')[1])
                elif "Journey post-arrival propellant consumed" in linecell[0]:
                    self.Journeys[self.ActiveJourney].journey_post_arrival_propellant_used = float(linecell[1].split(' ')[1])

                elif linecell[0] == "Journey entry/landing interface point (BCI frame)":
                    nextline = next(inputfile).split(' ')
                    self.Journeys[self.ActiveJourney].journey_entry_landing_interface_point = []
                    for entry in nextline:
                        self.Journeys[self.ActiveJourney].journey_entry_landing_interface_point.append(float(entry.strip('\n')))
                elif linecell[0] == "Journey entry/landing interface latitude (BCF spherical cow, degrees)":
                    self.Journeys[self.ActiveJourney].entry_latitude = float(linecell[1].split(' ')[1])
                elif linecell[0] == "Journey entry/landing interface longitude (BCF spherical cow, degrees)":
                    self.Journeys[self.ActiveJourney].entry_longitude = float(linecell[1].split(' ')[1])
                elif linecell[0] == "Journey entry/landing sun angle (degrees)":
                    self.Journeys[self.ActiveJourney].entry_sun_angle = float(linecell[1].split(' ')[1])
                elif linecell[0] == "Spacecraft-sun-Earth angle at arrival":
                    self.Journeys[self.ActiveJourney].journey_arrival_spacecraft_sun_Earth_angle = float(linecell[1].split(' ')[1])
                elif linecell[0] == "Journey electric propellant used":
                    self.Journeys[self.ActiveJourney].journey_electric_propellant_used = float(linecell[1].split(' ')[1])
                elif linecell[0] == "Journey chemical propellant used":
                    self.Journeys[self.ActiveJourney].journey_chemical_propellant_used = float(linecell[1].split(' ')[1])

                else:
                    #parse event lines
                    inputcell = line.split('|')
                    if inputcell[0].strip(' ').isdigit():
                        tempEvent = MissionEvent.MissionEvent( self.Journeys[self.ActiveJourney], inputcell)
                        self.Journeys[self.ActiveJourney].missionevents.append(tempEvent)
                        del tempEvent
                    del inputcell

            elif linecell[0] == "Total deterministic deltav (km/s)":
                self.total_deterministic_deltav = float(linecell[1])

            elif linecell[0] == "Total statistical deltav (km/s)":
                self.total_statistical_deltav = float(linecell[1])

            elif linecell[0] == "Flight time (y)":
                self.total_flight_time_years = float(linecell[1])
                
            elif linecell[0] == "First thrusting in control step":
                self.first_thrust_event = float(linecell[1])        

            elif "Worst constraint is" in linecell[0]:
                self.worst_constraint = line.split('is ')[1].lstrip(' ').strip('\n')

            elif "with violation" in linecell[0]:
                self.worst_violation = float(line.split("violation ")[1])

            elif linecell[0] == "Spacecraft":
                if linecell[1].lstrip(" ") == "Dry mass (kg)":
                    self.spacecraft_dry_mass = float(linecell[2])

                elif linecell[1].lstrip(" ") == "Final mass including propellant margin (kg)":
                    self.final_mass_including_propellant_margin = float(linecell[2]) 

                elif linecell[1].lstrip(" ") == "Total chemical fuel (kg)":
                    self.chemical_fuel_with_margin = float(linecell[2]) 
                    
                elif linecell[1].lstrip(" ") == "Total chemical oxidizer (kg)":
                    self.chemical_oxidizer_with_margin = float(linecell[2])
                    
                elif linecell[1].lstrip(" ") == "Chemical fuel used (kg)":
                    self.chemical_fuel_used = float(linecell[2]) 
                    
                elif linecell[1].lstrip(" ") == "Chemical oxidizer used (kg)":
                    self.chemical_oxidizer_used = float(linecell[2])
                    
                elif linecell[1].lstrip(" ") == "Total electric propellant (kg)":
                    self.total_electric_propellant_including_margin = float(linecell[2]) 
                    
                elif linecell[1].lstrip(" ") == "Electric propellant used (kg)":
                    self.total_electric_propellant_used = float(linecell[2]) 
                    
                elif linecell[1].lstrip(" ") == "Dry mass margin":
                    self.spacecraft_dry_mass_margin = float(linecell[1])
                
            elif line[0:4] == 'J = ':
                self.objective_value = float(line.replace('J = ',''))
        
        inputfile.close()

    def Generate1Dto2DInitialGuess(self,throttle_table_file,command_var):
        DecVec = []
        decvec_idx = 0
            
        thruster_throttle_table = ThrottleTable.ThrottleTable(throttle_table_file)
    
        if command_var == 0:
            max_mdot = thruster_throttle_table.max_mdot
        elif command_var == 1:
            max_voltage = thruster_throttle_table.max_voltage           
    
        for j in self.Journeys:
            for e in j.missionevents:
                if e.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'end_spiral','PSFBthrust']:
                    Rmat = np.array([[1.0,0.0,0.0],[0.0,0.91748206,-0.397777156],[0.0,0.39777716,0.9174820621]])
                
                    for i_control in range(0,3):
                        DecVec.append(self.DecisionVector[decvec_idx])
                        decvec_idx += 1
                
                    if e.Number_of_Active_Engines:
                        ThrottleLevel = thruster_throttle_table.find_nearest_throttle_setting_2D(e.MassFlowRate / e.Number_of_Active_Engines * 1e6 / e.DVmagorThrottle, e.AvailableThrust / e.Number_of_Active_Engines * 1e3)
                                                                                    
                        if command_var == 0:
                            DecVec.append((thruster_throttle_table.get_mdot(ThrottleLevel) + .05)/(max_mdot*1.1))
                        if command_var == 1:
                            DecVec.append((thruster_throttle_table.get_voltage(ThrottleLevel) + 5)/(max_voltage*1.1))
                    else:
                        DecVec.append(0.0)
            
                elif e.EventType in ['coast', 'nav-coast']:

                    for i_control in range(0,3):
                        DecVec.append(self.DecisionVector[decvec_idx])
                        decvec_idx += 1
                        
                    DecVec.append(0.0)
                        
                elif e.EventType == 'launch' or e.EventType == 'launch' or e.EventType == 'pwr_flyby' or e.EventType == 'upwr_flyby' or e.EventType == 'departure':
                    
                    # add header
                    while "u_" not in self.Xdescriptions[decvec_idx]:
                        DecVec.append(self.DecisionVector[decvec_idx])
                        decvec_idx += 1                        
    
        return DecVec

    def Generate2Dto1DInitialGuess(self):

        DecVec = []
        
        for i in range(0,len(self.DecisionVector)):
            if 'mdot' not in self.Xdescriptions[i] and 'voltage' not in self.Xdescriptions[i]:
                DecVec.append(self.DecisionVector[i])
                
        return DecVec

    def PlotMission(self, PlotOptions, ShowFigure=True, CustomLabels=None):
        import matplotlib
        import matplotlib.pyplot
        matplotlib.use('wxagg')
        self.MissionFigure = matplotlib.pyplot.figure()
        matplotlib.rcParams.update({'font.size': PlotOptions.FontSize})
        self.MissionFigure.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
        self.MissionAxes = self.MissionFigure.add_axes([0,0,1,1], projection='3d')
        self.MissionAxes.xaxis.pane.fill = False
        self.MissionAxes.yaxis.pane.fill = False
        self.MissionAxes.zaxis.pane.fill = False
        self.MissionAxes.xaxis._axinfo['tick']['inward_factor'] = 0
        self.MissionAxes.xaxis._axinfo['tick']['outward_factor'] = 0.2
        self.MissionAxes.yaxis._axinfo['tick']['inward_factor'] = 0
        self.MissionAxes.yaxis._axinfo['tick']['outward_factor'] = 0.2
        self.MissionAxes.zaxis._axinfo['tick']['inward_factor'] = 0
        self.MissionAxes.zaxis._axinfo['tick']['outward_factor'] = 0.2
        #self.MissionAxes.set_aspect('equal') # commented out becausee not supported by matplotlib as of Feb 2019
        self.MissionAxes.set_xlabel('\nx (km)')
        self.MissionAxes.set_ylabel('\ny (km)')
        # self.MissionAxes.set_zlabel('\nz (km)')
        self.MissionAxes.w_zaxis.line.set_lw(0.)
        self.MissionAxes.set_zticks([])
        self.MissionAxes.autoscale_view(tight=True, scalex=True, scaley=True, scalez=True)
        self.MissionAxes.view_init(elev=90, azim=-90)
        self.MissionAxes.grid(b=False)
        

        #reset the event label counter
        PlotOptions.EventCounter = 1

        #load DEfile if we have to
        mySpiceHandler = []
        if PlotOptions.PlotCentralBody != 'Journey central body':
            from sys import path
            from os.path import dirname, join
            from pathlib import Path

            PyEMTG_path = str(Path().absolute()).replace('\\','/')
            path.append(PyEMTG_path + '/SpiceyPy_Utilities')
            import SpiceyPy_Utilities
            mySpiceHandler = SpiceyPy_Utilities.SpiceHandler(PlotOptions.SPICEpath)
            mySpiceHandler.loadSpiceFiles()

        #plot journey boundary orbits
        if PlotOptions.ShowBoundaryOrbits:
            for ActiveJourney in self.ChosenJourneys:
                self.ActiveJourney = ActiveJourney
                if self.ActiveJourney <= len(self.Journeys) - 1:
                    self.Journeys[self.ActiveJourney].PlotJourneyBoundaryOrbits(self.MissionAxes, PlotOptions)
                else:
                    for CurrentJourney in self.Journeys:
                        CurrentJourney.PlotJourneyBoundaryOrbits(self.MissionAxes, PlotOptions)

        #plot journey central body orbits
        if PlotOptions.ShowCentralBodyOrbits:
            for ActiveJourney in self.ChosenJourneys:
                self.ActiveJourney = ActiveJourney
                if self.ActiveJourney <= len(self.Journeys) - 1:
                    self.Journeys[self.ActiveJourney].PlotJourneyCentralBodyOrbits(self.MissionAxes, PlotOptions)
                else:
                    for CurrentJourney in self.Journeys:
                        CurrentJourney.PlotJourneyCentralBodyOrbits(self.MissionAxes, PlotOptions)

        #plot the trajectory
        for ActiveJourney in self.ChosenJourneys:
            self.ActiveJourney = ActiveJourney
            if self.ActiveJourney <= len(self.Journeys) - 1:
                self.Journeys[self.ActiveJourney].PlotJourney(self.MissionAxes, PlotOptions, CustomLabels)
            else:
                for CurrentJourney in self.Journeys:
                    CurrentJourney.PlotJourney(self.MissionAxes, PlotOptions, CustomLabels)

        #unload DEfile if we have to
        if PlotOptions.PlotCentralBody != 'Journey central body':
            mySpiceHandler.unloadSpiceFiles()

        #plot the central body
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        central_body_radius = 1.0
        if self.ActiveJourney <= len(self.Journeys) - 1:
            central_body_radius = self.Journeys[self.ActiveJourney].central_body_radius
        else:
            central_body_radius = self.Journeys[0].central_body_radius
        x = central_body_radius * np.outer(np.cos(u), np.sin(v))
        y = central_body_radius * np.outer(np.sin(u), np.sin(v))
        z = central_body_radius * np.outer(np.ones(np.size(u)), np.cos(v))
        self.MissionAxes.plot_surface(x, y, z,  rstride=10, cstride=10, color='DarkOrange', linewidth=0.5, edgecolors='k')

        X = self.MissionAxes.get_xlim()
        Y = self.MissionAxes.get_ylim()
        Z = self.MissionAxes.get_zlim()
        
        Xrange = X[1] - X[0]
        Yrange = Y[1] - Y[0]
        Zrange = Z[1] - Z[0]
        
        # Create cubic bounding box to simulate equal aspect ratio
        max_range = np.array([Xrange, Yrange, Zrange]).max()
        Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(Xrange)
        Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Yrange)
        Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Zrange)
        # Comment or uncomment following both lines to test the fake bounding box:
        for xb, yb, zb in zip(Xb, Yb, Zb):
           self.MissionAxes.plot([xb], [yb], [zb], 'w')
        
        if PlotOptions.LabelEvents or CustomLabels != None:
            self.MissionFigure.canvas.mpl_connect('button_release_event', self.UpdateLabelPositionsEvent)
            
        if PlotOptions.ShowPlotAxes == False:
            matplotlib.pyplot.axis('off')
            matplotlib.pyplot.margins(0, 0, 0)
            matplotlib.pyplot.gca().xaxis.set_major_locator(matplotlib.pyplot.NullLocator())
            matplotlib.pyplot.gca().yaxis.set_major_locator(matplotlib.pyplot.NullLocator())

        self.MissionAxes.autoscale_view(tight=False, scalex=True, scaley=True, scalez=True)
        self.MissionAxes.autoscale(enable=False,axis='both')  #you will need this line to change the Z-axis
        self.MissionAxes.set_xbound(-max_range, max_range)
        self.MissionAxes.set_ybound(-max_range, max_range)
        self.MissionAxes.set_zbound(-max_range, max_range)
        
        if ShowFigure:
            self.MissionFigure.show()

        #these lines have to be run AFTER MissionFigure.show() is called
        if (PlotOptions.LabelEvents or CustomLabels != None) and ShowFigure == True:
            self.UpdateLabelPositions()

        return self.MissionFigure

    def UpdateLabelPositionsEvent(self, e):
        self.UpdateLabelPositions()

    def UpdateLabelPositions(self):    
        for ActiveJourney in self.ChosenJourneys:
            self.ActiveJourney = ActiveJourney
            if self.ActiveJourney <= len(self.Journeys) - 1:
                self.Journeys[self.ActiveJourney].UpdateLabelPosition(self.MissionFigure, self.MissionAxes)
            else:
                for CurrentJourney in self.Journeys:
                    CurrentJourney.UpdateLabelPosition(self.MissionFigure, self.MissionAxes)

    def GenerateDataPlot(self, PlotOptions):
        import matplotlib.pyplot
        if PlotOptions.PlotR or PlotOptions.PlotV or PlotOptions.PlotThrust or PlotOptions.PlotIsp or PlotOptions.PlotMdot or PlotOptions.PlotEfficiency or PlotOptions.PlotThrottle or PlotOptions.PlotPower or PlotOptions.PlotGamma or PlotOptions.PlotDelta or PlotOptions.PlotArray_Thrust_Angle or PlotOptions.PlotMass or PlotOptions.PlotNumberOfEngines or PlotOptions.PlotActivePower or PlotOptions.PlotWasteHeat or PlotOptions.PlotEarthDistance or PlotOptions.PlotSunSpacecraftEarthAngle or PlotOptions.PlotSpacecraftViewingAngle or PlotOptions.PlotThrottleLevel:
            self.DataFigure = matplotlib.pyplot.figure()
            self.DataAxesLeft = self.DataFigure.add_axes([0.1, 0.1, 0.8, 0.8])
            self.DataAxesRight = self.DataAxesLeft.twinx()
            self.MouseAxes = self.DataAxesLeft.figure.add_axes(self.DataAxesLeft.get_position(True), sharex=self.DataAxesLeft, sharey=self.DataAxesLeft, frameon=False)
            self.MouseAxes.xaxis.set_visible(False)
            self.MouseAxes.yaxis.set_visible(False)
            matplotlib.rcParams.update({'font.size': PlotOptions.FontSize})
            matplotlib.rcParams.update({'font.family': 'Times New Roman'})
            
            font = {'family': 'Times New Roman',
                    'color':  'black',
                    'weight': 'normal',
                    'size': PlotOptions.FontSize}
            
            needLegend = True
            for ActiveJourney in self.ChosenJourneys:
                self.ActiveJourney = ActiveJourney
                if self.ActiveJourney <= len(self.Journeys) - 1:
                    self.Journeys[self.ActiveJourney].GenerateJourneyDataPlot(self.DataAxesLeft, self.DataAxesRight, PlotOptions, needLegend)
                    needLegend = False

                else:
                    for journey in self.Journeys:
                        journey.GenerateJourneyDataPlot(self.DataAxesLeft, self.DataAxesRight, PlotOptions, needLegend)
                        needLegend = False
            
            needLegend = True
            if PlotOptions.PlotCriticalEvents:
                #determine the Y limits of the plot
                YboundsLeft = self.DataAxesLeft.get_ylim()
                YboundsRight = self.DataAxesRight.get_ylim()
                Ybounds = (min(YboundsLeft[0], YboundsRight[0]), max(YboundsLeft[1], YboundsRight[1]))
                
                for ActiveJourney in self.ChosenJourneys:
                    self.ActiveJourney = ActiveJourney
                    if self.ActiveJourney <= len(self.Journeys) - 1:
                        self.Journeys[self.ActiveJourney].PlotPhaseBoundariesOnDataPlot(self.DataAxesLeft, PlotOptions, needLegend, Ybounds)
                        needLegend = False
                    else:
                        for journey in self.Journeys:
                            journey.PlotPhaseBoundariesOnDataPlot(self.DataAxesLeft, PlotOptions, needLegend, Ybounds)
                            needLegend = False

                     
            self.DataAxesLeft.set_xlabel('Epoch', fontdict=font)
            def format_date(x, pos=None):
                import pylab
                return pylab.num2date(x).strftime('%m-%d-%Y')

            
            self.DataFigure.autofmt_xdate()

            h1 = []
            h2 = []
            l1 = []
            l2 = []
            
            
            if PlotOptions.PlotGamma or PlotOptions.PlotDelta or PlotOptions.PlotArray_Thrust_Angle or PlotOptions.PlotSunSpacecraftEarthAngle or PlotOptions.PlotCB_thrust_angle or PlotOptions.PlotSunBoresightAngle:
                h2, l2 = self.DataAxesRight.get_legend_handles_labels()
                self.DataAxesRight.set_ylabel('Angle Metric', fontdict=font)
                self.DataAxesRight.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(format_date))
                try:
                    for label in self.DataAxesRight.get_xticklabels():
                        label.set_rotation(30)
                except:
                    print("Failed to rotate x tick labels.")
                #self.DataAxesRight.grid(b=True, ls='-')
            else:
                yticks = self.DataAxesRight.yaxis.get_major_ticks()
                for tick in yticks:
                    tick.set_visible(False)

            if PlotOptions.PlotR or PlotOptions.PlotV or PlotOptions.PlotThrust or PlotOptions.PlotIsp or PlotOptions.PlotMdot or PlotOptions.PlotEfficiency or PlotOptions.PlotThrottle or PlotOptions.PlotPower or PlotOptions.PlotMass or PlotOptions.PlotNumberOfEngines or PlotOptions.PlotActivePower or PlotOptions.PlotWasteHeat or PlotOptions.PlotEarthDistance or PlotOptions.PlotSpacecraftViewingAngle or PlotOptions.PlotThrottleLevel:
                h1, l1 = self.DataAxesLeft.get_legend_handles_labels()
                self.DataAxesLeft.set_ylabel('Scalar Metric', fontdict=font)
                self.DataAxesLeft.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(format_date))
                try:
                    for label in self.DataAxesLeft.get_xticklabels():
                        label.set_rotation(30)
                except:
                    print("Failed to rotate x tick labels.")
                #self.DataAxesLeft.grid(b=True, ls='-')
            else:
                yticks = self.DataAxesLeft.yaxis.get_major_ticks()
                for tick in yticks:
                    tick.set_visible(False)

            leg = self.DataAxesLeft.legend(h1+h2, l1+l2, loc='upper center', fancybox=True, fontsize=PlotOptions.FontSize)    
            #leg = self.DataAxesRight.legend(h1+h2, l1+l2, loc='upper left', fancybox=True)
            leg.get_frame().set_alpha(0.5)
            leg.set_draggable(True, use_blit=True)
            self.DataFigure.show()

    def WriteDataReport(self, PlotOptions, reportfilename):
        if PlotOptions.PlotR or PlotOptions.PlotV or PlotOptions.PlotThrust or PlotOptions.PlotIsp or PlotOptions.PlotMdot or PlotOptions.PlotEfficiency or PlotOptions.PlotThrottle or PlotOptions.PlotPower or PlotOptions.PlotGamma or PlotOptions.PlotDelta or PlotOptions.PlotArray_Thrust_Angle or PlotOptions.PlotMass or PlotOptions.PlotNumberOfEngines or PlotOptions.PlotActivePower or PlotOptions.PlotWasteHeat or PlotOptions.PlotEarthDistance or PlotOptions.PlotSunSpacecraftEarthAngle or PlotOptions.PlotSpacecraftViewingAngle or PlotOptions.PlotThrottleLevel:
            
            reportfile = open(reportfilename, 'w')
            reportfile.write('#EMTG systems report file\n')
            reportfile.write('\n')

            
            for ActiveJourney in self.ChosenJourneys:
                self.ActiveJourney = ActiveJourney
                if self.ActiveJourney <= len(self.Journeys) - 1:
                    reportfile.write('Journey ' + str(self.ActiveJourney) + ': ' + self.Journeys[self.ActiveJourney].journey_name + '\n')
                    self.Journeys[self.ActiveJourney].WriteDateReport(PlotOptions, reportfile)

                else:
                    for j in range(0, len(self.Journeys)):
                        reportfile.write('Journey ' + str(j) + ': ' + self.Journeys[j].journey_name + '\n')
                        self.Journeys[j].WriteDateReport(PlotOptions, reportfile)
                        reportfile.write('\n')

    def BubbleSearch(self, BubbleOptions, outputfilename, OutputWindow = []):
        import astropy
        import SmallBodyList
        
        try:
            import spiceypy as spice
        except:
            print("spiceypy not available")
        
        if BubbleOptions.ifSpice:
            loaded_SPICE_files = []
            BubbleOptions.SPICE_IDs = []
            for IDfile in BubbleOptions.SPICEID_file:
                file = open(IDfile,'r')
                for line in file:
                    BubbleOptions.SPICE_IDs.append(int(line))
                file.close()

            for file in BubbleOptions.spiceFiles:
                if file not in loaded_SPICE_files:
                    loaded_SPICE_files.append(file)
                    spice.furnsh(file)


        bubblefile = open(outputfilename, mode = 'w')
        bubblefile.write('#Bubble file for ' + self.mission_name + '\n')
        bubblefile.write('#\n')
        bubblefile.write('#\n')
        SmallBodyDatabase = SmallBodyList.SmallBodyList(BubbleOptions.smallbodyfile, BubbleOptions.mu, BubbleOptions.LU)
        
        
        for ActiveJourney in self.ChosenJourneys:
            self.ActiveJourney = ActiveJourney
            if self.ActiveJourney <= len(self.Journeys) - 1:
                bubblefile.write('Journey ' + str(self.ActiveJourney+1) + ': ' + self.Journeys[self.ActiveJourney].journey_name + '\n')
                bubblefile.write('Julian Date, Gregorian Date, Body Name, SPICE ID, Body SMA, Tholen spectral type, SMASSII spectral type, Absolute Magnitude, Diameter, Distance (km), Distance (LU), Relative Velocity (km/s), Guesstimated intercept delta-v (km/s)\n')

                #the first thing we need to do is identify the major events in the journey so that, for each missionevent, we know which flybys and/or maneuvers bracket it
                major_event_list = []
                for event in self.Journeys[self.ActiveJourney].missionevents:
                    if not event.EventType in ['coast', 'force-coast','SFthrust','PSBIthrust','FBLTthrust','PSFBthrust','nav-coast']:
                        major_event_list.append((event.EventType, event.JulianDate))

                for event in self.Journeys[self.ActiveJourney].missionevents:
                    # print event.JulianDate
                    if hasattr(OutputWindow,'WriteText'):
                        OutputWindow.WriteText("Searching for targets of opportunity on JD " + str(event.JulianDate) + "\n")
                
                    event_body_list = SmallBodyDatabase.find_bubble_targets(event.SpacecraftState, event.JulianDate, BubbleOptions.RelativePositionFilterMagnitude, BubbleOptions.RelativeVelocityFilterMagnitude, BubbleOptions.MaximumMagnitude, BubbleOptions.ifSpice, BubbleOptions.spiceFiles, BubbleOptions.SPICE_IDs)

                    #determine which events bracket us
                    event_bracket = [0.0,0.0]
                    for eventIndex in range(0, len(major_event_list)-1):
                        if major_event_list[eventIndex][1] <= event.JulianDate and major_event_list[eventIndex + 1][1] >= event.JulianDate:
                            event_bracket[0] = major_event_list[eventIndex][1]
                            event_bracket[1] = major_event_list[eventIndex + 1][1]

                    for body in event_body_list:
                        if body.Tholen == '-1':
                            Tholen = ''
                        else:
                            Tholen = body.Tholen

                        if body.SMASSII == '-1':
                            SMASSII = ''
                        else:
                            SMASSII = body.SMASSII

                        if body.Diameter == -1:
                            Diameter = ''
                        else:
                            Diameter = str(body.Diameter)

                        #compute the Sutter-method guesstimated intercept delta-v relative to the forward and backward bracketing event and take the worst case
                        if abs(event.JulianDate - event_bracket[0]) > 1.0e-10:
                            forward_guesstimated_deltav = body.RelativePositionMagnitude / (event.JulianDate - event_bracket[0]) / 86400.0
                        else:
                            forward_guesstimated_deltav = 1.0e+100
                        if abs(event_bracket[1] - event.JulianDate) > 1.0e-10:
                            backward_guesstimated_deltav = body.RelativePositionMagnitude / (event_bracket[1] - event.JulianDate) / 86400.0
                        else:
                            backward_guesstimated_deltav = 1.0e+100

                        bubblefile.write(str(event.JulianDate) + ',' + event.GregorianDate + ',' + body.name + ',' + str(body.SPICEID) + ',' + str(body.SMA) + ',' + Tholen + ',' + SMASSII + ',' + str(body.H) + ',' + Diameter + ',' + str(body.RelativePositionMagnitude) + ',' + str(body.RelativePositionMagnitude / BubbleOptions.LU) + ',' + str(body.RelativeVelocityMagnitude) + ',' + str(max(forward_guesstimated_deltav, backward_guesstimated_deltav)) + '\n')
                
                    if hasattr(OutputWindow,'WriteText'):
                        OutputWindow.WriteText("Found " + str(len(event_body_list)) + " candidates\n")
                    bubblefile.flush()

                #do we want to keep going after the end of the mission?
                if (self.ActiveJourney == len(self.Journeys) - 1 and BubbleOptions.CheckForEncountersAfterMissionEnd):
                    bubblefile.write('Post-mission encounters\n')
                    bubblefile.write('Julian Date, Gregorian Date, Body Name, SPICE ID, Body SMA, Tholen spectral type, SMASSII spectral type, Absolute Magnitude, Diameter, Distance (km), Distance (LU), Relative Velocity (km/s), Guesstimated intercept delta-v (km/s)\n')

                    PostMissionState = self.Journeys[-1].missionevents[-1].SpacecraftState
                    PostMissionCOE = kepler.cart2kep(PostMissionState[0:3],PostMissionState[3:6],BubbleOptions.mu)
                    MissionEndEpoch = self.Journeys[-1].missionevents[-1].JulianDate
                    PostMissionPropagationStep = BubbleOptions.PostMissionCheckDuration / BubbleOptions.PostMissionCheckSteps

                    for step in range(0, BubbleOptions.PostMissionCheckSteps):
                        CurrentEpoch = MissionEndEpoch + (step + 1) * PostMissionPropagationStep
                        #propagate the spacecraft to this epoch
                        r, v = kepler.kepler(PostMissionCOE[0], PostMissionCOE[1], PostMissionCOE[2], PostMissionCOE[3], PostMissionCOE[4], PostMissionCOE[5], MissionEndEpoch, CurrentEpoch, BubbleOptions.mu)
                        state = [r[0], r[1], r[2], v[0], v[1], v[2]]

                        #check for bodies
                        OutputWindow.WriteText("Searching for targets of opportunity on JD " + str(CurrentEpoch) + "\n")
                        event_body_list = SmallBodyDatabase.find_bubble_targets(state, CurrentEpoch, BubbleOptions.RelativePositionFilterMagnitude, BubbleOptions.RelativeVelocityFilterMagnitude, BubbleOptions.MaximumMagnitude, BubbleOptions.ifSpice, BubbleOptions.spiceFiles, BubbleOptions.SPICE_IDs)

                        #write fields
                        for body in event_body_list:
                            if body.Tholen == '-1':
                                Tholen = ''
                            else:
                                Tholen = body.Tholen

                            if body.SMASSII == '-1':
                                SMASSII = ''
                            else:
                                SMASSII = body.SMASSII

                            if body.Diameter == -1:
                                Diameter = ''
                            else:
                                Diameter = str(body.Diameter)

                            #compute the Sutter-method guesstimated intercept delta-v relative to the forward and backward bracketing event and take the worst case
                        
                            forward_guesstimated_deltav = body.RelativePositionMagnitude / (CurrentEpoch - MissionEndEpoch + 1.0) / 86400.0

                            try:
                                myDate = astropy.time.Time(CurrentEpoch, format='jd').datetime.strftime('%m/%d/%Y')
                            except ValueError:
                                myDate = 'idate'

                            bubblefile.write(str(CurrentEpoch) + ',' + myDate + ',' + body.name + ',' + str(body.SPICEID) + ',' + str(body.SMA) + ',' + Tholen + ',' + SMASSII + ',' + str(body.H) + ',' + Diameter + ',' + str(body.RelativePositionMagnitude) + ',' + str(body.RelativePositionMagnitude / BubbleOptions.LU) + ',' + str(body.RelativeVelocityMagnitude) + ',' + str(forward_guesstimated_deltav) + '\n')
                    
                        OutputWindow.WriteText("Found " + str(len(event_body_list)) + " candidates\n")
                        bubblefile.flush()

                    bubblefile.write('\n')
                    

            else:
                for j in range(0, len(self.Journeys)):
                    bubblefile.write('Journey ' + str(j+1) + ': ' + self.Journeys[j].journey_name + '\n')
                    bubblefile.write('Julian Date, Gregorian Date, Body Name, SPICE ID, Body SMA, Tholen spectral type, SMASSII spectral type, Absolute Magnitude, Diameter, Distance (km), Distance (LU), Relative Velocity (km/s), Guesstimated intercept delta-v (km/s)\n')

                    #the first thing we need to do is identify the major events in the journey so that, for each missionevent, we know which flybys and/or maneuvers bracket it
                    major_event_list = []
                    for event in self.Journeys[j].missionevents:
                        if not event.EventType in ['coast', 'force-coast','SFthrust','PSBIthrust','FBLTthrust','PSFBthrust','nav-coast']:
                            major_event_list.append((event.EventType, event.JulianDate))

                    for event in self.Journeys[j].missionevents:
                        # print event.JulianDate
                        if hasattr(OutputWindow,'WriteText'):
                            OutputWindow.WriteText("Searching for targets of opportunity on JD " + str(event.JulianDate) + "\n")
                        event_body_list = SmallBodyDatabase.find_bubble_targets(event.SpacecraftState, event.JulianDate, BubbleOptions.RelativePositionFilterMagnitude, BubbleOptions.RelativeVelocityFilterMagnitude, BubbleOptions.MaximumMagnitude, BubbleOptions.ifSpice, BubbleOptions.spiceFiles, BubbleOptions.SPICE_IDs)

                        #determine which events bracket us
                        event_bracket = [0.0,0.0]
                        for eventIndex in range(0, len(major_event_list)-1):
                            if major_event_list[eventIndex][1] <= event.JulianDate and major_event_list[eventIndex + 1][1] >= event.JulianDate:
                                event_bracket[0] = major_event_list[eventIndex][1]
                                event_bracket[1] = major_event_list[eventIndex + 1][1]
                                break
                        if hasattr(OutputWindow,'WriteText'):
                            OutputWindow.WriteText("bracketing events are at JD [" + str(event_bracket[0]) + ', ' + str(event_bracket[1]) + '\n')

                        for body in event_body_list:
                            if body.Tholen == '-1':
                                Tholen = ''
                            else:
                                Tholen = body.Tholen

                            if body.SMASSII == '-1':
                                SMASSII = ''
                            else:
                                SMASSII = body.SMASSII

                            if body.Diameter == -1:
                                Diameter = ''
                            else:
                                Diameter = str(body.Diameter)

                            #compute the Sutter-method guesstimated intercept delta-v relative to the forward and backward bracketing event and take the worst case
                        
                            forward_guesstimated_deltav = body.RelativePositionMagnitude / (event.JulianDate - event_bracket[0] + 1.0) / 86400.0
                            backward_guesstimated_deltav = body.RelativePositionMagnitude / (event_bracket[1] - event.JulianDate + 1.0) / 86400.0

                            bubblefile.write(str(event.JulianDate) + ',' + event.GregorianDate + ',' + body.name + ',' + str(body.SPICEID) + ',' + str(body.SMA) + ',' + Tholen + ',' + SMASSII + ',' + str(body.H) + ',' + Diameter + ',' + str(body.RelativePositionMagnitude) + ',' + str(body.RelativePositionMagnitude / BubbleOptions.LU) + ',' + str(body.RelativeVelocityMagnitude) + ',' + str(max(forward_guesstimated_deltav, backward_guesstimated_deltav)) + '\n')
                    
                        if hasattr(OutputWindow,'WriteText'):
                            OutputWindow.WriteText("Found " + str(len(event_body_list)) + " candidates\n")
                        bubblefile.flush()

                    bubblefile.write('\n')
                    bubblefile.flush()

        if hasattr(OutputWindow,'WriteText'):
            OutputWindow.WriteText("Bubble search outputted to " + bubblefile.name + "\n")

            #do we want to keep going after the end of the mission?
            if (BubbleOptions.CheckForEncountersAfterMissionEnd):
                bubblefile.write('Post-mission encounters\n')
                bubblefile.write('Julian Date, Gregorian Date, Body Name, SPICE ID, Body SMA, Tholen spectral type, SMASSII spectral type, Absolute Magnitude, Diameter, Distance (km), Distance (LU), Relative Velocity (km/s), Guesstimated intercept delta-v (km/s)\n')

                PostMissionState = self.Journeys[-1].missionevents[-1].SpacecraftState
                PostMissionCOE = kepler.cart2kep(PostMissionState[0:3],PostMissionState[3:6],BubbleOptions.mu)
                MissionEndEpoch = self.Journeys[-1].missionevents[-1].JulianDate
                PostMissionPropagationStep = BubbleOptions.PostMissionCheckDuration / BubbleOptions.PostMissionCheckSteps

                for step in range(0, BubbleOptions.PostMissionCheckSteps):
                    CurrentEpoch = MissionEndEpoch + (step + 1) * PostMissionPropagationStep
                    #propagate the spacecraft to this epoch
                    r, v = kepler.kepler(PostMissionCOE[0], PostMissionCOE[1], PostMissionCOE[2], PostMissionCOE[3], PostMissionCOE[4], PostMissionCOE[5], MissionEndEpoch, CurrentEpoch, BubbleOptions.mu)
                    state = [r[0], r[1], r[2], v[0], v[1], v[2]]

                    #check for bodies
                    OutputWindow.WriteText("Searching for targets of opportunity on JD " + str(CurrentEpoch) + "\n")
                    event_body_list = SmallBodyDatabase.find_bubble_targets(state, CurrentEpoch, BubbleOptions.RelativePositionFilterMagnitude, BubbleOptions.RelativeVelocityFilterMagnitude, BubbleOptions.MaximumMagnitude)

                    #write fields
                    for body in event_body_list:
                        if body.Tholen == '-1':
                            Tholen = ''
                        else:
                            Tholen = body.Tholen

                        if body.SMASSII == '-1':
                            SMASSII = ''
                        else:
                            SMASSII = body.SMASSII

                        if body.Diameter == -1:
                            Diameter = ''
                        else:
                            Diameter = str(body.Diameter)

                        #compute the Sutter-method guesstimated intercept delta-v relative to the forward and backward bracketing event and take the worst case
                        
                        forward_guesstimated_deltav = body.RelativePositionMagnitude / (CurrentEpoch - MissionEndEpoch + 1.0) / 86400.0
                        try:
                            myDate = astropy.time.Time(CurrentEpoch, format='jd').datetime.strftime('%m/%d/%Y')
                        except ValueError:
                            myDate = 'idate'
                        bubblefile.write(str(CurrentEpoch) + ',' + myDate + ',' + body.name + ',' + str(body.SPICEID) + ',' + str(body.SMA) +',' + Tholen + ',' + SMASSII + ',' + str(body.H) + ',' + Diameter + ',' + str(body.RelativePositionMagnitude) + ',' + str(body.RelativePositionMagnitude / BubbleOptions.LU) + ',' + str(body.RelativeVelocityMagnitude) + ',' + str(forward_guesstimated_deltav) + '\n')
                    
                    OutputWindow.WriteText("Found " + str(len(event_body_list)) + " candidates\n")
                    bubblefile.flush()

                bubblefile.write('\n')
        OutputWindow.WriteText("Bubble search outputted to " + bubblefile.name + "\n")

        bubblefile.close()

        if BubbleOptions.ifSpice:
            for file in loaded_SPICE_files:
                    spice.unload(file)

    #method to generate a throttle table report
    def CalculateThrottleTableHistory(self):

        #find the throttle table setting for each time step
        requested_throttle_table = []
        mission_throttle_table = []
        time_step_length_array = []
        mass_per_step_array = [] # in kg
        mass_per_step_per_thruster_array = []
        hours_per_step_array = []
        hours_per_step_per_thruster_array = []
        
        for ActiveJourney in self.ChosenJourneys:
            self.ActiveJourney = ActiveJourney
            if self.ActiveJourney <= len(self.Journeys) - 1: #if one journey is selected
                for event in self.Journeys[self.ActiveJourney].missionevents:
                    if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'PSFBthrust', 'LT_spiral'] and event.Number_of_Active_Engines > 0:
                        
                        mission_throttle_table.append(event.ThrottleLevel)

                        time_step_length_array.append(event.TimestepLength * 86400) # in seconds
                        
                        mass_per_step_array.append(event.MassFlowRate * time_step_length_array[-1] * self.Journeys[self.ActiveJourney].thruster_duty_cycle)
                        mass_per_step_per_thruster_array.append(mass_per_step_array[-1] / event.Number_of_Active_Engines)
                        
                        hours_per_step_array.append(event.TimestepLength * 24.0 * self.Journeys[self.ActiveJourney].thruster_duty_cycle)
                        hours_per_step_per_thruster_array.append(hours_per_step_array[-1] / event.Number_of_Active_Engines)
            else:
                for j in range(0, len(self.Journeys)): #for all journeys
                    for event in self.Journeys[j].missionevents:
                        if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSBIthrust', 'end_spiral'] and event.Number_of_Active_Engines > 0:
                            
                            mission_throttle_table.append(event.ThrottleLevel)

                            time_step_length_array.append(event.TimestepLength * 86400) # in seconds
                           
                            mass_per_step_array.append(event.MassFlowRate * time_step_length_array[-1] * self.Journeys[j].thruster_duty_cycle)
                            mass_per_step_per_thruster_array.append(mass_per_step_array[-1] / event.Number_of_Active_Engines)
                            
                            hours_per_step_array.append(event.TimestepLength * 24.0 * self.Journeys[j].thruster_duty_cycle)
                            hours_per_step_per_thruster_array.append(hours_per_step_array[-1] / event.Number_of_Active_Engines)

        #bin the throttle table
        throttle_table_bins = [0.0] * max(mission_throttle_table)
        throttle_table_bins_per_thruster = [0.0] * max(mission_throttle_table)
        throttle_table_hours_bins = [0.0] * max(mission_throttle_table)
        throttle_table_hours_per_thruster = [0.0] * max(mission_throttle_table)
        for i in range(0, len(mission_throttle_table)):
            throttle_table_bins[int(mission_throttle_table[i] - 1)] += mass_per_step_array[i]
            throttle_table_hours_bins[int(mission_throttle_table[i] - 1)] += hours_per_step_array[i]
            throttle_table_bins_per_thruster[int(mission_throttle_table[i] - 1)] += mass_per_step_per_thruster_array[i]
            throttle_table_hours_per_thruster[int(mission_throttle_table[i] - 1)] += hours_per_step_per_thruster_array[i]

        return throttle_table_bins, mission_throttle_table, throttle_table_hours_bins, throttle_table_bins_per_thruster, throttle_table_hours_per_thruster

    def GenerateThrottleReport(self, reportfilename):
        import astropy
        
        #get the throttle table settings history and throttle table binning arrays
        throttle_table_bin_array, throttle_table_history_array, throttle_operation_hours_array, throttle_table_bin_per_thruster, throttle_operations_hours_per_thruster = self.CalculateThrottleTableHistory()

        #write out the throttle report for each time step
        reportfile = open(reportfilename, mode = 'w')
        reportfile.write('Throttle table report for ' + self.mission_name + '\n')
        reportfile.write('Segment start date, Segment width (days), Throttle Setting, Throughtput (kg), Total power in (kW), Number of active thrusters, Total thrust (N), Isp, Total mass flow rate (mg/s)\n')
        counter = 0
        
        for ActiveJourney in self.ChosenJourneys:
            self.ActiveJourney = ActiveJourney
            if self.ActiveJourney <= len(self.Journeys) - 1: #if one journey is selected
                for event in self.Journeys[self.ActiveJourney].missionevents:
                    if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSFBthrust', "PSBIthrust"]:
                        reportfile.write(astropy.time.Time(event.JulianDate - event.TimestepLength / 2.0, format='jd', scale='tdb', out_subfmt='date').utc.iso + ',' + 
                                         str(event.TimestepLength) + ',' +
                                         str(event.ThrottleLevel) + ',' +
                                         str(event.MassFlowRate * event.TimestepLength * 86400.0 * self.Journeys[self.ActiveJourney].thruster_duty_cycle) + ',' + 
                                         str(event.ActivePower) + ',' +
                                         str(event.Number_of_Active_Engines) + ',' +
                                         str(event.AvailableThrust) + ',' +
                                         str(event.Isp) + ',' +
                                         str(event.MassFlowRate * 1.0e+6) + '\n')
                        counter += 1
                    elif event.EventType == 'end_spiral':
                        reportfile.write(astropy.time.Time(event.JulianDate - event.TimestepLength, format='jd', scale='tdb', out_subfmt='date').utc.iso + ',' + 
                                         str(event.TimestepLength) + ',' +
                                         str(event.ThrottleLevel) + ',' +
                                         str(event.MassFlowRate * event.TimestepLength * 86400.0 * self.Journeys[self.ActiveJourney].thruster_duty_cycle) + ',' + 
                                         str(event.ActivePower) + ',' +
                                         str(event.Number_of_Active_Engines) + ',' +
                                         str(event.AvailableThrust) + ',' +
                                         str(event.Isp) + ',' +
                                         str(event.MassFlowRate * 1.0e+6) + '\n')
                        counter += 1

            else:
                for j in range(0, len(self.Journeys)): #for all journeys
                    for event in self.Journeys[j].missionevents:
                        if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', 'PSFBthrust', "PSBIthrust"] and event.Number_of_Active_Engines > 0:
                            reportfile.write(astropy.time.Time(event.JulianDate - event.TimestepLength / 2.0, format='jd', scale='tdb', out_subfmt='date').utc.iso + ',' + 
                                         str(event.TimestepLength) + ',' +
                                         str(event.ThrottleLevel) + ',' +
                                         str(event.MassFlowRate * event.TimestepLength * 86400.0 * self.Journeys[j].thruster_duty_cycle) + ',' + 
                                         str(event.ActivePower) + ',' +
                                         str(event.Number_of_Active_Engines) + ',' +
                                         str(event.AvailableThrust) + ',' +
                                         str(event.Isp) + ',' +
                                         str(event.MassFlowRate * 1.0e+6) + '\n')
                            counter += 1
                        elif event.EventType == 'end_spiral':
                            reportfile.write(astropy.time.Time(event.JulianDate - event.TimestepLength, format='jd', scale='tdb', out_subfmt='date').utc.iso + ',' + 
                                         str(event.TimestepLength) + ',' +
                                         str(event.ThrottleLevel) + ',' +
                                         str(event.MassFlowRate * event.TimestepLength * 86400.0 * self.Journeys[j].thruster_duty_cycle) + ',' + 
                                         str(event.ActivePower) + ',' +
                                         str(event.Number_of_Active_Engines) + ',' +
                                         str(event.AvailableThrust) + ',' +
                                         str(event.Isp) + ',' +
                                         str(event.MassFlowRate * 1.0e+6) + '\n')
                            counter += 1

        reportfile.write('\n')
        reportfile.write('\n')

        #write out a throttle binning report
        reportfile.write('Throttle table binning report for ' + self.mission_name + '\n')
        reportfile.write('Throttle setting, Throughput (kg), operating hours, Throughput per thruster (kg), operating hours per thruster\n')
        for i in range(0, len(throttle_table_bin_array)):
            reportfile.write(str(i+1) + ',' + str(throttle_table_bin_array[i]) + ',' + str(throttle_operation_hours_array[i]) + ','\
               + str(throttle_table_bin_per_thruster[i]) + ',' + str(throttle_operations_hours_per_thruster[i]) + '\n')

        reportfile.write('\n\n')
        reportfile.write('Worst case for a single thruster\n')
        reportfile.write('computed by always being on in multi-thruster operations and handling all of the single-thrust operations)\n')
        single_thruster_worst_case_throughput = 0.0
        single_thruster_worst_case_hours = 0.0
        for i in range(0, len(throttle_table_bin_array)):
            single_thruster_worst_case_throughput += throttle_table_bin_per_thruster[i]
            single_thruster_worst_case_hours += throttle_operations_hours_per_thruster[i]
        reportfile.write('throughput (kg): ' + str(single_thruster_worst_case_throughput) + '\n')
        reportfile.write('operating hours: ' + str(single_thruster_worst_case_hours) + '\n')

        reportfile.close()
        
    def GenerateThrottlePlot(self, throttletablefile, reportfilename):
        import matplotlib.pyplot as plt
        
        trajectory_mdot = []
        trajectory_thrust = []
        trajectory_power = []
        
        counter = 0
        for ActiveJourney in self.ChosenJourneys:
            self.ActiveJourney = ActiveJourney
            if self.ActiveJourney <= len(self.Journeys) - 1: #if one journey is selected
                for event in self.Journeys[self.ActiveJourney].missionevents:
                    if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', "PSBIthrust"]:
                        trajectory_mdot.append(event.MassFlowRate * 1.0e+6 / event.Number_of_Active_Engines / event.DVmagorThrottle)
                        trajectory_thrust.append(event.AvailableThrust * 1.0e+3 / event.Number_of_Active_Engines)
                        trajectory_power.append(event.ActivePower / event.Number_of_Active_Engines)
                        counter += 1
                    elif event.EventType == 'end_spiral':
                        trajectory_mdot.append(event.MassFlowRate * 1.0e+6 / event.Number_of_Active_Engines / event.DVmagorThrottle)
                        trajectory_thrust.append(event.AvailableThrust * 1.0e+3 / event.Number_of_Active_Engines)
                        trajectory_power.append(event.ActivePower / event.Number_of_Active_Engines)
                        counter += 1
            else:
                for j in range(0, len(self.Journeys)): #for all journeys
                    for event in self.Journeys[j].missionevents:
                        if event.EventType in ['SFthrust', 'SSFthrust', 'FBLTSthrust', 'FBLTthrust', "PSBIthrust"] and event.Number_of_Active_Engines > 0:
                            trajectory_mdot.append(event.MassFlowRate * 1.0e+6 / event.Number_of_Active_Engines / event.DVmagorThrottle)
                            trajectory_thrust.append(event.AvailableThrust * 1.0e+3 / event.Number_of_Active_Engines)
                            trajectory_power.append(event.ActivePower / event.Number_of_Active_Engines)
                            counter += 1
                        elif event.EventType == 'end_spiral':
                            trajectory_mdot.append(event.MassFlowRate * 1.0e+6 / event.Number_of_Active_Engines / event.DVmagorThrottle)
                            trajectory_thrust.append(event.AvailableThrust * 1.0e+3 / event.Number_of_Active_Engines)
                            trajectory_power.append(event.ActivePower / event.Number_of_Active_Engines)
                            counter += 1
                 

        thruster_throttle_table = ThrottleTable.ThrottleTable(throttletablefile)
                                
        thruster_throttle_table.plot_throttle_table(False, [], [], False,trajectory_mdot,trajectory_thrust,trajectory_power)

    def GenerateThrottleHistogram(self, PlotOptions):
        
        import matplotlib.pyplot
        #get the throttle table settings history and throttle table binning arrays
        throttle_table_bin_array, throttle_table_history_array, throttle_operation_hours_array, throttle_table_bin_per_thruster, throttle_operations_hours_per_thruster = self.CalculateThrottleTableHistory()

        #generate a histogram
        HistogramFigure = matplotlib.pyplot.figure()
        HistogramAxes = HistogramFigure.add_axes([0.1, 0.1, 0.8, 0.8])
        matplotlib.rcParams.update({'font.size': PlotOptions.FontSize})

        bins = range(1, len(throttle_table_bin_array)+1)
        HistogramAxes.bar(bins, throttle_table_bin_array)
        HistogramAxes.set_xlabel('Throttle level')
        HistogramAxes.set_ylabel('Throughput (kg)')
        HistogramAxes.grid(b=True)
        HistogramAxes.set_xlim(1, len(throttle_table_bin_array) + 1)
        HistogramAxes.set_xticks(bins)

        # Hide major tick labels
        HistogramAxes.xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())

        # Customize minor tick labels
        minorlabels = []
        minorlocators = []
        for b in bins:
            minorlabels.append(str(b))
            minorlocators.append(b + 0.5)
        HistogramAxes.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(minorlocators))
        HistogramAxes.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(minorlabels))

        HistogramFigure.show()

    def AutoTableMission(self, TableFileName, PlotOptions):
        #reset the event label counter
        PlotOptions.EventCounter = 1
        skipNext = False

        #do we need the delta-v and flyby altitude columns?
        for journey in self.Journeys:
            for event in journey.missionevents:
                if event.EventType in ['chem_burn','pwr_flyby','insertion','rendezvous','TCM']:
                    PlotOptions.AutoTableDeltavColumn = True
                elif event.EventType == 'departure' and event.DVmagorThrottle > 1.0e-6:
                    PlotOptions.AutoTableDeltavColumn = True
                if event.EventType in ['pwr_flyby','upwr_flyby','periapse']:
                    PlotOptions.AutoTableAltitudeColumn = True

        TableFile = open(TableFileName, mode = 'w')

        TableFile.write('%auto-generated table from PyEMTG\n')
        TableFile.write('\\documentclass{article}\n')
        TableFile.write('\\usepackage[margin=1in, landscape]{geometry}\n')
        TableFile.write('\\usepackage{amsmath}\n')
        TableFile.write('\\usepackage{float}\n')
        TableFile.write('\n')
        TableFile.write('\\begin{document}\n')
        TableFile.write('\n')
        TableFile.write('\\begin{table}\n')
        TableFile.write('   \centering\n')
        TableFile.write('   \caption{Itinerary for ' + self.mission_name.replace('_','-').strip(' ') + '}\n')
        TableFile.write('   \\begin{tabular}{lllllll')
        if PlotOptions.AutoTableAltitudeColumn:
            TableFile.write('l')
        if PlotOptions.AutoTableDeltavColumn:
            TableFile.write('ll')
        if PlotOptions.DisplayEventMass:
            TableFile.write('l')
        if PlotOptions.DisplayArrivalPhaseAngle:
            TableFile.write('l')
        TableFile.write('}\n')
        TableFile.write('       \hline\hline\n')
        TableFile.write('       Event \# & Date & Event & Location & $v_\\infty \\left(km/s\\right)$')
        if PlotOptions.AutoTableDeltavColumn:
           TableFile.write(' & $\Delta v \\left(km/s\\right)$')
           TableFile.write(' & $I_{sp}$')
        if PlotOptions.AutoTableAltitudeColumn:
           TableFile.write(' & Flyby altitude (km)')
        if PlotOptions.DisplayEventMass:
            TableFile.write(' & Mass (kg)')
        if PlotOptions.DisplayArrivalPhaseAngle:
            TableFile.write(' & ' + r'$\beta$' + '(' + r'$^\circ$' + ')')
        TableFile.write('\\\\\n')
        TableFile.write('       \hline\n')

        for journey in self.Journeys:
            skipNext = journey.AutoTableJourney(TableFile, PlotOptions, skipNext)

        TableFile.write('       \hline\hline\n')
        TableFile.write('   \end{tabular}\n')
        TableFile.write('   \label{tab:' + self.mission_name.replace('_','-').strip(' ') + '}\n')
        TableFile.write('\end{table}\n')
        TableFile.write('\n')
        TableFile.write('Total flight time: ' + str(self.total_flight_time_years) + ' years\n')
        TableFile.write('\n')
        TableFile.write('Total deterministic $\Delta v$: ' + str(self.total_deterministic_deltav) + ' km/s \n')
        TableFile.write('\n')
        if (self.total_statistical_deltav > 1.0e-8):
            TableFile.write('Total statistical $\Delta v$: ' + str(self.total_statistical_deltav) + ' km/s \n')
        TableFile.write('\n')
        TableFile.write('\end{document}\n')

        TableFile.close()

    #************************************************************************************getJourneyIndex()
    def getJourneyIndex(self, journeyNameString):
        for journeyIndex in range(0, len(self.Journeys)):
            if self.Journeys[journeyIndex].journey_name == journeyNameString:
                return journeyIndex
        
        #if you get this far, something went wrong                                                                 
        raise Exception("Journey '" + journeyNameString + "' not found.")  
        
    def Comparatron(self,baseline_path,csv_file_name=None,full_output=False,tolerance_dict={},default_tolerance=1e-15):
        '''
        NOTE FOR JACOB:
        The lines commented out near the end of the script (starting at line 1250)
        I left in because they may be useful if you decide you do not want the 
        entire dataframe showing up in the failFile. If you want to change this, 
        remove the "output.to_csv(outputdir + '/failed_tests.csv', mode='a', index=False)" 
        line from testatron.py (line ___) and go back to writing the separate 
        output file in this script. Whatever aspects of comparison you want to 
        write to the failFile can be returned to the driver as comp_return in this script.
        '''
        #Compare just constraint and decision vectors
        import pandas as pd
        #baseline_path is the file location of the emtg case to compare against. This is the only required argument.
        #csv_file_name can specify the output file name. A csv file is only output if there are mismatches between cases. By default, it uses the mission name as the name for the csv output file
        #full_output determines whether all cases that disagree between self and the baseline will be saved to csv. If True, all cases that are not in exact agreement will be saved. If False, only cases that do not meet the tolerance will be saved to the csv output.
        #tolerance_dict is a dictionary of tolerences that, if provided, will override the default value. Provide a dictionary with the attribute name as the key and the alternate tolerance as the value (i.e. {'total_statistical_deltav':1e-6,'Declination':1e-8})
        #default_tolerance is the tolerance value used for all attributes that do not appear in the tolerance_dict
        
        #Create mission object for baseline case to compare against, create empty dataframe to store comparison
        baseline = Mission(baseline_path)
        comparison = pd.DataFrame(columns=['Output Name','Baseline Value','New Value','Error','Match'])
        #This comparison df is what is written to the csv output file. The Output Name corresponds to the attribute name within the Mission, Journey, or Mission Event objects. Error is the absolute error 
        
        #can I test the current mission? did it actually load anything? if not, return
        if not hasattr(self, 'Journeys'):
            return False, comparison

        if csv_file_name==None:
            csv_file_name = self.mission_name + '_comparison.csv'
        
        #Make sure the number of journeys and mission events are the same between both cases
        if len(baseline.Journeys) != len(self.Journeys):
            pd.DataFrame(columns=['Journey Mismatch!']).to_csv(csv_file_name,index=False)
            return False, comparison
        for i in range(len(baseline.Journeys)): 
            if len(baseline.Journeys[i].missionevents) != len(self.Journeys[i].missionevents):
                pd.DataFrame(columns=['Journey ' + str(i) + ' Mission Events Mismatch']).to_csv(csv_file_name,index=False)
                return False, comparison
        
        #Check that all attributes agree between self and the baseline case. Any disagreements are added to the end of the csv output file
        all_baseline_attributes = ['Mission.'+a for a in dir(baseline)] + ['Journey.'+a for a in dir(baseline.Journeys[0])] + ['MissionEvent.'+a for a in dir(baseline.Journeys[0].missionevents[0])]
        all_new_attributes = ['Mission.'+a for a in dir(self)] + ['Journey.'+a for a in dir(self.Journeys[0])] + ['MissionEvent.'+a for a in dir(self.Journeys[0].missionevents[0])]
        attributes_only_in_baseline = list(set(all_baseline_attributes).difference(all_new_attributes))
        baseline_mission_attr_to_ignore = [a.strip('Mission.') for a in attributes_only_in_baseline if a[7]=='.']
        baseline_journey_attr_to_ignore = [a.strip('Journey.') for a in attributes_only_in_baseline if a[0]=='J']
        baseline_mevent_attr_to_ignore = [a.strip('MissionEvent.') for a in attributes_only_in_baseline if a[7]=='E']
        attributes_only_in_new = list(set(all_new_attributes).difference(all_baseline_attributes))
        new_mission_attr_to_ignore = [a.strip('Mission.') for a in attributes_only_in_new if a[7]=='.']
        new_journey_attr_to_ignore = [a.strip('Journey.') for a in attributes_only_in_new if a[0]=='J']
        new_mevent_attr_to_ignore = [a.strip('MissionEvent.') for a in attributes_only_in_new if a[7]=='E']
        attr_check = len(attributes_only_in_baseline) + len(attributes_only_in_new) #If this is greater than zero than attributes are not in complete agreement between baseline and new        
        
        #The following function takes the dataframes provided by subtracting the self and baseline dataframes and removes all nan/blank values
        def getDiff(df,correct_difference):
            finaldf = df.loc[(df != correct_difference).any(axis=1)]
            finaldf = finaldf.transpose()
            finaldf = finaldf.loc[(finaldf != correct_difference).any(axis=1)]
            return finaldf
        
        #The following function unpacks lists of lists
        def unpackLists(df):
            finaldf=df.copy()
            for col in finaldf.columns:
                if isinstance(finaldf[col][0], list): #checks if the seris is populated by lists
                    for i in range(0, finaldf[col].count()):
                        finaldf[col+str(i)]=pd.Series(finaldf[col][i])
                    finaldf=finaldf.drop([col],axis=1)
            return finaldf
        
        #Store all relevant mission values and find the values that are not in agreement
        baseline_mission_attributes = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in vars(baseline).items() if k not in ['user_data','Journeys'] + baseline_mission_attr_to_ignore]))
        new_mission_attributes = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in vars(self).items() if k not in ['user_data','Journeys'] + new_mission_attr_to_ignore]))
        
        temp_num = new_mission_attributes.select_dtypes(include=[np.number]).subtract(baseline_mission_attributes.select_dtypes(include=[np.number])).fillna(0)
        temp_str = new_mission_attributes.select_dtypes(include=['object']).fillna('') == baseline_mission_attributes.select_dtypes(include=['object']).fillna('')
        mission_attr_diff = getDiff(temp_num,0)
        mission_str_diff = getDiff(temp_str,True)
        
        #Store all disagreements into comparison dataframe
        for ix in mission_attr_diff.index:
            for col in mission_attr_diff.columns:
                if not np.isnan(mission_attr_diff.at[ix,col]):
                    if ix not in tolerance_dict.keys():
                        tolerance = default_tolerance
                    else:
                        tolerance = tolerance_dict[ix]
                    if abs(mission_attr_diff.at[ix,col]) <= tolerance:
                        match = True
                    else:
                        match = False
                    comparison = comparison.append({'Output Name':ix+'['+str(col)+']','Baseline Value':baseline_mission_attributes.at[col,ix],'New Value':new_mission_attributes.at[col,ix],'Error':mission_attr_diff.at[ix,col],'Tolerance':tolerance,'Match':match},ignore_index=True)
        for ix in mission_str_diff.index:
            for col in mission_str_diff.columns:
                if mission_str_diff.at[ix,col] != '':
                    comparison = comparison.append({'Output Name':ix+'['+str(col)+']','Baseline Value':baseline_mission_attributes.at[col,ix],'New Value':new_mission_attributes.at[col,ix],'Match':mission_str_diff.at[ix,col]},ignore_index=True)
        
        #Repeat for Journey objects
        for i in range(len(baseline.Journeys)):
            baseline_journey_attributes = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in vars(baseline.Journeys[i]).items() if k not in ['parent','missionevents']+baseline_journey_attr_to_ignore]))
            new_journey_attributes = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in vars(self.Journeys[i]).items() if k not in ['parent','missionevents']+new_journey_attr_to_ignore]))
            
            temp_num = new_journey_attributes.select_dtypes(include=[np.number]).subtract(baseline_journey_attributes.select_dtypes(include=[np.number])).fillna(0)
            temp_str = new_journey_attributes.select_dtypes(include=['object']).fillna('') == baseline_journey_attributes.select_dtypes(include=['object']).fillna('')
            journey_attr_diff = getDiff(temp_num,0)
            journey_str_diff = getDiff(temp_str,True)
            for ix in journey_attr_diff.index:
                for col in journey_attr_diff.columns:
                    if not np.isnan(journey_attr_diff.at[ix,col]):
                        if ix not in tolerance_dict.keys():
                            tolerance = default_tolerance
                        else:
                            tolerance = tolerance_dict[ix]
                        if abs(journey_attr_diff.at[ix,col]) <= tolerance:
                            match = True
                        else:
                            match = False
                        comparison = comparison.append({'Output Name':'Journey '+str(i)+' '+ix+'['+str(col)+']','Baseline Value':baseline_journey_attributes.at[col,ix],'New Value':new_journey_attributes.at[col,ix],'Error':journey_attr_diff.at[ix,col],'Tolerance':tolerance,'Match':match},ignore_index=True)
            for ix in journey_str_diff.index:
                for col in journey_str_diff.columns:
                    if journey_str_diff.at[ix,col] != '':
                        comparison = comparison.append({'Output Name':'Journey '+str(i)+' '+ix+'['+str(col)+']','Baseline Value':baseline_journey_attributes.at[col,ix],'New Value':new_journey_attributes.at[col,ix],'Match':journey_str_diff.at[ix,col]},ignore_index=True)
            
            #Repeat for Mission Events
            for j in range(len(baseline.Journeys[i].missionevents)):
                baseline_missionevent_attributes = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in vars(baseline.Journeys[i].missionevents[j]).items() if k not in ['parent']+baseline_mevent_attr_to_ignore]))
                baseline_missionevent_attributes['ThrottleLevel'] = baseline_missionevent_attributes['ThrottleLevel'].astype(str).replace('0.0','-')
                new_missionevent_attributes = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in vars(self.Journeys[i].missionevents[j]).items() if k not in ['parent']+new_mevent_attr_to_ignore]))
                new_missionevent_attributes['ThrottleLevel'] = new_missionevent_attributes['ThrottleLevel'].astype(str).replace('0.0','-')                
                
                temp_num = new_missionevent_attributes.select_dtypes(include=[np.number]).subtract(baseline_missionevent_attributes.select_dtypes(include=[np.number])).fillna(0)
                temp_str = new_missionevent_attributes.select_dtypes(include=['object']).fillna('')==baseline_missionevent_attributes.select_dtypes(include=['object']).fillna('')
                mevent_attr_diff = getDiff(temp_num,0)
                mevent_str_diff = getDiff(temp_str,True)
                for ix in mevent_attr_diff.index:
                    for col in mevent_attr_diff.columns:
                        if not np.isnan(mevent_attr_diff.at[ix,col]):
                            if ix not in tolerance_dict.keys():
                                tolerance = default_tolerance
                            else:
                                tolerance = tolerance_dict[ix]
                            if abs(mevent_attr_diff.at[ix,col]) <= tolerance:
                                match = True
                            else:
                                match = False
                            comparison = comparison.append({'Output Name':'Journey '+str(i)+' Mission Event '+str(j)+' '+ix+'['+str(col)+']','Baseline Value':baseline_missionevent_attributes.at[col,ix],'New Value':new_missionevent_attributes.at[col,ix],'Error':mevent_attr_diff.at[ix,col],'Tolerance':tolerance,'Match':match},ignore_index=True)
                for ix in mevent_str_diff.index:
                    for col in mevent_str_diff.columns:
                        if mevent_str_diff.at[ix,col] != '':
                            comparison = comparison.append({'Output Name':'Journey '+str(i)+' Mission Event '+str(j)+' '+ix+'['+str(col)+']','Baseline Value':baseline_missionevent_attributes.at[col,ix],'New Value':new_missionevent_attributes.at[col,ix],'Match':mevent_str_diff.at[ix,col]},ignore_index=True)
                
        #If there are any attributes in the attributes_only_in_xxx lists then append them to comparison df
        if len(attributes_only_in_baseline) != 0:
            comparison = comparison.append({'Output Name':'Attributes only in Baseline','Baseline Value':attributes_only_in_baseline,'Match':False},ignore_index=True)
        if len(attributes_only_in_new) != 0:
            comparison = comparison.append({'Output Name':'Attributes only in New','New Value':attributes_only_in_new,'Match':False},ignore_index=True)
        
        #pdb.set_trace()
        
        #Return a value of True if the missions are in complete agreement
        if len(comparison.loc[comparison['Match']==False]) == 0 and attr_check == 0:
            return True, comparison;
        else:
#            comp_return=list(comparison.loc[comparison['Match']==False]['Output Name'])
            if full_output: #if full_output=True, save all values present in comparison to a csv
                comparison = comparison[['Output Name','Baseline Value','New Value','Error','Tolerance','Match']]
                comparison.to_csv(csv_file_name,index=False)
                return False, comparison;
            else: #if full_output=False, save only values that are not within the tolerance to a csv
                comparison = comparison[['Output Name','Baseline Value','New Value','Error','Tolerance','Match']]
                comparison.loc[comparison['Match']==False].to_csv(csv_file_name,index=False)
                return False, comparison;