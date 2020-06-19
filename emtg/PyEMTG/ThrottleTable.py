#class for reading and parsing electric propulsion throttle tables

import os
import math
import operator
import numpy
import copy

class ThrottleSetting(object):
    #constructor for parsing a string from a throttle file
    def __init__(self):
        self.clear()

    def initialize_from_string(self, throttlestring):
        self.clear()
        self.parse_throttle_line(throttlestring)

    #constructor when the propulsion information is known
    def initialize_from_input_data(self, Mdot, Thrust, Isp, efficiency, ThrusterPower, BeamCurrent=None, BeamVoltage=None):
        self.TL = 'TL0'
        self.Mdot = Mdot #mg/s
        self.BeamCurrent = BeamCurrent #A
        self.BeamVoltage = BeamVoltage #V
        self.Thrust = Thrust #mN
        self.Isp = Isp # s
        self.efficiency = efficiency
        self.ThrusterPower = ThrusterPower #kW

    def clear(self):
        self.TL = [] 
        self.Mdot = [] #mg/s
        self.BeamCurrent = [] #A
        self.BeamVoltage = [] #V
        self.Thrust = [] #mN
        self.Isp = [] # s
        self.efficiency = []
        self.ThrusterPower = [] #kW

    #parse an input line
    def parse_throttle_line(self, throttlestring):
        throttlecell = throttlestring.split(',')
        self.TL = throttlecell[0]
        self.Mdot = float(throttlecell[1]) #mg/s
        self.BeamCurrent = float(throttlecell[2]) #A
        self.BeamVoltage = float(throttlecell[3]) #V
        self.Thrust = float(throttlecell[4]) #mN
        self.Isp = float(throttlecell[5]) # s
        self.efficiency = float(throttlecell[6])
        self.ThrusterPower = float(throttlecell[7].strip('\n')) #kW

    #function to compare this throttle setting to another
    def compare(self, otherThrottleSetting):
        diff_ThrusterPower = self.ThrusterPower - otherThrottleSetting.ThrusterPower
        diff_Thrust = self.Thrust - otherThrottleSetting.Thrust
        diff_Mdot = self.Mdot - otherThrottleSetting.Mdot
        diff_Isp = self.Isp - otherThrottleSetting.Isp
        diff_efficiency = self.efficiency - otherThrottleSetting.Isp

        if hasattr(self, 'Voltage'):
            diff_Voltage = self.Voltage - otherThrottleSetting.Voltage
        else:
            diff_Voltage = 0
        return diff_ThrusterPower, diff_Thrust, diff_Mdot, diff_Isp, diff_efficiency, diff_Voltage

    #function to find the four nearest neighbors in P-V space to a given throttle point
    def find_neighbors(self, ThrottleTable):
        #the nearest neighbors are the non-dominated set in +P, -P, +V, -V space
        neighbors = []
        for Setting in ThrottleTable.ThrottleSettings:
            if not Setting.TL == self.TL:
                isDominated = False
                for OtherSetting in ThrottleTable.ThrottleSettings:
                
                    if not Setting.TL == OtherSetting.TL:
                        #check by same voltage and different power, then check those that are left
                        if Setting.BeamVoltage == OtherSetting.BeamVoltage == self.BeamVoltage:
                            if Setting.ThrusterPower > self.ThrusterPower and OtherSetting.ThrusterPower > self.ThrusterPower and OtherSetting.ThrusterPower < Setting.ThrusterPower:
                                print(Setting.TL, ' is dominated because it has the same BeamVoltage as ', self.TL, ' but ', OtherSetting.TL, ' is closer to ', self.TL, ' in ThrusterPower')
                                isDominated = True
                                break
                            elif Setting.ThrusterPower < self.ThrusterPower and OtherSetting.ThrusterPower < self.ThrusterPower and OtherSetting.ThrusterPower > Setting.ThrusterPower:
                                print(Setting.TL, ' is dominated because it has the same BeamVoltage as ', self.TL, ' but ', OtherSetting.TL, ' is closer to ', self.TL, ' in ThrusterPower')
                                isDominated = True
                                break
                        #now the ones that don't have the same voltage
                        
                        else:
                            if Setting.BeamVoltage > self.BeamVoltage and OtherSetting.BeamVoltage > self.BeamVoltage:
                                if Setting.ThrusterPower > self.ThrusterPower and OtherSetting.ThrusterPower > self.ThrusterPower and OtherSetting.ThrusterPower < Setting.ThrusterPower:
                                    print(Setting.TL, ' is dominated because it has a greater BeamVoltage than ', self.TL, ' but ', OtherSetting.TL, ' is closer to ', self.TL, ' in ThrusterPower')
                                    #print Setting.TL, Setting.BeamVoltage, Setting.ThrusterPower, OtherSetting.TL, OtherSetting.BeamVoltage, OtherSetting.ThrusterPower
                                    isDominated = True
                                    break
                                elif Setting.ThrusterPower < self.ThrusterPower and OtherSetting.ThrusterPower < self.ThrusterPower and OtherSetting.ThrusterPower > Setting.ThrusterPower:
                                    print(Setting.TL, ' is dominated because it has a greater BeamVoltage than ', self.TL, ' but ', OtherSetting.TL, ' is closer to ', self.TL, ' in ThrusterPower')
                                    isDominated = True
                                    break
                            elif Setting.BeamVoltage < self.BeamVoltage and OtherSetting.BeamVoltage < self.BeamVoltage:
                                if Setting.ThrusterPower > self.ThrusterPower and OtherSetting.ThrusterPower > self.ThrusterPower and OtherSetting.ThrusterPower < Setting.ThrusterPower:
                                    print(Setting.TL, ' is dominated because it has a lesser BeamVoltage than ', self.TL, ' but ', OtherSetting.TL, ' is closer to ', self.TL, ' in ThrusterPower')
                                    #print Setting.TL, Setting.BeamVoltage, Setting.ThrusterPower, OtherSetting.TL, OtherSetting.BeamVoltage, OtherSetting.ThrusterPower
                                    isDominated = True
                                    break
                                elif Setting.ThrusterPower < self.ThrusterPower and OtherSetting.ThrusterPower < self.ThrusterPower and OtherSetting.ThrusterPower > Setting.ThrusterPower:
                                    print(Setting.TL, ' is dominated because it has a lesser BeamVoltage than ', self.TL, ' but ', OtherSetting.TL, ' is closer to ', self.TL, ' in ThrusterPower')
                                    isDominated = True
                                    break
                        
                print(isDominated)
                if not isDominated:
                    neighbors.append(Setting)

        return neighbors

    #def function to print the throttle setting
    def print_throttle_setting(self, outputfile):
        outputfile.write(str(self.TL) + ',' + str(self.Mdot) + ',' + str(self.BeamCurrent) + ',' + str(self.BeamVoltage) + ',' + str(self.Thrust) + ',' + str(self.Isp) + ',' + str(self.efficiency) + ',' + str(self.ThrusterPower) + '\n')

class ThrottleSet(object):
    def __init__(self, Mdot):
        self.clear()
        self.Mdot = Mdot

    def clear(self):
        self.ThrottleSettings = []

    def add_setting(self, setting):
        self.ThrottleSettings.append(setting)

    def sort(self):
        self.ThrottleSettings = sorted(self.ThrottleSettings, key=lambda setting: setting.ThrusterPower)
        

class ThrottleTable(object):
    def __init__(self, filename="none"):
        self.clear()
        if filename != "none":
            self.parse_throttle_file(filename)
            self.create_constant_Mdot_sets()

    def clear(self):
        self.ThrottleSettings = []
        self.PPUefficiency = 1.0
        self.PPUminpower = 0.5
        self.PPUmaxpower = 10.0
        self.min_Mdot_set = []
        self.min_Mdot_polynomial = []
        self.max_Mdot_set = []
        self.max_Mdot_polynomial = []
        self.high_thrust_Thrust_coeff = numpy.zeros(5)
        self.high_thrust_Mdot_coeff = numpy.zeros(5)
        self.high_Isp_Thrust_coeff = numpy.zeros(5)
        self.high_Isp_Mdot_coeff = numpy.zeros(5)
        self.poly2D_coefficients = numpy.zeros(25).reshape(5, 5)
        self.FileName = ''
        self.max_mdot = 0.0
        self.max_voltage = 0.0
        self.max_thrust = 0.0

    #function to parse a throttle file
    def parse_throttle_file(self, filename):
        self.clear()

        #read the throttle file
        if os.path.isfile(filename):
            inputfile = open(filename, "r")
            self.FileName = filename
        else:
            print("File ", inputfile, " does not exist!")
            return

        for line in inputfile:
            if not line[0] == '#':
                if line[0:3] == 'PPU':
                    linecell = line.split(',')
                    if linecell[0] == 'PPU efficiency':
                        self.PPUefficiency = float(linecell[1])
                    elif linecell[0] == 'PPU min power (kW)':
                        self.PPUminpower = float(linecell[1])
                    elif linecell[0] == 'PPU max power (kW)':
                        self.PPUmaxpower = float(linecell[1])

                elif line[0:2] == 'TL' or line[0:3] == 'ETL':
                    NewThrottleSetting = ThrottleSetting()
                    NewThrottleSetting.initialize_from_string(line)
                    self.ThrottleSettings.append(NewThrottleSetting)
                    
                    if self.ThrottleSettings[-1].Mdot > self.max_mdot:
                        self.max_mdot = self.ThrottleSettings[-1].Mdot
                    
                    if self.ThrottleSettings[-1].BeamVoltage > self.max_voltage:
                        self.max_voltage = self.ThrottleSettings[-1].BeamVoltage
                        
                    
                    if self.ThrottleSettings[-1].Thrust > self.max_thrust:
                        self.max_thrust = self.ThrottleSettings[-1].Thrust
                        
                else:   
                    linecell = line.split(',')
                    if (linecell[0] == "high_thrust_Thrust"):
                        for k in range(0, 5):
                            self.high_thrust_Thrust_coeff[k] = float(linecell[1 + k])
                    elif (linecell[0] == "high_thrust_Mdot"):
                        for k in range(0, 5):
                            self.high_thrust_Mdot_coeff[k] = float(linecell[1 + k])
                    elif (linecell[0] == "high_Isp_Thrust"):
                        for k in range(0, 5):
                            self.high_Isp_Thrust_coeff[k] = float(linecell[1 + k])
                    elif (linecell[0] == "high_Isp_Mdot"):
                        for k in range(0, 5):
                            self.high_Isp_Mdot_coeff[k] = float(linecell[1 + k])
                    elif (linecell[0] == "2Dpolyrow1"):
                        for k in range(0, 5):
                            self.poly2D_coefficients[0][k] = float(linecell[1 + k])
                    elif (linecell[0] == "2Dpolyrow2"):
                        for k in range(0, 5):
                            self.poly2D_coefficients[1][k] = float(linecell[1 + k])
                    elif (linecell[0] == "2Dpolyrow3"):
                        for k in range(0, 5):
                            self.poly2D_coefficients[2][k] = float(linecell[1 + k])
                    elif (linecell[0] == "2Dpolyrow4"):
                        for k in range(0, 5):
                            self.poly2D_coefficients[3][k] = float(linecell[1 + k])
                    elif (linecell[0] == "2Dpolyrow5"):
                        for k in range(0, 5):
                            self.poly2D_coefficients[4][k] = float(linecell[1 + k])


    def get_voltage(self,throttle_level):
        return self.ThrottleSettings[throttle_level].BeamVoltage
        
    def get_mdot(self,throttle_level):
        return self.ThrottleSettings[throttle_level].Mdot

    def find_nearest_throttle_setting_2D(self,mdot,thrust):
        import math
        least_dist = math.sqrt((mdot-self.ThrottleSettings[0].Mdot)*(mdot-self.ThrottleSettings[0].Mdot)/self.max_mdot/self.max_mdot+(thrust-self.ThrottleSettings[0].Thrust)*(thrust-self.ThrottleSettings[0].Thrust)/self.max_thrust/self.max_thrust)
        least_dist_idx = 0
        
        for idx in range(1,len(self.ThrottleSettings)):
            new_dist = math.sqrt((mdot-self.ThrottleSettings[idx].Mdot)*(mdot-self.ThrottleSettings[idx].Mdot)/self.max_mdot/self.max_mdot+(thrust-self.ThrottleSettings[idx].Thrust)*(thrust-self.ThrottleSettings[idx].Thrust)/self.max_thrust/self.max_thrust)
            if new_dist < least_dist:
                least_dist_idx = idx
                least_dist = new_dist
        
        return least_dist_idx

    def CalculateThrusterPerformance1D_Polynomial(self, inputPower, PreferredThrottleSet='HighThrust'):
        P2 = inputPower * inputPower;
        P3 = P2 * inputPower;
        P4 = P3 * inputPower;

        if PreferredThrottleSet == 'HighThrust':
            self.polyThrust_1D = (self.high_thrust_Thrust_coeff[0]
                + self.high_thrust_Thrust_coeff[1] * inputPower
                + self.high_thrust_Thrust_coeff[2] * P2
                + self.high_thrust_Thrust_coeff[3] * P3
                + self.high_thrust_Thrust_coeff[4] * P4) * 1.0e+3

            self.poly_dThrust_dP_1D = (self.high_thrust_Thrust_coeff[1]
                + 2.0 * self.high_thrust_Thrust_coeff[2] * inputPower 
                + 3.0 * self.high_thrust_Thrust_coeff[3] * P2 
                + 4.0 * self.high_thrust_Thrust_coeff[4] * P3 ) * 1.0e+3

            self.polyMdot_1D = (self.high_thrust_Mdot_coeff[0]
                + self.high_thrust_Mdot_coeff[1] * inputPower
                + self.high_thrust_Mdot_coeff[2] * P2
                + self.high_thrust_Mdot_coeff[3] * P3
                + self.high_thrust_Mdot_coeff[4] * P4) * 1.0e+6

            self.poly_dMdotdP_1D = (self.high_thrust_Mdot_coeff[1]
                + 2.0 * self.high_thrust_Mdot_coeff[2] * inputPower 
                + 3.0 * self.high_thrust_Mdot_coeff[3] * P2 
                + 4.0 * self.high_thrust_Mdot_coeff[4] * P3 ) * 1.0e+6

            self.d2HighThrust_1D_Mdot_dP2 = (self.high_thrust_Mdot_coeff[2]
                + 2.0 * self.high_thrust_Mdot_coeff[3] * inputPower 
                + 6.0 * self.high_thrust_Mdot_coeff[4] * P2 ) * 1.0e+6
        else:
            self.polyThrust_1D = (self.high_Isp_Thrust_coeff[0]
                + self.high_Isp_Thrust_coeff[1] * inputPower
                + self.high_Isp_Thrust_coeff[2] * P2
                + self.high_Isp_Thrust_coeff[3] * P3
                + self.high_Isp_Thrust_coeff[4] * P4) * 1.0e+3

            self.poly_dThrust_dP_1D = (self.high_Isp_Thrust_coeff[1]
                + 2.0 * self.high_Isp_Thrust_coeff[2] * inputPower 
                + 3.0 * self.high_Isp_Thrust_coeff[3] * P2 
                + 4.0 * self.high_Isp_Thrust_coeff[4] * P3 ) * 1.0e+3

            self.polyMdot_1D = (self.high_Isp_Mdot_coeff[0]
                + self.high_Isp_Mdot_coeff[1] * inputPower
                + self.high_Isp_Mdot_coeff[2] * P2
                + self.high_Isp_Mdot_coeff[3] * P3
                + self.high_Isp_Mdot_coeff[4] * P4) * 1.0e+6

            self.poly_dMdotdP_1D = (self.high_Isp_Mdot_coeff[1]
                + 2.0 * self.high_Isp_Mdot_coeff[2] * inputPower 
                + 3.0 * self.high_Isp_Mdot_coeff[3] * P2 
                + 4.0 * self.high_Isp_Mdot_coeff[4] * P3 ) * 1.0e+6

            self.d2HighIsp_1D_Mdot_dP2 = (self.high_thrust_Mdot_coeff[2]
                + 2.0 * self.high_thrust_Mdot_coeff[3] * inputPower 
                + 6.0 * self.high_thrust_Mdot_coeff[4] * P2 ) * 1.0e+6

    #function to print a throttle file
    def print_throttle_file(self, filename):
        outputfile = open(filename, "w")

        outputfile.write('PPU efficiency,' + str(self.PPUefficiency) + '\n')
        outputfile.write('PPU min power (kW),' + str(self.PPUminpower) + '\n')
        outputfile.write('PPU max power (kW),' + str(self.PPUmaxpower) + '\n')
        outputfile.write('Throttle level,Mass flow rate (mg/s),Beam Current (A),Beam Voltage (V),Thrust (mN),Isp (s),Efficiency,Thruster input power (kW)\n')

        for setting in self.ThrottleSettings:
            setting.print_throttle_setting(outputfile)
        
        outputfile.close()

    #function to find the closest throttle setting to a reference performance point
    #this can be used in conjunction with an EMTG output file to find the throttle setting closest to what the optimizer is asking for at that time step
    def find_closest_throttle_setting(self, ReferenceThrottleSetting):
        #start by initializing the throttle level to zero
        TL = 'TL0'

        #we want the solution to the minimum-norm problem comparing the reference throttle setting to the throttle table
        #first we need to compute the normalized difference vector between the reference throttle setting and each throttle setting in the table
        normalized_difference_vectors = []
        normalized_distance_array = []

        for setting in self.ThrottleSettings:
            difference_vector = setting.compare(ReferenceThrottleSetting)
            normalized_difference_vectors.append([difference_vector[0] / ReferenceThrottleSetting.ThrusterPower,
                                                  difference_vector[1] / ReferenceThrottleSetting.Thrust,
                                                  difference_vector[2] / ReferenceThrottleSetting.Mdot,
                                                  difference_vector[3] / ReferenceThrottleSetting.Isp])
            normalized_distance_array.append(math.sqrt(normalized_difference_vectors[-1][0]**2 +
                                                       normalized_difference_vectors[-1][1]**2 +
                                                       normalized_difference_vectors[-1][2]**2 +
                                                       normalized_difference_vectors[-1][3]**2))
        
        
        return self.ThrottleSettings[normalized_distance_array.index(min(normalized_distance_array))].TL

    #function to make scatter plot from the throttle table
    def plot_throttle_table(self, create_fit_line=False, referencepoint=[], neighbors=[], MakeLabels=False,trajectory_mdot=[],trajectory_thrust=[],trajectory_power=[]):
        import matplotlib        
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d import proj3d
        PPUpower = []
        Isp = []
        Thrust = []
        Voltage = []
        Mdot = []
        PointLabels = []
        Efficiency = []
        Current = []
        ScaledCurrent = []

        for setting in self.ThrottleSettings:
            if not (referencepoint == setting or setting in neighbors):
                PPUpower.append(setting.ThrusterPower / self.PPUefficiency)
                Isp.append(setting.Isp)
                Thrust.append(setting.Thrust)
                Voltage.append(setting.BeamVoltage)
                Mdot.append(setting.Mdot)
                Efficiency.append(setting.efficiency)
            
                Current.append(setting.BeamCurrent)
                ScaledCurrent.append(Efficiency[-1]*Current[-1])
                PointLabels.append(setting.TL)

        ThrottleFigure = matplotlib.pyplot.figure()
        ThrottleFigure.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
        ThrottleAxes = ThrottleFigure.add_subplot(111, projection='3d')
        ThrottleAxes.view_init(elev=90, azim=-90)
        ThrottleAxes.scatter(numpy.array(PPUpower), numpy.array(Mdot), numpy.array(Thrust), color='blue', s=50, alpha=1.0, label='throttle settings')

        if MakeLabels:
            for i, txt in enumerate(PointLabels):
                x2D, y2D, _ = proj3d.proj_transform(PPUpower[i],Mdot[i],Thrust[i], ThrottleAxes.get_proj())
                ThrottleAxes.annotate(PointLabels[i], (x2D,y2D))

        if create_fit_line:
            Voltage_polynomial = numpy.poly1d(numpy.polyfit(PPUpower, Voltage, 4))
            Isp_polynomial = numpy.poly1d(numpy.polyfit(PPUpower, Isp, 4))
            Mdot_polynomial = numpy.poly1d(numpy.polyfit(PPUpower, Mdot, 4))
            Thrust_polynomial = numpy.poly1d(numpy.polyfit(PPUpower, Thrust, 4))
            Current_polynomial = numpy.poly1d(numpy.polyfit(PPUpower, Current, 4))
            powerlinepoints = numpy.linspace(min(PPUpower), max(PPUpower), 100)
            ThrottleAxes.plot(powerlinepoints, Mdot_polynomial(powerlinepoints), Thrust(powerlinepoints))

        if referencepoint:
            ThrottleAxes.scatter(referencepoint.ThrusterPower / self.PPUefficiency, referencepoint.Mdot, referencepoint.Thrust, color='crimson', s=50, alpha=1.0, label='reference throttle setting')
            if MakeLabels:
                x2D, y2D, _ = proj3d.proj_transform(referencepoint.ThrusterPower / self.PPUefficiency, referencepoint.Mdot, referencepoint.Thrust, ThrottleAxes.get_proj())
                ThrottleAxes.annotate(referencepoint.TL, (x2D, y2D))

        if neighbors:
            NPPUpower = []
            NVoltage = []
            NIsp = []
            NThrust = []
            NPointLabels = []
            NMdot = []
            for neighbor in neighbors:
                print(neighbor.TL)
                NPPUpower.append(neighbor.ThrusterPower / self.PPUefficiency)
                NMdot.append(neighbor.Mdot)
                NIsp.append(neighbor.Isp)
                NThrust.append(neighbor.Thrust)
                NVoltage.append(neighbor.BeamVoltage)
                NPointLabels.append(neighbor.TL)
            ThrottleAxes.scatter(numpy.array(NPPUpower), numpy.array(NMdot), numpy.array(NThrust), color='orange', s=50, alpha=1.0, label='neighbors of reference throttle setting')
            if MakeLabels:
                for i, txt in enumerate(NPointLabels):
                    x2D, y2D, _ = proj3d.proj_transform(NPPUpower[i],NMdot[i],NThrust[i], ThrottleAxes.get_proj())
                    ThrottleAxes.annotate(NPointLabels[i], (x2D,y2D))

        if len(trajectory_mdot) and len(trajectory_thrust) and len(trajectory_power):
            ThrottleAxes.scatter(numpy.array(trajectory_power),numpy.array(trajectory_mdot),numpy.array(trajectory_thrust), marker='x',c='black',s=100,alpha=1.0,label='Throttle outputs used in trajectory')
            ThrottleAxes.view_init(elev=0, azim=0)

        ThrottleAxes.set_xlabel('PPU input power (kW)')
        ThrottleAxes.set_ylabel('Mass flow rate (mg/s)')
        ThrottleAxes.set_zlabel('Thrust (mN)')
        ThrottleAxes.set_title(self.FileName)
        leg = ThrottleAxes.legend()
        leg.draggable()

        ThrottleFigure.show()
        
        if create_fit_line:
            return Voltage_polynomial, Thrust_polynomial, Mdot_polynomial, Isp_polynomial

    #function to enable 3D rotation of labels in throttle-surface scatter plot
    def UpdateLabelPositionsEvent(self, Figure, Axes, referencepoint=[], neighbors=[]):
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d import proj3d
        PPUpower = []
        Isp = []
        Thrust = []
        Voltage = []
        Mdot = []
        PointLabels = []
        Efficiency = []
        Current = []
        ScaledCurrent = []

        for setting in self.ThrottleSettings:
            if not (referencepoint == setting or setting in neighbors):
                PPUpower.append(setting.ThrusterPower / self.PPUefficiency)
                Isp.append(setting.Isp)
                Thrust.append(setting.Thrust)
                Voltage.append(setting.BeamVoltage)
                Mdot.append(setting.Mdot)
                Efficiency.append(setting.efficiency)
            
                Current.append(setting.BeamCurrent)
                ScaledCurrent.append(Efficiency[-1]*Current[-1])
                PointLabels.append(setting.TL)
        
        for i, txt in enumerate(PointLabels):
            x2D, y2D, _ = proj3d.proj_transform(PPUpower[i],Mdot[i],Thrust[i], Axes.get_proj())
            self.MainBoxLabels[i].xy = x2D,y2D
            self.MainBoxLabels[i].update_positions(Figure.canvas.renderer)

    def create_constant_Mdot_sets(self):
        self.ThrottleSets = []

        #first find all of the unique Mdot values in the ThrottleTable and make a set for each one
        Mdot = []
        for setting in self.ThrottleSettings:
            if setting.Mdot not in Mdot:
                Mdot.append(setting.Mdot)
        Mdot.sort()
        for mdot in Mdot:
            self.ThrottleSets.append(ThrottleSet(mdot))

        #then populate the sets
        for setting in self.ThrottleSettings:
            for Set in self.ThrottleSets:
                if setting.Mdot == Set.Mdot:
                    Set.add_setting(setting)
                    break;

        #sort the sets by power, then print them
        for Set in self.ThrottleSets:
            Set.sort()
            # print 'Set ', Set.Mdot
            # for setting in Set.ThrottleSettings:
                # print setting.TL
            # print '----'

    #function to create throttle tables for the high-Thrust non-dominated set and the high-Isp non-dominated set
    def create_non_dominated_sets(self):
        #start with the high-Thrust set
        self.high_thrust_set = ThrottleTable()
        self.high_thrust_set.PPUefficiency = self.PPUefficiency
        self.high_thrust_set.FileName = 'high_thrust_set'
        self.high_thrust_set.PPUmaxpower = self.PPUmaxpower
        self.high_thrust_set.PPUminpower = self.PPUminpower

        for setting in self.ThrottleSettings:
            isdominated = False
            for comparison_setting in self.ThrottleSettings:
                if setting.Thrust < comparison_setting.Thrust and setting.ThrusterPower >= comparison_setting.ThrusterPower:
                    isdominated = True
                    break
            if not isdominated:
                self.high_thrust_set.ThrottleSettings.append(copy.deepcopy(setting))

        #then the high-Isp set
        self.high_Isp_set = ThrottleTable()
        self.high_Isp_set.PPUefficiency = self.PPUefficiency
        self.high_Isp_set.FileName = 'high_Isp_set'
        self.high_Isp_set.PPUmaxpower = self.PPUmaxpower
        self.high_Isp_set.PPUminpower = self.PPUminpower
        

        for setting in self.ThrottleSettings:
            isdominated = False
            for comparison_setting in self.ThrottleSettings:
                if setting.Isp < comparison_setting.Isp and setting.ThrusterPower >= comparison_setting.ThrusterPower:
                    isdominated = True
                    break
            if not isdominated:
                self.high_Isp_set.ThrottleSettings.append(copy.deepcopy(setting))

        return self.high_thrust_set, self.high_Isp_set

    #function to create 1D polynomials
    def create_1D_polynomials(self):
        PPUpower = []
        Thrust = []
        Mdot = []
        for setting in self.high_thrust_set.ThrottleSettings:
            PPUpower.append(setting.ThrusterPower / self.PPUefficiency)
            Thrust.append(setting.Thrust)
            Mdot.append(setting.Mdot)

        
        self.High_Thrust_Mdot_polynomial = numpy.poly1d(numpy.polyfit(PPUpower, Mdot, 4))
        #self.high_thrust_thrust_rawpoly, self.high_thrust_thrust_rawerror, self.c = numpy.polyfit(PPUpower, Thrust, 4, full = True)
        self.High_Thrust_Thrust_polynomial = numpy.poly1d(numpy.polyfit(PPUpower, Thrust, 4))

        PPUpower = []
        Thrust = []
        Mdot = []
        for setting in self.high_Isp_set.ThrottleSettings:
            PPUpower.append(setting.ThrusterPower / self.PPUefficiency)
            Thrust.append(setting.Thrust)
            Mdot.append(setting.Mdot)

        self.High_Isp_Mdot_polynomial = numpy.poly1d(numpy.polyfit(PPUpower, Mdot, 4))
        self.High_Isp_Thrust_polynomial = numpy.poly1d(numpy.polyfit(PPUpower, Thrust, 4))

    #function to create 2D polynomial in P, mdot
    def create_2D_thrust_polynomial(self):
        PPUpower = []
        Thrust = []
        Mdot = []
        for setting in self.ThrottleSettings:
            PPUpower.append(setting.ThrusterPower / self.PPUefficiency)
            Thrust.append(setting.Thrust)
            Mdot.append(setting.Mdot)

        def polyfit2d(x, y, f, deg):
            from numpy.polynomial import polynomial
            import numpy as np
            x = np.asarray(x)
            y = np.asarray(y)
            f = np.asarray(f)
            deg = np.asarray(deg)
            vander = polynomial.polyvander2d(x, y, deg)
            vander = vander.reshape((-1,vander.shape[-1]))
            f = f.reshape((vander.shape[0],))
            c = np.linalg.lstsq(vander, f)[0]
            return c.reshape(deg+1)

        self.thrust_2D_fit = polyfit2db.polyfit2d(numpy.array(PPUpower), numpy.array(Mdot), numpy.array(Thrust), [4, 4])

    def find_nearest_throttle_setting_1D(self, PPUPower):
        #Step 1: convert PPU power to thruster power
        ThrusterPower = PPUPower * self.PPUefficiency

        min_delta = abs(self.ThrottleSettings[0].ThrusterPower - ThrusterPower)
        best_setting = self.ThrottleSettings[0]

        #Step 2: loop through the throttle settings
        for setting in self.ThrottleSettings:
            delta = abs(setting.ThrusterPower - ThrusterPower)
            if delta < min_delta + 1.0e-2:
                min_delta = delta
                best_setting = setting

        #Step 3: return
        return best_setting