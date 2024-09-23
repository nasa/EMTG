#tool for reading NSGA-II population and archive files
#for use with EMTG-NSGAII outer-loop by Vavrina and Englander
#Python interface by Jacob Englander begun 3-16-2014

import Universe

import os
import platform
import numpy as np
from scipy.integrate import ode
import matplotlib
import matplotlib.dates as dates
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
import copy
import wx
import datetime
from matplotlib.backends.backend_agg import FigureCanvas as FigureCanvas

class NSGAII_outerloop_solution(object):
    
    #constructor
    def __init__(self, input_line, column_headers):
        self.initialize()
        self.parse_input_line(input_line, column_headers)

    #clear function
    def initialize(self):
        self.Xouter = [] #outer-loop decision vector
        self.Xinner = [] #inner-loop decision vector
        self.objective_values = [] #vector of objective function values
        self.power_system_size = [] #if applicable, power level of the array in kW
        self.thruster = []
        self.number_of_thrusters = []
        self.launch_vehicle = []
        self.launch_date = []
        self.description = '' # case name, can be parsed for data
        self.outputfilename = ''
        self.mission_sequence = []
        self.generation_found = [] #what generation was this solution found?
        self.timestamp = [] #at what time, in seconds from program start, was this solution found?
        self.Legal_Solution = False

    #line parser
    def parse_input_line(self, input_line, column_headers):
        #strip off the newline character
        input_line = input_line.strip("\n")

        #strip the input line by commas
        input_cell = input_line.split(',')

        #declare arrays of launch vehicle and thruster names
        LV_names = ['AV401','AV411','AV421','AV431','AV501','AV511','AV521','AV531','AV541','AV551',
                    'F910','F911','AV551s48','F9H','D4H','SLSb1']
        LV_preference_rank = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 13, 14, 15]

        thruster_names = ['CHEM','NSTAR','XIPS25','BPT4000HIsp','BPT4000Hthrust','BPT4000XHIsp',
                          'NEXTHIspv9','VASIMRargon','VSIxenonhall','NEXTHIspv10','NEXTHthrustv10',
                          'BPT4000MALTO','NEXIS','H6MS','BHT20K','HiVHAc','13kWSTMDHallHisp','13kWSTMDHallHthrust',
                          'NEXT-TT11-Hisp','NEXT-TT11-Hthrust','NEXT-TT11-expanded',
                          '13kWSTMDHall-10-1-2014-Hisp','13kWSTMDHall-10-1-2014-Mthrust','13kWSTMDHall-10-1-2014-Hthrust']

        for column_index in range(0, len(column_headers)):
            if column_index < len(input_cell):
                if column_headers[column_index] == 'Generation found':
                    self.generation_found = int(input_cell[column_index])

                elif column_headers[column_index] == 'File name':
                    self.outputfilename = input_cell[column_index]

                elif column_headers[column_index] == 'timestamp':
                    self.timestamp = int(input_cell[column_index])
            
                elif column_headers[column_index] == 'Description':
                    self.description = input_cell[column_index]

                    #find the mission sequence descriptor
                    left_parenthesis_index = self.description.find('(')
                    self.mission_sequence = self.description[left_parenthesis_index+1:].strip(')')
                
                    #reconstruct the full mission description from the case name
                    descriptioncell = self.description.split('_')

                    for descriptionitem in descriptioncell:
                        if descriptionitem.rfind('kW') == len(descriptionitem) - 2: #this entry encodes power system size
                            self.power_system_size = float(descriptionitem.strip('kW'))

                        if descriptionitem.find('nTh') > 0:
                            self.number_of_thrusters = descriptionitem.strip('nTh')

                        for LV_name in LV_names:
                            if descriptionitem == LV_name:
                                self.launch_vehicle = descriptionitem

                        for thruster_name in thruster_names:
                            if descriptionitem == thruster_name:
                                self.thruster = thruster_name

                                if self.thruster == 'CHEM':
                                    self.power_system_size = 0.0

                elif column_headers[column_index] == 'BOL power at 1 AU (kW)' \
                    or column_headers[column_index] == 'Launch epoch (MJD)' \
                    or column_headers[column_index] == 'Flight time (days)' \
                    or column_headers[column_index] == 'First journey departure C3 (km^2/s^2)' \
                    or column_headers[column_index] == 'Final journey arrival C3 (km^2/s^2)' \
                    or column_headers[column_index] == 'Final journey arrival declination(deg)' \
                    or column_headers[column_index] == 'Total deterministic delta-v (km/s)' \
                    or column_headers[column_index] == 'Total propellant mass including margin (kg)' \
                    or column_headers[column_index] == 'Final journey interface velocity (km/s)' \
                    or column_headers[column_index] == 'Thruster duty cycle' \
                    or column_headers[column_index] == 'Normalized aggregate control':
                    #this entry is an objective function value
                    self.objective_values.append(float(input_cell[column_index]))
                elif column_headers[column_index] == 'Thruster preference' \
                    or column_headers[column_index] == 'Number of thrusters' \
                    or column_headers[column_index] == 'Launch vehicle preference':
                    if float(input_cell[column_index]) < 1.0e+50:
                        if column_headers[column_index] == 'Launch vehicle preference':
                            self.objective_values.append(LV_preference_rank[int(input_cell[column_index]) - 1])
                        else:
                            self.objective_values.append(abs(int(input_cell[column_index])))
                    else:
                        self.objective_values.append(float(input_cell[column_index]))
                elif column_headers[column_index] == 'Delivered mass to final target (kg)' \
                    or column_headers[column_index] == 'Final journey mass increment (for maximizing sample return)' \
                    or column_headers[column_index] == 'Dry mass margin' \
                    or column_headers[column_index] == 'bus power (kW)':
                    self.objective_values.append(-float(input_cell[column_index]))
                elif column_headers[column_index].find('Gene ') > 0:
                    self.Xouter.append(int(input_cell[column_index]))
                elif  column_headers[column_index] == 'Number of journeys' \
                    or column_headers[column_index] == 'Point-group value':
                    if float(input_cell[column_index]) < 1.0e+50:
                        self.objective_values.append(abs(int(input_cell[column_index])))
                    else:
                        self.objective_values.append(float(input_cell[column_index]))

                elif input_cell[column_index] != '': #this entry is a member of the inner-loop decision vector
                    self.Xinner.append(float(input_cell[column_index]))

    def plot_solution(self, PopulationAxes, PopulationFigure, ordered_list_of_objectives, colorbar, lowerbounds, upperbounds):
        if len(ordered_list_of_objectives) == 2: #2D
            self.point = PopulationAxes.scatter(self.objective_values[ordered_list_of_objectives[0]], self.objective_values[ordered_list_of_objectives[1]], s=50, c='b', marker='o', lw=0, picker=1)
        elif len(ordered_list_of_objectives) == 3: #3D
                self.point = PopulationAxes.scatter(self.objective_values[ordered_list_of_objectives[0]], self.objective_values[ordered_list_of_objectives[1]], self.objective_values[ordered_list_of_objectives[2]], s=50, c='b', marker='o', lw=0, picker=1)
        else: #4D
            if self.colorbar is None:
                self.point = PopulationAxes.scatter(self.objective_values[ordered_list_of_objectives[0]], self.objective_values[ordered_list_of_objectives[1]], self.objective_values[ordered_list_of_objectives[2]], s=50, c=self.objective_values[ordered_list_of_objectives[4]], marker='o', lw=0, picker=1)
                self.point.set_clim([lowerbounds[-1],upperbounds[-1]])
                colorbar = PopulationFigure.colorbar(self.point)
            else:
                self.point = PopulationAxes.scatter(self.objective_values[ordered_list_of_objectives[0]], self.objective_values[ordered_list_of_objectives[1]], self.objective_values[ordered_list_of_objectives[2]], s=50, c=self.objective_values[ordered_list_of_objectives[4]], marker='o', lw=0, picker=1)
                self.point.set_clim([lowerbounds[-1],upperbounds[-1]])

        self.picker = self.point.figure.canvas.mpl_connect('pick_event', self.onpick)

    def onpick(self, event):
        #description = []
        #for objective_index in ordered_list_of_objectives:
        #    description = description + self.objective_column_headers[objective_index] + ': ' + str(
        ind = event.ind[0]
        x, y, z = event.artist._offsets3d
        print(self.description)
        #print self.description, x[ind], y[ind], z[ind]
        #print ind


#top-level container of NSGAII_outerloop_solution objects
class NSGAII_outerloop_population(object):
        
    #constructor
    def __init__(self, population_file_name):
        self.clear()
        self.parse_population_file(population_file_name)

    def clear(self):
        self.solutions = [] #vector of NSGAII_outerloop_solution objects
        self.legal_solutions = []
        self.global_column_headers = []
        self.gene_column_headers = []
        self.objective_column_headers = []
        self.number_of_feasible_solutions = []
        self.points_array = []
        self.ordered_list_of_objectives = []
        self.colorbar = []
        self.solution_points = []

    #method to read a population file
    def parse_population_file(self, population_file_name):
        #Step 1: attempt to open a population file
        if os.path.isfile(population_file_name):
            inputfile = open(population_file_name, "r")
            self.success = 1
        else:
            print("Unable to open", population_file_name, "EMTG Error")
            self.success = 0
            return

        #Step 2: scan through the file
        linenumber = 0
        for line in inputfile:
            #strip off the newline character
            line = line.replace("\n","")
            linenumber = linenumber + 1

            #the fourth line of the population file contains the column headers
            if linenumber == 4:
               self.global_column_headers = line.split(',')
               for header in self.global_column_headers:
                if header == 'BOL power at 1 AU (kW)' \
                    or header == 'Launch epoch (MJD)' \
                    or header == 'Flight time (days)' \
                    or header == 'Thruster preference' \
                    or header == 'Number of thrusters' \
                    or header == 'Launch vehicle preference' \
                    or header == 'Delivered mass to final target (kg)' \
                    or header == 'Final journey mass increment (for maximizing sample return)' \
                    or header == 'First journey departure C3 (km^2/s^2)' \
                    or header == 'Final journey arrival C3 (km^2/s^2)' \
                    or header == 'Final journey arrival declination(deg)' \
                    or header == 'Total deterministic delta-v (km/s)' \
                    or header == 'Inner-loop objective function' \
                    or header == 'Point-group value' \
                    or header == 'Total propellant mass including margin (kg)' \
                    or header == 'Number of journeys' \
                    or header == 'Dry mass margin' \
                    or header == 'Final journey interface velocity (km/s)' \
                    or header == 'Thruster duty cycle' \
                    or header == 'Normalized aggregate control'\
                    or header == 'bus power (kW)'\
                    or header == 'Final journey arrival declination (deg)':
                    self.objective_column_headers.append(header)

            #the fifth line of the population file contains the gene names
            elif linenumber == 5:
                self.gene_column_headers = line.split(',')

            #every line after the fifth is a solution line
            elif linenumber > 5:
                tempSolution = NSGAII_outerloop_solution(line, self.global_column_headers)
                self.solutions.append(tempSolution)
        inputfile.close()

    #method to plot the population
    #input is an ordered list of objectives, [x, y, z, color]. If there are two objectives, a monochrome 2D plot will be shown. If there are three objectives, a monochrome 3D plot will be shown.
    #if there are four, a colored 3D plot will be shown. If there are more than four there will be an error message.
    def plot_population(self, ordered_list_of_objectives, LowerBounds = [], UpperBounds = [], TimeUnit = 1, EpochUnit = 1, FontSize = 10, BaseMarkerSize = 20, OutputWindow = [], PopulationFigure = [], PopulationAxes = [], Filter = []):
        self.ordered_list_of_objectives = ordered_list_of_objectives
        self.LowerBounds = LowerBounds
        self.UpperBounds = UpperBounds
        self.TimeUnit = TimeUnit
        self.EpochUnit = EpochUnit
        self.BaseMarkerSize = BaseMarkerSize
        self.OutputWindow = OutputWindow
        self.Filter = Filter

        #first check to see if the correct number of objective function indices were supplied
        if len(self.ordered_list_of_objectives) < 2 or len(self.ordered_list_of_objectives) > 5:
            if not self.OutputWindow == []:
                self.OutputWindow.WriteText("NSGAII_outerloop_population::plot_population ERROR. You must specify between two and five objective functions to plot.")
            return

        NeedtoShow = False
        if PopulationFigure == [] or PopulationAxes == []:
            
            matplotlib.rcParams.update({'font.size': 20})
            self.PopulationFigure = matplotlib.pyplot.figure()
            #self.PopulationFigure.rcParams.update({'axes.titlesize': 'x-large'})


            #plt.rcParams.update({'axes.titlesize': 'small'})

            self.PopulationFigure.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
            if len(ordered_list_of_objectives) == 2:
                self.PopulationAxes = self.PopulationFigure.add_subplot(111)
            else:
                self.PopulationAxes = self.PopulationFigure.add_subplot(111, projection='3d')
            NeedtoShow = True
        else:
            self.PopulationFigure = PopulationFigure
            self.PopulationAxes = PopulationAxes

        #build up a list of objective values to be plotted
        self.objective_values_matrix = []
        self.legal_solutions = []

        def process_solution(solution):
            if max(solution.objective_values) < 1.0e+99:
                solution.Legal_Solution = True
                if not (self.LowerBounds == [] or self.UpperBounds == []):
                    #if bounds were supplied, check to see if the solution fits inside the bounds
                    for obj in range(0, len(self.ordered_list_of_objectives)):
                        if solution.objective_values[ordered_list_of_objectives[obj]] < self.LowerBounds[obj] or solution.objective_values[ordered_list_of_objectives[obj]] > self.UpperBounds[obj]:
                            solution.Legal_Solution = False
                            break
                if solution.Legal_Solution:
                    self.legal_solutions.append(solution)

                    if self.objective_column_headers[self.ordered_list_of_objectives[objective_index]] == 'Flight time (days)' and self.TimeUnit == 0:
                        objective_values_vector.append(solution.objective_values[ordered_list_of_objectives[objective_index]] / 365.25)
                    elif self.objective_column_headers[self.ordered_list_of_objectives[objective_index]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                        objective_values_vector.append(dates.date2num(datetime.datetime.fromtimestamp(wx.DateTimeFromJDN(solution.objective_values[ordered_list_of_objectives[objective_index]] + 2400000.5).GetTicks())))
                    else:
                        objective_values_vector.append(copy.deepcopy(float(solution.objective_values[ordered_list_of_objectives[objective_index]])))

        for objective_index in range(0, len(self.ordered_list_of_objectives)):
            objective_values_vector = []
            for solution in self.solutions:
                #we are only interested in solutions which contain the bodies in our filter list
                #if there is no filter list, show all of the bodies
                if self.Filter == []:
                    process_solution(solution)
                    
                else:
                    if self.Filter.filter_list == []:
                        FoundFlag = False
                        for FilterBody in self.Filter.object_list:
                            if FilterBody in solution.description:
                                process_solution(solution)
                                FoundFlag = True
                        if not FoundFlag:
                            solution.Legal_Solution = False
                    else:
                        FoundFlag = False
                        for FilterBody in self.Filter.filter_list:
                            if FilterBody in solution.description:
                                process_solution(solution)
                                FoundFlag = True
                        if not FoundFlag:
                            solution.Legal_Solution = False
            self.objective_values_matrix.append(np.array(objective_values_vector))

        #determine upper and lower bounds on each objective
        self.upperbounds = []
        self.lowerbounds = []
        for objective_index in range(0, len(self.ordered_list_of_objectives)):
            if (len(self.objective_values_matrix[objective_index]) > 0):
                self.upperbounds.append(self.objective_values_matrix[objective_index].max())
                self.lowerbounds.append(self.objective_values_matrix[objective_index].min())
            else:
                self.lowerbounds.append(1.0e-10)
                self.upperbounds.append(-1.0e-10)

        #plot each solution
        self.plot_solution_points()

        #set the axes labels
        if self.objective_column_headers[self.ordered_list_of_objectives[0]] == 'Flight time (days)' and self.TimeUnit == 0:
            self.PopulationAxes.set_xlabel('Flight time (years)')
        elif self.objective_column_headers[self.ordered_list_of_objectives[0]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
            self.PopulationAxes.set_xlabel('Launch Epoch (TDB Gregorian)')
            self.PopulationAxes.w_xaxis.set_major_formatter(ticker.FuncFormatter(self.format_date))
        else:
            self.PopulationAxes.set_xlabel(self.objective_column_headers[self.ordered_list_of_objectives[0]])

        if self.objective_column_headers[self.ordered_list_of_objectives[1]] == 'Flight time (days)' and self.TimeUnit == 0:
            self.PopulationAxes.set_ylabel('Flight time (years)')
        elif self.objective_column_headers[self.ordered_list_of_objectives[1]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
            self.PopulationAxes.set_ylabel('Launch Epoch (TDB Gregorian)')
            self.PopulationAxes.w_yaxis.set_major_formatter(ticker.FuncFormatter(self.format_date))
        else:
            self.PopulationAxes.set_ylabel(self.objective_column_headers[self.ordered_list_of_objectives[1]])

        if len(ordered_list_of_objectives) > 2:
            if self.objective_column_headers[self.ordered_list_of_objectives[2]] == 'Flight time (days)' and self.TimeUnit == 0:
                self.PopulationAxes.set_zlabel('Flight time (years)')
            elif self.objective_column_headers[self.ordered_list_of_objectives[2]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                self.PopulationAxes.set_zlabel('Launch Epoch (TDB Gregorian)')
                self.PopulationAxes.w_zaxis.set_major_formatter(ticker.FuncFormatter(self.format_date))
            else:
                self.PopulationAxes.set_zlabel(self.objective_column_headers[self.ordered_list_of_objectives[2]])
            self.PopulationAxes.autoscale_view(tight=True, scalex=True, scaley=True, scalez=True)
        else:
            self.PopulationAxes.autoscale_view(tight=True, scalex=True, scaley=True)
        self.PopulationAxes.grid(b=True)

        if NeedtoShow:
            self.PopulationFigure.show()

        for item in ([self.PopulationAxes.title, self.PopulationAxes.xaxis.label, self.PopulationAxes.yaxis.label, self.PopulationAxes.zaxis.label] +
             self.PopulationAxes.get_xticklabels() + self.PopulationAxes.get_yticklabels() + self.PopulationAxes.get_zticklabels()):
             item.set_fontsize(FontSize)

        self.colorbar.ax.tick_params(labelsize=FontSize)
        self.colorbar.ax.yaxis.label.set_size(FontSize)

    def plot_solution_points(self):
        self.solution_names = []
        X = []
        Y = []
        Z = []
        C = []
        S = []
        for solution in self.solutions:
            if solution.Legal_Solution:
                if self.objective_column_headers[self.ordered_list_of_objectives[0]] == 'Flight time (days)' and self.TimeUnit == 0:
                    X.append(solution.objective_values[self.ordered_list_of_objectives[0]] / 365.25)
                elif self.objective_column_headers[self.ordered_list_of_objectives[0]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                    X.append(dates.date2num(datetime.datetime.fromtimestamp(wx.DateTimeFromJDN(solution.objective_values[0] + 2400000.5).GetTicks())))
                else:
                    X.append(solution.objective_values[self.ordered_list_of_objectives[0]])

                if self.objective_column_headers[self.ordered_list_of_objectives[1]] == 'Flight time (days)' and self.TimeUnit == 0:
                    Y.append(solution.objective_values[self.ordered_list_of_objectives[1]] / 365.25)
                elif self.objective_column_headers[self.ordered_list_of_objectives[1]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                    Y.append(dates.date2num(datetime.datetime.fromtimestamp(wx.DateTimeFromJDN(solution.objective_values[1] + 2400000.5).GetTicks())))
                else:
                    Y.append(solution.objective_values[self.ordered_list_of_objectives[1]])

                if len(self.ordered_list_of_objectives) > 2:
                    if self.objective_column_headers[self.ordered_list_of_objectives[2]] == 'Flight time (days)' and self.TimeUnit == 0:
                        Z.append(solution.objective_values[self.ordered_list_of_objectives[2]] / 365.25)
                    elif self.objective_column_headers[self.ordered_list_of_objectives[2]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                        Z.append(dates.date2num(datetime.datetime.fromtimestamp(wx.DateTimeFromJDN(solution.objective_values[2] + 2400000.5).GetTicks())))
                    else:
                        Z.append(solution.objective_values[self.ordered_list_of_objectives[2]])
                else:
                    Z.append(0.0)

                if len(self.ordered_list_of_objectives) > 3:
                    if self.objective_column_headers[self.ordered_list_of_objectives[3]] == 'Flight time (days)' and self.TimeUnit == 0:
                        C.append(solution.objective_values[self.ordered_list_of_objectives[3]] / 365.25)
                    elif self.objective_column_headers[self.ordered_list_of_objectives[3]] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                        C.append(dates.date2num(datetime.datetime.fromtimestamp(wx.DateTimeFromJDN(solution.objective_values[3] + 2400000.5).GetTicks())))
                    else:
                        C.append(solution.objective_values[self.ordered_list_of_objectives[3]])
                else:
                    C.append(0.0)


                #if there is a fifth objective, size the markers to reflect it
                if len(self.ordered_list_of_objectives) > 4:
                    S.append(self.BaseMarkerSize * solution.objective_values[self.ordered_list_of_objectives[4]])# / (self.upperbounds[4] - self.lowerbounds[4]))
                else:
                    S.append(self.BaseMarkerSize)

        if not self.solution_points == []:
            self.solution_points.remove()
        self.solution_points = self.PopulationAxes.scatter(X, Y, Z, s=S, c=C, marker='o', edgecolors='none', picker=1)
        self.number_of_feasible_solutions = len(X)
        self.OutputWindow.WriteText(str(self.number_of_feasible_solutions) + ' solutions meet the applied filters.\n')
        
        self.picker = self.PopulationFigure.canvas.mpl_connect('pick_event', self.onpick)
        if len(self.ordered_list_of_objectives) > 3 and len(self.legal_solutions) > 0:
            if not self.colorbar == []:
                self.PopulationFigure.delaxes(self.PopulationFigure.axes[1])
                self.PopulationAxes.set_position(self.pre_colorbar_position)
            self.pre_colorbar_position = self.PopulationAxes.get_position()
            self.colorbar = self.PopulationFigure.colorbar(self.solution_points, label=self.objective_column_headers[self.ordered_list_of_objectives[3]], shrink=0.9)
    
    
    def force_update(self, event):
        for solution in self.solutions:
            if solution.Legal_Solution:
                solution.points.changed()

    def onpick(self, event):
        ind = event.ind[0]
        if not self.OutputWindow == []:
            if len(self.ordered_list_of_objectives) == 2: #2D
                self.OutputWindow.WriteText('2D picker not implemented\n')
            elif len(self.ordered_list_of_objectives) >= 3: #3D or 4D plot
                x, y, z = event.artist._offsets3d

                idx = np.where(self.objective_values_matrix[0] == x[ind])
                idy = np.where(self.objective_values_matrix[1] == y[ind])
                idz = np.where(self.objective_values_matrix[2] == z[ind])

                SolutionIndex = np.intersect1d(idx[0], np.intersect1d(idy[0], idz[0]))[0]
                ThisSolution = self.legal_solutions[SolutionIndex]

                self.OutputWindow.WriteText('Description: ' + ThisSolution.description + '\n')
                if ThisSolution.outputfilename != '':
                    
                    if platform.system() == 'Windows':
                        self.OutputWindow.WriteText('Output file: ' + ThisSolution.outputfilename.replace('//','\\') + '.emtg\n')
                    else:
                        self.OutputWindow.WriteText('Output file: ' + ThisSolution.outputfilename + '.emtg\n')

                for objIndex in range(0, len(self.objective_column_headers)):
                    if self.objective_column_headers[objIndex] == 'Flight time (days)' and self.TimeUnit == 0:
                        self.OutputWindow.WriteText('Flight time (years): ' + str(ThisSolution.objective_values[objIndex] / 365.25) + '\n')
                    elif self.objective_column_headers[objIndex] == 'Launch epoch (MJD)' and self.EpochUnit == 0:
                        dt = datetime.datetime.fromtimestamp(wx.DateTimeFromJDN(ThisSolution.objective_values[objIndex] + 2400000.5).GetTicks())
                        self.OutputWindow.WriteText('Launch Epoch (TDB Gregorian): '+ dt.strftime('%m/%d/%Y') + '\n')
                    elif self.objective_column_headers[objIndex] == 'Thruster preference':
                        self.OutputWindow.WriteText('Thruster type: ' + ThisSolution.thruster + '\n')
                    elif self.objective_column_headers[objIndex] == 'Launch vehicle preference':
                        self.OutputWindow.WriteText('Launch vehicle: ' + ThisSolution.launch_vehicle + '\n')
                    else:
                        self.OutputWindow.WriteText(str(self.objective_column_headers[objIndex]) + ': ' +  str(ThisSolution.objective_values[objIndex]) + '\n')

                self.OutputWindow.WriteText('---------------------------------------------------------------------------------------------\n')

    def format_date(self, x, pos=None):
        return dates.num2date(x).strftime('%Y-%m-%d')

    def generate_body_prevalence_report(self, Universe):
        #this method generates a list of tuples, (BodyName, NumberOfOccurrences)
        #pass in a Universe object, return a list of tuples
        
        prevalence_report = []
        for body in Universe.bodies:
            number_of_occurrences = 0
            for solution in self.solutions:
                #only book-keep feasible solutions

                if solution.objective_values[0] < 1.0e+99 and body.shortname in solution.description:
                    number_of_occurrences += 1
            prevalence_report.append((body.name, number_of_occurrences))

        return prevalence_report