import numpy as np
#from matplotlib import ticker
#import matplotlib.pyplot as plt
import os
import datetime

import astropy.time

#EMTG path item class
class PathItem(object):
    def __init__(self, step_index, time_index, number_of_steps_since_last_restart, seconds_since_last_restart, NLP_success, NLP_distance_from_filament, NLP_inform_code,  fitness):
        self.clear()

        if step_index >= 0:
            self.setvalues(step_index, time_index, number_of_steps_since_last_restart, seconds_since_last_restart, NLP_success, NLP_distance_from_filament, NLP_inform_code,  fitness)

    def clear(self):
        self.step_index = []
        self.time_index = []
        self.number_of_steps_since_last_restart = []
        self.seconds_since_last_restart = []
        self.NLP_success = []
        self.NLP_distance_from_filament = []
        self.NLP_inform_code = []
        self.fitness = []

    def setvalues(self, step_index, time_index, number_of_steps_since_last_restart, seconds_since_last_restart, NLP_success, NLP_distance_from_filament, NLP_inform_code, fitness):
        self.step_index = step_index
        self.time_index = time_index
        self.number_of_steps_since_last_restart = number_of_steps_since_last_restart
        self.seconds_since_last_restart = seconds_since_last_restart
        self.NLP_success = NLP_success
        self.NLP_distance_from_filament = NLP_distance_from_filament
        self.NLP_inform_code = NLP_inform_code
        self.fitness = fitness

    def write(self, outputfile):
        outputfile.write(str(int(self.step_index)) + ','
                         + str(self.time_index) + ','
                         + str(int(self.number_of_steps_since_last_restart)) + ','
                         + str(int(self.seconds_since_last_restart)) + ','
                         + str(int(self.NLP_success)) + ','
                         + str(self.NLP_distance_from_filament) + ','
                         + str(int(self.NLP_inform_code)) + ','
                         + str(self.fitness) + '\n')

#EMTG path class
class Path(object):
    def __init__(self, step_index, time_index, number_of_steps_since_last_restart, seconds_since_last_restart, NLP_success, NLP_distance_from_filament, NLP_inform_code, fitness):
        self.clear()

        if step_index >= 0:
            self.append(step_index, time_index, number_of_steps_since_last_restart, seconds_since_last_restart, NLP_success, NLP_distance_from_filament, NLP_inform_code, fitness)

    def clear(self):
        self.PathItems = []

    def append(self, step_index, time_index, number_of_steps_since_last_restart, seconds_since_last_restart, NLP_success, NLP_distance_from_filament, NLP_inform_code, fitness):
        self.PathItems.append(PathItem(step_index, time_index, number_of_steps_since_last_restart, seconds_since_last_restart, NLP_success, NLP_distance_from_filament, NLP_inform_code, fitness))

    def write(self, outputfile):
        for Item in self.PathItems:
            Item.write(outputfile)

#EMTG archive item class
class ArchiveItem(object):
    def __init__(self, namestring):
        self.clear()

        if namestring != []:
            self.setname(namestring)

    def clear(self):
        self.values = []
        self.name = []

    def setname(self, namestring):
        self.name = namestring

    def append(self, value):
        self.values.append(value)

    def pop(self, popvalue):
        self.values.pop(popvalue)

    def remove(self, removevalue):
        self.values.remove(removevalue)

#EMTG archive class
class Archive(object):
    def __init__(self, input_file_name):
        self.clear()
        if input_file_name != []:
            self.parse_archive(input_file_name)

    def clear(self):
        self.ArchiveItems = []
        self.Paths = []

    def parse_archive(self, input_file_name):
        #Step 1: open the file
        self.filename = input_file_name

        if os.path.isfile(self.filename):
            inputfile = open(input_file_name, "r")
            self.success = 1
        else:
            print("Unable to open", input_file_name, "EMTG Error")
            self.success = 0
            return

        #Step 2: scan through the file
        linenumber = 0
        for line in inputfile:
            #strip off the newline character
            line = line.replace("\n","")
            linenumber = linenumber + 1

            #the first line of the archive defines everything else
            if linenumber == 1:
               linecell = line.split(',')

               for entry in linecell:
                   self.ArchiveItems.append(ArchiveItem(entry))

            else:
                linecell = line.split(',')
                entryindex = 0
                for entry in linecell:
                    self.ArchiveItems[entryindex].append(float(entry))
                    entryindex += 1

        #Close the file
        inputfile.close()

    def write_archive(self, output_file_name):
        outputfile = open(output_file_name, "w")

        outputfile.write(self.ArchiveItems[0].name)
        for entry in self.ArchiveItems[1:]:
            outputfile.write(',' + entry.name)
        outputfile.write('\n')

        for linenumber in range(0,len(self.ArchiveItems[0].values)):
            outputfile.write(str(self.ArchiveItems[0].values[linenumber]))
            for item in self.ArchiveItems[1:]:
                outputfile.write(',' + str(item.values[linenumber]))
            outputfile.write('\n')

        outputfile.close()

    def plot_objective_vs_arrival_date(self):
        #first we need to create an array of arrival dates

        ArrivalDates = []
        for linenumber in range(0,len(self.ArchiveItems[0].values)):
            currentArrivalDate = 0.0
            for item in self.ArchiveItems:
                
                if "flight time" in item.name or "launch epoch" in item.name or "stay time" in item.name or "wait time" in item.name:
                    currentArrivalDate += item.values[linenumber]
            ArrivalDates.append(currentArrivalDate)

        date_vector = []
        for date in ArrivalDates:
            datestring = astropy.time.Time(date, format='mjd', scale='tdb', out_subfmt='date').iso
            date_vector.append(datetime.datetime.strptime(datestring,'%Y-%m-%d').date())
        
        self.DataFigure = plt.figure()
        self.DataAxes = self.DataFigure.add_axes([0.1, 0.1, 0.8, 0.8])
        self.DataAxes.scatter(date_vector, -np.array(self.ArchiveItems[-1].values))
        self.DataAxes.set_xlabel('Arrival Epoch')

        def format_date(x, pos=None):
            return datetime.datetime(x).strftime('%m-%d-%Y')

        self.DataAxes.xaxis.set_major_formatter(ticker.FuncFormatter(format_date))
        self.DataFigure.autofmt_xdate()
        plt.show()

    def assemble_paths(self):
        #step 1: set the current number of restarts to a negative number
        reset_number = -1
        reference_step_count = 0
        reference_timestamp = 0

        #step 2: iterate through the archive and detect path data
        for linenumber in range(0,len(self.ArchiveItems[0].values)):
            current_reset_count = []
            current_step_count = []
            current_fitness = []
            current_timestamp = []
            current_NLP_success = []
            current_NLP_distance_from_filament = []
            current_NLP_inform_code = []

            for item in self.ArchiveItems:
                if "reset count" in item.name:
                    current_reset_count = item.values[linenumber]
                if "step count" in item.name:
                    current_step_count = item.values[linenumber]
                if "solution timestamp" in item.name:
                    current_timestamp = item.values[linenumber]
                if "NLP success" in item.name:
                    current_NLP_success = item.values[linenumber]
                if "NLP distance from equality filament" in item.name:
                    current_NLP_distance_from_filament = item.values[linenumber]
                if "NLP inform code" in item.name:
                    current_NLP_inform_code = item.values[linenumber]
                if "Objective function" in item.name:
                    current_fitness = item.values[linenumber]
            
            #detect if we are on a new path
            if current_reset_count > reset_number:
                #create a new path
                print("New path at step ", int(current_step_count))
                reset_number = current_reset_count
                reference_step_count = current_step_count
                reference_timestamp = current_timestamp
                self.Paths.append(Path(current_step_count, current_timestamp, 0, 0, 0, 0, 0, current_fitness))
            else:
                #if we are not on a new path, append this step's results to the current path
                self.Paths[-1].append(current_step_count,
                                      current_timestamp,
                                      current_step_count - reference_step_count,
                                      current_timestamp - reference_timestamp, 
                                      current_NLP_success,
                                      current_NLP_distance_from_filament,
                                      current_NLP_inform_code,
                                      current_fitness)

    def write_paths(self, output_file_name):
        outputfile = open(output_file_name, "w")

        current_path = 1

        for Path in self.Paths:
            outputfile.write('Path #' + str(int(current_path)) + '\n')
            outputfile.write('Step #, Timestamp (s), steps since path start, seconds since path start, NLP success?, NLP exit distance from equality filament, NLP inform code, fitness\n')
            current_path += 1
            Path.write(outputfile)

        outputfile.close()