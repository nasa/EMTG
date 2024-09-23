import wx
import Archive

#class to post-process archive files
class ArchivePanel(wx.Panel):
    

    def __init__(self, parent, archive):
        wx.Panel.__init__(self, parent)

        self.archive = archive
        self.assemble_data_arrays()

        print('banana')

        #by calling methods of the Archive class,
        #search through header line for:
        #launch epoch
        #all flight times
        #mass at each encounter
        #objective function
        #C3
        #DLA

        #then compute:
        #total flight time
        #arrival epoch for each encounter
        #any wait times

        #assemble a set of radio button choices for each of these quantities to make them plottable
        #allow for easy plotting of any data set against any other data set in two, three, or four (color) dimensions
        #perhaps need to make an additional "ArchivePlot" class to handle that or make it part of the Archive class?

    def assemble_data_arrays(self):
        #create a vector of launch epochs - this is always the first column of the archive
        self.launch_epoch = self.archive.ArchiveItems[0]
        
        #create vectors of flight times
        self.flight_times = []
        for item in self.archive.ArchiveItems:
            if 'phase flight time' in item.name:
                self.flight_times.append(item)
        
        #create vectors of wait times
        self.wait_times = []
        for item in self.archive.ArchiveItems:
            if 'wait time' in item.name:
                self.wait_times.append(item)
        
        #create vector of mass at each encounter
        self.arrival_masses = []
        for item in self.archive.ArchiveItems:
            if 'arrival mass' in item.name:
                self.arrival_masses.append(item)

        #C3
        self.C3 = Archive.ArchiveItem('C3')
        for item in self.archive.ArchiveItems:
            if 'j0p0: magnitude of outgoing velocity asymptote' in item.name:
                for entry in item.values:
                    self.C3.values.append(entry ** 2)

        #DLA
        self.DLA = Archive.ArchiveItem('DLA')
        for item in self.archive.ArchiveItems:
            if 'j0p0: DEC of departure asymptote' in item.name:
                for entry in item.values:
                    self.DLA.values.append(entry ** 2)


        #solution time stamp
        self.solution_timestamps = self.archive.ArchiveItems[-2]

        #step count
        self.solution_step_count = self.archive.ArchiveItems[-3]

        #objective function - this is always the last entry
        self.objective_function = self.archive.ArchiveItems[-1]

        #compute total flight time
        self.total_flight_time = Archive.ArchiveItem('Total flight time')
        for i in range(0, len(self.archive.ArchiveItems[0].values) - 1):
            temp = 0.0
            for flight_time in self.flight_times:
                temp += flight_time.values[i]

            self.total_flight_time.values.append(temp)