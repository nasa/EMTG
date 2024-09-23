#source file for NSGAIIFilter class
#Jacob Englander 7-20-2015

import os

class NSGAIIFilter:
    def __init__(self, filename):
        self.clear()
        self.read_input_file(filename)
    
    def clear(self):
        self.object_list = []
        self.filter_list = []
    
    def read_input_file(self, filename):
        #read the small body list file
        if os.path.isfile(filename):
            inputfile = open(filename, "r")
        else:
            print("File ", inputfile, " does not exist!")
            return
    
        for line in inputfile:
            if not line[0] == '#':
                self.object_list.append(line.strip('\n'))