# EMTG: Evolutionary Mission Trajectory Generator
# An open-source global optimization tool for preliminary mission design
# Provided by NASA Goddard Space Flight Center
#
# Copyright (c) 2013 - 2020 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Other Rights Reserved.
#
# Licensed under the NASA Open Source License (the "License"); 
# You may not use this file except in compliance with the License. 
# You may obtain a copy of the License at:
# https://opensource.org/licenses/NASA-1.3
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
# express or implied.   See the License for the specific language
# governing permissions and limitations under the License.

class Body(object):

    def __init__(self, linestring = [], J2flag = False):
        #class to contain all data for a body
        self.name = 'ChinchillaLand'
        self.shortname = 'Fuzzy'
        self.number = 0
        self.SPICE_ID = 10000000000
        self.minimum_flyby_altitude = -1
        self.mu = 1.0
        self.G = 6.674280000000000367e-20
        self.J2 = 0.0010826265
        self.AbsoluteMagnitude = 10.0
        self.albedo = 0.367
        self.radius = 1000.0
        self.ephemeris_epoch = 51544
        self.alpha0 = 0.0
        self.alphadot = 0.0
        self.delta0 = 0.0
        self.deltadot = 0.0
        self.W = 0.0
        self.Wdot = 0.0
        self.SMA = 1.5e+9
        self.ECC = 0.0
        self.INC = 1.0e-6
        self.RAAN = 1.0e-6
        self.AOP = 1.0e-6
        self.MA = 0.0
        if linestring != []:
            self.parse_input_line(linestring, J2flag)

        self.mass = self.mu / self.G

    def parse_input_line(self, linestring, J2flag):
        linecell = linestring.split(' ')
        if J2flag:
            self.name = linecell[0]
            self.shortname = linecell[1]
            self.number = eval(linecell[2])
            self.SPICE_ID = eval(linecell[3])
            self.minimum_flyby_altitude = float(linecell[4])
            self.mu = float(linecell[5])
            self.radius = float(linecell[6])
            self.J2 = float(linecell[7])
            self.AbsoluteMagnitude = float(linecell[8])
            self.albedo = float(linecell[9])
            self.ephemeris_epoch = float(linecell[10])
            self.alpha0 = float(linecell[11])
            self.alphadot = float(linecell[12])
            self.delta0 = float(linecell[13])
            self.deltadot = float(linecell[14])
            self.W = float(linecell[15])
            self.Wdot = float(linecell[16])
            self.SMA = float(linecell[17])
            self.ECC = float(linecell[18])
            self.INC = float(linecell[19])
            self.RAAN = float(linecell[20])
            self.AOP = float(linecell[21])
            self.MA = float(linecell[22])
        else:
            self.name = linecell[0]
            self.shortname = linecell[1]
            self.number = eval(linecell[2])
            self.SPICE_ID = eval(linecell[3])
            self.minimum_flyby_altitude = float(linecell[4])
            self.mu = float(linecell[5])
            self.radius = float(linecell[6])
            self.ephemeris_epoch = float(linecell[7])
            self.alpha0 = float(linecell[8])
            self.alphadot = float(linecell[9])
            self.delta0 = float(linecell[10])
            self.deltadot = float(linecell[11])
            self.W = float(linecell[12])
            self.Wdot = float(linecell[13])
            self.SMA = float(linecell[14])
            self.ECC = float(linecell[15])
            self.INC = float(linecell[16])
            self.RAAN = float(linecell[17])
            self.AOP = float(linecell[18])
            self.MA = float(linecell[19])

    def body_line(self):

        body_line = ""
        body_line += str(self.name) + " "
        body_line += str(self.shortname) + " "
        body_line += str(self.number) + " "
        body_line += str(self.SPICE_ID) + " "
        body_line += str(self.minimum_flyby_altitude) + " "
        body_line += str(self.mu) + " "
        body_line += str(self.radius) + " "
        body_line += str(self.J2) + " "
        body_line += str(self.AbsoluteMagnitude) + " "
        body_line += str(self.albedo) + " "
        body_line += str(self.ephemeris_epoch) + " "
        body_line += str(self.alpha0) + " "
        body_line += str(self.alphadot) + " "
        body_line += str(self.delta0) + " "
        body_line += str(self.deltadot) + " "
        body_line += str(self.W) + " "
        body_line += str(self.Wdot) + " "
        body_line += str(self.SMA) + " "
        body_line += str(self.ECC) + " "
        body_line += str(self.INC) + " "
        body_line += str(self.RAAN) + " "
        body_line += str(self.AOP) + " "
        body_line += str(self.MA)

        return body_line