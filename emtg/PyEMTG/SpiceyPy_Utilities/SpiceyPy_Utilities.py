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

import spiceypy
try:
    import spiceypy as spice
    import spiceypy.utils.support_types as stypes
except:
    print("spiceypy not available")

class SpiceHandler(object):
    def __init__(self, ephem_directory_in):
        self.SPICE_ephem_directory = ephem_directory_in
        self.loadedSpiceFiles = []

    def loadSpiceFiles(self):
        from os import listdir
        from os.path import isfile, join

        # load SPICE ephemeris files
        spiceBinaryFiles = [f for f in listdir(self.SPICE_ephem_directory) if isfile(join(self.SPICE_ephem_directory, f))]
        for spiceFile in spiceBinaryFiles:
            try:
                if spiceFile.endswith('.bsp') or spiceFile.endswith('.tls') or spiceFile.endswith('.tpc'):
                    spice.furnsh(self.SPICE_ephem_directory + spiceFile)
                    self.loadedSpiceFiles.append(spiceFile)
                    print("loaded " + spiceFile)
            except:
                print("Could not furnish ", spiceFile)

    def unloadSpiceFiles(self):
        # unload spice files
        for spiceFile in self.loadedSpiceFiles:
            spice.unload(spiceFile)
        spice.kclear()

def greg2Greg(gregorian_date, time_system_string):
    return spice.timout(spice.str2et(gregorian_date), "YYYY Mon DD ::" + time_system_string + " HR:MN:SC.######")

def gregDateOffsetCalculator(gregorian_date, delta_days, time_system_string):

    # returns the time_system_string Gregorian date, offset from the input date string by delta_days
    # gregorian_date: SPICE compatible Gregorian date string
    # delta_days: number of days that you want to advance the Gregorian date string by
    # time_system_string: string representation of the output requested for the gregorian_date's time system (e.g. "TDB", "UTC" etc.)

    return spice.timout(spice.str2et(str(spice.str2et(gregorian_date) / 86400.0 + 2451545.0 + delta_days) + " JD TDB"), "YYYY Mon DD ::" + time_system_string + " HR:MN:SC.######")

def greg2Julian(gregorian_date):

    # returns the TDB Julian date corresponding to the input Gregorian date
    # gregorian_date: SPICE compatible Gregorian date string

    return spice.str2et(gregorian_date) / 86400.0 + 2451545.0

def greg2GMATJulian(gregorian_date):

    # returns the TDB Julian date corresponding to the input Gregorian date
    # gregorian_date: SPICE compatible Gregorian date string

    return greg2Julian(gregorian_date) - 2430000.0

def julian2Greg(julian_date):
    # returns the TDB Gregorian date
    # julian_date: TDB julian epoch
    # delta_days: number of days that you want to advance the Gregorian date string by
    # time_system_string: string representation of the input gregorian_date's time system (e.g. "TDB", "UTC" etc.)

    return spice.timout((julian_date - 2451545.0) * 86400.0, "YYYY Mon DD ::TDB HR:MN:SC.###### TDB")

def secSinceJ20002Greg(secSinceJ2000):
    # returns the TDB Gregorian date
    # secSinceJ2000: TDB julian seconds since J2000
    # delta_days: number of days that you want to advance the Gregorian date string by
    # time_system_string: string representation of the input gregorian_date's time system (e.g. "TDB", "UTC" etc.)

    return spice.timout(secSinceJ2000, "YYYY Mon DD ::TDB HR:MN:SC.###### TDB")