#body distance report generator
#ingests .ephemeris files
#creates a file of distance from (your favorite place in the Spiceiverse) with the same number of rows as the .ephemeris
#Jacob Englander 10-18-2018

import SpiceyPy_Utilities as SpiceyUtil
import ConOpsPeriod as ConOps

class BodyDistanceReport(object):
    def __init__(self, ephemerisfilename = None, optionsfilename = None):
        self.optionsfilename = optionsfilename
        import MissionOptions
        self.myOptions = MissionOptions.MissionOptions(optionsfilename)

        # SPICE
        spice_handler = SpiceyUtil.SpiceHandler(self.myOptions.universe_folder + '/ephemeris_files/')
        spice_handler.loadSpiceFiles()

        self.ephemeris_reader = ConOps.EphemerisFileReader()        

        if not ephemerisfilename == None:
            self.ephemeris_file_data = self.ephemeris_reader.parseEMTGephemerisFile(ephemerisfilename)
        
        spice_handler.unloadSpiceFiles()


    def parseEphemerisFile(self, ephemerisfilename, optionsfilename):
        import EMTG_ephemeris_reader
        import sys
        sys.path.append('C:/Users/delliso2/NASA/EMTG')

        leapsecondspath = self.myOptions.universe_folder + '/ephemeris_files/' + self.myOptions.SPICE_leap_seconds_kernel

        self.myReader = EMTG_ephemeris_reader.EMTG_ephemeris_reader(ephemerisfilename, leapsecondspath)

    def printReport(self, SPICEbody = 399, reportfilename = 'default_body_distance_report.csv', units = 'AU'):
        #we need to load all of the SPICEness
        import spiceypy
        import os

        spiceypy.furnsh(self.myOptions.universe_folder + '/ephemeris_files/' + self.myOptions.SPICE_leap_seconds_kernel)
        spiceypy.furnsh(self.myOptions.universe_folder + '/ephemeris_files/' + self.myOptions.SPICE_reference_frame_kernel)

        for dirpath, dirnames, filenames in os.walk(self.myOptions.universe_folder + '/ephemeris_files/'):
            for file in filenames:
                sourcefile = os.path.join(dirpath, file)

                if sourcefile.endswith('.bsp'):
                    spiceypy.furnsh(sourcefile)


        with open(reportfilename, 'w') as file:
            #all operations
            file.write('Gregorian date (ET/TDB), Julian date (ET/TDB), Spacecraft distance from ' + spiceypy.bodc2n(SPICEbody) + ' (' + str(SPICEbody) + ') ' + ' (' + units + ')\n')

            for recordIndex in range(0, len(self.ephemeris_file_data)):
                myRecord = self.ephemeris_file_data[recordIndex]

                epoch = myRecord.julian_date

                seconds_since_J2000 = spiceypy.str2et(str(epoch) + " JD TDB")
                x_spacecraft = myRecord.spacecraft_position_x
                y_spacecraft = myRecord.spacecraft_position_y
                z_spacecraft = myRecord.spacecraft_position_z

                try:
                    [body_state, light_time] = spiceypy.spkez(SPICEbody, seconds_since_J2000, 'J2000', 'NONE', 10)
                except:
                    print ('I don\'t have enough information to tell you the position of ' + str(SPICEbody) + ' with respect to the Sun on ' + epoch + '\n')

                x_body = body_state[0]
                y_body = body_state[1]
                z_body = body_state[2]

                x_relative = x_spacecraft - x_body
                y_relative = y_spacecraft - y_body
                z_relative = z_spacecraft - z_body

                r_relative = (x_relative**2 + y_relative**2 + z_relative**2)**0.5

                if units == 'AU':
                    r_relative /= 149597870.691
                elif units != 'km':
                    raise('I don\'t know what a ' + units + ' is!')

                file.write(myRecord.gregorian_date + ',' + str(myRecord.julian_date) + ',' + str(r_relative) + '\n')

        #now we can unload all the SPICE kernels
        spiceypy.unload(self.myOptions.universe_folder + '/ephemeris_files/' + self.myOptions.SPICE_leap_seconds_kernel)
        spiceypy.unload(self.myOptions.universe_folder + '/ephemeris_files/' + self.myOptions.SPICE_reference_frame_kernel)

        for dirpath, dirnames, filenames in os.walk(self.myOptions.universe_folder + '/ephemeris_files/'):
            for file in filenames:
                sourcefile = os.path.join(dirpath, file)

                if sourcefile.endswith('.bsp'):
                    spiceypy.unload(sourcefile)



if __name__ == '__main__':

    ephemeris_filename = "C:/Discovery/MAGIC/high_fidelity/Launch_to_Callisto_SOI_HighFidelity_July_2025_AtlasV401/Launch_to_Callisto_SOI_HighFidelity_July_2025_AtlasV401_10172018_161254/Launch_to_Callisto_SOI_HighFidelity_July_2025_AtlasV401.ephemeris"
    options_filename = "C:/Discovery/MAGIC/high_fidelity/Launch_to_Callisto_SOI_HighFidelity_July_2025_AtlasV401/Launch_to_Callisto_SOI_HighFidelity_July_2025_AtlasV401_10172018_161254/Launch_to_Callisto_SOI_HighFidelity_July_2025_AtlasV401.emtgopt"
    myBodyDistanceReport = BodyDistanceReport(ephemeris_filename, options_filename)

    report_filename = "C:/Discovery/MAGIC/high_fidelity/Launch_to_Callisto_SOI_HighFidelity_July_2025_AtlasV401/Launch_to_Callisto_SOI_HighFidelity_July_2025_AtlasV401_10172018_161254/distance_from_the_Sun.csv"
    myBodyDistanceReport.printReport(10, report_filename, "AU")
