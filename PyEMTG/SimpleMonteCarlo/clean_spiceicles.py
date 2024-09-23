#command line utility to clean spicicles from an ephemeris file

def do_the_stuff(option_list):

    import ConOpsPeriod
    import SpiceyPy_Utilities as SpiceyUtil
    try:
        import spiceypy as spice
    except:
        print("spiceypy not available")
        
    working_directory = option_list[0]
    SPICE_ephem_directory = option_list[2]
    spice_handler = SpiceyUtil.SpiceHandler(SPICE_ephem_directory)
    spice_handler.loadSpiceFiles()
    
    my_ephemeris_reader = ConOpsPeriod.EphemerisFileReader()
    ephemeris_file_data = my_ephemeris_reader.parseEMTGephemerisFile(working_directory + '/' + option_list[1])

    my_ephemeris_reader.generateCleanEphemerisFileForBSP()
    
    spice_handler.unloadSpiceFiles()

if __name__ == '__main__':
    import sys
    import os

    if len(sys.argv) != 4:
        raise Exception("Unknown number of command line options!\n\
                         Syntax: python clean_spiceicles.py working_directory input_ephemeris_file SPICE_folder")
    
    thisdir = os.path.dirname(os.path.realpath(option_list[0]))
    sys.path.append(thisdir)
    sys.path.append(thisdir + '/../')
    sys.path.append(thisdir + '/../SpiceyPy_Utilities')
    
    do_the_stuff(sys.argv[1:])



