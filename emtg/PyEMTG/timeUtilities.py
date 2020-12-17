#time conversion utilities
#leverages the rest of PyEMTG and SpiceyPy

def stringToJD(dateString, universe_folder):
    import sys
    import inspect
    import os
    currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    sys.path.append(currentdir + "./" + 'SpiceyPy_Utilities')

    import SpiceyPy_Utilities

    JDate = 0
    try:
        JDate = eval(dateString)
    except:
        #can't eval a string, so it must not be a Julian date. Try converting from Gregorian
        import SpiceyPy_Utilities

        mySpiceHandler = SpiceyPy_Utilities.SpiceHandler(universe_folder + '/ephemeris_files/', only_tls=True)
        mySpiceHandler.loadSpiceFiles()

        JDate = SpiceyPy_Utilities.greg2Julian(dateString.replace('GMT','UTC'))

        mySpiceHandler.unloadSpiceFiles()

        
    #convert from GMAT Julian date to real Julian Date
    if JDate < 30000:
        JDate += 2430000

    #convert from JD to MJD if applicable
    if JDate > 2400000.5:
        JDate -= 2400000.5

    return JDate