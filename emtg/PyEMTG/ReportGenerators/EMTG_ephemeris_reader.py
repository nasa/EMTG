#class to ingest an EMTG ephemeris reader
#Jacob Englander 10-10-2018

class EMTG_ephemeris_reader(object):
    def __init__(self, filename = 'default.ephemeris', leapsecondspath = None):
        self.filename = filename
        self.__data__ = []

        self.parse(filename)

        if not leapsecondspath == None:
            self.computeTimeWidths(leapsecondspath)

    def parse(self, filename):
        import csv
        with open(filename, 'r') as file:
            reader = csv.DictReader(self.cleanfile(file))

            for line in reader:
                self.__data__.append(line)

    def cleanfile(self, file):
        for row in file:
            raw = row.replace('#','').replace(', ',',') #strip comment character and leading space

            if raw: yield raw

    def getRecord(self, recordIndex = 0):
        try:
            return self.__data__[recordIndex]
        except:
            print('invalid recordIndex, [' + str(recordIndex) + "]\n")

    def getValue(self, recordIndex = 0, key = 'epoch'):
        try:
            return self.__data__[recordIndex][key]
        except:
            print('invalid recordIndex or key, [' + str(recordIndex) + '][' + key + ']\n')

    def getLength(self):
        return len(self.__data__)

    def computeTimeWidths(self, leapsecondspath = None):
        #note we need SPICE to compute time widths because we have to convert from string ET to seconds-past-J2000
        if not leapsecondspath == None:
            try:
                import spiceypy
            except:
                print('spiceypy not available')
                return

            spiceypy.furnsh(leapsecondspath)
            self.__data__[-1]['timeWidth'] = 0.0
            for recordIndex in range(0, self.getLength() - 1):
                self.__data__[recordIndex]['timeWidth'] = spiceypy.str2et(self.getValue(recordIndex + 1, 'epoch')) - spiceypy.str2et(self.getValue(recordIndex, 'epoch'))

            spiceypy.unload(leapsecondspath)