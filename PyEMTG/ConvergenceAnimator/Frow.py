#container and plotter class for read an F row of an XF file
#Jacob Englander 7/19/2018

class Frow(object):
    def __init__(self, textstring):
        textcell = textstring.split(',')
        self.Findex = int(textcell[0])
        self.description = textcell[1]
        self.lowerbound = float(textcell[2])
        self.upperbound = float(textcell[3])
        self.value = float(textcell[4])

        self.violation = 0.0
        if self.value > self.upperbound:
            self.violation = self.value - self.upperbound
        elif self.value < self.lowerbound:
            self.violation = self.value - self.lowerbound

    def getViolation(self):
        return self.violation

    def getLogAbsViolation(self):
        from math import log10, fabs
        return log10(fabs(self.violation + 1.0e-13))

    def getLowerBound(self):
        return self.lowerbound

    def getUpperBound(self):
        return self.upperbound

    def getValue(self):
        return self.value

    def getIndex(self):
        return self.Findex