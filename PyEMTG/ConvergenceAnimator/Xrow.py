#container and plotter class for read an X row of an XF file
#Jacob Englander 7/19/2018

class Xrow(object):
    def __init__(self, textstring):
        textcell = textstring.split(',')
        self.Xindex = int(textcell[0])
        self.description = textcell[1]
        self.lowerbound = float(textcell[2])
        self.upperbound = float(textcell[3])
        self.ScaleFactor = float(textcell[4])
        self.value = float(textcell[5])

        self.violation = 0.0
        if self.value > self.upperbound:
            self.violation = self.value - self.upperbound
        elif self.value < self.lowerbound:
            self.violation = self.value - self.lowerbound

    def getViolation(self):
        return self.violation

    def getLogAbsViolation(self):
        from math import log10, fabs
        return log10(fabs(self.violation))

    def getLowerBound(self):
        return self.lowerbound

    def getUpperBound(self):
        return self.upperbound

    def getValue(self):
        return self.value

    def getIndex(self):
        return self.Xindex