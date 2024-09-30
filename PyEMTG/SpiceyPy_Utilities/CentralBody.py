class CentralBody(object):
    def __init__(self):
        self.name = ""
        self.spiceId = 0
        self.equatorialRadius = 1. # km
        self.mu = 1. # km3/s2
        self.flatteningCoef = 0. # unitless
        self.stateFrame = ""
        self.BplaneFrame = ""
        self.parentBody = None
        self.frameHandler = None
        return
        
    def setName(self, name):
        self.name = name
        import spiceypy
        self.spiceId = spiceypy.bodn2c(self.name.lower())
        return