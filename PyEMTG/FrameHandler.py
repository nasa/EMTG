import numpy as np


class FrameHandler(object):
    
    def __init__(self):
        
        self.alpha0 = 0.0
        self.alphadot = 0.0
        self.delta0 = 0.0
        self.deltadot = 0.0
        self.W0 = 0.0
        self.Wdot = 0.0
        
        self.alpha = 0.0
        self.delta = 0.0
        self.W = 0.0
        
        self.R_J2000BCI_to_ICRF = np.zeros([3,3])
        self.R_J2000BCI_to_J2000BCF = np.zeros([3,3])
        self.R_J2000BCI_to_BCI = np.zeros([3,3])
        self.R_J2000BCI_to_BCF = np.zeros([3,3])
        
        self.R_J2000BCF_to_J2000BCI = np.zeros([3,3])
        self.R_J2000BCF_to_ICRF = np.zeros([3,3])
        self.R_J2000BCF_to_BCI = np.zeros([3,3])
        self.R_J2000BCF_to_BCF = np.zeros([3,3])
        
        self.R_BCI_to_ICRF = np.zeros([3,3])
        self.R_BCI_to_BCF = np.zeros([3,3])
        self.R_BCI_to_J2000BCI = np.zeros([3,3])
        self.R_BCI_to_J2000BCF = np.zeros([3,3])
        
        self.R_BCF_to_BCI = np.zeros([3,3])
        self.R_BCF_to_ICRF = np.zeros([3,3])
        self.R_BCF_to_J2000BCI = np.zeros([3,3])
        self.R_BCF_to_J2000BCF = np.zeros([3,3])
        
        self.R_ICRF_to_J2000BCI = np.zeros([3,3])
        self.R_ICRF_to_J2000BCF = np.zeros([3,3])
        self.R_ICRF_to_BCI = np.zeros([3,3])
        self.R_ICRF_to_BCF = np.zeros([3,3])

        
    def getR(self, from_string, to_string):
        
        
        from_string = from_string.upper()
        to_string = to_string.upper()
        
        if from_string == to_string:
            return np.identity(3)
        
        if   from_string == "J2000BCI" and to_string == "ICRF":
            return self.R_J2000BCI_to_ICRF
        elif from_string == "J2000BCI" and to_string == "J2000BCF":
            return self.R_J2000BCI_to_J2000BCF
        elif from_string == "J2000BCI" and to_string == "BCI":
            return self.R_J2000BCI_to_BCI
        elif from_string == "J2000BCI" and to_string == "BCF":
            return self.R_J2000BCI_to_BCF
        
        elif from_string == "J2000BCF" and to_string == "J2000BCI":
            return self.R_J2000BCF_to_J2000BCI
        elif from_string == "J2000BCF" and to_string == "ICRF":
            return self.R_J2000BCF_to_ICRF
        elif from_string == "J2000BCF" and to_string == "BCI":
            return self.R_J2000BCF_to_BCI
        elif from_string == "J2000BCF" and to_string == "BCF":
            return self.R_J2000BCF_to_BCF
        
        elif from_string == "BCI" and to_string == "ICRF":
            return self.R_BCI_to_ICRF
        elif from_string == "BCI" and to_string == "BCF":
            return self.R_BCI_to_BCF
        elif from_string == "BCI" and to_string == "J2000BCI":
            return self.R_BCI_to_J2000BCI
        elif from_string == "BCI" and to_string == "J2000BCF":
            return self.R_BCI_to_J2000BCF
            
        elif from_string == "BCF" and to_string == "BCI":
            return self.R_BCF_to_BCI
        elif from_string == "BCF" and to_string == "ICRF":
            return self.R_BCF_to_ICRF
        elif from_string == "BCF" and to_string == "J2000BCI":
            return self.R_BCF_to_J2000BCI
        elif from_string == "BCF" and to_string == "J2000BCF":
            return self.R_BCF_to_J2000BCF
        
        elif from_string == "ICRF" and to_string == "J2000BCI":
            return self.R_ICRF_to_J2000BCI
        elif from_string == "ICRF" and to_string == "J2000BCF":
            return self.R_ICRF_to_J2000BCF
        elif from_string == "ICRF" and to_string == "BCI":
            return self.R_ICRF_to_BCI
        elif from_string == "ICRF" and to_string == "BCF":
            return self.R_ICRF_to_BCF    

        else:
            print("The conversion from " + from_string + " to " + to_string + " is not supported")
            return
            
    def rotateVector(self, vec, from_string, to_string):   
        return np.dot(self.getR(from_string, to_string), vec)
            
    def Rx(self, angle):
        return np.array([[1.0, 0.0          ,  0.0           ],
                            [0.0, np.cos(angle), -np.sin(angle)],
                            [0.0, np.sin(angle),  np.cos(angle) ]])
                        
    def Ry(self, angle):
        return np.array([[ np.cos(angle), 0.0, np.sin(angle)],
                            [ 0.0          , 1.0, 0.0          ],
                            [-np.sin(angle), 0.0, np.cos(angle)]])
                        
    def Rz(self, angle):
        return np.array([[np.cos(angle), -np.sin(angle), 0.0],
                            [np.sin(angle),  np.cos(angle), 0.0],
                            [0.0          ,  0.0          , 1.0]])
        
    def initJ2000frames(self, alpha0, alphadot, delta0, deltadot, W0, Wdot):
        self.alpha0 = alpha0 * np.pi / 180.0
        self.alphadot = alphadot * np.pi / 180.0
        self.delta0 = delta0 * np.pi / 180.0
        self.deltadot = deltadot * np.pi / 180.0
        self.W0 = W0 * np.pi / 180.0
        self.Wdot = Wdot * np.pi / 180.0
        self.R_J2000BCI_to_ICRF = np.dot(self.Rz(np.pi / 2.0 + self.alpha0), self.Rx(np.pi / 2.0 - self.delta0))
        self.R_ICRF_to_J2000BCI = np.transpose(self.R_J2000BCI_to_ICRF)
        
    def construct_rotation_matrices(self, ET_sec_past_J2000):
        days_since_reference_epoch = ET_sec_past_J2000 / 86400.0
        centuries_since_reference_epoch = days_since_reference_epoch / 36525.0
        
        self.alpha = self.alpha0 + self.alphadot * centuries_since_reference_epoch
        self.delta = self.delta0 + self.deltadot * centuries_since_reference_epoch
        self.W = self.W0 + self.Wdot * days_since_reference_epoch
        
        self.R_BCI_to_ICRF = np.dot(self.Rz(np.pi / 2.0 + self.alpha), self.Rx(np.pi / 2.0 - self.delta))
        self.R_ICRF_to_BCI = np.transpose(self.R_BCI_to_ICRF)
        self.R_BCF_to_BCI = self.Rz(self.W)
        self.R_BCI_to_BCF = np.transpose(self.R_BCF_to_BCI)
        
        self.R_J2000BCF_to_J2000BCI = self.Rz(self.W)
        self.R_J2000BCI_to_J2000BCF = np.transpose(self.R_J2000BCF_to_J2000BCI)
        
        self.R_BCF_to_ICRF = np.dot(self.R_BCI_to_ICRF, self.R_BCF_to_BCI)
        self.R_ICRF_to_BCF = np.transpose(self.R_BCF_to_ICRF)
        
        self.R_J2000BCF_to_ICRF = np.dot(self.R_J2000BCI_to_ICRF, self.R_J2000BCF_to_J2000BCI)
        self.R_ICRF_to_J2000BCF = np.transpose(self.R_J2000BCF_to_ICRF)
        
        self.R_J2000BCI_to_BCI = np.dot(self.R_ICRF_to_BCI, self.R_J2000BCI_to_ICRF)
        self.R_J2000BCI_to_BCF = np.dot(self.R_ICRF_to_BCF, self.R_J2000BCI_to_ICRF)
        self.R_J2000BCF_to_BCI = np.dot(self.R_ICRF_to_BCI, self.R_J2000BCF_to_ICRF)
        self.R_J2000BCF_to_BCF = np.dot(self.R_ICRF_to_BCF, self.R_J2000BCF_to_ICRF)
        self.R_BCI_to_J2000BCI = np.transpose(self.R_J2000BCI_to_BCI)
        self.R_BCI_to_J2000BCF = np.transpose(self.R_J2000BCF_to_BCI)
        self.R_BCF_to_J2000BCI = np.transpose(self.R_J2000BCI_to_BCF)
        self.R_BCF_to_J2000BCF = np.transpose(self.R_J2000BCF_to_BCF)


if __name__ == '__main__':

    import spiceypy as spice
    import StateConverter
    import posVel2BPlane
    import Universe

    DE_file            = "C:/emtg/missions/Mission1/universe/ephemeris_files/de430.bsp"
    leap_second_kernel = "C:/emtg/missions/Mission1/universe/ephemeris_files/naif0012.tls"
    PCK_kernel         = "C:/emtg/missions/Mission1/universe/ephemeris_files/pck00010.tpc"
    launch_open_bus_SPK_file =     "C:/C:/emtg/missions/Mission1/output/bus.bsp"
    spice.furnsh(DE_file)
    spice.furnsh(leap_second_kernel)
    spice.furnsh(PCK_kernel)
    spice.furnsh(launch_open_bus_SPK_file)

    venus_universe_file = "C:/emtg/missions/Mission1/universe/Venus.emtg_universe"
    venus_universe = Universe.Universe(venus_universe_file)

    epoch = "3-SEP-2027 14:37:27.347905159"
    state, light_time = spice.spkez(-123, spice.str2et(epoch), "J2000", 'NONE', 299)
    state = np.array([[state[0]],[state[1]],[state[2]],[state[3]],[state[4]],[state[5]]])

    state_converter = StateConverter.StateConverter()
    converter = posVel2BPlane.posVel2BPlane()

    rot_M = spice.tisbod("J2000", 299, spice.str2et(epoch))
    print(rot_M)

    
    frame_handler = FrameHandler()
    frame_handler.initJ2000frames(float(venus_universe.reference_angles[0]),
                                  float(venus_universe.reference_angles[1]),
                                  float(venus_universe.reference_angles[2]),
                                  float(venus_universe.reference_angles[3]),
                                  float(venus_universe.reference_angles[4]),
                                  float(venus_universe.reference_angles[5]))

    frame_handler.construct_rotation_matrices(spice.str2et(epoch))

    rot_M = frame_handler.getR("ICRF", "BCF")

    state_BCF = np.concatenate((np.dot(rot_M, state[0:3]), np.dot(rot_M, state[3:6])))
    state_BCF = state_BCF.flatten()
    
    BdotR = converter.bDotR(state_BCF[0:3], state_BCF[3:6], venus_universe.mu)
    BdotT = converter.bDotT(state_BCF[0:3], state_BCF[3:6], venus_universe.mu)
    Btheta = converter.bTheta(state_BCF[0:3], state_BCF[3:6], venus_universe.mu) * 180.0 / np.pi
    print(BdotR, BdotT, Btheta)

    Bplane_state = state_converter.CartesiantoIncomingBplane(state_BCF, venus_universe.mu)
    print(Bplane_state)

    rot_M = frame_handler.getR("ICRF", "BCI")

    state_BCI = np.concatenate((np.dot(rot_M, state[0:3]), np.dot(rot_M, state[3:6])))
    state_BCI = state_BCI.flatten()
    
    BdotR = converter.bDotR(state_BCI[0:3], state_BCI[3:6], venus_universe.mu)
    BdotT = converter.bDotT(state_BCI[0:3], state_BCI[3:6], venus_universe.mu)
    Btheta = converter.bTheta(state_BCI[0:3], state_BCI[3:6], venus_universe.mu) * 180.0 / np.pi
    print(BdotR, BdotT, Btheta)

    Bplane_state = state_converter.CartesiantoIncomingBplane(state_BCI, venus_universe.mu)
    print(Bplane_state)


    spice.unload(DE_file)
    spice.unload(leap_second_kernel)
    spice.unload(PCK_kernel)
    spice.unload(launch_open_bus_SPK_file)