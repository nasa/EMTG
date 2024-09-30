# maneuverExecutionError.py
#
# Take in inertial control vector and variances on non-inertial control parameterization.
# output randomly perturbed inertial control vector.
#
# Revision history:
#
# 10/01/2018; Noble Hatten; initial revision

import numpy as np
from math import atan2, sin, cos

class maneuverExecutionError(object):
    def RotX(self,vec,angle):
        out_x = vec[0]
        out_y = cos(angle) * vec[1]  + sin(angle) * vec[2]
        out_z = -sin(angle) * vec[1] + cos(angle) * vec[2]
        return [out_x,out_y,out_z]
    
    def RotY(self,vec,angle):
        out_x = cos(angle) * vec[0] - sin(angle) * vec[2]
        out_y = vec[1]
        out_z = sin(angle) * vec[0] + cos(angle) * vec[2]
        return [out_x,out_y,out_z]

    def RotZ(self,vec,angle):
        out_x = cos(angle) * vec[0]  + sin(angle) * vec[1]
        out_y = -sin(angle) * vec[0] + cos(angle) * vec[1]
        out_z = vec[2]
        return [out_x,out_y,out_z]
    
    """Super-simple way of doing this stuff"""
    def __init__(self, uNomInertial, sigmaUMag, muAngle1, sigmaAngle1,muAngle2,sigmaAngle2):
        """Constructor"""
        # Just do everything here. You can then reorganize how you see fit
        # Random variables are assumed to be Gaussian
        # uNomInertial: 3-element numpy array of inertial nominal control vector
        # UMag is magnitude of control
        # Angle 1 is off-nominal control angle about body 1 axis
        # Angle 2 is off-nonimal control angle about body 2 axis
        # mu is mean (mu of uMag is assumed to be the magnitude of uNomInertial)
        # sigma is standard deviation (for uMag, taken to be as a fraction of actual control magnitude)
        # Returns:
        # self uRandomInertial: 3-element numpy array; control vector, perturbed from nominal control vector, in inertial frame

        uNomMag = np.linalg.norm(uNomInertial) # magnitude of nominal control
        uMag = uNomMag + sigmaUMag * np.random.randn() # realization of uMag
        angle_about_x = np.pi * np.random.randn() # realization of angle 1
        angle1 = muAngle1 + sigmaAngle1 * (np.random.randn()+1.0)/2.0 # realization of angle 2
        angle2 = muAngle2 + sigmaAngle2 * (np.random.randn()+1.0)/2.0 # realization of angle 2
        
        angle_from_u_nom_to_positive_x =  atan2(uNomInertial[1],uNomInertial[0])
        
        u_xz = self.RotZ(uNomInertial,angle_from_u_nom_to_positive_x)
        
        angle_from_u_xz_to_zplane = atan2(u_xz[2],u_xz[0])
        
        u_x = self.RotY(u_xz,-angle_from_u_xz_to_zplane)
        
        if abs(u_x[1]) > 1e-10 or u_x[2] > 1e-10:
            raise Exception("Rotation to x axis went wrong")
        
        new_u_x = self.RotZ(self.RotY([uMag,0,0],angle1),angle2)
        
        new_u_xz = self.RotY(new_u_x,angle_from_u_xz_to_zplane)
        
        self.uRandomInertial = self.RotZ(new_u_xz,-angle_from_u_nom_to_positive_x)

if __name__ == '__main__':
    import math
    # example driver
    # set values
    uInertial = [100,10,-10]
    
    uNomInertial = np.array([uInertial[0]/np.linalg.norm(uInertial), uInertial[1]/np.linalg.norm(uInertial), uInertial[2]/np.linalg.norm(uInertial)]) # nominal control vector in inertial frame
    sigmaUMag = 0.0023333333 # standard deviation of magnitude of control, expressed as fraction of nominal control magnitude
    muAngle1 = 0.0 # mean of angle1 (probably assume 0?)
    sigmaAngle1 = np.deg2rad(.059) # standard deviation of angle1 (assume angle1 and angle2 errors are evenly distributed)
    muAngle2 = 0.0 # mean of angle1 (probably assume 0?)
    sigmaAngle2 = np.deg2rad(.059) # standard deviation of angle1 (assume angle1 and angle2 errors are evenly distributed)

    print(sigmaAngle1)
    print(sigmaAngle2)

    # create the object
    
    error = []
    
    magError = []
    
    for i in range(100000):
        objTest = maneuverExecutionError(uNomInertial, sigmaUMag, muAngle1, sigmaAngle1, muAngle2, sigmaAngle2)
        
        newU = objTest.uRandomInertial
        
        magError 
        
        angle = math.acos(np.dot(newU,uNomInertial)/(np.linalg.norm(newU)*np.linalg.norm(uNomInertial)))

        error.append(angle)
        magError.append(np.linalg.norm(newU))
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    import matplotlib.pyplot as plt
    
    plt.figure(1)
    plt.hist(error,bins=100)
    
    plt.figure(2)
    plt.hist(magError,bins=100)
    plt.show()

    # check result
    # print('Original control:')
    # print(uNomInertial, np.linalg.norm(uNomInertial))
    # print('Perturbed control:')
    # print(objTest.uRandomInertial, np.linalg.norm(objTest.uRandomInertial))