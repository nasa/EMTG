# math_utilities.py

# Purpose
# Contains basic math utility functions in the class math_utilities

# Revision history
# Noble Hatten; 08/22/2018; Began
# Noble Hatten; 08/27/2018; Added unit vector getters
# Noble Hatten; 09/04/2018; Added angle_between_2_vectors

import numpy as np

class math_utilities_class:
    """Basic support math"""

    def i_unit(self):
        """Returns i unit columnn vector (3x1)"""
        v = np.array([[1.0], [0.0], [0.0]])
        return v

    def j_unit(self):
        """Returns j unit columnn vector (3x1)"""
        v = np.array([[0.0], [1.0], [0.0]])
        return v

    def k_unit(self):
        """Returns k unit columnn vector (3x1)"""
        v = np.array([[0.0], [0.0], [1.0]])
        return v

    def transform1(self, th):
        """Returns transformation matrix about 1-axis by angle th (rad)"""
        mat = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(th), np.sin(th)], [0.0, -np.sin(th), np.cos(th)]])
        return mat

    def transform2(self, th):
        """Returns transformation matrix about 2-axis by angle th (rad)"""
        mat = np.array([[np.cos(th), 0.0, -np.sin(th)], [0.0, 1.0, 0.0], [np.sin(th), 0.0, np.cos(th)]])
        return mat

    def transform3(self, th):
        """Returns transformation matrix about 3-axis by angle th (rad)"""
        mat = np.array([[np.cos(th), np.sin(th), 0.0], [-np.sin(th), np.cos(th), 0.0], [0.0, 0.0, 1.0]])
        return mat

    def atan2(self, y, x):
        """Calculate atan2 manually. Doing this so that we can use complex-step derivatives because np.arctan2 does not support complex arguments."""
        if x > 0.0:
            z = np.arctan(y/x)
        elif x < 0.0 and y >= 0.0:
            z = np.arctan(y/x) + np.pi
        elif x < 0.0 and y < 0.0:
            z = np.arctan(y/x) - np.pi
        elif x == 0.0 and y > 0.0:
            z = 0.5 * np.pi
        elif x == 0.0 and y < 0.0:
            z = -0.5 * np.pi
        else:
            z = 1.0e30 # undefined
        return z

    def d_atan2(self, y, x):
        """Calculate first derivatives of atan2 w.r.t. its arguments. First element is d/dy, second element is d/dx"""
        denom = x**2 + y**2
        if denom == 0.0:
            datan2_darg = np.array([[1.0e30], [1.0e30]]) # divide by zero (atan is undefined)
        datan2_darg = np.array([[x/denom], [-y/denom]])
        return datan2_darg

    def crossmat(self, x):
        """Calculate skew-symmetric cross matrix of 3x1 array x"""
        y = np.array([[0.0, -np.asscalar(x[2]), np.asscalar(x[1])], [np.asscalar(x[2]), 0.0, -np.asscalar(x[0])], [-np.asscalar(x[1]), np.asscalar(x[0]), 0.0]])
        return y

    def angle_between_2_vectors(self, x, y):
        """Calculate angle between 2 3D vectors in rad"""

        c = np.asscalar(np.dot(np.transpose(x), y)) # cosine is the dot product
        xcy = np.dot(self.crossmat(x), y) # sine is the norm of the cross product
        s = (xcy[0]*xcy[0] + xcy[1]*xcy[1] + xcy[2]*xcy[2])**0.5 # sine is the norm of the cross product

        if (s==0. and c==0.):
            angle = 0. # degenerate case
        else:
            angle = self.atan2(s, c)

        return angle

    def explicit_norm2(self, x):
        """Calculate 2 norm of vector explicitly rather than using numpy's built-in utilities"""
        dims = np.shape(x)
        dim = np.amax(dims)
        norm = 0.
        for i in range(dim):
            norm = norm + x[i]**2
        norm = norm**0.5
        return norm




