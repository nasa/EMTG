# mathUtilities.py

# Revision history
# 2018-12-18; Noble Hatten; began (adapted from math_utilities_autograd.py)

#import numpy as np
import autograd.numpy as np # use the autograd version of np

class mathUtilities:
    """Basic support math"""
    def __init__(self):
        return

    def i_unit(self):
        """Returns i unit vector"""
        v = np.array([1.0, 0.0, 0.0])
        return v

    def j_unit(self):
        """Returns j unit vector"""
        v = np.array([0.0, 1.0, 0.0])
        return v

    def k_unit(self):
        """Returns k unit vector"""
        v = np.array([0.0, 0.0, 1.0])
        return v

    def transform1(self, th):
        """Returns transformation matrix about 1-axis by angle th (rad)"""
        mat = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(th), np.sin(th)], [0.0, -np.sin(th), np.cos(th)]])
        return mat

    def d_transform1(self, th):
        """Returns derivative of transformation matrix about 1-axis by angle th (rad) wrt th"""
        mat = np.array([[0.0, 0.0, 0.0], [0.0, -np.sin(th), np.cos(th)], [0.0, -np.cos(th), -np.sin(th)]])
        return mat

    def transform2(self, th):
        """Returns transformation matrix about 2-axis by angle th (rad)"""
        mat = np.array([[np.cos(th), 0.0, -np.sin(th)], [0.0, 1.0, 0.0], [np.sin(th), 0.0, np.cos(th)]])
        return mat

    def d_transform2(self, th):
        """Returns derivative of transformation matrix about 2-axis by angle th (rad) wrt th"""
        mat = np.array([[-np.sin(th), 0.0, -np.cos(th)], [0.0, 0.0, 0.0], [np.cos(th), 0.0, -np.sin(th)]])
        return mat

    def transform3(self, th):
        """Returns transformation matrix about 3-axis by angle th (rad)"""
        mat = np.array([[np.cos(th), np.sin(th), 0.0], [-np.sin(th), np.cos(th), 0.0], [0.0, 0.0, 1.0]])
        return mat

    def d_transform3(self, th):
        """Returns derivative of transformation matrix about 3-axis by angle th (rad) wrt th"""
        mat = np.array([[-np.sin(th), np.cos(th), 0.0], [-np.cos(th), -np.sin(th), 0.0], [0.0, 0.0, 0.0]])
        return mat

    def d_atan2(self, y, x):
        """Calculate first derivatives of atan2 w.r.t. its arguments. First element is d/dy, second element is d/dx"""
        denom = x**2 + y**2
        if denom == 0.0:
            datan2_darg = np.array([1.0e30, 1.0e30]) # divide by zero (atan is undefined)
        datan2_darg = np.array([x/denom, -y/denom])
        return datan2_darg

    def crossmat(self, x):
        """Calculate skew-symmetric cross matrix of 3 array x"""
        y = np.array([[0.0, -x[2], x[1]], [x[2], 0.0, -x[0]], [-x[1], x[0], 0.0]])
        return y

    def d_crossproduct(self, x, y):
        """
        derivative of cross product w.r.t. its arguments (x \times y)
        """
        dcp_dx = -self.crossmat(y)
        dcp_dy = self.crossmat(x)

        return dcp_dx, dcp_dy

    def angle_between_2_vectors(self, x, y):
        """Calculate angle between 2 3D vectors in rad"""

        c = np.dot(np.transpose(x), y) # cosine is the dot product
        xcy = np.dot(self.crossmat(x), y) # sine is the norm of the cross product
        s = self.column_vector_norm2(xcy) # sine is the norm of the cross product

        if (s==0. and c==0.):
            angle = 0. # degenerate case
        else:
            angle = np.arctan2(s, c)

        return angle

    def column_vector_norm2(self, x):
        """Calculate 2 norm of column vector explicitly rather than using numpy's built-in utilities"""
        dims = np.shape(x)
        dim = np.amax(dims)
        norm = 0.
        for i in range(dim):
            norm = norm + x[i]**2
        norm = norm**0.5
        return norm

    def column_vector_norm2_deriv(self, x):
        """Calculate derivative of norm2 of column vector w.r.t. the vector itself (x)"""
        normx = self.column_vector_norm2(x) # get the norm
        d_norm_d_x = (1.0/normx) * x
        return d_norm_d_x

    def tensor_bullet2_vector(self, t, v):
        """Calculate product of m x n x p tensor and n vector to produce m x p matrix"""
        dims = np.shape(v)
        n = np.amax(dims)
        dimsTensor = np.shape(t)
        m = dimsTensor[0]
        p = dimsTensor[2]

        mat = np.zeros((m,p))

        for i in range(n):
            mat[:,i] = np.dot(t[:,:,i], v[:])
        return mat

    def unit_vector_deriv(self, x):
        """Calculate derivative of unit vector x_unit w.r.t. its non-unit vector x"""
        dims = np.shape(x)
        n = np.amax(dims)
        xmag = self.column_vector_norm2(x)

        d_x_unit_d_x = (1./xmag) * (np.identity(n) - (1./xmag**2) * np.outer(x, x)) 
        return d_x_unit_d_x






