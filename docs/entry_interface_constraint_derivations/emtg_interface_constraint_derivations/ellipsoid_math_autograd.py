# ellipsoid_math_autograd.py

# Purpose
# ellipsoid interface math use with autograd
# just use functions rather than classes to make it better for autograd

# Revision history
# Noble Hatten; 09/13/2018; Began


# standard imports
#import numpy as np
import autograd.numpy as np # use the autograd version of np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# custom imports
import math_utilities_autograd

class ellipsoid_math_autograd_class:
    def calc_ellipsmat(self, a, b, c):
        """Turn geometric values that define ellipse into matrix form"""
        a2 = a*a
        b2 = b*b
        c2 = c*c
        ellipsmat = np.diagflat([1./a2, 1./b2, 1./c2])
        return ellipsmat

    def calc_thetas_linear(self, t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad):
        """Calculate Euler angles relating PA frame and BCF frame as functions of time, assuming they all change linearly and t0=0"""
        theta1rad = theta1_0rad + theta1dot_0rad * t
        theta2rad = theta2_0rad + theta2dot_0rad * t
        theta3rad = theta3_0rad + theta3dot_0rad * t
        return theta1rad, theta2rad, theta3rad

    def calc_thetas_linear_derivs(self, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad):
        """Calculate derivatives of Euler angles relating PA frame and BCF frame as functions of time, assuming they all change linearly and t0=0"""

        # only the time derivatives are nonzero
        d_theta1rad_d_t = theta1dot_0rad
        d_theta2rad_d_t = theta2dot_0rad
        d_theta3rad_d_t = theta3dot_0rad
        return d_theta1rad_d_t, d_theta2rad_d_t, d_theta3rad_d_t

    def calc_w(self, t, w0rad, w0dotrad):
        """Calculate angle w between x axes of BCI and BCF frames (rad)"""
        wrad = w0rad + w0dotrad * t
        return wrad

    def calc_omega(self, w0dotrad):
        """Calculate angular velocity of BCF frame w.r.t. BCI frame (rad/s)"""
        omega = np.array([[0.0], [0.0], [w0dotrad]])
        return omega
    
    def calc_R_bcf2pa(self, theta1rad, theta2rad, theta3rad):
        """Calculate transformation matrix s.t. r_PA = R_bcf2pa * r_BCF"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        R_bcf2fprime = m_u.transform3(theta1rad)
        R_fprime2fdprime = m_u.transform1(theta2rad)
        R_fdprime2pa = m_u.transform3(theta3rad)
        R_bcf2pa = np.dot(R_fdprime2pa, np.dot(R_fprime2fdprime, R_bcf2fprime))
        R_pa2bcf = np.transpose(R_bcf2pa)
        return R_bcf2pa

    def calc_d_R_bcf2pa_d_t(self, theta1rad, theta2rad, theta3rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad):
        """Calculate time derivatives of transformation matrix s.t. r_PA = R_bcf2pa * r_BCF, assuming linear rates of change of the thetas"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()
        R_bcf2fprime = m_u.transform3(theta1rad)
        R_fprime2fdprime = m_u.transform1(theta2rad)
        R_fdprime2pa = m_u.transform3(theta3rad)
        d_R_d_theta1 = np.dot(R_fdprime2pa, np.dot(R_fprime2fdprime, m_u.d_transform3(theta1rad)))
        d_R_d_theta2 = np.dot(R_fdprime2pa, np.dot(m_u.d_transform1(theta2rad), R_bcf2fprime))
        d_R_d_theta3 = np.dot(m_u.d_transform3(theta3rad), np.dot(R_fprime2fdprime, R_bcf2fprime))
        d_R_bcf2pa_d_t = d_R_d_theta1 * theta1dot_0rad + d_R_d_theta2 * theta2dot_0rad + d_R_d_theta3 * theta3dot_0rad
        return d_R_bcf2pa_d_t

    def calc_rbcf(self, t, r_pa, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad):
        """Calculate position vector in BCF frame given position vector in PA frame"""
        rbcf = np.zeros((3,1))
        theta1rad, theta2rad, theta3rad = self.calc_thetas_linear(t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad)
        R_bcf2pa = self.calc_R_bcf2pa(theta1rad, theta2rad, theta3rad)
        rbcf = np.dot(np.transpose(R_bcf2pa), r_pa)
        return rbcf

    def calc_vbcf_mag(self, vbcf):
        """Calculate magnitude of velocity wrt BCF frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()
        vbcf_mag = m_u.column_vector_norm2(vbcf)
        return vbcf_mag

    def calc_vbcf_mag_derivs(self, vbcf):
        """Calculate derivatives of velocity wrt BCF frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()
        d_vbcf_mag_d_vbcf = m_u.column_vector_norm2_deriv(vbcf)

        # others zero
        d_vbcf_mag_d_t = 0.0
        d_vbcf_mag_d_rbcf = np.zeros((1,3))
        return d_vbcf_mag_d_t, d_vbcf_mag_d_rbcf, d_vbcf_mag_d_vbcf

    def calc_vbci_in_bcf(self, rbcf, vbcf, w0dotrad):
        """Calculate velocity wrt BCI frame expressed in BCF frame"""
        omega = self.calc_omega(w0dotrad)
        vbci_in_bcf = vbcf + np.cross(omega, rbcf, axis=0)
        return vbci_in_bcf

    def calc_vbci_in_bcf_derivs(self, w0dotrad):
        """Calculate derivatives of velocity wrt BCI frame expressed in BCF frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()
        omega = self.calc_omega(w0dotrad)
        d_vbci_in_bcf_d_t = np.zeros((3,1)) # this assumes that d [omega] / d [t] is zero!
        d_vbci_in_bcf_d_rbcf = m_u.crossmat(omega)
        d_vbci_in_bcf_d_vbcf = np.identity(3)
        return d_vbci_in_bcf_d_t, d_vbci_in_bcf_d_rbcf, d_vbci_in_bcf_d_vbcf

    def calc_vbci_mag(self, rbcf, vbcf, w0dotrad):
        """Calculate magnitude of velocity wrt BCI frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()
        vbci_in_bcf = self.calc_vbci_in_bcf(rbcf, vbcf, w0dotrad)
        vbci_mag = m_u.column_vector_norm2(vbci_in_bcf)
        return vbci_mag

    def calc_vbci_mag_derivs(self, rbcf, vbcf, w0dotrad):
        """Calculate derivatives of velocity magnitude wrt BCI frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()
        d_vbci_in_bcf_d_t, d_vbci_in_bcf_d_rbcf, d_vbci_in_bcf_d_vbcf = self.calc_vbci_in_bcf_derivs(w0dotrad) # derivatives of the vector vbci_in_bcf
        vbci_in_bcf = self.calc_vbci_in_bcf(rbcf, vbcf, w0dotrad) # vbci_in_bcf
        d_vbci_mag_d_vbci_in_bcf = m_u.column_vector_norm2_deriv(vbci_in_bcf) # derivative of scalar wrt vector
        d_vbci_mag_d_t = np.zeros((3,1)) # assumes that d [omega] / [t] is zero! 
        d_vbci_mag_d_rbcf = np.dot(d_vbci_mag_d_vbci_in_bcf, d_vbci_in_bcf_d_rbcf)
        d_vbci_mag_d_vbcf = np.dot(d_vbci_mag_d_vbci_in_bcf, d_vbci_in_bcf_d_vbcf)

        return d_vbci_mag_d_t, d_vbci_mag_d_rbcf, d_vbci_mag_d_vbcf

    def calc_longitude_bodycentric(self, rbcf):
        """Calculate bodycentric longitude in radians"""
        longitude_bodycentric = np.arctan2(rbcf[1,0], rbcf[0,0])
        return longitude_bodycentric

    def calc_longitude_bodycentric_derivs(self, rbcf):
        """Calculate bodycentric longitude derivatives"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # position derivatives
        d_atan2 = m_u.d_atan2(rbcf[1,0], rbcf[0,0])
        d_lon_d_rbcf = np.array([d_atan2[1], d_atan2[0], 0.0])

        # only depends on position, so:
        d_lon_d_t = 0.0 # d/d t
        d_lon_d_vbcf = np.zeros((1,3)) # d/d vbcf
        return d_lon_d_t, d_lon_d_rbcf, d_lon_d_vbcf

    def calc_longitude_bodydetic(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate bodydetic longitude in radians"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # calculate normal vector of project of ellipsoid onto BCF xy plane:
        longitude_bodycentric = self.calc_longitude_bodycentric(rbcf) # need bodycentric longitude

        rproj_bcf_mag = m_u.column_vector_norm2(rbcf[0:2]) # magnitude of projection into xy plane
        rproj_bcf = rproj_bcf_mag * np.array([[np.cos(longitude_bodycentric)], [np.sin(longitude_bodycentric)], [0.]]) # magnitude + direction
        theta1rad, theta2rad, theta3rad = self.calc_thetas_linear(t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad) # get the thetas
        R_bcf2pa = self.calc_R_bcf2pa(theta1rad, theta2rad, theta3rad) # convert thetas to rotation matrix
        R_pa2bcf = np.transpose(R_bcf2pa)
        rproj_pa = np.dot(R_bcf2pa, rproj_bcf) # transform to PA frame
        ellipsmat = self.calc_ellipsmat(a, b, c)
        nproj_pa = 2. * np.dot(ellipsmat, rproj_pa) # get normal vector of rproj in PA frame
        nproj_bcf = np.dot(R_pa2bcf, nproj_pa) # transform to BCF frame
        arg1 = nproj_bcf
        i_bcf = m_u.i_unit()
        j_bcf = m_u.j_unit()
        arg2 = i_bcf
        temp = m_u.angle_between_2_vectors(arg1, arg2) # bodydetic longitude is angle between nproj and i_bcf

        # quadrant check
        if (np.dot(np.transpose(nproj_bcf), j_bcf) >= 0.):
            longitude_bodydetic = temp
        else:
            longitude_bodydetic = 2.*np.pi - temp

        return longitude_bodydetic

    def calc_longitude_bodydetic_derivs(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of bodydetic longitude in radians"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()
        i_bcf = m_u.i_unit()
        j_bcf = m_u.j_unit()

        longitude_bodycentric = self.calc_longitude_bodycentric(rbcf) # need bodycentric longitude
        rproj_bcf_mag = m_u.column_vector_norm2(rbcf[0:2]) # magnitude of projection into xy plane
        rproj_bcf = rproj_bcf_mag * np.array([[np.cos(longitude_bodycentric)], [np.sin(longitude_bodycentric)], [0.]]) # magnitude + direction
        theta1rad, theta2rad, theta3rad = self.calc_thetas_linear(t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad) # get the thetas
        R_bcf2pa = self.calc_R_bcf2pa(theta1rad, theta2rad, theta3rad) # convert thetas to rotation matrix
        R_pa2bcf = np.transpose(R_bcf2pa)
        rproj_pa = np.dot(R_bcf2pa, rproj_bcf) # transform to PA frame
        ellipsmat = self.calc_ellipsmat(a, b, c)
        nproj_pa = 2. * np.dot(ellipsmat, rproj_pa) # get normal vector of rproj in PA frame
        nproj_bcf = np.dot(R_pa2bcf, nproj_pa) # transform to BCF frame

        nproj_dot_i = np.dot(np.transpose(nproj_bcf[:,0]), i_bcf[:,0]) # dot product
        nproj_bcf_cross = m_u.crossmat(nproj_bcf) # cross matrix
        icross = m_u.crossmat(i_bcf) # cross matrix
        nproj_cross_i = np.dot(nproj_bcf_cross, i_bcf) # cross product
        nproj_cross_i_mag = m_u.column_vector_norm2(nproj_cross_i) # mag of cross product

        d_lon_d_cross = nproj_dot_i / (nproj_cross_i_mag**2 + nproj_dot_i**2) # derivative w.r.t. nproj cross i mag
        d_lon_d_dot = (-nproj_cross_i_mag) / (nproj_cross_i_mag**2 + nproj_dot_i**2) # derivative w.r.t. nproj dot i mag

        # time
        d_R_bcf2pa_d_t = self.calc_d_R_bcf2pa_d_t(theta1rad, theta2rad, theta3rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad)
        d_R_pa2bcf_d_t = np.transpose(d_R_bcf2pa_d_t)
        d_rproj_bcf_d_t = np.zeros((3,1))
        d_nproj_d_t = 2. * (np.dot(d_R_pa2bcf_d_t, np.dot(ellipsmat, np.dot(R_bcf2pa, rproj_bcf))) +
                            np.dot(R_pa2bcf, np.dot(ellipsmat, np.dot(d_R_bcf2pa_d_t, rproj_bcf))))
        d_cross_d_t = -(1./nproj_cross_i_mag) * np.dot(np.transpose(nproj_cross_i), np.dot(icross, d_nproj_d_t))
        d_dot_d_t = np.dot(np.transpose(i_bcf), d_nproj_d_t)
        d_lon_d_t = d_lon_d_cross * d_cross_d_t + d_lon_d_dot * d_dot_d_t

        # position
        d_lon_centric_d_t, d_lon_centric_d_rbcf, d_lon_centric_d_vbcf = self.calc_longitude_bodycentric_derivs(rbcf)
        d_rproj_bcf_d_rbcf = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
        d_nproj_d_rbcf = 2. * np.dot(R_pa2bcf, np.dot(ellipsmat, np.dot(R_bcf2pa, d_rproj_bcf_d_rbcf)))
        d_cross_d_rbcf = -(1./nproj_cross_i_mag) * np.dot(np.transpose(nproj_cross_i), np.dot(icross, d_nproj_d_rbcf))
        d_dot_d_rbcf = np.dot(np.transpose(i_bcf), d_nproj_d_rbcf)
        d_lon_d_rbcf = (np.dot(d_lon_d_cross, d_cross_d_rbcf) + np.dot(d_lon_d_dot, d_dot_d_rbcf))

        # quadrant check
        if (np.dot(np.transpose(nproj_bcf), j_bcf) < 0.):
            d_lon_d_t = -d_lon_d_t
            d_lon_d_rbcf = -d_lon_d_rbcf

        # no velocity dependence
        d_lon_d_vbcf = np.zeros((1,3)) # d/d vbcf
        return d_lon_d_t, d_lon_d_rbcf, d_lon_d_vbcf

    def calc_latitude_bodycentric(self, rbcf):
        """Calculate bodycentric latitude in radians"""
        longitude_bodycentric = self.calc_longitude_bodycentric(rbcf)
        latitude_bodycentric = np.arctan2(rbcf[2,0], (rbcf[0,0]*np.cos(longitude_bodycentric) + rbcf[1,0]*np.sin(longitude_bodycentric)))
        return latitude_bodycentric

    def calc_latitude_bodycentric_derivs(self, rbcf):
        """Calculate bodycentric latitude derivatives"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # position derivative
        longitude_bodycentric = self.calc_longitude_bodycentric(rbcf) # bodycentric longitude
        d_lon_d_t, d_lon_d_rbcf, d_lon_d_vbcf = self.calc_longitude_bodycentric_derivs(rbcf) # bodycentric longitude derivatives

        d_atan2 = m_u.d_atan2(rbcf[2,0], rbcf[0,0]*np.cos(longitude_bodycentric) + rbcf[1,0]*np.sin(longitude_bodycentric)) # derivatives of atan2
        d_coslon_d_rbcf = -np.sin(longitude_bodycentric) *d_lon_d_rbcf
        d_sinlon_d_rbcf = np.cos(longitude_bodycentric) * d_lon_d_rbcf
        d_lat_x_d_rbcf = np.array([np.cos(longitude_bodycentric) + rbcf[0,0] * d_coslon_d_rbcf[0] + rbcf[1,0] * d_sinlon_d_rbcf[0],
                                   rbcf[0,0] * d_coslon_d_rbcf[1] + np.sin(longitude_bodycentric) + rbcf[1,0] * d_sinlon_d_rbcf[1],
                                   0.0])
        d_lat_y_d_rbcf = np.array([0.0, 0.0, 1.0])
        d_lat_d_rbcf = d_lat_y_d_rbcf * d_atan2[0] + d_lat_x_d_rbcf * d_atan2[1] # final chain rule

        # only depends on position, so:
        d_lat_d_t = 0.0 # d/d t
        d_lat_d_vbcf = np.zeros((1,3)) # d/d vbcf
        return d_lat_d_t, d_lat_d_rbcf, d_lat_d_vbcf

    def calc_latitude_bodydetic(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate bodydetic latitude in rad"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()
        arg1 = m_u.k_unit()
        u_bcf = self.calc_u_bcf(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # outward normal in bcf
        temp = m_u.angle_between_2_vectors(arg1, u_bcf)
        latitude_bodydetic = np.pi*0.5 - temp
        return latitude_bodydetic

    def calc_latitude_bodydetic_derivs(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of bodydetic latitude in rad"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()
        k_bcf = m_u.k_unit()
        kcross = m_u.crossmat(k_bcf)
        u_bcf = self.calc_u_bcf(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        k_dot_u = np.dot(np.transpose(k_bcf[:,0]), u_bcf[:,0]) # dot product
        k_cross_u = np.dot(kcross, u_bcf) # cross product
        k_cross_u_mag = m_u.column_vector_norm2(k_cross_u)
        d_lat_d_cross = k_dot_u / (k_cross_u_mag**2 + k_dot_u**2) # derivative w.r.t. k cross u_bcf
        d_lat_d_dot = (-k_cross_u_mag) / (k_cross_u_mag**2 + k_dot_u**2) # derivative w.r.t. k dot u_bcf mag
        d_u_bcf_d_t, d_u_bcf_d_rbcf, d_u_bcf_d_vbcf = self.calc_u_bcf_derivs(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # derivatives of u_bcf

        # time
        d_cross_d_t = (1./k_cross_u_mag) * np.dot(np.transpose(k_cross_u), np.dot(kcross, d_u_bcf_d_t))
        d_dot_d_t = np.dot(np.transpose(k_bcf), d_u_bcf_d_t)
        d_lat_d_t = -(d_lat_d_cross * d_cross_d_t + d_lat_d_dot * d_dot_d_t)

        # position
        d_cross_d_rbcf = (1./k_cross_u_mag) * np.dot(np.transpose(k_cross_u), np.dot(kcross, d_u_bcf_d_rbcf))
        d_dot_d_rbcf = np.dot(np.transpose(k_bcf), d_u_bcf_d_rbcf)
        d_lat_d_rbcf = -(np.dot(d_lat_d_cross, d_cross_d_rbcf) + np.dot(d_lat_d_dot, d_dot_d_rbcf))

        # no velocity dependence
        d_lat_d_vbcf = np.zeros((1,3)) # d/d vbcf
        return d_lat_d_t, d_lat_d_rbcf, d_lat_d_vbcf

    def calc_n_pa(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calcuate outward normal in PA coordinates"""
        theta1rad, theta2rad, theta3rad = self.calc_thetas_linear(t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad) # euler angles
        R_bcf2pa = self.calc_R_bcf2pa(theta1rad, theta2rad, theta3rad) # transformation matrix
        r_pa = np.dot(R_bcf2pa, rbcf) # position in pa
        ellipsmat = self.calc_ellipsmat(a, b, c)
        n_pa = 2.0 * np.dot(ellipsmat, r_pa) # outward normal in pa
        return n_pa

    def calc_u_bcf(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate outward normal vector in bcf coordinates"""
        theta1rad, theta2rad, theta3rad = self.calc_thetas_linear(t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad) # euler angles
        R_bcf2pa = self.calc_R_bcf2pa(theta1rad, theta2rad, theta3rad) # transformation matrix
        n_pa = self.calc_n_pa(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # outward normal in pa
        R_pa2bcf = np.transpose(R_bcf2pa) # transformation matrix
        u_bcf = np.dot(R_pa2bcf, n_pa) # outward normal in bcf
        return u_bcf

    def calc_u_bcf_derivs(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of outward normal vector in bcf coordinates"""

        # time derivative
        ellipsmat = self.calc_ellipsmat(a, b, c)
        theta1rad, theta2rad, theta3rad = self.calc_thetas_linear(t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad)
        R_bcf2pa = self.calc_R_bcf2pa(theta1rad, theta2rad, theta3rad)
        R_pa2bcf = np.transpose(R_bcf2pa)
        d_R_bcf2pa_d_t = self.calc_d_R_bcf2pa_d_t(theta1rad, theta2rad, theta3rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad) # d R^{BCF->PA} / d t
        d_R_pa2bcf_d_t = np.transpose(d_R_bcf2pa_d_t) # d R^{PA->BCF} / d t
        n_pa = self.calc_n_pa(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # outward normal in pa
        d_u_bcf_d_t = np.dot(d_R_pa2bcf_d_t, n_pa) + 2. * np.dot(R_pa2bcf, np.dot(ellipsmat, np.dot(d_R_bcf2pa_d_t, rbcf))) # d/dt

        # rbcf derivative
        d_n_pa_d_r_pa = 2.0 * ellipsmat
        d_u_bcf_d_r_pa = np.dot(R_pa2bcf, d_n_pa_d_r_pa)
        d_r_pa_d_rbcf = R_bcf2pa
        d_u_bcf_d_rbcf = np.dot(d_u_bcf_d_r_pa, d_r_pa_d_rbcf)

        # no velocity dependence
        d_u_bcf_d_vbcf = np.zeros((3,3)) # d/d vbcf
        return d_u_bcf_d_t, d_u_bcf_d_rbcf, d_u_bcf_d_vbcf

    def calc_heading_bcf_topocentric(self, t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate heading angle using velocity w.r.t. BCF frame and topocentric SEU coordinates. Due south in topocentric SEU is zero heading, with positive toward east."""
        R_bcf2seu_topocentric = self.calc_R_bcf2seu_topocentric(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        vbcf_expressed_in_seu_topocentric = np.dot(R_bcf2seu_topocentric, vbcf)
        heading_bcf_topocentric = np.arctan2(vbcf_expressed_in_seu_topocentric[1,0], vbcf_expressed_in_seu_topocentric[0,0])
        return heading_bcf_topocentric

    def calc_heading_bcf_topocentric_derivs(self, t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of heading angle using velocity w.r.t. BCF frame and topocentric SEU coordinates. Due south in topocentric SEU is zero heading, with positive toward east."""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        R_bcf2seu_topocentric = self.calc_R_bcf2seu_topocentric(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        d_R_bcf2seu_d_t, d_R_bcf2seu_d_rbcf, d_R_bcf2seu_d_vbcf = self.calc_R_bcf2seu_topocentric_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        vbcf_expressed_in_seu_topocentric = np.dot(R_bcf2seu_topocentric, vbcf)
        d_head_d_vx_topo = (-vbcf_expressed_in_seu_topocentric[1,0]) / (vbcf_expressed_in_seu_topocentric[0,0]**2 + vbcf_expressed_in_seu_topocentric[1,0]**2)
        d_head_d_vy_topo = (vbcf_expressed_in_seu_topocentric[0,0]) / (vbcf_expressed_in_seu_topocentric[0,0]**2 + vbcf_expressed_in_seu_topocentric[1,0]**2)

        # time
        d_vbcf_expressed_in_seu_topocentric_d_t = np.dot(d_R_bcf2seu_d_t, vbcf)
        d_vx_topo_d_t = d_vbcf_expressed_in_seu_topocentric_d_t[0,0]
        d_vy_topo_d_t = d_vbcf_expressed_in_seu_topocentric_d_t[1,0]
        d_head_d_t = d_head_d_vy_topo * d_vy_topo_d_t + d_head_d_vx_topo * d_vx_topo_d_t # final chain rule for d/d t

        # position
        d_vbcf_expressed_in_seu_topocentric_d_rbcf = m_u.tensor_bullet2_vector(d_R_bcf2seu_d_rbcf, vbcf)
        d_vx_topo_d_rbcf = d_vbcf_expressed_in_seu_topocentric_d_rbcf[0,:]
        d_vy_topo_d_rbcf = d_vbcf_expressed_in_seu_topocentric_d_rbcf[1,:]
        d_head_d_rbcf = d_head_d_vy_topo * d_vy_topo_d_rbcf + d_head_d_vx_topo * d_vx_topo_d_rbcf # final chain rule for d/d rbcf


        # velocity
        d_vbcf_expressed_in_seu_topocentric_d_vbcf = m_u.tensor_bullet2_vector(d_R_bcf2seu_d_vbcf, vbcf) + R_bcf2seu_topocentric
        d_vx_topo_d_vbcf = d_vbcf_expressed_in_seu_topocentric_d_vbcf[0,:]
        d_vy_topo_d_vbcf = d_vbcf_expressed_in_seu_topocentric_d_vbcf[1,:]
        d_head_d_vbcf = d_head_d_vy_topo * d_vy_topo_d_vbcf + d_head_d_vx_topo * d_vx_topo_d_vbcf # final chain rule for d/d vbcf
        
        return d_head_d_t, d_head_d_rbcf, d_head_d_vbcf

    def calc_heading_bci_topocentric(self, t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate heading angle using velocity w.r.t. BCI frame and topocentric SEU coordinates. Due south in topocentric SEU is zero heading, with positive toward east."""
        # convert vbcf to vbci_expressed_in_bcf
        vbci_expressed_in_bcf = self.calc_vbci_in_bcf(rbcf, vbcf, w0dotrad)

        # express in seu topocentric
        R_bcf2seu_topocentric = self.calc_R_bcf2seu_topocentric(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        vbci_expressed_in_seu_topocentric = np.dot(R_bcf2seu_topocentric, vbci_expressed_in_bcf)
        heading_bci_topocentric = np.arctan2(vbci_expressed_in_seu_topocentric[1,0], vbci_expressed_in_seu_topocentric[0,0])
        return heading_bci_topocentric

    def calc_heading_bci_topocentric_derivs(self, t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of heading angle using velocity w.r.t. BCI frame and topocentric SEU coordinates. Due south in topocentric SEU is zero heading, with positive toward east."""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # BCI velocity
        vbci_in_bcf = self.calc_vbci_in_bcf(rbcf, vbcf, w0dotrad)

        R_bcf2seu_topocentric = self.calc_R_bcf2seu_topocentric(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        d_R_bcf2seu_d_t, d_R_bcf2seu_d_rbcf, d_R_bcf2seu_d_vbcf = self.calc_R_bcf2seu_topocentric_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        vbci_expressed_in_seu_topocentric = np.dot(R_bcf2seu_topocentric, vbci_in_bcf) # only one transformation because vbci is expressed in bcf to start with
        d_head_d_vx_topo = (-vbci_expressed_in_seu_topocentric[1,0]) / (vbci_expressed_in_seu_topocentric[0,0]**2 + vbci_expressed_in_seu_topocentric[1,0]**2)
        d_head_d_vy_topo = (vbci_expressed_in_seu_topocentric[0,0]) / (vbci_expressed_in_seu_topocentric[0,0]**2 + vbci_expressed_in_seu_topocentric[1,0]**2)

        # time
        d_vbci_expressed_in_seu_topocentric_d_t = np.dot(d_R_bcf2seu_d_t, vbci_in_bcf)
        d_vx_topo_d_t = d_vbci_expressed_in_seu_topocentric_d_t[0,0]
        d_vy_topo_d_t = d_vbci_expressed_in_seu_topocentric_d_t[1,0]
        d_head_d_t = d_head_d_vy_topo * d_vy_topo_d_t + d_head_d_vx_topo * d_vx_topo_d_t # final chain rule for d/d t

        # position
        omega = self.calc_omega(w0dotrad)
        omega_cross = m_u.crossmat(omega)
        d_vbci_expressed_in_seu_topocentric_d_rbcf = m_u.tensor_bullet2_vector(d_R_bcf2seu_d_rbcf, vbci_in_bcf) + np.dot(R_bcf2seu_topocentric, omega_cross)
        d_vx_topo_d_rbcf = d_vbci_expressed_in_seu_topocentric_d_rbcf[0,:]
        d_vy_topo_d_rbcf = d_vbci_expressed_in_seu_topocentric_d_rbcf[1,:]
        d_head_d_rbcf = d_head_d_vy_topo * d_vy_topo_d_rbcf + d_head_d_vx_topo * d_vx_topo_d_rbcf # final chain rule for d/d rbcf


        # velocity
        d_vbci_expressed_in_seu_topocentric_d_vbcf = m_u.tensor_bullet2_vector(d_R_bcf2seu_d_vbcf, vbci_in_bcf) + R_bcf2seu_topocentric
        d_vx_topo_d_vbcf = d_vbci_expressed_in_seu_topocentric_d_vbcf[0,:]
        d_vy_topo_d_vbcf = d_vbci_expressed_in_seu_topocentric_d_vbcf[1,:]
        d_head_d_vbcf = d_head_d_vy_topo * d_vy_topo_d_vbcf + d_head_d_vx_topo * d_vx_topo_d_vbcf # final chain rule for d/d vbcf
        
        return d_head_d_t, d_head_d_rbcf, d_head_d_vbcf

    def calc_heading_bcf_polar(self, t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate heading angle using velocity w.r.t. BCF frame and polar SEU coordinates. Due south in polar SEU is zero heading, with positive toward east."""
        R_bcf2seu_polar = self.calc_R_bcf2seu_polar(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        vbcf_expressed_in_seu_polar = np.dot(R_bcf2seu_polar, vbcf)
        heading_bcf_polar = np.arctan2(vbcf_expressed_in_seu_polar[1,0], vbcf_expressed_in_seu_polar[0,0])
        return heading_bcf_polar

    def calc_heading_bcf_polar_derivs(self, t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of heading angle using velocity w.r.t. BCF frame and polar SEU coordinates. Due south in polar SEU is zero heading, with positive toward east."""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        R_bcf2seu_polar = self.calc_R_bcf2seu_polar(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        d_R_bcf2seu_d_t, d_R_bcf2seu_d_rbcf, d_R_bcf2seu_d_vbcf = self.calc_R_bcf2seu_polar_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        vbcf_expressed_in_seu_polar = np.dot(R_bcf2seu_polar, vbcf)
        d_head_d_vx_polar = (-vbcf_expressed_in_seu_polar[1,0]) / (vbcf_expressed_in_seu_polar[0,0]**2 + vbcf_expressed_in_seu_polar[1,0]**2)
        d_head_d_vy_polar = (vbcf_expressed_in_seu_polar[0,0]) / (vbcf_expressed_in_seu_polar[0,0]**2 + vbcf_expressed_in_seu_polar[1,0]**2)

        # time
        d_vbcf_expressed_in_seu_polar_d_t = np.dot(d_R_bcf2seu_d_t, vbcf)
        d_vx_polar_d_t = d_vbcf_expressed_in_seu_polar_d_t[0,0]
        d_vy_polar_d_t = d_vbcf_expressed_in_seu_polar_d_t[1,0]
        d_head_d_t = d_head_d_vy_polar * d_vy_polar_d_t + d_head_d_vx_polar * d_vx_polar_d_t # final chain rule for d/d t

        # position
        d_vbcf_expressed_in_seu_polar_d_rbcf = m_u.tensor_bullet2_vector(d_R_bcf2seu_d_rbcf, vbcf)
        d_vx_polar_d_rbcf = d_vbcf_expressed_in_seu_polar_d_rbcf[0,:]
        d_vy_polar_d_rbcf = d_vbcf_expressed_in_seu_polar_d_rbcf[1,:]
        d_head_d_rbcf = d_head_d_vy_polar * d_vy_polar_d_rbcf + d_head_d_vx_polar * d_vx_polar_d_rbcf # final chain rule for d/d rbcf


        # velocity
        d_vbcf_expressed_in_seu_polar_d_vbcf = m_u.tensor_bullet2_vector(d_R_bcf2seu_d_vbcf, vbcf) + R_bcf2seu_polar
        d_vx_polar_d_vbcf = d_vbcf_expressed_in_seu_polar_d_vbcf[0,:]
        d_vy_polar_d_vbcf = d_vbcf_expressed_in_seu_polar_d_vbcf[1,:]
        d_head_d_vbcf = d_head_d_vy_polar * d_vy_polar_d_vbcf + d_head_d_vx_polar * d_vx_polar_d_vbcf # final chain rule for d/d vbcf
        
        return d_head_d_t, d_head_d_rbcf, d_head_d_vbcf

    def calc_heading_bci_polar(self, t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate heading angle using velocity w.r.t. BCI frame and polar SEU coordinates. Due south in topocentric SEU is zero heading, with positive toward east."""
        # convert vbcf to vbci_expressed_in_bcf
        vbci_expressed_in_bcf = self.calc_vbci_in_bcf(rbcf, vbcf, w0dotrad)

        # express in seu polar
        R_bcf2seu_polar = self.calc_R_bcf2seu_polar(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        vbci_expressed_in_seu_polar = np.dot(R_bcf2seu_polar, vbci_expressed_in_bcf)
        heading_bci_polar = np.arctan2(vbci_expressed_in_seu_polar[1,0], vbci_expressed_in_seu_polar[0,0])
        return heading_bci_polar

    def calc_heading_bci_polar_derivs(self, t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of heading angle using velocity w.r.t. BCI frame and polar SEU coordinates. Due south in polar SEU is zero heading, with positive toward east."""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # BCI velocity
        vbci_in_bcf = self.calc_vbci_in_bcf(rbcf, vbcf, w0dotrad)

        R_bcf2seu_polar = self.calc_R_bcf2seu_polar(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        d_R_bcf2seu_d_t, d_R_bcf2seu_d_rbcf, d_R_bcf2seu_d_vbcf = self.calc_R_bcf2seu_polar_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        vbci_expressed_in_seu_polar = np.dot(R_bcf2seu_polar, vbci_in_bcf) # only one transformation because vbci is expressed in bcf to start with
        d_head_d_vx_polar = (-vbci_expressed_in_seu_polar[1,0]) / (vbci_expressed_in_seu_polar[0,0]**2 + vbci_expressed_in_seu_polar[1,0]**2)
        d_head_d_vy_polar = (vbci_expressed_in_seu_polar[0,0]) / (vbci_expressed_in_seu_polar[0,0]**2 + vbci_expressed_in_seu_polar[1,0]**2)

        # time
        d_vbci_expressed_in_seu_polar_d_t = np.dot(d_R_bcf2seu_d_t, vbci_in_bcf)
        d_vx_polar_d_t = d_vbci_expressed_in_seu_polar_d_t[0,0]
        d_vy_polar_d_t = d_vbci_expressed_in_seu_polar_d_t[1,0]
        d_head_d_t = d_head_d_vy_polar * d_vy_polar_d_t + d_head_d_vx_polar * d_vx_polar_d_t # final chain rule for d/d t

        # position
        omega = self.calc_omega(w0dotrad)
        omega_cross = m_u.crossmat(omega)
        d_vbci_expressed_in_seu_polar_d_rbcf = m_u.tensor_bullet2_vector(d_R_bcf2seu_d_rbcf, vbci_in_bcf) + np.dot(R_bcf2seu_polar, omega_cross)
        d_vx_polar_d_rbcf = d_vbci_expressed_in_seu_polar_d_rbcf[0,:]
        d_vy_polar_d_rbcf = d_vbci_expressed_in_seu_polar_d_rbcf[1,:]
        d_head_d_rbcf = d_head_d_vy_polar * d_vy_polar_d_rbcf + d_head_d_vx_polar * d_vx_polar_d_rbcf # final chain rule for d/d rbcf


        # velocity
        d_vbci_expressed_in_seu_polar_d_vbcf = m_u.tensor_bullet2_vector(d_R_bcf2seu_d_vbcf, vbci_in_bcf) + R_bcf2seu_polar
        d_vx_polar_d_vbcf = d_vbci_expressed_in_seu_polar_d_vbcf[0,:]
        d_vy_polar_d_vbcf = d_vbci_expressed_in_seu_polar_d_vbcf[1,:]
        d_head_d_vbcf = d_head_d_vy_polar * d_vy_polar_d_vbcf + d_head_d_vx_polar * d_vx_polar_d_vbcf # final chain rule for d/d vbcf
        
        return d_head_d_t, d_head_d_rbcf, d_head_d_vbcf

    def calc_seu_topocentric(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate basis vectors of topocentric SEU frame, expressed in BCF frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # u (up or normal) vector in BCF frame
        u_bcf = self.calc_u_bcf(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        u_bcf_unit = u_bcf / m_u.column_vector_norm2(u_bcf)

        # e (east) vector in BCF frame
        k_bcf = m_u.k_unit()
        arg1 = m_u.crossmat(k_bcf)
        theta1rad, theta2rad, theta3rad = self.calc_thetas_linear(t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad) # euler angles
        R_bcf2pa = self.calc_R_bcf2pa(theta1rad, theta2rad, theta3rad) # transformation matrix
        R_pa2bcf = np.transpose(R_bcf2pa) # transformation matrix
        ellipsmat = self.calc_ellipsmat(a, b, c)
        arg2 = 2.0 * np.dot(R_pa2bcf, np.dot(ellipsmat, np.dot(R_bcf2pa, rbcf)))
        e_bcf = np.dot(arg1, arg2)
        e_bcf_topocentric_unit = e_bcf / m_u.column_vector_norm2(e_bcf)

        # pseudo-s (south) vector in BCF frame
        e_bcf_cross = m_u.crossmat(e_bcf)
        s_bcf = np.dot(e_bcf_cross, u_bcf)
        s_bcf_topocentric_unit = s_bcf / m_u.column_vector_norm2(s_bcf)

        return s_bcf_topocentric_unit, e_bcf_topocentric_unit, u_bcf_unit

    def calc_seu_polar(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate basis vectors of polar SEU frame, expressed in BCF frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # u (up or normal) vector in BCF frame
        u_bcf = self.calc_u_bcf(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        u_bcf_unit = u_bcf / m_u.column_vector_norm2(u_bcf)

        # pseudo-east vector in BCF frame
        k_bcf = m_u.k_unit()
        kcross = m_u.crossmat(k_bcf)
        e_tilde_bcf = np.dot(kcross, rbcf)
        e_tilde_bcf_polar_unit = e_tilde_bcf / m_u.column_vector_norm2(e_tilde_bcf)

        # south vector in BCF frame
        e_tilde_bcf_cross = m_u.crossmat(e_tilde_bcf)
        s_bcf = np.dot(e_tilde_bcf_cross, u_bcf)
        s_bcf_polar_unit = s_bcf / m_u.column_vector_norm2(s_bcf)

        # east vector in BCF frame (still a sort of pseudo-east)
        u_bcf_cross = m_u.crossmat(u_bcf)
        e_bcf = np.dot(u_bcf_cross, s_bcf)
        e_bcf_polar_unit = e_bcf / m_u.column_vector_norm2(e_bcf)

        return s_bcf_polar_unit, e_bcf_polar_unit, u_bcf_unit

    def calc_seu_topocentric_derivs(self, t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of basis vectors of topocentric SEU frame, expressed in BCF frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        ellipsmat = self.calc_ellipsmat(a, b, c)
        theta1rad, theta2rad, theta3rad = self.calc_thetas_linear(t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad)
        R_bcf2pa = self.calc_R_bcf2pa(theta1rad, theta2rad, theta3rad)
        R_pa2bcf = np.transpose(R_bcf2pa)
        d_R_bcf2pa_d_t = self.calc_d_R_bcf2pa_d_t(theta1rad, theta2rad, theta3rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad)
        d_R_pa2bcf_d_t = np.transpose(d_R_bcf2pa_d_t)
        k_bcf = m_u.k_unit()
        kcross = m_u.crossmat(k_bcf)

        # up
        u_bcf = self.calc_u_bcf(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)

        # up vector derivatives already known from other function
        d_u_bcf_d_t, d_u_bcf_d_rbcf, d_u_bcf_d_vbcf = self.calc_u_bcf_derivs(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)

        # convert to unit vector derivatives
        d_u_bcf_unit_d_u_bcf = m_u.unit_vector_deriv(u_bcf)
        d_u_bcf_unit_d_t = np.dot(d_u_bcf_unit_d_u_bcf, d_u_bcf_d_t) 
        d_u_bcf_unit_d_rbcf = np.dot(d_u_bcf_unit_d_u_bcf, d_u_bcf_d_rbcf)
        d_u_bcf_unit_d_vbcf = np.dot(d_u_bcf_unit_d_u_bcf, d_u_bcf_d_vbcf)
        
        # east
        arg2 = 2.0 * np.dot(R_pa2bcf, np.dot(ellipsmat, np.dot(R_bcf2pa, rbcf)))
        e_bcf = np.dot(kcross, arg2)
        d_e_bcf_unit_d_e_bcf = m_u.unit_vector_deriv(e_bcf)

        # time
        d_e_bcf_d_t = 2. * np.dot(kcross, np.dot(np.dot(d_R_pa2bcf_d_t, np.dot(ellipsmat, R_bcf2pa)) + np.dot(R_pa2bcf, np.dot(ellipsmat, d_R_bcf2pa_d_t)), rbcf))        
        d_e_bcf_unit_d_t = np.dot(d_e_bcf_unit_d_e_bcf, d_e_bcf_d_t) 

        # position
        d_e_bcf_d_rbcf = 2. * np.dot(kcross, np.dot(R_pa2bcf, np.dot(ellipsmat, np.dot(R_bcf2pa, np.identity(3)))))
        d_e_bcf_unit_d_rbcf = np.dot(d_e_bcf_unit_d_e_bcf, d_e_bcf_d_rbcf)

        # velocity
        d_e_bcf_unit_d_vbcf = np.zeros((3,3))

        # south
        e_bcf_cross = m_u.crossmat(e_bcf)
        s_bcf = np.dot(e_bcf_cross, u_bcf)
        d_s_bcf_unit_d_s_bcf = m_u.unit_vector_deriv(s_bcf)

        # time
        d_e_bcf_cross_d_t = np.array([[0., -d_e_bcf_d_t[2], d_e_bcf_d_t[1]], [d_e_bcf_d_t[2], 0.0, -d_e_bcf_d_t[0]], [-d_e_bcf_d_t[1], d_e_bcf_d_t[0], 0.0]])
        d_s_bcf_d_t = np.dot(d_e_bcf_cross_d_t, u_bcf) + np.dot(e_bcf_cross, d_u_bcf_d_t)
        d_s_bcf_unit_d_t = np.dot(d_s_bcf_unit_d_s_bcf, d_s_bcf_d_t) 

        # position
        d_e_bcf_cross_d_rbcf_x = np.array([[0., -d_e_bcf_d_rbcf[2,0], d_e_bcf_d_rbcf[1,0]], [d_e_bcf_d_rbcf[2,0], 0.0, -d_e_bcf_d_rbcf[0,0]], [-d_e_bcf_d_rbcf[1,0], d_e_bcf_d_rbcf[0,0], 0.0]])
        d_e_bcf_cross_d_rbcf_y = np.array([[0., -d_e_bcf_d_rbcf[2,1], d_e_bcf_d_rbcf[1,1]], [d_e_bcf_d_rbcf[2,1], 0.0, -d_e_bcf_d_rbcf[0,1]], [-d_e_bcf_d_rbcf[1,1], d_e_bcf_d_rbcf[0,1], 0.0]])
        d_e_bcf_cross_d_rbcf_z = np.array([[0., -d_e_bcf_d_rbcf[2,2], d_e_bcf_d_rbcf[1,2]], [d_e_bcf_d_rbcf[2,2], 0.0, -d_e_bcf_d_rbcf[0,2]], [-d_e_bcf_d_rbcf[1,2], d_e_bcf_d_rbcf[0,2], 0.0]])
        d_e_bcf_cross_d_rbcf = np.zeros((3,3,3))
        d_e_bcf_cross_d_rbcf[:,:,0] = d_e_bcf_cross_d_rbcf_x
        d_e_bcf_cross_d_rbcf[:,:,1] = d_e_bcf_cross_d_rbcf_y
        d_e_bcf_cross_d_rbcf[:,:,2] = d_e_bcf_cross_d_rbcf_z
        tens_times_vec = m_u.tensor_bullet2_vector(d_e_bcf_cross_d_rbcf, u_bcf)
        d_s_bcf_d_rbcf = tens_times_vec + np.dot(e_bcf_cross, d_u_bcf_d_rbcf)
        d_s_bcf_unit_d_rbcf = np.dot(d_s_bcf_unit_d_s_bcf, d_s_bcf_d_rbcf)

        # velocity
        d_s_bcf_d_vcbf = np.dot(e_bcf_cross, d_u_bcf_d_vbcf)
        d_s_bcf_unit_d_vbcf = np.dot(d_s_bcf_unit_d_s_bcf, d_s_bcf_d_vcbf)

        return d_s_bcf_unit_d_t, d_s_bcf_unit_d_rbcf, d_s_bcf_unit_d_vbcf, d_e_bcf_unit_d_t, d_e_bcf_unit_d_rbcf, d_e_bcf_unit_d_vbcf, d_u_bcf_unit_d_t, d_u_bcf_unit_d_rbcf, d_u_bcf_unit_d_vbcf

    def calc_seu_polar_derivs(self, t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of basis vectors of polar SEU frame, expressed in BCF frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # up
        u_bcf = self.calc_u_bcf(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)

        # up vector derivatives already known from other function
        d_u_bcf_d_t, d_u_bcf_d_rbcf, d_u_bcf_d_vbcf = self.calc_u_bcf_derivs(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)

        # convert to unit vector derivatives
        d_u_bcf_unit_d_u_bcf = m_u.unit_vector_deriv(u_bcf)
        d_u_bcf_unit_d_t = np.dot(d_u_bcf_unit_d_u_bcf, d_u_bcf_d_t) 
        d_u_bcf_unit_d_rbcf = np.dot(d_u_bcf_unit_d_u_bcf, d_u_bcf_d_rbcf)
        d_u_bcf_unit_d_vbcf = np.dot(d_u_bcf_unit_d_u_bcf, d_u_bcf_d_vbcf)

        # pseudo-east (e tilde)
        k_bcf = m_u.k_unit()
        kcross = m_u.crossmat(k_bcf)
        e_tilde_bcf = np.dot(kcross, rbcf)
        d_e_tilde_bcf_unit_d_e_tilde_bcf = m_u.unit_vector_deriv(e_tilde_bcf)

        # time
        d_e_tilde_bcf_d_t = np.zeros((3,1))
        d_e_tilde_bcf_unit_d_t = np.zeros((3,1))

        # position
        d_e_tilde_bcf_d_rbcf = kcross
        d_e_tilde_bcf_unit_d_rbcf = np.dot(d_e_tilde_bcf_unit_d_e_tilde_bcf, d_e_tilde_bcf_d_rbcf)

        # velocity
        d_e_tilde_bcf_d_vbcf = np.zeros((3,3))
        d_e_tilde_bcf_unit_d_vbcf = np.zeros((3,3))

        # south
        e_tilde_bcf_cross = m_u.crossmat(e_tilde_bcf)
        s_bcf = np.dot(e_tilde_bcf_cross, u_bcf)
        d_s_bcf_unit_d_s_bcf = m_u.unit_vector_deriv(s_bcf)

        # time
        d_s_bcf_d_t = np.dot(e_tilde_bcf_cross, d_u_bcf_d_t)
        d_s_bcf_unit_d_t = np.dot(d_s_bcf_unit_d_s_bcf, d_s_bcf_d_t)

        # position
        d_e_tilde_bcf_cross_d_rbcf_x = np.array([[0., -d_e_tilde_bcf_d_rbcf[2,0], d_e_tilde_bcf_d_rbcf[1,0]], [d_e_tilde_bcf_d_rbcf[2,0], 0.0, -d_e_tilde_bcf_d_rbcf[0,0]], [-d_e_tilde_bcf_d_rbcf[1,0], d_e_tilde_bcf_d_rbcf[0,0], 0.0]])
        d_e_tilde_bcf_cross_d_rbcf_y = np.array([[0., -d_e_tilde_bcf_d_rbcf[2,1], d_e_tilde_bcf_d_rbcf[1,1]], [d_e_tilde_bcf_d_rbcf[2,1], 0.0, -d_e_tilde_bcf_d_rbcf[0,1]], [-d_e_tilde_bcf_d_rbcf[1,1], d_e_tilde_bcf_d_rbcf[0,1], 0.0]])
        d_e_tilde_bcf_cross_d_rbcf_z = np.array([[0., -d_e_tilde_bcf_d_rbcf[2,2], d_e_tilde_bcf_d_rbcf[1,2]], [d_e_tilde_bcf_d_rbcf[2,2], 0.0, -d_e_tilde_bcf_d_rbcf[0,2]], [-d_e_tilde_bcf_d_rbcf[1,2], d_e_tilde_bcf_d_rbcf[0,2], 0.0]])
        d_e_tilde_bcf_cross_d_rbcf = np.zeros((3,3,3))
        d_e_tilde_bcf_cross_d_rbcf[:,:,0] = d_e_tilde_bcf_cross_d_rbcf_x
        d_e_tilde_bcf_cross_d_rbcf[:,:,1] = d_e_tilde_bcf_cross_d_rbcf_y
        d_e_tilde_bcf_cross_d_rbcf[:,:,2] = d_e_tilde_bcf_cross_d_rbcf_z
        tens_times_vec = m_u.tensor_bullet2_vector(d_e_tilde_bcf_cross_d_rbcf, u_bcf)
        d_s_bcf_d_rbcf = tens_times_vec + np.dot(e_tilde_bcf_cross, d_u_bcf_d_rbcf)
        d_s_bcf_unit_d_rbcf = np.dot(d_s_bcf_unit_d_s_bcf, d_s_bcf_d_rbcf)

        # velocity
        #d_s_bcf_d_vbcf = np.dot(e_tilde_bcf_cross, d_u_bcf_d_vbcf)
        #d_s_bcf_unit_d_vbcf = np.dot(d_s_bcf_unit_d_s_bcf, d_s_bcf_d_vbcf)
        d_s_bcf_d_vbcf = np.zeros((3,3))
        d_s_bcf_unit_d_vbcf = np.zeros((3,3))

        # east
        u_bcf_cross = m_u.crossmat(u_bcf)
        e_bcf = np.dot(u_bcf_cross, s_bcf)
        d_e_bcf_unit_d_e_bcf = m_u.unit_vector_deriv(e_bcf)

        # time
        d_u_bcf_cross_d_t = m_u.crossmat(d_u_bcf_d_t)
        d_e_bcf_d_t = np.dot(d_u_bcf_cross_d_t, s_bcf) + np.dot(u_bcf_cross, d_s_bcf_d_t)
        d_e_bcf_unit_d_t = np.dot(d_e_bcf_unit_d_e_bcf, d_e_bcf_d_t)

        # position
        d_u_bcf_cross_d_rbcf_x = np.array([[0., -d_u_bcf_d_rbcf[2,0], d_u_bcf_d_rbcf[1,0]], [d_u_bcf_d_rbcf[2,0], 0.0, -d_u_bcf_d_rbcf[0,0]], [-d_u_bcf_d_rbcf[1,0], d_u_bcf_d_rbcf[0,0], 0.0]])
        d_u_bcf_cross_d_rbcf_y = np.array([[0., -d_u_bcf_d_rbcf[2,1], d_u_bcf_d_rbcf[1,1]], [d_u_bcf_d_rbcf[2,1], 0.0, -d_u_bcf_d_rbcf[0,1]], [-d_u_bcf_d_rbcf[1,1], d_u_bcf_d_rbcf[0,1], 0.0]])
        d_u_bcf_cross_d_rbcf_z = np.array([[0., -d_u_bcf_d_rbcf[2,2], d_u_bcf_d_rbcf[1,2]], [d_u_bcf_d_rbcf[2,2], 0.0, -d_u_bcf_d_rbcf[0,2]], [-d_u_bcf_d_rbcf[1,2], d_u_bcf_d_rbcf[0,2], 0.0]])
        d_u_bcf_cross_d_rbcf = np.zeros((3,3,3))
        d_u_bcf_cross_d_rbcf[:,:,0] = d_u_bcf_cross_d_rbcf_x
        d_u_bcf_cross_d_rbcf[:,:,1] = d_u_bcf_cross_d_rbcf_y
        d_u_bcf_cross_d_rbcf[:,:,2] = d_u_bcf_cross_d_rbcf_z
        tens_times_vec = m_u.tensor_bullet2_vector(d_u_bcf_cross_d_rbcf, s_bcf)
        d_e_bcf_d_rbcf = tens_times_vec + np.dot(u_bcf_cross, d_s_bcf_d_rbcf)
        d_e_bcf_unit_d_rbcf = np.dot(d_e_bcf_unit_d_e_bcf, d_e_bcf_d_rbcf)

        # velocity
        #d_e_bcf_d_vbcf = np.dot(u_bcf_cross, d_s_bcf_d_vbcf)
        #d_e_bcf_unit_d_vbcf = np.dot(d_e_bcf_unit_d_e_bcf, d_e_bcf_d_vbcf)
        d_e_bcf_d_vbcf = np.zeros((3,3))
        d_e_bcf_unit_d_vbcf = np.zeros((3,3))

        return d_s_bcf_unit_d_t, d_s_bcf_unit_d_rbcf, d_s_bcf_unit_d_vbcf, d_e_bcf_unit_d_t, d_e_bcf_unit_d_rbcf, d_e_bcf_unit_d_vbcf, d_u_bcf_unit_d_t, d_u_bcf_unit_d_rbcf, d_u_bcf_unit_d_vbcf

    def calc_R_bcf2seu_topocentric(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate transformation matrix from BCF frame to topocentric SEU frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # get the unit vectors of the topocentric SEU frame in BCF frame
        s_bcf_topocentric_unit, e_bcf_topocentric_unit, u_bcf_unit = self.calc_seu_topocentric(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)

        # unit vectors of BCF frame in BCF frame
        i_bcf = m_u.i_unit()
        j_bcf = m_u.j_unit()
        k_bcf = m_u.k_unit()

        # transformation matrix
        R_bcf2seu = np.array([[np.dot(np.transpose(s_bcf_topocentric_unit[:,0]), i_bcf[:,0]), np.dot(np.transpose(s_bcf_topocentric_unit[:,0]), j_bcf[:,0]), np.dot(np.transpose(s_bcf_topocentric_unit[:,0]), k_bcf[:,0])],
                              [np.dot(np.transpose(e_bcf_topocentric_unit[:,0]), i_bcf[:,0]), np.dot(np.transpose(e_bcf_topocentric_unit[:,0]), j_bcf[:,0]), np.dot(np.transpose(e_bcf_topocentric_unit[:,0]), k_bcf[:,0])],
                              [np.dot(np.transpose(u_bcf_unit[:,0]), i_bcf[:,0]), np.dot(np.transpose(u_bcf_unit[:,0]), j_bcf[:,0]), np.dot(np.transpose(u_bcf_unit[:,0]), k_bcf[:,0])]])
 
        return R_bcf2seu

    def calc_R_bcf2seu_topocentric_derivs(self, t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of transformation matrix from BCF frame to topocentric SEU frame"""

        # get the derivatives of the unit vectors that make up the transformation matrix
        d_s_bcf_unit_d_t, d_s_bcf_unit_d_rbcf, d_s_bcf_unit_d_vbcf, d_e_bcf_unit_d_t, d_e_bcf_unit_d_rbcf, d_e_bcf_unit_d_vbcf, d_u_bcf_unit_d_t, d_u_bcf_unit_d_rbcf, d_u_bcf_unit_d_vbcf = self.calc_seu_topocentric_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)

        # time
        d_s_bcf_unit_d_t = np.array([[np.asscalar(d_s_bcf_unit_d_t[0,0])], [np.asscalar(d_s_bcf_unit_d_t[1,0])], [np.asscalar(d_s_bcf_unit_d_t[2,0])]]) # need to reorganize the array to make it cleaner
        d_R_d_t = np.zeros((3,3))
        d_R_d_t[0,:] = np.reshape(d_s_bcf_unit_d_t, (3))
        d_R_d_t[1,:] = np.reshape(d_e_bcf_unit_d_t, (3))
        d_R_d_t[2,:] = np.reshape(d_u_bcf_unit_d_t, (3))
        #d_R_d_t = np.array([[np.reshape(d_s_bcf_unit_d_t, (3))], [np.reshape(d_e_bcf_unit_d_t, (3))], [np.reshape(d_u_bcf_unit_d_t, (3))]])

        # position
        d_R_d_rbcf = np.zeros((3,3,3))
        d_R_d_rbcf[0,:,:] = np.reshape(d_s_bcf_unit_d_rbcf, (3,3))
        d_R_d_rbcf[1,:,:] = np.reshape(d_e_bcf_unit_d_rbcf, (3,3))
        d_R_d_rbcf[2,:,:] = np.reshape(d_u_bcf_unit_d_rbcf, (3,3))

        # velocity (no dependence)
        d_R_d_vbcf = np.zeros((3,3,3))
        return d_R_d_t, d_R_d_rbcf, d_R_d_vbcf

    def calc_R_bcf2seu_polar(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate basis vectors of topocentric SEU frame, expressed in BCF frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # get the unit vectors of the topocentric SEU frame in BCF frame
        s_bcf_polar_unit, e_bcf_polar_unit, u_bcf_unit = self.calc_seu_polar(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)

        # unit vectors of BCF frame in BCF frame
        i_bcf = m_u.i_unit()
        j_bcf = m_u.j_unit()
        k_bcf = m_u.k_unit()

        # transformation matrix
        R_bcf2seu = np.array([[np.dot(np.transpose(s_bcf_polar_unit[:,0]), i_bcf[:,0]), np.dot(np.transpose(s_bcf_polar_unit[:,0]), j_bcf[:,0]), np.dot(np.transpose(s_bcf_polar_unit[:,0]), k_bcf[:,0])],
                              [np.dot(np.transpose(e_bcf_polar_unit[:,0]), i_bcf[:,0]), np.dot(np.transpose(e_bcf_polar_unit[:,0]), j_bcf[:,0]), np.dot(np.transpose(e_bcf_polar_unit[:,0]), k_bcf[:,0])],
                              [np.dot(np.transpose(u_bcf_unit[:,0]), i_bcf[:,0]), np.dot(np.transpose(u_bcf_unit[:,0]), j_bcf[:,0]), np.dot(np.transpose(u_bcf_unit[:,0]), k_bcf[:,0])]])
        return R_bcf2seu

    def calc_R_bcf2seu_polar_derivs(self, t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of transformation matrix from BCF frame to polar SEU frame"""

        # get the derivatives of the unit vectors that make up the transformation matrix
        d_s_bcf_unit_d_t, d_s_bcf_unit_d_rbcf, d_s_bcf_unit_d_vbcf, d_e_bcf_unit_d_t, d_e_bcf_unit_d_rbcf, d_e_bcf_unit_d_vbcf, d_u_bcf_unit_d_t, d_u_bcf_unit_d_rbcf, d_u_bcf_unit_d_vbcf = self.calc_seu_polar_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)

        # time
        d_s_bcf_unit_d_t = np.array([[np.asscalar(d_s_bcf_unit_d_t[0,0])], [np.asscalar(d_s_bcf_unit_d_t[1,0])], [np.asscalar(d_s_bcf_unit_d_t[2,0])]]) # need to reorganize the array to make it cleaner
        d_R_d_t = np.zeros((3,3))
        d_R_d_t[0,:] = np.reshape(d_s_bcf_unit_d_t, (3))
        d_R_d_t[1,:] = np.reshape(d_e_bcf_unit_d_t, (3))
        d_R_d_t[2,:] = np.reshape(d_u_bcf_unit_d_t, (3))
        #d_R_d_t = np.array([[np.reshape(d_s_bcf_unit_d_t, (3))], [np.reshape(d_e_bcf_unit_d_t, (3))], [np.reshape(d_u_bcf_unit_d_t, (3))]])

        # position
        d_R_d_rbcf = np.zeros((3,3,3))
        d_R_d_rbcf[0,:,:] = np.reshape(d_s_bcf_unit_d_rbcf, (3,3))
        d_R_d_rbcf[1,:,:] = np.reshape(d_e_bcf_unit_d_rbcf, (3,3))
        d_R_d_rbcf[2,:,:] = np.reshape(d_u_bcf_unit_d_rbcf, (3,3))

        # velocity (no dependence)
        d_R_d_vbcf = np.zeros((3,3,3))
        return d_R_d_t, d_R_d_rbcf, d_R_d_vbcf

    def calc_flight_path_angle_bcf(self, t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate flight path angle using velocity wrt the BCF frame as the reference vector"""
        R_bcf2seu_topocentric = self.calc_R_bcf2seu_topocentric(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        vbcf_expressed_in_seu_topocentric = np.dot(R_bcf2seu_topocentric, vbcf)
        vxy = (vbcf_expressed_in_seu_topocentric[0,0]**2 + vbcf_expressed_in_seu_topocentric[1,0]**2)**0.5
        gamma_bcf = np.arctan2(vbcf_expressed_in_seu_topocentric[2,0], vxy)
        return gamma_bcf

    def calc_flight_path_angle_bcf_derivs(self, t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of flight path angle using velocity wrt the BCF frame as the reference vector"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        vbcf_mag = m_u.column_vector_norm2(vbcf)

        # in plane and out of plane refer to the plane normal to the ellipsoid. use seu topocentric coordinates for this
        R_bcf2seu_topocentric = self.calc_R_bcf2seu_topocentric(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        d_R_bcf2seu_d_t, d_R_bcf2seu_d_rbcf, d_R_bcf2seu_d_vbcf = self.calc_R_bcf2seu_topocentric_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        vbcf_expressed_in_seu_topocentric = np.dot(R_bcf2seu_topocentric, vbcf)
        v_in_plane_mag = (vbcf_expressed_in_seu_topocentric[0,0]**2 + vbcf_expressed_in_seu_topocentric[1,0]**2)**0.5
        v_out_of_plane = vbcf_expressed_in_seu_topocentric[2,0]
        d_gamma_d_v_out_of_plane = v_in_plane_mag / (vbcf_mag**2)
        d_gamma_d_v_in_plane_mag = -v_out_of_plane / (vbcf_mag**2)

        d_v_out_of_plane_d_v = np.array([0.0, 0.0, 1.0]) # both v's are seu topocentric expressed in bcf
        d_v_in_plane_mag_d_v = (1. / v_in_plane_mag) * np.array([vbcf_expressed_in_seu_topocentric[0,0], vbcf_expressed_in_seu_topocentric[1,0], 0.0])

        # time
        d_vbcf_expressed_in_seu_topocentric_d_t = np.dot(d_R_bcf2seu_d_t, vbcf)
        d_v_out_of_plane_d_t = np.asscalar(np.dot(np.reshape(d_v_out_of_plane_d_v, (1,3)), d_vbcf_expressed_in_seu_topocentric_d_t))
        d_v_in_plane_mag_d_t = np.asscalar(np.dot(np.reshape(d_v_in_plane_mag_d_v, (1,3)), d_vbcf_expressed_in_seu_topocentric_d_t))
        d_gamma_d_t = np.asscalar(d_gamma_d_v_out_of_plane * d_v_out_of_plane_d_t + d_gamma_d_v_in_plane_mag * d_v_in_plane_mag_d_t)

        # position
        d_vbcf_expressed_in_seu_topocentric_d_rbcf = m_u.tensor_bullet2_vector(d_R_bcf2seu_d_rbcf, vbcf)
        d_v_out_of_plane_d_rbcf = np.dot(np.reshape(d_v_out_of_plane_d_v, (1,3)), d_vbcf_expressed_in_seu_topocentric_d_rbcf)
        d_v_in_plane_mag_d_rbcf = np.dot(np.reshape(d_v_in_plane_mag_d_v, (1,3)), d_vbcf_expressed_in_seu_topocentric_d_rbcf)
        d_gamma_d_rbcf = d_gamma_d_v_out_of_plane * d_v_out_of_plane_d_rbcf + d_gamma_d_v_in_plane_mag * d_v_in_plane_mag_d_rbcf

        # velocity
        d_vbcf_expressed_in_seu_topocentric_d_vbcf = m_u.tensor_bullet2_vector(d_R_bcf2seu_d_vbcf, vbcf) + R_bcf2seu_topocentric
        d_v_out_of_plane_d_vbcf = np.dot(np.reshape(d_v_out_of_plane_d_v, (1,3)), d_vbcf_expressed_in_seu_topocentric_d_vbcf)
        d_v_in_plane_mag_d_vbcf = np.dot(np.reshape(d_v_in_plane_mag_d_v, (1,3)), d_vbcf_expressed_in_seu_topocentric_d_vbcf)
        d_gamma_d_vbcf = d_gamma_d_v_out_of_plane * d_v_out_of_plane_d_vbcf + d_gamma_d_v_in_plane_mag * d_v_in_plane_mag_d_vbcf
        return d_gamma_d_t, d_gamma_d_rbcf, d_gamma_d_vbcf

    def calc_flight_path_angle_bci(self, t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate flight path angle using velocity wrt the BCI frame as the reference vector"""

        # convert vbcf to vbci_expressed_in_bcf
        vbci_expressed_in_bcf = self.calc_vbci_in_bcf(rbcf, vbcf, w0dotrad)

        # express in seu topocentric
        R_bcf2seu_topocentric = self.calc_R_bcf2seu_topocentric(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        vbci_expressed_in_seu_topocentric = np.dot(R_bcf2seu_topocentric, vbci_expressed_in_bcf)
        vxy = (vbci_expressed_in_seu_topocentric[0,0]**2 + vbci_expressed_in_seu_topocentric[1,0]**2)**0.5
        gamma_bci = np.arctan2(vbci_expressed_in_seu_topocentric[2,0], vxy)
        return gamma_bci

    def calc_flight_path_angle_bci_derivs(self, t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate derivatives of flight path angle using velocity wrt the BCI frame as the reference vector"""

        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # BCI velocity
        vbci_in_bcf = self.calc_vbci_in_bcf(rbcf, vbcf, w0dotrad)

        vbci_mag = self.calc_vbci_mag(rbcf, vbcf, w0dotrad)

        # in plane and out of plane refer to the plane normal to the ellipsoid. use seu topocentric coordinates for this
        R_bcf2seu_topocentric = self.calc_R_bcf2seu_topocentric(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        d_R_bcf2seu_d_t, d_R_bcf2seu_d_rbcf, d_R_bcf2seu_d_vbcf = self.calc_R_bcf2seu_topocentric_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        vbci_expressed_in_seu_topocentric = np.dot(R_bcf2seu_topocentric, vbci_in_bcf) # only one transformation because vbci is expressed in bcf to start with
        v_in_plane_mag = (vbci_expressed_in_seu_topocentric[0,0]**2 + vbci_expressed_in_seu_topocentric[1,0]**2)**0.5
        v_out_of_plane = vbci_expressed_in_seu_topocentric[2,0]
        d_gamma_d_v_out_of_plane = v_in_plane_mag / (vbci_mag**2)
        d_gamma_d_v_in_plane_mag = -v_out_of_plane / (vbci_mag**2)

        d_v_out_of_plane_d_v = np.array([0.0, 0.0, 1.0]) # both v's are seu topocentric expressed in bcf
        d_v_in_plane_mag_d_v = (1. / v_in_plane_mag) * np.array([vbci_expressed_in_seu_topocentric[0,0], vbci_expressed_in_seu_topocentric[1,0], 0.0])

        # time
        d_vbci_expressed_in_seu_topocentric_d_t = np.dot(d_R_bcf2seu_d_t, vbci_in_bcf)
        d_v_out_of_plane_d_t = np.asscalar(np.dot(np.reshape(d_v_out_of_plane_d_v, (1,3)), d_vbci_expressed_in_seu_topocentric_d_t))
        d_v_in_plane_mag_d_t = np.asscalar(np.dot(np.reshape(d_v_in_plane_mag_d_v, (1,3)), d_vbci_expressed_in_seu_topocentric_d_t))
        d_gamma_d_t = np.asscalar(d_gamma_d_v_out_of_plane * d_v_out_of_plane_d_t + d_gamma_d_v_in_plane_mag * d_v_in_plane_mag_d_t)

        # position
        omega = self.calc_omega(w0dotrad)
        omega_cross = m_u.crossmat(omega)
        d_vbci_expressed_in_seu_topocentric_d_rbcf = m_u.tensor_bullet2_vector(d_R_bcf2seu_d_rbcf, vbci_in_bcf) + np.dot(R_bcf2seu_topocentric, omega_cross)
        d_v_out_of_plane_d_rbcf = np.dot(np.reshape(d_v_out_of_plane_d_v, (1,3)), d_vbci_expressed_in_seu_topocentric_d_rbcf)
        d_v_in_plane_mag_d_rbcf = np.dot(np.reshape(d_v_in_plane_mag_d_v, (1,3)), d_vbci_expressed_in_seu_topocentric_d_rbcf)
        d_gamma_d_rbcf = d_gamma_d_v_out_of_plane * d_v_out_of_plane_d_rbcf + d_gamma_d_v_in_plane_mag * d_v_in_plane_mag_d_rbcf

        # velocity
        d_vbci_expressed_in_seu_topocentric_d_vbcf = m_u.tensor_bullet2_vector(d_R_bcf2seu_d_vbcf, vbci_in_bcf) + R_bcf2seu_topocentric
        d_v_out_of_plane_d_vbcf = np.dot(np.reshape(d_v_out_of_plane_d_v, (1,3)), d_vbci_expressed_in_seu_topocentric_d_vbcf)
        d_v_in_plane_mag_d_vbcf = np.dot(np.reshape(d_v_in_plane_mag_d_v, (1,3)), d_vbci_expressed_in_seu_topocentric_d_vbcf)
        d_gamma_d_vbcf = d_gamma_d_v_out_of_plane * d_v_out_of_plane_d_vbcf + d_gamma_d_v_in_plane_mag * d_v_in_plane_mag_d_vbcf
        return d_gamma_d_t, d_gamma_d_rbcf, d_gamma_d_vbcf

    def calc_seu_polar_s(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate basis vectors of polar SEU frame, expressed in BCF frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # u (up or normal) vector in BCF frame
        u_bcf = self.calc_u_bcf(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        u_bcf_unit = u_bcf / m_u.column_vector_norm2(u_bcf)

        # pseudo-east vector in BCF frame
        k_bcf = m_u.k_unit()
        kcross = m_u.crossmat(k_bcf)
        e_tilde_bcf = np.dot(kcross, rbcf)
        e_tilde_bcf_polar_unit = e_tilde_bcf / m_u.column_vector_norm2(e_tilde_bcf)

        # south vector in BCF frame
        e_tilde_bcf_cross = m_u.crossmat(e_tilde_bcf)
        s_bcf = np.dot(e_tilde_bcf_cross, u_bcf)
        s_bcf_polar_unit = s_bcf / m_u.column_vector_norm2(s_bcf)

        # east vector in BCF frame (still a sort of pseudo-east)
        u_bcf_cross = m_u.crossmat(u_bcf)
        e_bcf = np.dot(u_bcf_cross, s_bcf)
        e_bcf_polar_unit = e_bcf / m_u.column_vector_norm2(e_bcf)

        return s_bcf_polar_unit

    def calc_seu_polar_e(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate basis vectors of polar SEU frame, expressed in BCF frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # u (up or normal) vector in BCF frame
        u_bcf = self.calc_u_bcf(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        u_bcf_unit = u_bcf / m_u.column_vector_norm2(u_bcf)

        # pseudo-east vector in BCF frame
        k_bcf = m_u.k_unit()
        kcross = m_u.crossmat(k_bcf)
        e_tilde_bcf = np.dot(kcross, rbcf)
        e_tilde_bcf_polar_unit = e_tilde_bcf / m_u.column_vector_norm2(e_tilde_bcf)

        # south vector in BCF frame
        e_tilde_bcf_cross = m_u.crossmat(e_tilde_bcf)
        s_bcf = np.dot(e_tilde_bcf_cross, u_bcf)
        s_bcf_polar_unit = s_bcf / m_u.column_vector_norm2(s_bcf)

        # east vector in BCF frame (still a sort of pseudo-east)
        u_bcf_cross = m_u.crossmat(u_bcf)
        e_bcf = np.dot(u_bcf_cross, s_bcf)
        e_bcf_polar_unit = e_bcf / m_u.column_vector_norm2(e_bcf)

        return e_bcf_polar_unit

    def calc_seu_polar_u(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
        """Calculate basis vectors of polar SEU frame, expressed in BCF frame"""
        m_u = math_utilities_autograd.math_utilities_autograd_class()

        # u (up or normal) vector in BCF frame
        u_bcf = self.calc_u_bcf(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
        u_bcf_unit = u_bcf / m_u.column_vector_norm2(u_bcf)

        # pseudo-east vector in BCF frame
        k_bcf = m_u.k_unit()
        kcross = m_u.crossmat(k_bcf)
        e_tilde_bcf = np.dot(kcross, rbcf)
        e_tilde_bcf_polar_unit = e_tilde_bcf / m_u.column_vector_norm2(e_tilde_bcf)

        # south vector in BCF frame
        e_tilde_bcf_cross = m_u.crossmat(e_tilde_bcf)
        s_bcf = np.dot(e_tilde_bcf_cross, u_bcf)
        s_bcf_polar_unit = s_bcf / m_u.column_vector_norm2(s_bcf)

        # east vector in BCF frame (still a sort of pseudo-east)
        u_bcf_cross = m_u.crossmat(u_bcf)
        e_bcf = np.dot(u_bcf_cross, s_bcf)
        e_bcf_polar_unit = e_bcf / m_u.column_vector_norm2(e_bcf)

        return u_bcf_unit



    #def calc_seu_topocentric_s(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
    #    """Calculate basis vectors of topocentric SEU frame, expressed in BCF frame"""
    #    m_u = math_utilities_autograd.math_utilities_autograd_class()

    #    # u (up or normal) vector in BCF frame
    #    u_bcf = self.calc_u_bcf(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
    #    u_bcf_unit = u_bcf / m_u.column_vector_norm2(u_bcf)

    #    # e (east) vector in BCF frame
    #    k_bcf = m_u.k_unit()
    #    arg1 = m_u.crossmat(k_bcf)
    #    theta1rad, theta2rad, theta3rad = self.calc_thetas_linear(t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad) # euler angles
    #    R_bcf2pa = self.calc_R_bcf2pa(theta1rad, theta2rad, theta3rad) # transformation matrix
    #    R_pa2bcf = np.transpose(R_bcf2pa) # transformation matrix
    #    ellipsmat = self.calc_ellipsmat(a, b, c)
    #    arg2 = 2.0 * np.dot(R_pa2bcf, np.dot(ellipsmat, np.dot(R_bcf2pa, rbcf)))
    #    e_bcf = np.dot(arg1, arg2)
    #    e_bcf_topocentric_unit = e_bcf / m_u.column_vector_norm2(e_bcf)

    #    # pseudo-s (south) vector in BCF frame
    #    e_bcf_cross = m_u.crossmat(e_bcf)
    #    s_bcf = np.dot(e_bcf_cross, u_bcf)
    #    s_bcf_topocentric_unit = s_bcf / m_u.column_vector_norm2(s_bcf)

    #    return s_bcf_topocentric_unit

    #def calc_seu_topocentric_e(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
    #    """Calculate basis vectors of topocentric SEU frame, expressed in BCF frame"""
    #    m_u = math_utilities_autograd.math_utilities_autograd_class()

    #    # u (up or normal) vector in BCF frame
    #    u_bcf = self.calc_u_bcf(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
    #    u_bcf_unit = u_bcf / m_u.column_vector_norm2(u_bcf)

    #    # e (east) vector in BCF frame
    #    k_bcf = m_u.k_unit()
    #    arg1 = m_u.crossmat(k_bcf)
    #    theta1rad, theta2rad, theta3rad = self.calc_thetas_linear(t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad) # euler angles
    #    R_bcf2pa = self.calc_R_bcf2pa(theta1rad, theta2rad, theta3rad) # transformation matrix
    #    R_pa2bcf = np.transpose(R_bcf2pa) # transformation matrix
    #    ellipsmat = self.calc_ellipsmat(a, b, c)
    #    arg2 = 2.0 * np.dot(R_pa2bcf, np.dot(ellipsmat, np.dot(R_bcf2pa, rbcf)))
    #    e_bcf = np.dot(arg1, arg2)
    #    e_bcf_topocentric_unit = e_bcf / m_u.column_vector_norm2(e_bcf)

    #    # pseudo-s (south) vector in BCF frame
    #    e_bcf_cross = m_u.crossmat(e_bcf)
    #    s_bcf = np.dot(e_bcf_cross, u_bcf)
    #    s_bcf_topocentric_unit = s_bcf / m_u.column_vector_norm2(s_bcf)

    #    return e_bcf_topocentric_unit

    #def calc_seu_topocentric_u(self, t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c):
    #    """Calculate basis vectors of topocentric SEU frame, expressed in BCF frame"""
    #    m_u = math_utilities_autograd.math_utilities_autograd_class()

    #    # u (up or normal) vector in BCF frame
    #    u_bcf = self.calc_u_bcf(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
    #    u_bcf_unit = u_bcf / m_u.column_vector_norm2(u_bcf)

    #    # e (east) vector in BCF frame
    #    k_bcf = m_u.k_unit()
    #    arg1 = m_u.crossmat(k_bcf)
    #    theta1rad, theta2rad, theta3rad = self.calc_thetas_linear(t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad) # euler angles
    #    R_bcf2pa = self.calc_R_bcf2pa(theta1rad, theta2rad, theta3rad) # transformation matrix
    #    R_pa2bcf = np.transpose(R_bcf2pa) # transformation matrix
    #    ellipsmat = self.calc_ellipsmat(a, b, c)
    #    arg2 = 2.0 * np.dot(R_pa2bcf, np.dot(ellipsmat, np.dot(R_bcf2pa, rbcf)))
    #    e_bcf = np.dot(arg1, arg2)
    #    e_bcf_topocentric_unit = e_bcf / m_u.column_vector_norm2(e_bcf)

    #    # pseudo-s (south) vector in BCF frame
    #    e_bcf_cross = m_u.crossmat(e_bcf)
    #    s_bcf = np.dot(e_bcf_cross, u_bcf)
    #    s_bcf_topocentric_unit = s_bcf / m_u.column_vector_norm2(s_bcf)

    #    return u_bcf_unit