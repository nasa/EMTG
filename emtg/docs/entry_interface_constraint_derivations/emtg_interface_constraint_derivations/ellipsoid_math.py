# ellipsoid_math.py

# Purpose
# Contains basic math utility functions in the class math_utilities

# Revision history
# Noble Hatten; 08/31/2018; Began
# Noble Hatten; 09/04/2018; Added bodydetic latitude and longitude calculations
# Noble Hatten; 09/07/2018; Changed definition of polar frame s.t. East is in tangent plane but does not necessarily have zero z component
                            # Added topocentric and polar versions of heading and flight path angle

# standard imports
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# custom imports
import math_utilities

class ellipsoid_math_class:
    """Contains routines for calculating ellipsoid interface quantities"""
    # Assumes a>b>c or a=b>c or a=b=c
    
    def __init__(self, a=1.0, b=1.0, c=1.0, t=0.0, r_pa=[1.0, 0.0, 0.0], vbcf=[0.0, 0.0, 0.0], theta1=0.0, theta2=0.0, theta3=0.0, w0=0.0, w0dot=0.0):
        """Constructor"""
        # Note that we take in r_pa rather than r_bcf so that we can be sure that the position is actually on the ellipse
        self.m_u = math_utilities.math_utilities_class()
        self.a = a
        self.b = b
        self.c = c
        self.ex = (a**2-c**2)/a**2
        self.ee = (a**2-b**2)/a**2
        self.calc_ellipsmat()
        # BCF frame unit vectors, expressed in BCF frame
        self.i_bcf = self.m_u.i_unit()
        self.j_bcf = self.m_u.j_unit()
        self.k_bcf = self.m_u.k_unit()
        self.set_state(t, r_pa, vbcf, theta1, theta2, theta3, w0, w0dot)
        return
    
    def calc_ellipsmat(self):
        """Turn geometric values that define ellipse into matrix form"""
        a2 = self.a*self.a
        b2 = self.b*self.b
        c2 = self.c*self.c
        self.ellipsmat = np.diagflat([1./a2, 1./b2, 1./c2])
        return

    def set_state(self, t, r_pa, vbcf, theta1, theta2, theta3, w0, w0dot):
        """Set everything that is just extracted directly from the state of the spacecraft and the body"""
        self.t = t
        self.r_pa = r_pa
        self.vbcf = np.reshape(vbcf, (3,1))
        self.theta1 = theta1
        self.theta2 = theta2
        self.theta3 = theta3
        self.theta1rad = np.deg2rad(theta1)
        self.theta2rad = np.deg2rad(theta2)
        self.theta3rad = np.deg2rad(theta3)
        self.w0 = w0
        self.w0dot = w0dot
        self.w0rad = np.deg2rad(w0)
        self.w0dotrad = np.deg2rad(w0dot)

        # derived quantities
        self.calc_R_bcf2pa()

        # state in BCF frame
        self.calc_rbcf()
        self.xbcf = np.concatenate((self.rbcf, self.vbcf), axis=0)

        # BCF velocity expressed in PA frame
        self.v_bcf_expressed_in_pa = np.dot(self.R_bcf2pa, self.vbcf)

        # angular velocity of BCF frame w.r.t. BCI frame
        self.w = self.w0 + self.w0dot * self.t
        self.omega = np.zeros((3,1))
        self.omega[2] = self.w0dotrad # in radians

        # vector normal to ellipsoid at point of interface
        self.n_pa = 2.0 * np.dot(self.ellipsmat, self.r_pa) # in principal axis frame
        n_pa_mag = (self.n_pa[0]**2 + self.n_pa[1]**2 + self.n_pa[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
        self.n_pa_unit = self.n_pa / n_pa_mag # unit vector

        # u (up) vector in BCF frame, i.e., n_bcf
        self.u_bcf = np.dot(self.R_pa2bcf, self.n_pa)
        u_bcf_mag = (self.u_bcf[0]**2 + self.u_bcf[1]**2 + self.u_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
        self.u_bcf_unit = self.u_bcf / u_bcf_mag # unit vector
        return

    def calc_R_bcf2pa(self):
        """Calculate transformation matrix s.t. r_PA = R_bcf2pa * r_BCF"""
        R_bcf2fprime = self.m_u.transform3(self.theta1rad)
        R_fprime2fdprime = self.m_u.transform1(self.theta2rad)
        R_fdprime2pa = self.m_u.transform3(self.theta3rad)
        self.R_bcf2pa = np.dot(R_fdprime2pa, np.dot(R_fprime2fdprime, R_bcf2fprime))
        self.R_pa2bcf = np.transpose(self.R_bcf2pa)
        return
    
    def calc_rbcf(self):
        """Calculate position vector in BCF frame given position vector in PA frame"""
        self.rbcf = np.zeros((3,1))
        self.rbcf = np.dot(np.transpose(self.R_bcf2pa), self.r_pa)
        return

    def calc_longitude_bodycentric(self):
        """Calculate bodycentric longitude in radians"""
        self.longitude_bodycentric = self.m_u.atan2(self.rbcf[1], self.rbcf[0])
        return

    def calc_longitude_bodydetic(self):
        """Calculate bodydetic longitude in radians"""

        # calculate normal vector of project of ellipsoid onto BCF xy plane:
        self.calc_longitude_bodycentric() # need bodycentric longitude

        #rproj_bcf_mag = (self.rbcf[0]**2 + self.rbcf[1]**2)**0.5
        rproj_bcf_mag = self.m_u.explicit_norm2(self.rbcf[0:2]) # magnitude of projection into xy plane
        rproj_bcf = rproj_bcf_mag * np.array([[np.cos(self.longitude_bodycentric)], [np.sin(self.longitude_bodycentric)], [0.]]) # magnitude + direction
        rproj_pa = np.dot(self.R_bcf2pa, rproj_bcf) # transform to PA frame
        nproj_pa = 2. * np.dot(self.ellipsmat, rproj_pa) # get normal vector of rproj in PA frame
        nproj_bcf = np.dot(self.R_pa2bcf, nproj_pa) # transform to BCF frame
        arg1 = nproj_bcf
        arg2 = self.i_bcf
        temp = self.m_u.angle_between_2_vectors(arg1, arg2) # bodydetic longitude is angle between nproj and i_bcf

        # quadrant check
        if (np.dot(np.transpose(nproj_bcf), self.j_bcf) >= 0.):
            self.longitude_bodydetic = temp
        else:
            self.longitude_bodydetic = 2.*np.pi - temp

        return

    def calc_longitude_bodydetic_bekta(self):
        """Calculate bodydetic longitude in rad"""
        # From Bekta, S., "Geodetic Computations on Triaxial Ellipsoid"
        # This form only works for point on surface of ellipsoid
        # Currently calculates longitude with respect to x axis defined by the PA axes, not the BCF axes and in the xy plane of the PA frame, not the BCF frame
        ee2 = self.ee**2
        argy = self.r_pa[1]
        argx = (1.-ee2)*self.r_pa[0]
        self.longitude_bodydetic_bekta = self.m_u.atan2(argy, argx)

        # quadrant check
        if (self.r_pa[0] < 0 and self.r_pa[1] >= 0):
            self.longitude_bodydetic_bekta = self.longitude_bodydetic_bekta + np.pi
        elif (self.r_pa[0] < 0 and self.r_pa[1] < 0):
            self.longitude_bodydetic_bekta = self.longitude_bodydetic_bekta - np.pi

        return

    def calc_latitude_bodycentric(self):
        """Calculate bodycentric latitude in radians"""
        self.calc_longitude_bodycentric()
        self.latitude_bodycentric = self.m_u.atan2(self.rbcf[2], (self.rbcf[0]*np.cos(self.longitude_bodycentric) + self.rbcf[1]*np.sin(self.longitude_bodycentric)))
        return

    def calc_latitude_bodydetic(self):
        """Calculate bodydetic latitude in rad"""
        # use angle between k and normal vector
        arg1 = self.k_bcf
        arg2 = self.u_bcf_unit
        temp = self.m_u.angle_between_2_vectors(arg1, arg2)
        complement = np.pi*0.5 - temp
        self.latitude_bodydetic = complement # bodydetic latitude is the complement of temp
        return

    def calc_latitude_bodydetic_bekta(self):
        """Calculate bodydetic latitude in rad"""
        # From Bekta, S., "Geodetic Computations on Triaxial Ellipsoid"
        # This form only works for point on surface of ellipsoid
        # Currently calculates latitude with respect to equator defined by the PA axes, not the BCF axes
        ee2 = self.ee**2
        ex2 = self.ex**2
        argy = (1. - ee2) * self.r_pa[2]
        argx = (1. - ex2) * np.sqrt((1. - ee2)**2 * self.r_pa[0]**2 + self.r_pa[1]**2)
        self.latitude_bodydetic_bekta = self.m_u.atan2(argy, argx)
        return

    def calc_latitude_bodydetic_oblate(self):
        """Calculate bodydetic latitude in rad."""
        # Only valid for oblate body (a=b), not general triaxial body
        # From https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=5&ved=2ahUKEwiaio-08KPdAhVKq1MKHWSiCv8QFjAEegQICBAC&url=http%3A%2F%2Fccar.colorado.edu%2Fasen5070%2Fhandouts%2Fgeodeticgeocentric.doc&usg=AOvVaw3qp43xImCsrvdcIHRI_arP
        self.calc_latitude_bodycentric()
        f = (self.a - self.c) / self.a # flattening parameter
        self.latitude_bodydetic_oblate = self.m_u.atan2(np.tan(self.latitude_bodycentric), (1.-f)**2)
        return

    def calc_vbcf_mag(self):
        """Calculate magnitude of velocity wrt BCF frame"""
        self.vbcf_mag = (self.vbcf[0]**2 + self.vbcf[1]**2 + self.vbcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
        return

    def calc_vbci_expressed_in_bcf(self):
        """Calculate velocity wrt BCI frame, expressed in BCF coordinates"""
        self.vbci_expressed_in_bcf = self.vbcf + np.cross(self.omega, self.rbcf, axis=0)
        return

    def calc_vbci_mag(self):
        """Calculate magnitude of velocity wrt BCI frame"""
        self.calc_vbci_expressed_in_bcf()
        self.vbci_mag = (self.vbci_expressed_in_bcf[0]**2 + self.vbci_expressed_in_bcf[1]**2 + self.vbci_expressed_in_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
        return

    def calc_seu_topocentric(self):
        """Calculate a south-east-up coordinate system in which the east vector points along a line of constant z and is in the plane tangent to the ellipsoid"""
        # In STK terminology, this would be called a topocentric frame

        # e (east) vector in BCF frame
        firstarg = self.m_u.crossmat(self.k_bcf)
        secondarg = 2.0 * np.dot(self.R_pa2bcf, np.dot(self.ellipsmat, np.dot(self.R_bcf2pa, self.rbcf)))
        e_bcf = np.dot(firstarg, secondarg)
        e_bcf_mag = (e_bcf[0]**2 + e_bcf[1]**2 + e_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
        self.e_bcf_topocentric_unit = e_bcf / e_bcf_mag # unit vector

        # pseudo-s (south) vector in BCF frame
        e_bcf_cross = self.m_u.crossmat(e_bcf)
        s_bcf = np.dot(e_bcf_cross, self.u_bcf_unit)
        s_bcf_mag = (s_bcf[0]**2 + s_bcf[1]**2 + s_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
        self.s_bcf_topocentric_unit = s_bcf / s_bcf_mag # unit vector

        # transformation matrix
        R_bcf2seu = np.zeros((3,3))

        R_bcf2seu[0,0] = np.dot(np.transpose(self.s_bcf_topocentric_unit), self.i_bcf)
        R_bcf2seu[0,1] = np.dot(np.transpose(self.s_bcf_topocentric_unit), self.j_bcf)
        R_bcf2seu[0,2] = np.dot(np.transpose(self.s_bcf_topocentric_unit), self.k_bcf)
        R_bcf2seu[1,0] = np.dot(np.transpose(self.e_bcf_topocentric_unit), self.i_bcf)
        R_bcf2seu[1,1] = np.dot(np.transpose(self.e_bcf_topocentric_unit), self.j_bcf)
        R_bcf2seu[1,2] = np.dot(np.transpose(self.e_bcf_topocentric_unit), self.k_bcf)
        R_bcf2seu[2,0] = np.dot(np.transpose(self.u_bcf_unit), self.i_bcf)
        R_bcf2seu[2,1] = np.dot(np.transpose(self.u_bcf_unit), self.j_bcf)
        R_bcf2seu[2,2] = np.dot(np.transpose(self.u_bcf_unit), self.k_bcf)

        self.R_bcf2seu_topocentric = R_bcf2seu
        return

    def calc_seu_polar(self):
        """Calculate a south-east-up coordinate system in which the north vector points toward the pole"""
        # east vector is tangent to ellipsoid but does not necessarily have zero z component
        firstarg = self.m_u.crossmat(self.k_bcf)
        secondarg = self.rbcf
        e_pseudo_bcf = np.dot(firstarg, secondarg) # pseudo-east: zero z component, but not necessarily in tangent plane
        e_pseudo_bcf_mag = (e_pseudo_bcf[0]**2 + e_pseudo_bcf[1]**2 + e_pseudo_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
        self.e_pseudo_bcf_unit = e_pseudo_bcf / e_pseudo_bcf_mag # unit vector
        e_pseudo_bcf_cross = self.m_u.crossmat(e_pseudo_bcf)
        s_polar_bcf = np.dot(e_pseudo_bcf_cross, self.u_bcf_unit)
        s_polar_bcf_mag = (s_polar_bcf[0]**2 + s_polar_bcf[1]**2 + s_polar_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
        self.s_polar_bcf_unit = s_polar_bcf / s_polar_bcf_mag # unit vector
        u_bcf_cross = self.m_u.crossmat(self.u_bcf_unit)
        e_polar_bcf = np.dot(u_bcf_cross, self.s_polar_bcf_unit)
        e_polar_bcf_mag = (e_polar_bcf[0]**2 + e_polar_bcf[1]**2 + e_polar_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
        self.e_polar_bcf_unit = e_polar_bcf / e_polar_bcf_mag # unit vector

        # transformation matrix
        R_bcf2seu = np.zeros((3,3))

        R_bcf2seu[0,0] = np.dot(np.transpose(self.s_polar_bcf_unit), self.i_bcf)
        R_bcf2seu[0,1] = np.dot(np.transpose(self.s_polar_bcf_unit), self.j_bcf)
        R_bcf2seu[0,2] = np.dot(np.transpose(self.s_polar_bcf_unit), self.k_bcf)
        R_bcf2seu[1,0] = np.dot(np.transpose(self.e_polar_bcf_unit), self.i_bcf)
        R_bcf2seu[1,1] = np.dot(np.transpose(self.e_polar_bcf_unit), self.j_bcf)
        R_bcf2seu[1,2] = np.dot(np.transpose(self.e_polar_bcf_unit), self.k_bcf)
        R_bcf2seu[2,0] = np.dot(np.transpose(self.u_bcf_unit), self.i_bcf)
        R_bcf2seu[2,1] = np.dot(np.transpose(self.u_bcf_unit), self.j_bcf)
        R_bcf2seu[2,2] = np.dot(np.transpose(self.u_bcf_unit), self.k_bcf)

        self.R_bcf2seu_polar = R_bcf2seu
        return

    def calc_seu_polar_rotmat(self):
        """Calculate a south-east-up coordinate system in which the north vector points toward the pole using rotation matrices rather than cross products"""
        ############
        # DEPRECATED
        ############

        # rotation matrix from BCF frame to Fprime frame
        self.calc_longitude_bodycentric()
        R_bcf2fprime = self.m_u.transform3(self.longitude_bodycentric)
        R_fprime2bcf = np.transpose(R_bcf2fprime)
        
        # rotation matrix from Fprime frame to Fdprime frame
        self.calc_latitude_bodydetic()
        R_fprime2fdprime = self.m_u.transform2(-self.latitude_bodydetic)
        R_fdprime2fprime = np.transpose(R_fprime2fdprime)

        # combine the rotation matrices
        R_fdprime2bcf = np.dot(R_fprime2bcf, R_fdprime2fprime)

        # north vector
        self.n_polar_bcf_unit_rotmat = np.dot(R_fdprime2bcf, self.k_bcf)

        # east
        self.e_pseudo_bcf_unit_rotmat = np.dot(R_fdprime2bcf, self.j_bcf)

        # up
        self.u_bcf_unit_rotmat = np.dot(R_fdprime2bcf, self.i_bcf)

        # this should be rproj in BCF

        return

    def calc_heading_bcf_polar(self):
        """Calculate heading angle (rad)"""
        # angle between v_bcf projection in plane tangent to ellipsoid and some reference direction
        # measured from s_polar_bcf_unit toward e_polar_bcf_unit
        self.calc_seu_polar()
        vbcf_expressed_in_seu_polar = np.dot(self.R_bcf2seu_polar, self.vbcf)
        self.heading_bcf_polar = self.m_u.atan2(vbcf_expressed_in_seu_polar[1], vbcf_expressed_in_seu_polar[0])
        return

    def calc_heading_bcf_topocentric(self):
        """Calculate heading angle (rad)"""
        # angle between v_bcf projection in plane tangent to ellipsoid and some reference direction
        # measured from s_topocentric_bcf_unit toward e_topocentric_bcf_unit
        self.calc_seu_topocentric()
        vbcf_expressed_in_seu_topocentric = np.dot(self.R_bcf2seu_topocentric, self.vbcf)
        self.heading_bcf_topocentric = self.m_u.atan2(vbcf_expressed_in_seu_topocentric[1], vbcf_expressed_in_seu_topocentric[0])
        return

    def calc_heading_bci_polar(self):
        """Calculate heading angle (rad)"""
        # angle between v_bci projection in plane tangent to ellipsoid and some reference direction
        # measured from s_polar_bcf_unit toward e_polar_bcf_unit
        self.calc_vbci_expressed_in_bcf()
        self.calc_seu_polar()
        vbci_expressed_in_seu_polar = np.dot(self.R_bcf2seu_polar, self.vbci_expressed_in_bcf)
        self.heading_bci_polar = self.m_u.atan2(vbci_expressed_in_seu_polar[1], vbci_expressed_in_seu_polar[0])
        return

    def calc_heading_bci_topocentric(self):
        """Calculate heading angle (rad)"""
        # angle between v_bci projection in plane tangent to ellipsoid and some reference direction
        # measured from s_topocentric_bcf_unit toward e_topocentric_bcf_unit
        self.calc_vbci_expressed_in_bcf()
        self.calc_seu_topocentric()
        vbci_expressed_in_seu_topocentric = np.dot(self.R_bcf2seu_topocentric, self.vbci_expressed_in_bcf)
        self.heading_bci_topocentric = self.m_u.atan2(vbci_expressed_in_seu_topocentric[1], vbci_expressed_in_seu_topocentric[0])
        return

    def calc_flight_path_angle_bcf(self):
        """Calculate flight path angle (rad)"""
        # angle between v_bcf and plane tangent to ellipsoid
        self.calc_vbcf_mag()
        self.calc_seu_polar()
        vbcf_expressed_in_seu_polar = np.dot(self.R_bcf2seu_polar, self.vbcf)
        self.gamma_bcf = np.arcsin(vbcf_expressed_in_seu_polar[2]/self.vbcf_mag) # same as atan
        self.gamma_bcf_atan = self.m_u.atan2(vbcf_expressed_in_seu_polar[2], (vbcf_expressed_in_seu_polar[0]**2 + vbcf_expressed_in_seu_polar[1]**2)**0.5) # same as asin
        return

    def calc_flight_path_angle_bci(self):
        """Calculate flight path angle (rad)"""
        # angle between v_bci and plane tangent to ellipsoid
        self.calc_vbci_mag()
        self.calc_seu_polar()
        vbci_expressed_in_seu_polar = np.dot(self.R_bcf2seu_polar, self.vbci_expressed_in_bcf)
        self.gamma_bci = np.arcsin(vbci_expressed_in_seu_polar[2]/self.vbci_mag)
        return

    def calc_kyle(self):
        """Calculate stuff using Kyle's methods"""
        xi1 = self.m_u.atan2(self.rbcf[1,0], self.rbcf[0,0]) # xi1
        s1 = np.sin(xi1)
        c1 = np.cos(xi1)
        xi2 = self.m_u.atan2(self.u_bcf_unit[0,0] * c1 + self.u_bcf_unit[1,0] * s1, self.u_bcf_unit[2,0]) # xi2
        s2 = np.sin(xi2)
        c2 = np.cos(xi2)
        xi3 = self.m_u.atan2(self.u_bcf_unit[0,0] * c1 * s2 + self.u_bcf_unit[1,0] * s1 * s2 + self.u_bcf_unit[2,0] * c2, -self.u_bcf_unit[0,0] * s1 + self.u_bcf_unit[1,0] * c1)
        s3 = np.sin(xi3)
        c3 = np.cos(xi3)
        self.south_kyle = np.array([[c1 * c2], [s1 * c2], [-s2]])
        self.westerly_kyle = np.array([[c1 * s2 * c3 + s1 * s3], [s1 * s2 * c3 - c1 * s3], [c2 * c3]])

        return

    def calc_everything(self):
        """Calculate everything"""
        self.calc_longitude_bodycentric()
        self.calc_longitude_bodydetic()
        self.calc_latitude_bodycentric()
        self.calc_latitude_bodydetic()
        self.calc_vbcf_mag()
        self.calc_vbci_mag()
        self.calc_seu_topocentric()
        self.calc_seu_polar()
        #self.calc_seu_polar_rotmat() # Deprecated
        self.calc_heading_bcf_polar()
        self.calc_heading_bcf_topocentric()
        self.calc_heading_bci_polar()
        self.calc_heading_bci_topocentric()
        self.calc_flight_path_angle_bcf()
        self.calc_flight_path_angle_bci()
        self.calc_kyle()

    def plot(self):
        """Plot stuff"""
        self.calc_everything() # need to calculate everything to make sure it exists

        fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
        ax = fig.add_subplot(111, projection='3d')

        # Radii corresponding to the coefficients of the ellipsoid
        rx, ry, rz = self.a, self.b, self.c

        # loop through all spherical angles:
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        # Cartesian coordinates that correspond to the spherical angles:
        # (this is the equation of an ellipsoid):
        x = rx * np.outer(np.cos(u), np.sin(v))
        y = ry * np.outer(np.sin(u), np.sin(v))
        z = rz * np.outer(np.ones_like(u), np.cos(v))

        # create the surface plot of the ellipsoid
        ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='r')

        # Adjustment of the axes, so that they all have the same span:
        max_radius = max(rx, ry, rz)
        for axis in 'xyz':
            getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))

        # plot relevant planes
        xx, yy = np.meshgrid(np.arange(-self.a,self.a,0.1), np.arange(-self.b,self.b,0.1))

        # plane tangent to surface at r
        d = -np.dot(np.transpose(self.rbcf), self.u_bcf_unit)
        zz = (-self.u_bcf_unit[0] * xx - self.u_bcf_unit[1] * yy - d) * 1. / self.u_bcf_unit[2]
        ax.plot_surface(xx, yy, zz, alpha=0.2)

        # plane that passes through origin and north vector
        xx, zz = np.meshgrid(np.arange(-self.a,self.a,0.1), np.arange(-self.c*1.5,self.c*1.5,0.1))
        yy = (-self.e_pseudo_bcf_unit[0] * xx) / self.e_pseudo_bcf_unit[1]
        ax.plot_surface(xx, yy, zz, alpha=0.2)

        # now plot vectors
        n_quiver = 4 # number of vectors to show
        quiver_length = self.a # length of arrows
        quiver_x = np.zeros(n_quiver)
        quiver_y = np.zeros(n_quiver)
        quiver_z = np.zeros(n_quiver)

        # north vector starting at the north pole
        quiver_x[0] = 0.
        quiver_y[0] = 0.
        quiver_z[0] = self.c

        # all the other vectors start at the position
        for i in range(1,n_quiver):
            quiver_x[i] = self.rbcf[0]
            quiver_y[i] = self.rbcf[1]
            quiver_z[i] = self.rbcf[2]

        # set the directions of the vectors
        cp_vectors_x = [0., self.n_pa[0], self.e_bcf_topocentric_unit[0], -self.s_polar_bcf_unit[0]]
        cp_vectors_y = [0., self.n_pa[1], self.e_bcf_topocentric_unit[1], -self.s_polar_bcf_unit[1]]
        cp_vectors_z = [1., self.n_pa[2], self.e_bcf_topocentric_unit[2], -self.s_polar_bcf_unit[2]]

        # create the quiver
        q = ax.quiver(quiver_x, quiver_y, quiver_z, cp_vectors_x, cp_vectors_y, cp_vectors_z, length=quiver_length, normalize=False)

        # velocity vector can be either into the ellipse or out of the ellipse depending on v dot n

        # velocity w.r.t. BCF frame
        if np.dot(np.transpose(self.vbcf), self.u_bcf) >= 0: # outgoing velocity
            q_vbcf = ax.quiver(quiver_x[1], quiver_y[1], quiver_z[1], self.vbcf[0], self.vbcf[1], self.vbcf[2], length=quiver_length, normalize=True)
        else:
            q_vbcf = ax.quiver(quiver_x[1], quiver_y[1], quiver_z[1], self.vbcf[0], self.vbcf[1], self.vbcf[2], length=quiver_length, normalize=True, pivot='tip')

        # label the axes
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        ax.axis('equal')

        # show the plot
        #plt.show()




