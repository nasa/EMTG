# emtg_ellipsoid_math_tester.py

# Purpose
# Driver program for numerically testing derivations related to ellipsoidal interface

# Revision history
# Noble Hatten; 08/22/2018; Began
# Noble Hatten; 08/27/2018; Testing definition of south-east-up frame
# Noble Hatten; 08/30/2018; Added quiver element just pointing North through the z axis

# imports
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# personal imports
import math_utilities

def get_R_bcf2pa(theta_1, theta_2, theta_3):
    m_u = math_utilities.math_utilities_class()

    R_bcf2fprime = m_u.transform3(theta_1)
    R_fprime2fdprime = m_u.transform1(theta_2)
    R_fdprime2pa = m_u.transform3(theta_3)
    R_bcf2pa = np.dot(R_fdprime2pa, np.dot(R_fprime2fdprime, R_bcf2fprime))
    return R_bcf2pa

def get_r_pa(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    m_u = math_utilities.math_utilities_class()

    R_bcf2pa = get_R_bcf2pa(theta_1, theta_2, theta_3)
    r_pa = np.dot(R_bcf2pa, rbcf)
    return r_pa

def get_n_pa(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    r_pa = get_r_pa(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
    n_pa = 2.0 * np.dot(ellipsmat, r_pa) # normal vector in PA frame (km)
    return n_pa

def get_longitude(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    m_u = math_utilities.math_utilities_class()

    # BCF longitude (rad)
    longitude = m_u.atan2(np.asscalar(rbcf[1,0]), np.asscalar(rbcf[0,0]))

    # derivatives
    d_atan2 = m_u.d_atan2(rbcf[1,0], rbcf[0,0])
    d_long_drbcf = np.zeros((1,3),dtype=np.complex128)
    d_long_drbcf[0,0:3] = np.array([np.asscalar(d_atan2[1]), np.asscalar(d_atan2[0]), 0.0])
    d_long_dvbcf = np.zeros((1,3))
    d_long_dvbcf[0,0:3] = np.array([0.0, 0.0, 0.0])
    d_long_dt = np.array([0.0])
    d_long_dx = np.zeros((1,6))
    d_long_dx = np.concatenate((d_long_drbcf, d_long_dvbcf), axis=1)
    #d_long_dx = np.concatenate((d_long_dx, d_long_dt), axis=0)
    return longitude, d_long_dx

def get_latitude(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    m_u = math_utilities.math_utilities_class()

    longitude, d_long_dx = get_longitude(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
    # BCF bodycentric latitude (rad)
    phi = m_u.atan2(rbcf[2,0], (rbcf[0,0]*np.cos(longitude) + rbcf[1,0]*np.sin(longitude)))

    # derivatives
    d_phi_d_rbcf = np.zeros((1,3), dtype=np.complex128)
    d_phi_d_arg = m_u.d_atan2(rbcf[2,0], rbcf[0,0]*np.cos(longitude) + rbcf[1,0]*np.sin(longitude)) 
    d_argx_d_rbcf = np.zeros((1,3), dtype=np.complex128)
    d_coslong_d_rbcf = np.array([-np.sin(longitude) * d_long_dx[0,0], -np.sin(longitude) * d_long_dx[0,1], 0.0])
    d_sinlong_d_rbcf = np.array([np.cos(longitude) * d_long_dx[0,0], np.cos(longitude) * d_long_dx[0,1], 0.0])
    d_argx_d_rbcf[0,0] = np.asscalar(np.cos(longitude) + rbcf[0,0] * d_coslong_d_rbcf[0] + rbcf[1,0] * d_sinlong_d_rbcf[0])
    d_argx_d_rbcf[0,1] = rbcf[0,0] * d_coslong_d_rbcf[1] + np.sin(longitude) + rbcf[1,0] * d_sinlong_d_rbcf[1]
    d_argx_d_rbcf[0,2] = 0.0
    d_argy_d_rbcf = np.array([0.0, 0.0, 1.0],dtype=np.complex128)
    d_phi_d_rbcf = np.asscalar(d_phi_d_arg[1]) * d_argx_d_rbcf + np.asscalar(d_phi_d_arg[0]) * d_argy_d_rbcf
    d_phi_d_vbcf = np.zeros((1,3),dtype=np.complex128)
    d_phi_dt = np.array([0.0],dtype=np.complex128)
    d_phi_dx = np.concatenate((d_phi_d_rbcf, d_phi_d_vbcf), axis=1)
    #d_phi_dx = np.concatenate((d_phi_dx, d_phi_dt), axis=1)
    return phi, d_phi_dx

def get_vbcf_mag(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    # magnitude of velocity vector w.r.t. BCF frame (km/s)
    vbcf_mag = (vbcf[0]**2 + vbcf[1]**2 + vbcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives

    # derivatives
    d_vbcf_mag_dx = np.zeros((1,6), dtype=np.complex128)
    d_vbcf_mag_dx[0,3:6] = (1.0/vbcf_mag) * np.transpose(vbcf)
    return vbcf_mag, d_vbcf_mag_dx

def get_vbci_in_bcf(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    m_u = math_utilities.math_utilities_class()

    omega = get_omega(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
    # velocity vector w.r.t. BCI frame (km/s). vector is still EXPRESSED in the BCF frame
    vbci_in_bcf = vbcf + np.cross(omega, rbcf, axis=0)

    # derivatives
    d_vbci_in_bcf_dx = np.zeros((3,6), dtype=np.complex128)

    d_vbci_in_bcf_dx[0:3,0:3] = m_u.crossmat(omega) # d/drbcf
    d_vbci_in_bcf_dx[0:3,3:6] = np.diag([1.0, 1.0, 1.0]) # d/dvbcf
    return vbci_in_bcf, d_vbci_in_bcf_dx

def get_vbci_in_bcf_mag(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    vbci_in_bcf, d_vbci_in_bcf_dx = get_vbci_in_bcf(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
    vbci_in_bcf_mag = (vbci_in_bcf[0]**2 + vbci_in_bcf[1]**2 + vbci_in_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives

    # derivatives
    d_vbci_in_bcf_mag_dx = np.zeros((1,6), dtype=np.complex128)
    d_vbci_in_bcf_mag_d_vbci_in_bcf = (1.0/vbci_in_bcf_mag) * np.transpose(vbci_in_bcf)
    d_vbci_in_bcf_mag_dx = np.dot(d_vbci_in_bcf_mag_d_vbci_in_bcf, d_vbci_in_bcf_dx)

    return vbci_in_bcf_mag, d_vbci_in_bcf_mag_dx

def get_w(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    w = w0 + wdot * t # current value of offset (rad)
    return w

def get_omega(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    omega = np.array([[0.0], [0.0], [wdot]]) # time derivative of w (rad/s) (for now)
    return omega

def get_ellipsmat(a, b, c):
    a2 = a*a
    b2 = b*b
    c2 = c*c
    ellipsmat = np.diagflat([1./a2, 1./b2, 1./c2])
    return ellipsmat

def get_heading_bcf(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    m_u = math_utilities.math_utilities_class()

    xi1, d_long_dx = get_longitude(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # aux. angle xi1 (rad)

    # transformation matrix from BCF frame to F1 frame
    R_bcf2f1 = m_u.transform3(xi1)

    n_pa = get_n_pa(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
    R_bcf2pa = get_R_bcf2pa(theta_1, theta_2, theta_3)
    R_pa2bcf = np.transpose(R_bcf2pa)

    n_f1 = np.dot(R_bcf2f1, np.dot(R_pa2bcf, n_pa)) # normal vector expressed in F1 frame

    xi2 = np.asscalar(m_u.atan2(n_f1[0], n_f1[2])) # aux. angle xi2 (rad)

    # transformation matrix from F1 frame to F2 frame
    R_f12f2 = m_u.transform2(xi2)

    vbcf_in_f2 = np.dot(R_f12f2, np.dot(R_bcf2f1, vbcf)) # velocity w.r.t. BCF frame EXPRESSED in F2 frame (km/s)

    # heading angle using velocity with respect to the BCF frame (rad)
    heading_bcf = np.asscalar(m_u.atan2(vbcf_in_f2[1], vbcf_in_f2[0]))
    return heading_bcf

def get_heading_bci(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    m_u = math_utilities.math_utilities_class()

    xi1, d_long_dx = get_longitude(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # aux. angle xi1 (rad)

    # transformation matrix from BCF frame to F1 frame
    R_bcf2f1 = m_u.transform3(xi1)

    n_pa = get_n_pa(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
    R_bcf2pa = get_R_bcf2pa(theta_1, theta_2, theta_3)
    R_pa2bcf = np.transpose(R_bcf2pa)

    n_f1 = np.dot(R_bcf2f1, np.dot(R_pa2bcf, n_pa)) # normal vector expressed in F1 frame

    xi2 = np.asscalar(m_u.atan2(n_f1[0], n_f1[2])) # aux. angle xi2 (rad)

    # transformation matrix from F1 frame to F2 frame
    R_f12f2 = m_u.transform2(xi2)

    vbci_in_f2 = np.dot(R_f12f2, np.dot(R_bcf2f1, vbci_in_bcf)) # velocity w.r.t. BCI frame EXPRESSED in F2 frame (km/s)

    # heading angle using velocity with respect to the BCI frame (rad)
    heading_bci = np.asscalar(m_u.atan2(vbci_in_f2[1], vbci_in_f2[0]))
    return heading_bci

def get_gamma_bcf(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    m_u = math_utilities.math_utilities_class()

    xi1, d_long_dx = get_longitude(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # aux. angle xi1 (rad)

    # transformation matrix from BCF frame to F1 frame
    R_bcf2f1 = m_u.transform3(xi1)

    n_pa = get_n_pa(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
    R_bcf2pa = get_R_bcf2pa(theta_1, theta_2, theta_3)
    R_pa2bcf = np.transpose(R_bcf2pa)

    n_f1 = np.dot(R_bcf2f1, np.dot(R_pa2bcf, n_pa)) # normal vector expressed in F1 frame

    xi2 = np.asscalar(m_u.atan2(n_f1[0], n_f1[2])) # aux. angle xi2 (rad)

    # transformation matrix from F1 frame to F2 frame
    R_f12f2 = m_u.transform2(xi2)

    vbcf_in_f2 = np.dot(R_f12f2, np.dot(R_bcf2f1, vbcf)) # velocity w.r.t. BCF frame EXPRESSED in F2 frame (km/s)

    # heading angle using velocity with respect to the BCF frame (rad)
    heading_bcf = np.asscalar(m_u.atan2(vbcf_in_f2[1], vbcf_in_f2[0]))

    vbcf_in_f2 = np.dot(R_f12f2, np.dot(R_bcf2f1, vbcf)) # velocity w.r.t. BCF frame EXPRESSED in F2 frame (km/s)

    # transformation matrix from F2 frame to F3 frame, using BCF heading angle
    R_f22f3_bcf = m_u.transform3(heading_bcf)

    # velocity w.r.t. BCF frame EXPRESSED in F3 frame (km/s)
    vbcf_in_f3 = np.dot(R_f22f3_bcf, vbcf_in_f2)

    # flight path angle using velocty w.r.t. BCF frame (rad)
    gamma_bcf = np.asscalar(m_u.atan2(vbcf_in_f3[2], vbcf_in_f3[0]))
    return gamma_bcf

def get_gamma_bci(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    m_u = math_utilities.math_utilities_class()

    xi1, d_long_dx = get_longitude(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # aux. angle xi1 (rad)

    # transformation matrix from BCF frame to F1 frame
    R_bcf2f1 = m_u.transform3(xi1)

    n_pa = get_n_pa(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
    R_bcf2pa = get_R_bcf2pa(theta_1, theta_2, theta_3)
    R_pa2bcf = np.transpose(R_bcf2pa)

    n_f1 = np.dot(R_bcf2f1, np.dot(R_pa2bcf, n_pa)) # normal vector expressed in F1 frame

    xi2 = np.asscalar(m_u.atan2(n_f1[0], n_f1[2])) # aux. angle xi2 (rad)

    # transformation matrix from F1 frame to F2 frame
    R_f12f2 = m_u.transform2(xi2)

    vbci_in_f2 = np.dot(R_f12f2, np.dot(R_bcf2f1, vbci_in_bcf)) # velocity w.r.t. BCI frame EXPRESSED in F2 frame (km/s)

    # heading angle using velocity with respect to the BCI frame (rad)
    heading_bci = np.asscalar(m_u.atan2(vbci_in_f2[1], vbci_in_f2[0]))

    vbci_in_f2 = np.dot(R_f12f2, np.dot(R_bcf2f1, vbci_in_bcf)) # velocity w.r.t. BCI frame EXPRESSED in F2 frame (km/s)

    # transformation matrix from F2 frame to F3 frame, using BCI heading angle
    R_f22f3_bci = m_u.transform3(heading_bci)

    # velocity w.r.t. BCI frame EXPRESSED in F3 frame (km/s)
    vbci_in_f3 = np.dot(R_f22f3_bci, vbci_in_f2)

    # flight path angle using velocty w.r.t. BCI frame (rad)
    gamma_bci = np.asscalar(m_u.atan2(vbci_in_f3[2], vbci_in_f3[0]))
    return gamma_bci

def get_bcf2seu_mat_crossp(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    """calculate the transformation matrix from BCF frame to SEU frame using cross products to define orthogonal vectors"""
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    m_u = math_utilities.math_utilities_class()

    # BCF frame unit vectors, expressed in BCF frame
    i_bcf = m_u.i_unit()
    j_bcf = m_u.j_unit()
    k_bcf = m_u.k_unit()

    # ellipsoid normal vector in PA frame
    n_pa = get_n_pa(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
    n_pa_mag = (n_pa[0]**2 + n_pa[1]**2 + n_pa[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
    n_pa_unit = n_pa / n_pa_mag # unit vector

    # transformation matrix between BCF and PA frames
    R_bcf2pa = get_R_bcf2pa(theta_1, theta_2, theta_3)
    R_pa2bcf = np.transpose(R_bcf2pa)

    # u (up) vector in BCF frame, i.e., n_bcf
    u_bcf = np.dot(R_pa2bcf, n_pa)
    u_bcf_mag = (u_bcf[0]**2 + u_bcf[1]**2 + u_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
    u_bcf_unit = u_bcf / u_bcf_mag # unit vector

    # e (east) vector in BCF frame
    firstarg = m_u.crossmat(k_bcf)
    secondarg = 2.0 * np.dot(R_pa2bcf, np.dot(ellipsmat, np.dot(R_bcf2pa, rbcf)))
    e_bcf = np.dot(firstarg, secondarg)
    e_bcf_mag = (e_bcf[0]**2 + e_bcf[1]**2 + e_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
    e_bcf_unit = e_bcf / e_bcf_mag # unit vector

    # s (south) vector in BCF frame
    e_bcf_cross = m_u.crossmat(e_bcf)
    s_bcf = np.dot(e_bcf_cross, u_bcf)
    s_bcf_mag = (s_bcf[0]**2 + s_bcf[1]**2 + s_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
    s_bcf_unit = s_bcf / s_bcf_mag # unit vect

    # s (south) vector in BCF frame that actually points toward the pole; not necessarily perpendicular to east as defined in this routine
    firstarg = m_u.crossmat(k_bcf)
    secondarg = rbcf
    eps_bcf = np.dot(firstarg, secondarg)
    eps_bcf_mag = (eps_bcf[0]**2 + eps_bcf[1]**2 + eps_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
    eps_bcf_unit = eps_bcf / eps_bcf_mag # unit vector
    eps_bcf_cross = m_u.crossmat(eps_bcf)
    s_pole_bcf = np.dot(eps_bcf_cross, u_bcf)
    s_pole_bcf_mag = (s_pole_bcf[0]**2 + s_pole_bcf[1]**2 + s_pole_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
    s_pole_bcf_unit = s_pole_bcf / s_pole_bcf_mag # unit vect

    # transformation matrix
    R_bcf2seu = np.zeros((3,3))

    R_bcf2seu[0,0] = np.dot(np.transpose(s_bcf_unit), i_bcf)
    R_bcf2seu[0,1] = np.dot(np.transpose(s_bcf_unit), j_bcf)
    R_bcf2seu[0,2] = np.dot(np.transpose(s_bcf_unit), k_bcf)
    R_bcf2seu[1,0] = np.dot(np.transpose(e_bcf_unit), i_bcf)
    R_bcf2seu[1,1] = np.dot(np.transpose(e_bcf_unit), j_bcf)
    R_bcf2seu[1,2] = np.dot(np.transpose(e_bcf_unit), k_bcf)
    R_bcf2seu[2,0] = np.dot(np.transpose(u_bcf_unit), i_bcf)
    R_bcf2seu[2,1] = np.dot(np.transpose(u_bcf_unit), j_bcf)
    R_bcf2seu[2,2] = np.dot(np.transpose(u_bcf_unit), k_bcf)

    print("Unit vectors from cross product method")
    print('s: ', s_bcf_unit)
    print('s_pole: ', s_pole_bcf_unit)
    print('e: ', e_bcf_unit)
    print('eps: ', eps_bcf_unit)
    print('u: ', u_bcf_unit)
    print('')

    return R_bcf2seu, e_bcf_unit, s_bcf_unit, s_pole_bcf_unit

def get_bcf2seu_mat_rotations(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot):
    """calculate the transformation matrix from BCF frame to SEU frame using rotations to define orthogonal vectors"""
    rbcf = x_bcf[0:3]
    vbcf = x_bcf[3:6]
    m_u = math_utilities.math_utilities_class()

    # ellipsoid normal vector in PA frame
    n_pa = get_n_pa(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
    n_pa_mag = (n_pa[0]**2 + n_pa[1]**2 + n_pa[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
    n_pa_unit = n_pa / n_pa_mag # unit vector

    # transformation matrix between BCF and PA frames
    R_bcf2pa = get_R_bcf2pa(theta_1, theta_2, theta_3)
    R_pa2bcf = np.transpose(R_bcf2pa)

    # u (up) vector in BCF frame, i.e., n_bcf
    u_bcf = np.dot(R_pa2bcf, n_pa)
    u_bcf_mag = (u_bcf[0]**2 + u_bcf[1]**2 + u_bcf[2]**2)**0.5 # using this instead of norm() for the benefit of CX derivatives
    u_bcf_unit = u_bcf / u_bcf_mag # unit vector

    # get longitude (= xi1)
    longitude, d_long_dx = get_longitude(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)

    # transformation matrix from BCF frame to F1 frame
    R_bcf2f1 = m_u.transform3(longitude)

    # up/normal vector in F1 frame
    u_f1 = np.dot(R_bcf2f1, u_bcf)

    # angle xi2 between F1 and F2 frames
    xi2 = m_u.atan2(np.asscalar(u_f1[0]), np.asscalar(u_f1[2]))

    # transformation matrix from F1 frame to F2 frame
    R_f12f2 = m_u.transform2(xi2)

    # up/normal vector in F2 frame
    u_f2 = np.dot(R_f12f2, u_f1)

    # angle xi3 between F2 and F3 frames
    xi3 = m_u.atan2(np.asscalar(u_f2[1]), np.asscalar(u_f2[2]))

    # transformation matrix from F2 frame to F3 frame
    R_f22f3 = m_u.transform1(xi3)

    # combine rotation matrices
    R_bcf2seu = np.dot(R_f22f3, np.dot(R_f12f2, R_bcf2f1))

    # unit vectors of SEU in BCF coordinates just for testing

    # s is i unit vector in F3
    s_bcf_unit = np.dot(R_bcf2seu, m_u.i_unit())

    # e is j unit vector in F3
    e_bcf_unit = np.dot(R_bcf2seu, m_u.j_unit())

    #print("Unit vectors from rotation matrix method")
    #print('s: ', s_bcf_unit)
    #print('e: ', e_bcf_unit)
    #print('u: ', u_bcf_unit)
    #print('')

    return R_bcf2seu



# class instantiation
m_u = math_utilities.math_utilities_class()

# Parameters defining ellipsoid (km)
#a = 6378.14 # semimajor axis (Earth)
#b = 6378.14 # semiminaor axis (Earth)
#c = 6356.75 # semiintermediate axis (Earth)
a = 1.7
b = 1.1
c = 0.8

ellipsmat = get_ellipsmat(a, b, c) # get matrix of ellipsoid parameters

# Euler angles between BCF frame and PA frame in rad (angles assumed to be constant)
#theta_1 = np.deg2rad(40.)
#theta_2 = np.deg2rad(70.)
#theta_3 = np.deg2rad(10.)
theta_1 = np.deg2rad(0.)
theta_2 = np.deg2rad(0.)
theta_3 = np.deg2rad(0.)

R_bcf2pa = get_R_bcf2pa(theta_1, theta_2, theta_3) # transformation matrix between BCF and PA frames

r_pa = np.zeros((3,1))
r_pa[0] = 0.1*a
#r_pa[1] = 0.0 #0.4*b
r_pa[2] = 0.8*c
r_pa[1] = np.sqrt(b**2 * (1.0 - r_pa[0]**2/a**2 - r_pa[2]**2/c**2))
#r_pa[2] = np.sqrt(c**2 * (1.0 - r_pa[0]**2/a**2 - r_pa[1]**2/b**2))

# complex step
h = 1.0e-16

# time (s)
t = 0.0

# Position vector in BCF frame (km)
#rbcf = np.array([[-1431.734007], [1486.260530], [6014.813181]]) # (Earth)
#rbcf = np.array([[400.0], [100.0], [-50.0]])
rbcf = np.dot(np.transpose(R_bcf2pa), r_pa)
print('r_bcf: ', rbcf)

# Velocity vector w.r.t. BCF frame, expressed in BCF frame (km/s)
vbcf = np.array([[0.4], [-1.1], [-0.6]])

x_bcf = np.concatenate((rbcf, vbcf), axis=0)

# Definition of rotation of BCF frame w.r.t. BCI frame
# Assumed to be rotation about z axis only (BCI z = BCF z)
w0 = np.deg2rad(15.) # initial offset (rad)
wdot = np.deg2rad(1.e-2) # linear rate of change of offset (rad/s)
w = w0 + wdot * t # current value of offset (rad)
omega = np.array([[0.0], [0.0], [wdot]]) # time derivative of w (rad/s) (for now)
r_pa = get_r_pa(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # transform position vector to PA frame (km)
omega = get_omega(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # get angular velocity of BCF frame w.r.t. BCI frame (rad/s)
n_pa = get_n_pa(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # ellipsoid normal vector in PA frame (km)
longitude, d_long_dx = get_longitude(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # BCF longitude (rad)
phi, d_phi_dx = get_latitude(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # BCF bodycentric latitude (rad)
vbcf_mag, d_vbcf_mag_dx = get_vbcf_mag(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # magnitude of velocity vector w.r.t. BCF frame (km/s)
vbci_in_bcf, d_vbci_in_bcf_dx = get_vbci_in_bcf(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # velocity vector w.r.t. BCI frame (km/s). vector is still EXPRESSED in the BCF frame
vbci_in_bcf_mag, d_vbci_in_bcf_mag_dx = get_vbci_in_bcf_mag(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # magnitude of velocity vector w.r.t. BCI frame (km/s)

R_bcf2seu_crossp, e_bcf_cp, s_bcf_cp, s_pole_bcf_cp = get_bcf2seu_mat_crossp(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
R_bcf2seu_rotations = get_bcf2seu_mat_rotations(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)

############
# plot stuff
############
fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
ax = fig.add_subplot(111, projection='3d')

# Radii corresponding to the coefficients of the ellipsoid
rx, ry, rz = a, b, c

# loop through all spherical angles:
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

# Cartesian coordinates that correspond to the spherical angles:
# (this is the equation of an ellipsoid):
x = rx * np.outer(np.cos(u), np.sin(v))
y = ry * np.outer(np.sin(u), np.sin(v))
z = rz * np.outer(np.ones_like(u), np.cos(v))

# rotate coordinate system so that the axes of the plot are BCF (if no rotate, then the axes are axes of PA frame)

# create the surface plot of the ellipsoid
ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='r')

# Adjustment of the axes, so that they all have the same span:
max_radius = max(rx, ry, rz)
for axis in 'xyz':
    getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))

n_quiver = 4 # number of vectors to show
quiver_x = np.zeros(n_quiver)
quiver_y = np.zeros(n_quiver)
quiver_z = np.zeros(n_quiver)

# north vector starting at the north pole
quiver_x[0] = 0.
quiver_y[0] = 0.
quiver_z[0] = c

# all the other vectors start at the position
for i in range(1,n_quiver):
    quiver_x[i] = rbcf[0]
    quiver_y[i] = rbcf[1]
    quiver_z[i] = rbcf[2]

# set the directions of the vectors
cp_vectors_x = [0., n_pa[0], e_bcf_cp[0], -s_pole_bcf_cp[0]]
cp_vectors_y = [0., n_pa[1], e_bcf_cp[1], -s_pole_bcf_cp[1]]
cp_vectors_z = [1., n_pa[2], e_bcf_cp[2], -s_pole_bcf_cp[2]]

# create the quiver
q = ax.quiver(quiver_x, quiver_y, quiver_z, cp_vectors_x, cp_vectors_y, cp_vectors_z, length=a, normalize=True)
#q_pole = ax.quiver(quiver_x[1], quiver_y[1], quiver_z[1], -s_pole_bcf_cp[0], -s_pole_bcf_cp[1], -s_pole_bcf_cp[2], length=a, normalize=True)

# label the axes
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# show the plot
plt.show()

###################################################
# complex-step derivative of longitude and latitude
###################################################
#d_long_dx_cx = np.zeros((1,6))
#d_phi_dx_cx = np.zeros((1,6))
#d_vbcf_mag_dx_cx = np.zeros((1,6))
#d_vbci_in_bcf_mag_dx_cx = np.zeros((1,6))
#x_pert = np.zeros((1,6),dtype=np.complex128)
#for i in range(0,6):
#    x_pert = np.complex128(x_bcf)
#    x_pert[i] = x_pert[i] + h*1j

#    # longitude
#    pert, scratch = get_longitude(x_pert, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
#    d_long_dx_cx[0,i] = np.imag(pert)/h

#    # latitude
#    pert, scratch = get_latitude(x_pert, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
#    d_phi_dx_cx[0,i] = np.imag(pert)/h

#    # velocity magnitude w.r.t body-fixed frame
#    pert, scratch = get_vbcf_mag(x_pert, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
#    d_vbcf_mag_dx_cx[0,i] = np.imag(pert)/h

#    # velocity magnitude w.r.t body-inertial frame, expressed in BCF coordinates
#    pert, scratch = get_vbci_in_bcf_mag(x_pert, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot)
#    d_vbci_in_bcf_mag_dx_cx[0,i] = np.imag(pert)/h


#heading_bcf = get_heading_bcf(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # heading angle using velocity with respect to the BCF frame (rad)
#heading_bci = get_heading_bci(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # heading angle using velocity with respect to the BCI frame (rad)
#gamma_bcf = get_gamma_bcf(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # flight path angle using velocty w.r.t. BCF frame (rad)
#gamma_bci = get_gamma_bci(x_bcf, t, theta_1, theta_2, theta_3, ellipsmat, w0, wdot) # flight path angle using velocty w.r.t. BCI frame (rad)

