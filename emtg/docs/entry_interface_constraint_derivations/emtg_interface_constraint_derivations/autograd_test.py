# autograd_test.py

# Purpose
# Testing the autograd algorithmic differentiation package

# Revision history
# Noble Hatten; 09/13/2018; Began


import autograd
import autograd.numpy as np
from autograd import grad
from autograd import jacobian
import ellipsoid_math_autograd


# test ellipsoid stuff

# set parameters of ellipsoid
a = 10.0
b = 8.0
c = 6.0

# set Euler angles relating BCF frame to PA frame (3-1-3)
# Note on comparisons with STK:
# The Body_Axes_Body_Fixed system was created as a "Fixed in Axes" type set of axes.
# The reference axes is principal axis frame.
# Theta1 corresponds to Euler Angle A
# Theta2 corresponds to Euler Angle B
# Theta3 corresponds to Euler Angle C
# However, all the signs of the angles need to be reversed because the rotations are going the other direction.
theta1_0 = -21.0 # deg
theta2_0 = -9.0 # deg
theta3_0 = -14.0 # deg
theta1_0rad = np.deg2rad(theta1_0) # rad
theta2_0rad = np.deg2rad(theta2_0)
theta3_0rad = np.deg2rad(theta3_0)

# assume the theta's can vary linearly in time with constant rates
theta1dot_0 = 0.2 # deg/s
theta2dot_0 = 0.1 # deg/s
theta3dot_0 = 0.4 # deg/s
#theta1dot_0 = 0. # deg/s
#theta2dot_0 = 0. # deg/s
#theta3dot_0 = 0. # deg/s
theta1dot_0rad = np.deg2rad(theta1dot_0) # rad/s
theta2dot_0rad = np.deg2rad(theta2dot_0)
theta3dot_0rad = np.deg2rad(theta3dot_0)

# calculate a point on that ellipsoid
r_pa = np.zeros((3,1))
r_pa[0] = 0.2*a
r_pa[1] = -0.4*b
#r_pa[2] = 0.5*c
#r_pa[0] = np.sqrt(a**2 * (1.0 - r_pa[1]**2/b**2 - r_pa[2]**2/c**2))
#r_pa[1] = np.sqrt(b**2 * (1.0 - r_pa[0]**2/a**2 - r_pa[2]**2/c**2))
r_pa[2] = np.sqrt(c**2 * (1.0 - r_pa[0]**2/a**2 - r_pa[1]**2/b**2))

# set velocity in BCF frame
vbcf = np.array([[0.2], [-0.1], [-0.4]])

# set rotation of BCF w.r.t. BCI
w0 = 18.0 # deg
w0dot = 1.5 # deg/s
w0rad = np.deg2rad(w0) # rad
w0dotrad = np.deg2rad(w0dot) # rad

# set the time
t = 0.4 # s

# create ellipsoid object
ellipse = ellipsoid_math_autograd.ellipsoid_math_autograd_class()

# calculate r_bcf from r_pa
rbcf = ellipse.calc_rbcf(t, r_pa, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad) # calculate a function from the object

# calculate bodycentric longitude
longitude_bodycentric = ellipse.calc_longitude_bodycentric(rbcf) # value
d_lon_centric_d_rbcf_fun = jacobian(ellipse.calc_longitude_bodycentric, 0) # jacobian function w.r.t. rbcf
d_lon_centric_d_rbcf_ad = d_lon_centric_d_rbcf_fun(rbcf) # d/drbcf AD
d_lon_centric_d_t, d_lon_centric_d_rbcf, d_lon_centric_d_vbcf = ellipse.calc_longitude_bodycentric_derivs(rbcf) # d/dx analytical

# calculate bodydetic longitude
longitude_bodydetic = ellipse.calc_longitude_bodydetic(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # value
d_lon_detic_d_t_fun = jacobian(ellipse.calc_longitude_bodydetic, 0) # jacobian function w.r.t. t
d_lon_detic_d_rbcf_fun = jacobian(ellipse.calc_longitude_bodydetic, 1) # jacobian function w.r.t. rbcf
d_lon_detic_d_t_ad = d_lon_detic_d_t_fun(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dt
d_lon_detic_d_rbcf_ad = d_lon_detic_d_rbcf_fun(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/drbcf
d_lon_detic_d_t, d_lon_detic_d_rbcf, d_lon_detic_d_vbcf = ellipse.calc_longitude_bodydetic_derivs(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dx analytical


# calculate bodycentric latitude
latitude_bodycentric = ellipse.calc_latitude_bodycentric(rbcf) # value
d_lat_centric_d_rbcf_fun = jacobian(ellipse.calc_latitude_bodycentric, 0) # jacobian function w.r.t. rbcf
d_lat_centric_d_rbcf_ad = d_lat_centric_d_rbcf_fun(rbcf) # d/drbcf AD
d_lat_centric_d_t, d_lat_centric_d_rbcf, d_lat_centric_d_vbcf = ellipse.calc_latitude_bodycentric_derivs(rbcf) # d/dx analytical

## test derivatives of R_bcf2pa
#theta1rad, theta2rad, theta3rad = ellipse.calc_thetas_linear(t, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad) # the thetas
#R_bcf2pa = ellipse.calc_R_bcf2pa(theta1rad, theta2rad, theta3rad) # value
#d_Rbcf2pa_d_theta1_fun = jacobian(ellipse.calc_R_bcf2pa, 0)
#d_Rbcf2pa_d_theta2_fun = jacobian(ellipse.calc_R_bcf2pa, 1)
#d_Rbcf2pa_d_theta3_fun = jacobian(ellipse.calc_R_bcf2pa, 2)
#d_Rbcf2pa_d_theta1_ad = d_Rbcf2pa_d_theta1_fun(theta1rad, theta2rad, theta3rad)
#d_Rbcf2pa_d_theta2_ad = d_Rbcf2pa_d_theta2_fun(theta1rad, theta2rad, theta3rad)
#d_Rbcf2pa_d_theta3_ad = d_Rbcf2pa_d_theta3_fun(theta1rad, theta2rad, theta3rad)
#d_Rbcf2pa_d_t_ad = d_Rbcf2pa_d_theta1_ad * theta1dot_0rad + d_Rbcf2pa_d_theta2_ad * theta2dot_0rad + d_Rbcf2pa_d_theta3_ad * theta3dot_0rad
#d_R_bcf2pa_d_t = ellipse.calc_d_R_bcf2pa_d_t(theta1rad, theta2rad, theta3rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad) # analytical derivatives

## test derivatives of u_bcf
#d_u_bcf_d_t, d_u_bcf_d_rbcf, d_u_bcf_d_vbcf = ellipse.calc_u_bcf_derivs(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # value
#du_dt_fun = jacobian(ellipse.calc_u_bcf, 0)
#du_dr_fun = jacobian(ellipse.calc_u_bcf, 1)
#du_dt_ad = du_dt_fun(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
#du_dr_ad = du_dr_fun(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)


# calculate bodydetic latitude
latitude_bodydetic = ellipse.calc_latitude_bodydetic(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # value
d_lat_detic_d_t_fun = jacobian(ellipse.calc_latitude_bodydetic, 0) # jacobian function w.r.t. t
d_lat_detic_d_rbcf_fun = jacobian(ellipse.calc_latitude_bodydetic, 1) # jacobian function w.r.t. rbcf
d_lat_detic_d_t_ad = d_lat_detic_d_t_fun(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dt
d_lat_detic_d_rbcf_ad = d_lat_detic_d_rbcf_fun(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/rbcf
d_lat_detic_d_t, d_lat_detic_d_rbcf, d_lat_detic_d_vbcf = ellipse.calc_latitude_bodydetic_derivs(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dx analytical

# calculate velocity magitude w.r.t. BCF frame
vbcf_mag = ellipse.calc_vbcf_mag(vbcf) # value
d_vbcf_mag_d_vbcf_fun = jacobian(ellipse.calc_vbcf_mag, 0) # jacobian function w.r.t. vbcf
d_vbcf_mag_d_vbcf_ad = d_vbcf_mag_d_vbcf_fun(vbcf) # d/dvbcf
d_vbcf_mag_d_t, d_vbcf_mag_d_rbcf, d_vbcf_mag_d_vbcf = ellipse.calc_vbcf_mag_derivs(vbcf) # d/dx analytical

# calculate velocity w.r.t. BCI frame expressed in BCF frame
vbci_in_bcf = ellipse.calc_vbci_in_bcf(rbcf, vbcf, w0dotrad) # value
d_vbci_in_bcf_d_rbcf_fun = jacobian(ellipse.calc_vbci_in_bcf, 0) # jacobian function w.r.t. rbcf
d_vbci_in_bcf_d_vbcf_fun = jacobian(ellipse.calc_vbci_in_bcf, 1) # jacobian function w.r.t. vbcf
d_vbci_in_bcf_d_rbcf_ad = d_vbci_in_bcf_d_rbcf_fun(rbcf, vbcf, w0dotrad) # d/drbcf AD
d_vbci_in_bcf_d_vbcf_ad = d_vbci_in_bcf_d_vbcf_fun(rbcf, vbcf, w0dotrad) # d/dvbcf AD
d_vbci_in_bcf_d_t, d_vbci_in_bcf_d_rbcf, d_vbci_in_bcf_d_vbcf = ellipse.calc_vbci_in_bcf_derivs(w0dotrad) # d/dx analytical

# calculate velocity magnitude w.r.t. BCI frame
vbci_mag = ellipse.calc_vbci_mag(rbcf, vbcf, w0dotrad) # value
d_vbci_mag_d_rbcf_fun = jacobian(ellipse.calc_vbci_mag, 0) # jacobian function w.r.t. rbcf
d_vbci_mag_d_vbcf_fun = jacobian(ellipse.calc_vbci_mag, 1) # jacobian function w.r.t. vbcf
d_vbci_mag_d_rbcf_ad = d_vbci_mag_d_rbcf_fun(rbcf, vbcf, w0dotrad) # d/drbcf AD
d_vbci_mag_d_vbcf_ad = d_vbci_mag_d_vbcf_fun(rbcf, vbcf, w0dotrad) # d/dvbcf AD
d_vbci_mag_d_t, d_vbci_mag_d_rbcf, d_vbci_mag_d_vbcf = ellipse.calc_vbci_mag_derivs(rbcf, vbcf, w0dotrad) # d/dx analytical

## test derivs of SEU topocentric unit vectors
#d_s_bcf_topocentric_unit_d_t_fun = jacobian(ellipse.calc_seu_topocentric_s, 0)
#d_s_bcf_topocentric_unit_d_r_fun = jacobian(ellipse.calc_seu_topocentric_s, 1)
#d_s_bcf_topocentric_unit_d_t_ad =  d_s_bcf_topocentric_unit_d_t_fun(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
#d_s_bcf_topocentric_unit_d_r_ad =  d_s_bcf_topocentric_unit_d_r_fun(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
#d_s_bcf_unit_d_t, d_s_bcf_unit_d_rbcf, d_s_bcf_unit_d_vbcf, d_e_bcf_unit_d_t, d_e_bcf_unit_d_rbcf, d_e_bcf_unit_d_vbcf, d_u_bcf_unit_d_t, d_u_bcf_unit_d_rbcf, d_u_bcf_unit_d_vbcf = ellipse.calc_seu_topocentric_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)

# test derivs of SEU topocentric transformation matrix
#d_R_d_t = ellipse.calc_R_bcf2seu_topocentric_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)

# test derivs of SEU polar unit vectors
#d_e_bcf_polar_unit_d_t_fun = jacobian(ellipse.calc_seu_polar_e, 0)
#d_e_bcf_polar_unit_d_r_fun = jacobian(ellipse.calc_seu_polar_e, 1)
#d_e_bcf_polar_unit_d_t_ad =  d_e_bcf_polar_unit_d_t_fun(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
#d_e_bcf_polar_unit_d_r_ad =  d_e_bcf_polar_unit_d_r_fun(t, rbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)
#d_s_bcf_unit_d_t, d_s_bcf_unit_d_rbcf, d_s_bcf_unit_d_vbcf, d_e_bcf_unit_d_t, d_e_bcf_unit_d_rbcf, d_e_bcf_unit_d_vbcf, d_u_bcf_unit_d_t, d_u_bcf_unit_d_rbcf, d_u_bcf_unit_d_vbcf = ellipse.calc_seu_polar_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c)


# calculate heading (topocentric, BCF)
heading_bcf_topocentric = ellipse.calc_heading_bcf_topocentric(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # value
d_heading_bcf_topocentric_d_t_fun = jacobian(ellipse.calc_heading_bcf_topocentric, 0) # jacobian function w.r.t. rbcf
d_heading_bcf_topocentric_d_rbcf_fun = jacobian(ellipse.calc_heading_bcf_topocentric, 1) # jacobian function w.r.t. rbcf
d_heading_bcf_topocentric_d_vbcf_fun = jacobian(ellipse.calc_heading_bcf_topocentric, 2) # jacobian function w.r.t. vbcf
d_heading_bcf_topocentric_d_t_ad = d_heading_bcf_topocentric_d_t_fun(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dt AD
d_heading_bcf_topocentric_d_rbcf_ad = d_heading_bcf_topocentric_d_rbcf_fun(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/drbcf AD
d_heading_bcf_topocentric_d_vbcf_ad = d_heading_bcf_topocentric_d_vbcf_fun(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dvbcf AD
d_heading_bcf_topocentric_d_t, d_heading_bcf_topocentric_d_rbcf, d_heading_bcf_topocentric_d_vbcf = ellipse.calc_heading_bcf_topocentric_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dx analytical

# calculate heading (polar, BCF)
heading_bcf_polar = ellipse.calc_heading_bcf_polar(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # value
d_heading_bcf_polar_d_t_fun = jacobian(ellipse.calc_heading_bcf_polar, 0) # jacobian function w.r.t. rbcf
d_heading_bcf_polar_d_rbcf_fun = jacobian(ellipse.calc_heading_bcf_polar, 1) # jacobian function w.r.t. rbcf
d_heading_bcf_polar_d_vbcf_fun = jacobian(ellipse.calc_heading_bcf_polar, 2) # jacobian function w.r.t. vbcf
d_heading_bcf_polar_d_t_ad = d_heading_bcf_polar_d_t_fun(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dt AD
d_heading_bcf_polar_d_rbcf_ad = d_heading_bcf_polar_d_rbcf_fun(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/drbcf AD
d_heading_bcf_polar_d_vbcf_ad = d_heading_bcf_polar_d_vbcf_fun(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dvbcf AD
d_heading_bcf_polar_d_t, d_heading_bcf_polar_d_rbcf, d_heading_bcf_polar_d_vbcf = ellipse.calc_heading_bcf_polar_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dx analytical

# calculate heading (topocentric BCI)
heading_bci_topocentric = ellipse.calc_heading_bci_topocentric(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # value
d_heading_bci_topocentric_d_t_fun = jacobian(ellipse.calc_heading_bci_topocentric, 0) # jacobian function w.r.t. rbcf
d_heading_bci_topocentric_d_rbcf_fun = jacobian(ellipse.calc_heading_bci_topocentric, 1) # jacobian function w.r.t. rbcf
d_heading_bci_topocentric_d_vbcf_fun = jacobian(ellipse.calc_heading_bci_topocentric, 2) # jacobian function w.r.t. vbcf
d_heading_bci_topocentric_d_t_ad = d_heading_bci_topocentric_d_t_fun(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dt AD
d_heading_bci_topocentric_d_rbcf_ad = d_heading_bci_topocentric_d_rbcf_fun(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/drbcf AD
d_heading_bci_topocentric_d_vbcf_ad = d_heading_bci_topocentric_d_vbcf_fun(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dvbcf AD
d_heading_bci_topocentric_d_t, d_heading_bci_topocentric_d_rbcf, d_heading_bci_topocentric_d_vbcf = ellipse.calc_heading_bci_topocentric_derivs(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dx analytical

# calculate heading (polar BCI)
heading_bci_polar = ellipse.calc_heading_bci_polar(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # value
d_heading_bci_polar_d_t_fun = jacobian(ellipse.calc_heading_bci_polar, 0) # jacobian function w.r.t. rbcf
d_heading_bci_polar_d_rbcf_fun = jacobian(ellipse.calc_heading_bci_polar, 1) # jacobian function w.r.t. rbcf
d_heading_bci_polar_d_vbcf_fun = jacobian(ellipse.calc_heading_bci_polar, 2) # jacobian function w.r.t. vbcf
d_heading_bci_polar_d_t_ad = d_heading_bci_polar_d_t_fun(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dt AD
d_heading_bci_polar_d_rbcf_ad = d_heading_bci_polar_d_rbcf_fun(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/drbcf AD
d_heading_bci_polar_d_vbcf_ad = d_heading_bci_polar_d_vbcf_fun(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dvbcf AD
d_heading_bci_polar_d_t, d_heading_bci_polar_d_rbcf, d_heading_bci_polar_d_vbcf = ellipse.calc_heading_bci_polar_derivs(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dx analytical

# calculate flight path angle (BCF)
gamma_bcf = ellipse.calc_flight_path_angle_bcf(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # value
d_gamma_bcf_d_t_fun = jacobian(ellipse.calc_flight_path_angle_bcf, 0) # jacobian function w.r.t. rbcf
d_gamma_bcf_d_rbcf_fun = jacobian(ellipse.calc_flight_path_angle_bcf, 1) # jacobian function w.r.t. rbcf
d_gamma_bcf_d_vbcf_fun = jacobian(ellipse.calc_flight_path_angle_bcf, 2) # jacobian function w.r.t. vbcf
d_gamma_bcf_d_t_ad = d_gamma_bcf_d_t_fun(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dt AD
d_gamma_bcf_d_rbcf_ad = d_gamma_bcf_d_rbcf_fun(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/drbcf AD
d_gamma_bcf_d_vbcf_ad = d_gamma_bcf_d_vbcf_fun(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dvbcf AD
d_gamma_bcf_d_t, d_gamma_bcf_d_rbcf, d_gamma_bcf_d_vbcf = ellipse.calc_flight_path_angle_bcf_derivs(t, rbcf, vbcf, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dx analytical

# calculate flight path angle (BCI)
gamma_bci = ellipse.calc_flight_path_angle_bci(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # value
d_gamma_bci_d_t_fun = jacobian(ellipse.calc_flight_path_angle_bci, 0) # jacobian function w.r.t. rbcf
d_gamma_bci_d_rbcf_fun = jacobian(ellipse.calc_flight_path_angle_bci, 1) # jacobian function w.r.t. rbcf
d_gamma_bci_d_vbcf_fun = jacobian(ellipse.calc_flight_path_angle_bci, 2) # jacobian function w.r.t. vbcf
d_gamma_bci_d_t_ad = d_gamma_bci_d_t_fun(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dt AD
d_gamma_bci_d_rbcf_ad = d_gamma_bci_d_rbcf_fun(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/drbcf AD
d_gamma_bci_d_vbcf_ad = d_gamma_bci_d_vbcf_fun(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dvbcf AD
d_gamma_bci_d_t, d_gamma_bci_d_rbcf, d_gamma_bci_d_vbcf = ellipse.calc_flight_path_angle_bci_derivs(t, rbcf, vbcf, w0dotrad, theta1_0rad, theta2_0rad, theta3_0rad, theta1dot_0rad, theta2dot_0rad, theta3dot_0rad, a, b, c) # d/dx analytical

# set what we want to print
print_all = 0
print_pos = 0
print_vel = 0
print_lat = 0
print_lon = 0
print_head = 1
print_fpa = 0

print('=========')
print('OUTPUT:::')
print('=========')
print('')

if (print_all == 1 or print_pos == 1):
    print('=========')
    print('Position:')
    print('=========')
    print('')
    print('r_pa = ', np.transpose(r_pa[:,0]))
    print('')
    print('r_bcf = ', np.transpose(rbcf[:,0]))
    print('')

if (print_all == 1 or print_vel == 1):
    print('=========')
    print('Velocity:')
    print('=========')
    print('')
    print('v_bcf in bcf = ', np.transpose(vbcf[:,0]))
    print('')
    print('v_bci in bci = ')
    print('')
    print('v_bci in bcf = ', np.transpose(vbci_in_bcf[:,0]))
    print('')
    print('d [v_bci in bcf] / d [r_bcf] analytical = ', d_vbci_in_bcf_d_rbcf)
    print('')
    print('d [v_bci in bcf] / d [v_bcf] analytical = ', d_vbci_in_bcf_d_vbcf)
    print('')
    print('d [v_bci in bcf] / d [r_bcf] AD = ', d_vbci_in_bcf_d_rbcf_ad)
    print('')
    print('d [v_bci in bcf] / d [v_bcf] AD = ', d_vbci_in_bcf_d_vbcf_ad)
    print('')
    print('v_bcf_mag = ', vbcf_mag)
    print('')
    print('d [v_bcf_mag] / d [v_bcf] analytical = ', d_vbcf_mag_d_vbcf)
    print('')
    print('d [v_bcf_mag] / d [v_bcf] AD = ', np.transpose(d_vbcf_mag_d_vbcf_ad[:,0]))
    print('')
    print('v_bci_mag = ', vbci_mag)
    print('')
    print('d [v_bci_mag] / d [r_bcf] analytical = ', d_vbci_mag_d_rbcf)
    print('')
    print('d [v_bci_mag] / d [v_bcf] analytical = ', d_vbci_mag_d_vbcf)
    print('')
    print('d [v_bci_mag] / d [r_bcf] AD = ', np.transpose(d_vbci_mag_d_rbcf_ad[:,0]))
    print('')
    print('d [v_bci_mag] / d [v_bcf] AD = ', np.transpose(d_vbci_mag_d_vbcf_ad[:,0]))
    print('')
    d_vbci_in_bcf_d_rbcf_ad
    

if (print_all == 1 or print_lat == 1):
    print('=========')
    print('Latitude:')
    print('=========')
    print('')
    print('bodycentric latitude (rad) = ', latitude_bodycentric)
    print('')
    print('d [bodycentric latitude] / d [rbcf] analytical = ', d_lat_centric_d_rbcf)
    print('')
    print('d [bodycentric latitude] / d [rbcf] AD = ', np.transpose(d_lat_centric_d_rbcf_ad[:,0]))
    print('')
    print('bodydetic latitude (rad) = ', latitude_bodydetic)
    print('')
    print('d [bodydetic latitude] / d [t] analytical = ', d_lat_detic_d_t[0,0])
    print('')
    print('d [bodydetic latitude] / d [rbcf] analytical = ', d_lat_detic_d_rbcf[0,:])
    print('')
    print('d [bodydetic latitude] / d [t] AD = ', d_lat_detic_d_t_ad)
    print('')
    print('d [bodydetic latitude] / d [rbcf] AD = ', np.transpose(d_lat_detic_d_rbcf_ad[:,0]))
    print('')

if (print_all == 1 or print_lon == 1):
    print('==========')
    print('Longitude:')
    print('==========')
    print('')
    print('bodycentric longitude (rad) = ', longitude_bodycentric)
    print('')
    print('d [bodycentric longitude] / d [rbcf] analytical = ', d_lon_centric_d_rbcf)
    print('')
    print('d [bodycentric longitude] / d [rbcf] AD = ', np.transpose(d_lon_centric_d_rbcf_ad[:,0]))
    print('')
    print('bodydetic longitude (rad) = ', longitude_bodydetic)
    print('')
    print('d [bodydetic longitude] / d [t] analytical = ', d_lon_detic_d_t[0,0])
    print('')
    print('d [bodydetic longitude] / d [rbcf] analytical = ', d_lon_detic_d_rbcf[0,:])
    print('')
    print('d [bodydetic longitude] / d [t] AD = ', d_lon_detic_d_t_ad)
    print('')
    print('d [bodydetic longitude] / d [rbcf] AD = ', np.transpose(d_lon_detic_d_rbcf_ad[:,0]))
    print('')
    d_lon_detic_d_t
    print('')

if (print_all == 1 or print_head == 1):
    print('========')
    print('Heading:')
    print('========')
    print('')
    print('BCF topocentric heading (rad) = ', heading_bcf_topocentric)
    print('')
    print('d [ BCF topocentric heading] / d [t] (rad) analytical = ', d_heading_bcf_topocentric_d_t)
    print('')
    print('d [ BCF topocentric heading] / d [rbcf] (rad) analytical = ', d_heading_bcf_topocentric_d_rbcf)
    print('')
    print('d [ BCF topocentric heading] / d [vbcf] (rad) analytical = ', d_heading_bcf_topocentric_d_vbcf)
    print('')
    print('d [ BCF topocentric heading] / d [t] (rad) AD = ', d_heading_bcf_topocentric_d_t_ad)
    print('')
    print('d [ BCF topocentric heading] / d [rbcf] (rad) AD = ', np.transpose(d_heading_bcf_topocentric_d_rbcf_ad[:,0]))
    print('')
    print('d [ BCF topocentric heading] / d [vbcf] (rad) AD = ', np.transpose(d_heading_bcf_topocentric_d_vbcf_ad[:,0]))
    print('')
    print('BCF polar heading (rad) = ', heading_bcf_polar)
    print('')
    print('d [ BCF polar heading] / d [t] (rad) analytical = ', d_heading_bcf_polar_d_t)
    print('')
    print('d [ BCF polar heading] / d [rbcf] (rad) analytical = ', d_heading_bcf_polar_d_rbcf)
    print('')
    print('d [ BCF polar heading] / d [vbcf] (rad) analytical = ', d_heading_bcf_polar_d_vbcf)
    print('')
    print('d [ BCF polar heading] / d [t] (rad) AD = ', d_heading_bcf_polar_d_t_ad)
    print('')
    print('d [ BCF polar heading] / d [rbcf] (rad) AD = ', np.transpose(d_heading_bcf_polar_d_rbcf_ad[:,0]))
    print('')
    print('d [ BCF polar heading] / d [vbcf] (rad) AD = ', np.transpose(d_heading_bcf_polar_d_vbcf_ad[:,0]))
    print('')
    print('BCI topocentric heading (rad) = ', heading_bci_topocentric)
    print('')
    print('d [ BCI topocentric heading] / d [t] (rad) = ', d_heading_bci_topocentric_d_t)
    print('')
    print('d [ BCI topocentric heading] / d [rbcf] (rad) = ', d_heading_bci_topocentric_d_rbcf)
    print('')
    print('d [ BCI topocentric heading] / d [vbcf] (rad) = ', d_heading_bci_topocentric_d_vbcf)
    print('')
    print('d [ BCI topocentric heading] / d [t] (rad) AD = ', d_heading_bci_topocentric_d_t_ad)
    print('')
    print('d [ BCI topocentric heading] / d [rbcf] (rad) AD = ', np.transpose(d_heading_bci_topocentric_d_rbcf_ad[:,0]))
    print('')
    print('d [ BCI topocentric heading] / d [vbcf] (rad) AD = ', np.transpose(d_heading_bci_topocentric_d_vbcf_ad[:,0]))
    print('')
    print('BCI polar heading (rad) = ', heading_bci_polar)
    print('')
    print('d [ BCI polar heading] / d [t] (rad) analytical = ', d_heading_bci_polar_d_t)
    print('')
    print('d [ BCI polar heading] / d [rbcf] (rad) analytical = ', d_heading_bci_polar_d_rbcf)
    print('')
    print('d [ BCI polar heading] / d [vbcf] (rad) analytical = ', d_heading_bci_polar_d_vbcf)
    print('')
    print('d [ BCI polar heading] / d [t] (rad) AD = ', d_heading_bci_polar_d_t_ad)
    print('')
    print('d [ BCI polar heading] / d [rbcf] (rad) AD = ', np.transpose(d_heading_bci_polar_d_rbcf_ad[:,0]))
    print('')
    print('d [ BCI polar heading] / d [vbcf] (rad) AD = ', np.transpose(d_heading_bci_polar_d_vbcf_ad[:,0]))
    print('')
    print('')

if (print_all == 1 or print_fpa == 1):
    print('==================')
    print('Flight Path Angle:')
    print('==================')
    print('')
    print('BCF flight path angle (rad) = ', gamma_bcf)
    print('')
    print('d [BCF flight path angle] / d [t] (rad) analytical = ', d_gamma_bcf_d_t)
    print('')
    print('d [BCF flight path angle] / d [rbcf] (rad) analytical = ', d_gamma_bcf_d_rbcf[0,:])
    print('')
    print('d [BCF flight path angle] / d [vbcf] (rad) analytical = ', d_gamma_bcf_d_vbcf[0,:])
    print('')
    print('d [BCF flight path angle] / d [t] (rad) AD = ', d_gamma_bcf_d_t_ad)
    print('')
    print('d [BCF flight path angle] / d [rbcf] (rad) AD = ', np.transpose(d_gamma_bcf_d_rbcf_ad[:,0]))
    print('')
    print('d [BCF flight path angle] / d [vbcf] (rad) AD = ', np.transpose(d_gamma_bcf_d_vbcf_ad[:,0]))
    print('')
    print('BCI flight path angle (rad) = ', gamma_bci)
    print('')
    print('d [BCI flight path angle] / d [t] (rad) analytical = ', d_gamma_bci_d_t)
    print('')
    print('d [BCI flight path angle] / d [rbcf] (rad) analytical = ', d_gamma_bci_d_rbcf[0,:])
    print('')
    print('d [BCI flight path angle] / d [vbcf] (rad) analytical = ', d_gamma_bci_d_vbcf[0,:])
    print('')
    print('d [BCI flight path angle] / d [t] (rad) AD = ', d_gamma_bci_d_t_ad)
    print('')
    print('d [BCI flight path angle] / d [rbcf] (rad) AD = ', np.transpose(d_gamma_bci_d_rbcf_ad[:,0]))
    print('')
    print('d [BCI flight path angle] / d [vbcf] (rad) AD = ', np.transpose(d_gamma_bci_d_vbcf_ad[:,0]))
    print('')

#d_rbcf_fun = grad(ellipse.calc_rbcf)

#d_rbcf_d_rpa = d_rbcf_fun(r_pa)

#print(d_rbcf_d_rpa)




## test function
#def get_longitude(rbcf):
#    longitude = np.arctan2(rbcf[1], rbcf[0])
#    return longitude

## test script
#rbcf = np.array([[0.2], [0.4], [0.1]])

#dlong_fun = grad(get_longitude)

#long = get_longitude(rbcf)
#print(long)

#dlong = dlong_fun(rbcf)
#print(dlong)