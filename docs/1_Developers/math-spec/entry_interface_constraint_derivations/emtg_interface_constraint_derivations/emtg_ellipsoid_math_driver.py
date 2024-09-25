import numpy as np
import ellipsoid_math

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
theta1 = -21.0
theta2 = -9.0
theta3 = -14.0

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
w0 = 0.0 # deg
w0dot = 1.5 # deg/s

ellipse = ellipsoid_math.ellipsoid_math_class(a=a, b=b, c=c, r_pa=r_pa, vbcf=vbcf, w0=w0, w0dot=w0dot, theta1=theta1, theta2=theta2, theta3=theta3)
ellipse.calc_everything()
ellipse.calc_longitude_bodydetic_bekta()
ellipse.calc_latitude_bodydetic_bekta()
#ellipse.calc_latitude_bodydetic_oblate()

print('========')
print('POSITION')
print('========')
print('')
print('r_pa = ', np.transpose(r_pa))
print('')
print('r_bcf = ', np.transpose(ellipse.rbcf))
print('')
print('========')
print('VELOCITY')
print('========')
print('')
print('v_bcf = ', np.transpose(vbcf))
print('')
print('v_bcf magnitude = ', ellipse.vbcf_mag)
print('')
print('v_bci in bcf = ', np.transpose(ellipse.vbci_expressed_in_bcf))
print('')
print('v_bci magnitude = ', ellipse.vbci_mag)
print('')
print('v_bcf in pa = ', np.transpose(ellipse.v_bcf_expressed_in_pa))
print('')
print('v_bci in bci = ')
print('')
print('=======================')
print('BODYCENTRIC COORDINATES')
print('=======================')
print('')
print('bodycentric longitude (deg) = ', np.rad2deg(ellipse.longitude_bodycentric))
print('')
print('bodycentric latitude (deg) = ', np.rad2deg(ellipse.latitude_bodycentric))
print('')
print('=====================')
print('BODYDETIC COORDINATES')
print('=====================')
print('')
print('bodydetic longitude (deg) = ', np.rad2deg(ellipse.longitude_bodydetic))
print('')
print('bodydetic latitude (deg) = ', np.rad2deg(ellipse.latitude_bodydetic))
print('')
print('========================')
print('TOPOCENTRIC UNIT VECTORS')
print('========================')
print('')
print('South vector (topocentric) = ', np.transpose(ellipse.s_bcf_topocentric_unit))
print('')
print('East vector (topocentric) = ', np.transpose(ellipse.e_bcf_topocentric_unit))
print('')
print('Normal/Up vector = ', np.transpose(ellipse.u_bcf_unit))
print('')
print('========================')
print('LOCAL POLAR UNIT VECTORS')
print('========================')
print('')
print('South vector (polar) = ', np.transpose(ellipse.s_polar_bcf_unit))
print('')
print('East vector (polar) = ', np.transpose(ellipse.e_polar_bcf_unit))
print('')
print('Normal/Up vector = ', np.transpose(ellipse.u_bcf_unit))
print('')
print('South vector (Kyle) = ', np.transpose(ellipse.south_kyle))
print('')
print('Easterly vector (Kyle) = ', np.transpose(-ellipse.westerly_kyle))
print('')
print('========================')
print('HEADING')
print('========================')
print('')
print('Heading, polar (BCF) (deg) = ', np.rad2deg(ellipse.heading_bcf_polar))
print('')
print('Heading, topocentric (BCF) (deg) = ', np.rad2deg(ellipse.heading_bcf_topocentric))
print('')
print('Heading, polar (BCI) (deg) = ', np.rad2deg(ellipse.heading_bci_polar))
print('')
print('Heading, topocentric (BCI) (deg) = ', np.rad2deg(ellipse.heading_bci_topocentric))
print('')
print('========================')
print('FLIGHT PATH ANGLE')
print('========================')
print('')
print('Flight path angle (BCF) (deg) = ', np.rad2deg(ellipse.gamma_bcf))
print('')
print('Flight path angle atan (BCF) (deg) = ', np.rad2deg(ellipse.gamma_bcf_atan)) 
print('')
print('Flight path angle (BCI) (deg) = ', np.rad2deg(ellipse.gamma_bci))

#ellipse.plot()