# ad_test.py

# Purpose
# Testing the ad algorithmic differentiation package

# Revision history
# Noble Hatten; 09/13/2018; Began

import numpy as np
import ad
import ad.admath

# test function
def get_longitude(rbcf):
    longitude = ad.admath.atan2(rbcf[1], rbcf[0])
    #d_longitude = longitude.gradient(rbcf)
    return longitude


# test script
x = ad.adnumber(0.2) # create a scalar number for algorithmic differentation
y = ad.adnumber(0.4)
z = ad.admath.atan2(y,x) # ad.admath has math functions that can be used with AD
# z is an ad object
# z.x is just the number
# z.d(x) is dz/dx

rbcf = ad.adnumber(np.array([[0.5], [1.2], [0.4]]))
print(rbcf)

get_longitude_d, get_longitude_dd = ad.gh(get_longitude) # gradient and hessian functions

longitude = get_longitude(rbcf)
print(longitude)
longitude_d = get_longitude_d(rbcf)



print(longitude_d)
