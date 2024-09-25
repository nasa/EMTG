import numpy
from math import tanh

g0 = 9.80665

P_sel = numpy.array([4375, 3752, 3460, 3008])

T_sel = numpy.array([246.0, 221.0, 184.0, 87.0])

c_sel = numpy.array([1790, 1617, 1099, 1579]) * g0

mdot_sel = numpy.array([14, 13.9, 17.1, 11.4])

m = 2400.0

lambda_m = 0.89

BtransposeLambda_norm = 88500

P_av = 4000

rho_b = 1.0e-5
rho_c = 1.0e-5
rho_m = 1.0e-5

Pw = numpy.zeros(4)
SF = numpy.zeros(4)
S_delta = numpy.zeros(4)

eta_m = numpy.zeros(4)

for i in range(0, 4):
    Pw[i] = 0.5 * (1.0 + tanh((P_av - P_sel[i]) / rho_c))
    
    SF[i] = c_sel[i] * BtransposeLambda_norm / m + lambda_m - 1.0 #switching function
    
    S_delta[i] = 0.5 * (1 + tanh(SF[i] / rho_b)) #mask
    
S_L = numpy.multiply(S_delta, Pw) #element multiply

S_H = mdot_sel + T_sel / m * BtransposeLambda_norm - lambda_m * mdot_sel

S_m = numpy.multiply(S_H, S_L)

for i in range(0, 4):