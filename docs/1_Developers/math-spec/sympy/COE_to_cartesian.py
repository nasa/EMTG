#constraint for reference-body-probe


import sympy
from sympy import cos, sin, sqrt

outputfile = open('COE_to_cartesian.out', 'w')

#actual angle
SMA, ECC, INC, RAAN, AOP, TA, mu = sympy.symbols('SMA ECC INC RAAN AOP TA mu')

theta = AOP + TA

r = SMA * (1 - ECC**2) / (1 + ECC * cos(TA))

v = sqrt(mu * (2 / r - 1 / SMA))

h = sqrt(mu * SMA * (1 - ECC**2))

rx = r * (cos(RAAN) * cos(theta) - sin(RAAN) * sin(theta) * cos(INC))
ry = r * (sin(RAAN) * cos(theta) + cos(RAAN) * sin(theta) * cos(INC))
rz = r * sin(theta) * sin(INC)

vx = -mu / h * (cos(RAAN) * (sin(theta) + ECC * sin(AOP)) + sin(RAAN) * (cos(theta) + ECC * cos(AOP)) * cos(INC))
vy = -mu / h * (sin(RAAN) * (sin(theta) + ECC * sin(AOP)) - cos(RAAN) * (cos(theta) + ECC * cos(AOP)) * cos(INC))
vz =  mu / h * (cos(theta) + ECC * cos(AOP)) * sin(INC)

outputfile.write('doubleType rx = ' + str(rx).replace('**2','2') + ';\n')
outputfile.write('doubleType ry = ' + str(ry).replace('**2','2') + ';\n')
outputfile.write('doubleType rz = ' + str(rz).replace('**2','2') + ';\n')
outputfile.write('doubleType vx = ' + str(vx).replace('**2','2') + ';\n')
outputfile.write('doubleType vy = ' + str(vy).replace('**2','2') + ';\n')
outputfile.write('doubleType vz = ' + str(vz).replace('**2','2') + ';\n')
outputfile.write('\n')
outputfile.write('double drx_dSMA = ('  + str(sympy.diff(rx, SMA)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double drx_dECC = ('  + str(sympy.diff(rx, ECC)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double drx_dINC = ('  + str(sympy.diff(rx, INC)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double drx_dRAAN = (' + str(sympy.diff(rx, RAAN)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double drx_dAOP = ('  + str(sympy.diff(rx, AOP)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double drx_dTA = ('   + str(sympy.diff(rx, TA)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('\n')
outputfile.write('double dry_dSMA = ('  + str(sympy.diff(ry, SMA)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dry_dECC = ('  + str(sympy.diff(ry, ECC)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dry_dINC = ('  + str(sympy.diff(ry, INC)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dry_dRAAN = (' + str(sympy.diff(ry, RAAN)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dry_dAOP = ('  + str(sympy.diff(ry, AOP)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dry_dTA = ('   + str(sympy.diff(ry, TA)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('\n')
outputfile.write('double drz_dSMA = ('  + str(sympy.diff(rz, SMA)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double drz_dECC = ('  + str(sympy.diff(rz, ECC)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double drz_dINC = ('  + str(sympy.diff(rz, INC)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double drz_dRAAN = (' + str(sympy.diff(rz, RAAN)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double drz_dAOP = ('  + str(sympy.diff(rz, AOP)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double drz_dTA = ('   + str(sympy.diff(rz, TA)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('\n')
outputfile.write('double dvx_dSMA = ('  + str(sympy.diff(vx, SMA)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvx_dECC = ('  + str(sympy.diff(vx, ECC)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvx_dINC = ('  + str(sympy.diff(vx, INC)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvx_dRAAN = (' + str(sympy.diff(vx, RAAN)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvx_dAOP = ('  + str(sympy.diff(vx, AOP)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvx_dTA = ('   + str(sympy.diff(vx, TA)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('\n')
outputfile.write('double dvy_dSMA = ('  + str(sympy.diff(vy, SMA)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvy_dECC = ('  + str(sympy.diff(vy, ECC)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvy_dINC = ('  + str(sympy.diff(vy, INC)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvy_dRAAN = (' + str(sympy.diff(vy, RAAN)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvy_dAOP = ('  + str(sympy.diff(vy, AOP)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvy_dTA = ('   + str(sympy.diff(vy, TA)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('\n')
outputfile.write('double dvz_dSMA = ('  + str(sympy.diff(vz, SMA)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvz_dECC = ('  + str(sympy.diff(vz, ECC)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvz_dINC = ('  + str(sympy.diff(vz, INC)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvz_dRAAN = (' + str(sympy.diff(vz, RAAN)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvz_dAOP = ('  + str(sympy.diff(vz, AOP)).replace('**2','2') + ')_GETVALUE;\n')
outputfile.write('double dvz_dTA = ('   + str(sympy.diff(vz, TA)).replace('**2','2') + ')_GETVALUE;\n')


outputfile.close()