#b-plane TA wrt state vector


import sympy
from sympy import cos, sin, sqrt, atan2
from sympy.printing.cxxcode import cxxcode

outputfile = open('bplaneTA.out', 'w')

#actual angle
ex, ey, ez, rx, ry, rz = sympy.symbols('ex ez ez rx ry rz', real='True')

R = sympy.Matrix([rx, ry, rz])
E = sympy.Matrix([ex, ey, ez])

A = (E.cross(R)).normalized()
B = E.transpose() * R

TA = atan2(A, B)

outputfile.write('TA = ' + str(TA) + '\n')
outputfile.write('\n')
outputfile.write('dA_dex = ' + str(sympy.simplify(sympy.diff(A, ex))) + '\n')
outputfile.write('dA_dey = ' + str(sympy.simplify(sympy.diff(A, ey))) + '\n')
outputfile.write('dA_dez = ' + str(sympy.simplify(sympy.diff(A, ez))) + '\n')
outputfile.write('\n')
outputfile.write('dA_drx = ' + str(sympy.simplify(sympy.diff(A, rx))) + '\n')
outputfile.write('dA_dry = ' + str(sympy.simplify(sympy.diff(A, ry))) + '\n')
outputfile.write('dA_drz = ' + str(sympy.simplify(sympy.diff(A, rz))) + '\n')
outputfile.write('\n')
outputfile.write('dB_dex = ' + str(sympy.simplify(sympy.diff(B, ex))) + '\n')
outputfile.write('dB_dey = ' + str(sympy.simplify(sympy.diff(B, ey))) + '\n')
outputfile.write('dB_dez = ' + str(sympy.simplify(sympy.diff(B, ez))) + '\n')
outputfile.write('\n')
outputfile.write('dB_drx = ' + str(sympy.simplify(sympy.diff(B, rx))) + '\n')
outputfile.write('dB_dry = ' + str(sympy.simplify(sympy.diff(B, ry))) + '\n')
outputfile.write('dB_drz = ' + str(sympy.simplify(sympy.diff(B, rz))) + '\n')
outputfile.write('\n')

outputfile.close()