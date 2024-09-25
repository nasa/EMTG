#constraint for reference-body-probe
#it's a lot easier to constrain cos(RBP) instead of RBP directly


import sympy

outputfile = open('RBPangle.out', 'w')

#actual angle
rx, ry, rz, vx, vy, vz = sympy.symbols('rx ry rz vx vy vz')

R = sympy.Matrix([rx, ry, rz])
V = sympy.Matrix([vx, vy, vz])

r = sympy.Function('r')(rx, ry, rz)
v = sympy.Function('v')(vx, vy, vz)

cosRBP = R.dot(V) / r / v

outputfile.write('cosRBP = ' + str(cosRBP).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v') + '\n')
outputfile.write('\n')
outputfile.write('dcosRBP_drx = ' + str(sympy.diff(cosRBP, rx)).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v').replace('(rx*vx + ry*vy + rz*vz)', 'rdotv') + '\n')
outputfile.write('dcosRBP_dry = ' + str(sympy.diff(cosRBP, ry)).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v').replace('(rx*vx + ry*vy + rz*vz)', 'rdotv') + '\n')
outputfile.write('dcosRBP_drz = ' + str(sympy.diff(cosRBP, rz)).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v').replace('(rx*vx + ry*vy + rz*vz)', 'rdotv') + '\n')
outputfile.write('dcosRBP_dvx = ' + str(sympy.diff(cosRBP, vx)).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v').replace('(rx*vx + ry*vy + rz*vz)', 'rdotv') + '\n')
outputfile.write('dcosRBP_dvy = ' + str(sympy.diff(cosRBP, vy)).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v').replace('(rx*vx + ry*vy + rz*vz)', 'rdotv') + '\n')
outputfile.write('dcosRBP_dvz = ' + str(sympy.diff(cosRBP, vz)).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v').replace('(rx*vx + ry*vy + rz*vz)', 'rdotv') + '\n')
outputfile.write('\n')


outputfile.close()