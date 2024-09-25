#constraint for vertical flight path angle
#it's a lot easier to constrain cos(VFPA) instead of VFPA or HFPA directly


import sympy

outputfile = open('VFPAangle.out', 'w')

#actual angle
rx, ry, rz, vx, vy, vz = sympy.symbols('rx ry rz vx vy vz')

R = sympy.Matrix([rx, ry, rz])
V = sympy.Matrix([vx, vy, vz])

r = sympy.Function('r')(rx, ry, rz)
v = sympy.Function('v')(vx, vy, vz)

cosVFPA = R.dot(V) / r / v

outputfile.write('cosVFPA = ' + str(cosVFPA).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v') + '\n')
outputfile.write('\n')
outputfile.write('dcosVFPA_drx = ' + str(sympy.diff(cosVFPA, rx)).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v').replace('(rx*vx + ry*vy + rz*vz)', 'rdotv') + '\n')
outputfile.write('dcosVFPA_dry = ' + str(sympy.diff(cosVFPA, ry)).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v').replace('(rx*vx + ry*vy + rz*vz)', 'rdotv') + '\n')
outputfile.write('dcosVFPA_drz = ' + str(sympy.diff(cosVFPA, rz)).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v').replace('(rx*vx + ry*vy + rz*vz)', 'rdotv') + '\n')
outputfile.write('dcosVFPA_dvx = ' + str(sympy.diff(cosVFPA, vx)).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v').replace('(rx*vx + ry*vy + rz*vz)', 'rdotv') + '\n')
outputfile.write('dcosVFPA_dvy = ' + str(sympy.diff(cosVFPA, vy)).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v').replace('(rx*vx + ry*vy + rz*vz)', 'rdotv') + '\n')
outputfile.write('dcosVFPA_dvz = ' + str(sympy.diff(cosVFPA, vz)).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v').replace('(rx*vx + ry*vy + rz*vz)', 'rdotv') + '\n')
outputfile.write('\n')


outputfile.close()