import sympy

outputfile = open('BPTangle.out', 'w')

#actual angle
rx, ry, rz, ux, uy, uz = sympy.symbols('rx ry rz ux uy uz')

R = sympy.Matrix([rx, ry, rz])
U = sympy.Matrix([ux, uy, uz])

r = sympy.Function('r')(rx, ry, rz)
u = sympy.Function('u')(ux, uy, uz)

cosAngleThrustVectorToSun = R.dot(U) / r / u

outputfile.write('cosAngleThrustVectorToSun = ' + str(cosAngleThrustVectorToSun) + '\n')
outputfile.write('\n')
outputfile.write('dcosAngleThrustVectorToSun_drx = ' + str(sympy.diff(cosAngleThrustVectorToSun, rx)).replace('(rx, ry ,rz)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('dcosAngleThrustVectorToSun_dry = ' + str(sympy.diff(cosAngleThrustVectorToSun, ry)).replace('(rx, ry ,rz)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('dcosAngleThrustVectorToSun_drz = ' + str(sympy.diff(cosAngleThrustVectorToSun, rz)).replace('(rx, ry ,rz)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('dcosAngleThrustVectorToSun_dux = ' + str(sympy.diff(cosAngleThrustVectorToSun, ux)).replace('(rx, ry ,rz)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('dcosAngleThrustVectorToSun_duy = ' + str(sympy.diff(cosAngleThrustVectorToSun, uy)).replace('(rx, ry ,rz)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('dcosAngleThrustVectorToSun_duz = ' + str(sympy.diff(cosAngleThrustVectorToSun, uz)).replace('(rx, ry ,rz)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('\n')

#heaviside bounds

Mid, Ceiling, Sharpness, rstar, Floor, Scale = sympy.symbols('Mid Ceiling Sharpness rstar Floor Scale')

#Floor = Mid - (Ceiling - Mid)
#Scale = Ceiling - Floor

H = Floor + Scale / (1 + sympy.exp(-Sharpness * (r - rstar)))

outputfile.write('H = ' + str(H).replace('(rx, ry ,rz)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('\n')
outputfile.write('dH_dr = ' + str(sympy.diff(H, r)).replace('(rx, ry ,rz)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('\n')
outputfile.write('dH_drx = ' + str(sympy.diff(H, rx)).replace('(rx, ry ,rz)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('dH_dry = ' + str(sympy.diff(H, ry)).replace('(rx, ry ,rz)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('dH_drz = ' + str(sympy.diff(H, rz)).replace('(rx, ry ,rz)','').replace('(ux, uy, uz)','') + '\n')


outputfile.close()