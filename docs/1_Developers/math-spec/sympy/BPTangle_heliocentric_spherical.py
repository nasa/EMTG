import sympy

outputfile = open('BPTangle_heliocentric_spherical.out', 'w')

#assumes a unit position vector
#actual angle
RA, DEC, ux, uy, uz = sympy.symbols('RA DEC ux uy uz')

R = sympy.Matrix([-sympy.cos(RA)*sympy.cos(DEC), -sympy.sin(RA)*sympy.cos(DEC), -sympy.sin(DEC)])
U = sympy.Matrix([ux, uy, uz])

r = sympy.Function('r')(RA, DEC)
u = sympy.Function('u')(ux, uy, uz)

cosAngleThrustVectorToSun = R.dot(U) / u

outputfile.write('cosAngleThrustVectorToSun = ' + str(cosAngleThrustVectorToSun) + '\n')
outputfile.write('\n')
outputfile.write('dcosAngleThrustVectorToSun_dRA = ' + str(sympy.diff(cosAngleThrustVectorToSun, RA)).replace('(RA, DEC)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('dcosAngleThrustVectorToSun_dDEC = ' + str(sympy.diff(cosAngleThrustVectorToSun, DEC)).replace('(RA, DEC)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('dcosAngleThrustVectorToSun_dux = ' + str(sympy.diff(cosAngleThrustVectorToSun, ux)).replace('(RA, DEC)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('dcosAngleThrustVectorToSun_duy = ' + str(sympy.diff(cosAngleThrustVectorToSun, uy)).replace('(RA, DEC)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('dcosAngleThrustVectorToSun_duz = ' + str(sympy.diff(cosAngleThrustVectorToSun, uz)).replace('(RA, DEC)','').replace('(ux, uy, uz)','') + '\n')
outputfile.write('\n')