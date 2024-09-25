#scale factor for fixed array
#simplified version, sun is the central body


import sympy

outputfile = open('fixed_array.out', 'w')

#actual angle
rx, ry, rz, ux, uy, uz = sympy.symbols('rx ry rz ux uy uz', real=True)

R = sympy.Matrix([rx, ry, rz])
U = sympy.Matrix([ux, uy, uz])

r = R.norm()
u = U.norm()

F = (R.cross(U)).norm() / r / u

outputfile.write('F = ' + str(F).replace('sqrt(rx**2 + ry**2 + rz**2)','r').replace('sqrt(ux**2 + uy**2 + uz**2)','u') + '\n')
outputfile.write('\n')
outputfile.write('dF_drx = ' + str(sympy.diff(F, rx)).replace('sqrt(rx**2 + ry**2 + rz**2)','r').replace('sqrt(ux**2 + uy**2 + uz**2)','u') + '\n')
outputfile.write('dF_dry = ' + str(sympy.diff(F, ry)).replace('sqrt(rx**2 + ry**2 + rz**2)','r').replace('sqrt(ux**2 + uy**2 + uz**2)','u') + '\n')
outputfile.write('dF_drz = ' + str(sympy.diff(F, rz)).replace('sqrt(rx**2 + ry**2 + rz**2)','r').replace('sqrt(ux**2 + uy**2 + uz**2)','u') + '\n')
outputfile.write('dF_dux = ' + str(sympy.diff(F, ux)).replace('sqrt(rx**2 + ry**2 + rz**2)','r').replace('sqrt(ux**2 + uy**2 + uz**2)','u') + '\n')
outputfile.write('dF_duy = ' + str(sympy.diff(F, uy)).replace('sqrt(rx**2 + ry**2 + rz**2)','r').replace('sqrt(ux**2 + uy**2 + uz**2)','u') + '\n')
outputfile.write('dF_duz = ' + str(sympy.diff(F, uz)).replace('sqrt(rx**2 + ry**2 + rz**2)','r').replace('sqrt(ux**2 + uy**2 + uz**2)','u') + '\n')
outputfile.write('\n')


outputfile.close()