from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False)


import sympy

outputfile = open('SAMz.out', 'w')

rx, ry, rz = sympy.symbols('rx ry rz', real=True)


R = sympy.Matrix([rx, ry, rz])
zhat = sympy.Matrix([0, 0, 1])

X = R
Y = zhat.cross(R)
Z = R.cross(Y)

Xhat = X.normalized()
Yhat = Y.normalized()
Zhat = Z.normalized()

T_J2000BCI_SAM = sympy.Matrix([Xhat.transpose(),
                               Yhat.transpose(),
                               Zhat.transpose()])

T_SAM_J2000BCI = T_J2000BCI_SAM.transpose()

outputfile.write('Sun Angular Momentum(SAM) is accessed relative to J2000BCI\n')
outputfile.write('The $x$ axis points from the central body to a reference body (\\textit{i.e.} the Sun).\n')
outputfile.write('The $y$ axis is the cross product of the \\texttt{J2000BCI} spin pole and the $x$ axis.\n')
outputfile.write('The $z$ axis complete the right-handed set. (not yet implemented)\n')
outputfile.write('\n')
outputfile.write('X = ' + str(X).replace('(t)','') + '\n')
outputfile.write('Y = ' + str(Y).replace('(t)','') + '\n')
outputfile.write('Z = ' + str(Z).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('T_J2000BCI_SAM = \n')
outputfile.write(sympy.pretty(T_J2000BCI_SAM).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('T_SAM_J2000BCI = \n')
outputfile.write(sympy.pretty(T_SAM_J2000BCI).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('derivatives with respect to position\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('\n')
outputfile.write('dT_J2000BCI_SAM_drx = \n')
outputfile.write(sympy.pretty(sympy.simplify(sympy.diff(T_J2000BCI_SAM, rx))).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('dT_J2000BCI_SAM_dry = \n')
outputfile.write(sympy.pretty(sympy.simplify(sympy.diff(T_J2000BCI_SAM, ry))).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('dT_J2000BCI_SAM_drz = \n')
outputfile.write(sympy.pretty(sympy.simplify(sympy.diff(T_J2000BCI_SAM, rz))).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('z components\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('dT_J2000BCI_SAM_drx[0,2] = ' + str(sympy.simplify(sympy.diff(T_J2000BCI_SAM, rx)[0,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('dT_J2000BCI_SAM_dry[0,2] = ' + str(sympy.simplify(sympy.diff(T_J2000BCI_SAM, ry)[0,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('dT_J2000BCI_SAM_drz[0,2] = ' + str(sympy.simplify(sympy.diff(T_J2000BCI_SAM, rz)[0,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('\n')
outputfile.write('dT_J2000BCI_SAM_drx[1,2] = ' + str(sympy.simplify(sympy.diff(T_J2000BCI_SAM, rx)[1,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('dT_J2000BCI_SAM_dry[1,2] = ' + str(sympy.simplify(sympy.diff(T_J2000BCI_SAM, ry)[1,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('dT_J2000BCI_SAM_drz[1,2] = ' + str(sympy.simplify(sympy.diff(T_J2000BCI_SAM, rz)[1,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('\n')
outputfile.write('dT_J2000BCI_SAM_drx[2,2] = ' + str(sympy.simplify(sympy.diff(T_J2000BCI_SAM, rx)[2,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('rx**2','rx2').replace('ry**2','ry2').replace('rz**2','rz2').replace('rx**4','rx4').replace('ry**4','ry4') + '\n')
outputfile.write('dT_J2000BCI_SAM_dry[2,2] = ' + str(sympy.simplify(sympy.diff(T_J2000BCI_SAM, ry)[2,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('rx**2','rx2').replace('ry**2','ry2').replace('rz**2','rz2').replace('rx**4','rx4').replace('ry**4','ry4') + '\n')
outputfile.write('dT_J2000BCI_SAM_drz[2,2] = ' + str(sympy.simplify(sympy.diff(T_J2000BCI_SAM, rz)[2,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('rx**2','rx2').replace('ry**2','ry2').replace('rz**2','rz2').replace('rx**4','rx4').replace('ry**4','ry4') + '\n')





outputfile.close()