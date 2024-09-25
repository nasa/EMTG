from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False)


import sympy

outputfile = open('ObjectReferencedz.out', 'w')

rx, ry, rz, vx, vy, vz = sympy.symbols('rx ry rz vx vy vz', real=True)


R = sympy.Matrix([rx, ry, rz])
V = sympy.Matrix([vx, vy, vz])
zhat = sympy.Matrix([0, 0, 1])

X = R
Y = V
Z = R.cross(V)

Xhat = X.normalized()
Yhat = Y.normalized()
Zhat = Z.normalized()

T_ICRF_ObjectReferenced = sympy.Matrix([Xhat.transpose(),
                                        Yhat.transpose(),
                                        Zhat.transpose()])

T_ObjectReferenced_ICRF = T_ICRF_ObjectReferenced.transpose()

outputfile.write('ObjectReferenced is accessed relative to ICRF\n')
outputfile.write('The $x$ axis points from the central body to a reference body.\n')
outputfile.write('The $y$ axis is the reference body\'s velocity vector wrt the reference body.\n')
outputfile.write('The $z$ axis is the reference body\'s angular momentum vector wrt the reference body.\n')
outputfile.write('\n')
outputfile.write('X = ' + str(X).replace('(t)','') + '\n')
outputfile.write('Y = ' + str(Y).replace('(t)','') + '\n')
outputfile.write('Z = ' + str(Z).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('T_ICRF_ObjectReferenced = \n')
outputfile.write(sympy.pretty(T_ICRF_ObjectReferenced).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('T_ObjectReferenced_ICRF = \n')
outputfile.write(sympy.pretty(T_ObjectReferenced_ICRF).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('derivatives with respect to position\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('\n')
outputfile.write('dT_ICRF_ObjectReferenced_drx = \n')
outputfile.write(sympy.pretty(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, rx))).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('dT_ICRF_ObjectReferenced_dry = \n')
outputfile.write(sympy.pretty(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, ry))).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('dT_ICRF_ObjectReferenced_drz = \n')
outputfile.write(sympy.pretty(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, rz))).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('z components\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('dT_ICRF_ObjectReferenced_drx[0,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, rx)[0,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('dT_ICRF_ObjectReferenced_dry[0,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, ry)[0,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('dT_ICRF_ObjectReferenced_drz[0,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, rz)[0,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('\n')
outputfile.write('dT_ICRF_ObjectReferenced_drx[1,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, rx)[1,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('dT_ICRF_ObjectReferenced_dry[1,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, ry)[1,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('dT_ICRF_ObjectReferenced_drz[1,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, rz)[1,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('\n')
outputfile.write('dT_ICRF_ObjectReferenced_drx[2,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, rx)[2,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('rx**2','rx2').replace('ry**2','ry2').replace('rz**2','rz2').replace('rx**4','rx4').replace('ry**4','ry4') + '\n')
outputfile.write('dT_ICRF_ObjectReferenced_dry[2,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, ry)[2,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('rx**2','rx2').replace('ry**2','ry2').replace('rz**2','rz2').replace('rx**4','rx4').replace('ry**4','ry4') + '\n')
outputfile.write('dT_ICRF_ObjectReferenced_drz[2,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, rz)[2,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('rx**2','rx2').replace('ry**2','ry2').replace('rz**2','rz2').replace('rx**4','rx4').replace('ry**4','ry4') + '\n')
outputfile.write('\n')
outputfile.write('dT_ICRF_ObjectReferenced_dvx[0,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, vx)[0,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('dT_ICRF_ObjectReferenced_dvy[0,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, vy)[0,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('dT_ICRF_ObjectReferenced_dvz[0,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, vz)[0,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('\n')
outputfile.write('dT_ICRF_ObjectReferenced_dvx[1,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, vx)[1,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('dT_ICRF_ObjectReferenced_dvy[1,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, vy)[1,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('dT_ICRF_ObjectReferenced_dvz[1,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, vz)[1,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('(rx**2 + ry**2)**(3/2)','rxy3') + '\n')
outputfile.write('\n')
outputfile.write('dT_ICRF_ObjectReferenced_dvx[2,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, vx)[2,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('rx**2','rx2').replace('ry**2','ry2').replace('rz**2','rz2').replace('rx**4','rx4').replace('ry**4','ry4') + '\n')
outputfile.write('dT_ICRF_ObjectReferenced_dvy[2,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, vy)[2,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('rx**2','rx2').replace('ry**2','ry2').replace('rz**2','rz2').replace('rx**4','rx4').replace('ry**4','ry4') + '\n')
outputfile.write('dT_ICRF_ObjectReferenced_dvz[2,2] = ' + str(sympy.simplify(sympy.diff(T_ICRF_ObjectReferenced, vz)[2,2])).replace('(rx**2 + ry**2 + rz**2)**(3/2)','r3').replace('rx**2','rx2').replace('ry**2','ry2').replace('rz**2','rz2').replace('rx**4','rx4').replace('ry**4','ry4') + '\n')

outputfile.close()