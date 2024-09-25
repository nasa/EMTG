from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False)


import sympy

outputfile = open('SAM.out', 'w')

t = sympy.symbols('t', real=True)
rx = sympy.Function('rx')(t)
ry = sympy.Function('ry')(t)
rz = sympy.Function('rz')(t)

R = sympy.Matrix([rx, ry, rz])
zhat = sympy.Matrix([0, 0, 1])

X = R
Y = zhat.cross(R)
Z = R.cross(Y)

Xhat = X / sympy.sqrt(X[0]**2 + X[1]**2 + X[2]**2)
Yhat = Y / sympy.sqrt(Y[0]**2 + Y[1]**2 + Y[2]**2)
Zhat = Z / sympy.sqrt(Z[0]**2 + Z[1]**2 + Z[2]**2)

T_J2000BCI_SAM = sympy.Matrix([Xhat.transpose(),
                               Yhat.transpose(),
                               Zhat.transpose()])

T_SAM_J2000BCI = T_J2000BCI_SAM.transpose()

outputfile.write('Sun Angular Momentum(SAM) is accessed relative to J2000BCI\n')
outputfile.write('The $x$ axis points from the central body to a reference body (\\textit{i.e.} the Sun).\n')
outputfile.write('The $y$ axis is the cross product of the \\texttt{J2000BCI} spin pole and the $x$ axis.\n')
outputfile.write('The $z$ axis complete the right-handed set. (not yet implemented)\n')
outputfile.write('\n')
outputfile.write('Xhat = ' + str(Xhat).replace('(t)','') + '\n')
outputfile.write('Yhat = ' + str(Yhat).replace('(t)','') + '\n')
outputfile.write('Zhat = ' + str(Zhat).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('T_J2000BCI_SAM = \n')
outputfile.write(sympy.pretty(T_J2000BCI_SAM) + '\n')
outputfile.write('\n')
outputfile.write('T_SAM_J2000BCI = \n')
outputfile.write(sympy.pretty(T_SAM_J2000BCI) + '\n')
outputfile.write('\n')
# outputfile.write('--------------------------------------------\n')
# outputfile.write('derivatives with respect to position\n')
# outputfile.write('--------------------------------------------\n')
# outputfile.write('\n')
# outputfile.write('T_SAM_J2000BCI = \n')
# outputfile.write(sympy.pretty(sympy.diff(T_SAM_J2000BCI, rx)).replace('(t)','') + '\n')
# outputfile.write('\n')
# outputfile.write('T_SAM_J2000BCI = \n')
# outputfile.write(sympy.pretty(sympy.diff(T_SAM_J2000BCI, ry)).replace('(t)','') + '\n')
# outputfile.write('\n')
# outputfile.write('T_SAM_J2000BCI = \n')
# outputfile.write(sympy.pretty(sympy.diff(T_SAM_J2000BCI, rz)).replace('(t)','') + '\n')
# outputfile.write('\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('derivatives with respect to time\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('\n')
outputfile.write('dT_J2000BCI_SAM_dt = \n')
outputfile.write(sympy.pretty(sympy.simplify(sympy.diff(T_J2000BCI_SAM, t))) + '\n')
outputfile.write('\n')
outputfile.write('dT_SAM_J2000BCI_dt = \n')
outputfile.write(sympy.pretty(sympy.simplify(sympy.diff(T_SAM_J2000BCI, t))) + '\n')
outputfile.write('\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('derivatives with respect to time without pretty-print\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('\n')
outputfile.write('dT_J2000BCI_SAM_dt = \n')
outputfile.write(str(sympy.simplify(sympy.diff(T_J2000BCI_SAM, t))).replace('(t)','').replace('Derivative(rx, t)','drx_dt').replace('Derivative(ry, t)','dry_dt').replace('Derivative(rz, t)','drz_dt').replace('Derivative(vx, t)','dvx_dt').replace('Derivative(vy, t)','dvy_dt').replace('Derivative(vz, t)','dvz_dt') + '\n')
outputfile.write('\n')
outputfile.write('dT_SAM_J2000BCI_dt = \n')
outputfile.write(str(sympy.simplify(sympy.diff(T_SAM_J2000BCI, t))).replace('(t)','').replace('Derivative(rx, t)','drx_dt').replace('Derivative(ry, t)','dry_dt').replace('Derivative(rz, t)','drz_dt').replace('Derivative(vx, t)','dvx_dt').replace('Derivative(vy, t)','dvy_dt').replace('Derivative(vz, t)','dvz_dt') + '\n')






outputfile.close()