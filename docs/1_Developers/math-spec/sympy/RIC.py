from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False)


import sympy

outputfile = open('RIC.out', 'w')

t = sympy.symbols('t')
rx = sympy.Function('rx')(t)
ry = sympy.Function('ry')(t)
rz = sympy.Function('rz')(t)
vx = sympy.Function('vx')(t)
vy = sympy.Function('vy')(t)
vz = sympy.Function('vz')(t)

#rx, ry, rz, vx, vy, vz = sympy.symbols('rx ry rz vx vy vz')


R = sympy.Matrix([rx, ry, rz])
V = sympy.Matrix([vx, vy, vz])
H = R.cross(V)

T_inertial_RIC = sympy.Matrix([R.transpose(),
                               H.cross(R).transpose(),
                               H.transpose()])

T_RIC_inertial = T_inertial_RIC.transpose()

outputfile.write('R and V are the position and velocity vectors of the central body relative to some user-supplied reference body\n')
outputfile.write('For example, an asteroid relative to the Sun\n')
outputfile.write('\n')
outputfile.write('R = ' + str(R).replace('(t)','') + '\n')
outputfile.write('V = ' + str(V).replace('(t)','') + '\n')
outputfile.write('H = ' + str(H).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('T_inertial_RIC = \n')
outputfile.write(sympy.pretty(T_inertial_RIC).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('T_RIC_inertial = \n')
outputfile.write(sympy.pretty(T_RIC_inertial).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('derivatives with respect to position\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('\n')
outputfile.write('dT_RIC_inertial_drx = \n')
outputfile.write(sympy.pretty(sympy.diff(T_RIC_inertial, rx)).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('dT_RIC_inertial_dry = \n')
outputfile.write(sympy.pretty(sympy.diff(T_RIC_inertial, ry)).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('dT_RIC_inertial_drz = \n')
outputfile.write(sympy.pretty(sympy.diff(T_RIC_inertial, rz)).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('derivatives with respect to velocity\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('\n')
outputfile.write('dT_RIC_inertial_dvx = \n')
outputfile.write(sympy.pretty(sympy.diff(T_RIC_inertial, vx)).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('dT_RIC_inertial_dvy = \n')
outputfile.write(sympy.pretty(sympy.diff(T_RIC_inertial, vy)).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('dT_RIC_inertial_dvz = \n')
outputfile.write(sympy.pretty(sympy.diff(T_RIC_inertial, vz)).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('derivatives with respect to time\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('\n')
outputfile.write('dT_RIC_inertial_dt = \n')
outputfile.write(sympy.pretty(sympy.simplify(sympy.diff(T_RIC_inertial, t))).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('dT_inertial_RIC_dt = \n')
outputfile.write(sympy.pretty(sympy.simplify(sympy.diff(T_inertial_RIC, t))).replace('(t)','') + '\n')
outputfile.write('\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('derivatives with respect to time without pretty-print\n')
outputfile.write('--------------------------------------------\n')
outputfile.write('\n')
outputfile.write('dT_RIC_inertial_dt = \n')
outputfile.write(str(sympy.simplify(sympy.diff(T_RIC_inertial, t))).replace('(t)','').replace('Derivative(rx, t)','drx_dt').replace('Derivative(ry, t)','dry_dt').replace('Derivative(rz, t)','drz_dt').replace('Derivative(vx, t)','dvx_dt').replace('Derivative(vy, t)','dvy_dt').replace('Derivative(vz, t)','dvz_dt') + '\n')
outputfile.write('\n')
outputfile.write('dT_inertial_RIC_dt = \n')
outputfile.write(str(sympy.simplify(sympy.diff(T_inertial_RIC, t))).replace('(t)','').replace('Derivative(rx, t)','drx_dt').replace('Derivative(ry, t)','dry_dt').replace('Derivative(rz, t)','drz_dt').replace('Derivative(vx, t)','dvx_dt').replace('Derivative(vy, t)','dvy_dt').replace('Derivative(vz, t)','dvz_dt') + '\n')






outputfile.close()