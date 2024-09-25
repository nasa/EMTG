#derivatives of Edelbaum spiral segment

import sympy

outputfile = open('spiral.out', 'w')

#actual angle
dv, t, g0 = sympy.symbols('dv t g0')

Isp = sympy.Function('Isp')(t)
mdot = sympy.Function('mdot')(t)

expfun = sympy.exp(-dv * 1000.0 / Isp / g0)

dexpfun_dt = sympy.diff(expfun, t)

outputfile.write('expfun = ' + str(expfun) + '\n')
outputfile.write('\n')
outputfile.write('dexpfun_dt = ' + str(dexpfun_dt) + '\n')
outputfile.write('\n')


outputfile.close()