#non-central body sun


import sympy

outputfile = open('vector_from_sun.out', 'w')

#actual angle
#rx, ry, rz, rcx, rcy, rcz = sympy.symbols('rx ry rz rcx rcy rcz', real=True)
t = sympy.symbols('t', real=True)
rx = sympy.Function('rx', real=True)(t)
ry = sympy.Function('ry', real=True)(t)
rz = sympy.Function('rz', real=True)(t)
rcx = sympy.Function('rcx', real=True)(t)
rcy = sympy.Function('rcy', real=True)(t)
rcz = sympy.Function('rcz', real=True)(t)

R_cb_to_sc = sympy.Matrix([rx, ry, rz])
R_sun_to_cb = sympy.Matrix([rcx, rcy, rcz])

r_cb_to_sc = R_cb_to_sc.norm()
r_sun_to_cb = R_sun_to_cb.norm()

R = R_cb_to_sc + R_sun_to_cb

r = R.norm()

outputfile.write('R[0] = ' + str(R[0]) + '\n')
outputfile.write('R[1] = ' + str(R[1]) + '\n')
outputfile.write('R[2] = ' + str(R[2]) + '\n')
outputfile.write('r = ' + str(r) + '\n')
outputfile.write('\n')
outputfile.write('dr_drcx = ' + str(sympy.diff(r, rcx)).replace('sqrt((rcx(t) + rx(t))**2 + (rcy(t) + ry(t))**2 + (rcz(t) + rz(t))**2)','r(t)') + '\n')
outputfile.write('dr_drcy = ' + str(sympy.diff(r, rcy)).replace('sqrt((rcx(t) + rx(t))**2 + (rcy(t) + ry(t))**2 + (rcz(t) + rz(t))**2)','r(t)') + '\n')
outputfile.write('dr_drcz = ' + str(sympy.diff(r, rcz)).replace('sqrt((rcx(t) + rx(t))**2 + (rcy(t) + ry(t))**2 + (rcz(t) + rz(t))**2)','r(t)') + '\n')
outputfile.write('dr_drx = ' + str(sympy.diff(r, rx)).replace('sqrt((rcx(t) + rx(t))**2 + (rcy(t) + ry(t))**2 + (rcz(t) + rz(t))**2)','r(t)') + '\n')
outputfile.write('dr_dry = ' + str(sympy.diff(r, ry)).replace('sqrt((rcx(t) + rx(t))**2 + (rcy(t) + ry(t))**2 + (rcz(t) + rz(t))**2)','r(t)') + '\n')
outputfile.write('dr_drz = ' + str(sympy.diff(r, rz)).replace('sqrt((rcx(t) + rx(t))**2 + (rcy(t) + ry(t))**2 + (rcz(t) + rz(t))**2)','r(t)') + '\n')
outputfile.write('dr_dt = ' + str(sympy.simplify(sympy.diff(r, t))).replace('Derivative','d').replace('sqrt((rcx(t) + rx(t))**2 + (rcy(t) + ry(t))**2 + (rcz(t) + rz(t))**2)','r(t)') + '\n')
outputfile.write('\n')


outputfile.close()