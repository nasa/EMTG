import sympy

outputfile = open('probe_separation_impulse.out', 'w')

vx, vy, vz, I, m = sympy.symbols('vx vy vz I m', real=True)
V = sympy.Matrix([vx, vy, vz])
v = sympy.sqrt(vx**2 + vy**2 + vz**2)
unitV = V / v

dV = unitV * I / m



outputfile.write('V[0] = ' + str(V[0]) + '\n')
outputfile.write('V[1] = ' + str(V[1]) + '\n')
outputfile.write('V[2] = ' + str(V[2]) + '\n')
outputfile.write('v = ' + str(v) + '\n')
outputfile.write('\n')
outputfile.write('dV[0] = ' + str(dV[0]) + '\n')
outputfile.write('dV[1] = ' + str(dV[1]) + '\n')
outputfile.write('dV[2] = ' + str(dV[2]) + '\n')
outputfile.write('\n')
outputfile.write('ddVx_dVx = ' + str(sympy.diff(dV[0], vx)).replace('(vx**2 + vy**2 + vz**2)**(3/2)','v^3').replace('sqrt(vx**2 + vy**2 + vz**2)','v') + '\n')
outputfile.write('ddVx_dVy = ' + str(sympy.diff(dV[0], vy)).replace('(vx**2 + vy**2 + vz**2)**(3/2)','v^3').replace('sqrt(vx**2 + vy**2 + vz**2)','v') + '\n')
outputfile.write('ddVx_dVz = ' + str(sympy.diff(dV[0], vz)).replace('(vx**2 + vy**2 + vz**2)**(3/2)','v^3').replace('sqrt(vx**2 + vy**2 + vz**2)','v') + '\n')
outputfile.write('ddVy_dVx = ' + str(sympy.diff(dV[1], vx)).replace('(vx**2 + vy**2 + vz**2)**(3/2)','v^3').replace('sqrt(vx**2 + vy**2 + vz**2)','v') + '\n')
outputfile.write('ddVy_dVy = ' + str(sympy.diff(dV[1], vy)).replace('(vx**2 + vy**2 + vz**2)**(3/2)','v^3').replace('sqrt(vx**2 + vy**2 + vz**2)','v') + '\n')
outputfile.write('ddVy_dVz = ' + str(sympy.diff(dV[1], vz)).replace('(vx**2 + vy**2 + vz**2)**(3/2)','v^3').replace('sqrt(vx**2 + vy**2 + vz**2)','v') + '\n')
outputfile.write('ddVz_dVx = ' + str(sympy.diff(dV[2], vx)).replace('(vx**2 + vy**2 + vz**2)**(3/2)','v^3').replace('sqrt(vx**2 + vy**2 + vz**2)','v') + '\n')
outputfile.write('ddVz_dVy = ' + str(sympy.diff(dV[2], vy)).replace('(vx**2 + vy**2 + vz**2)**(3/2)','v^3').replace('sqrt(vx**2 + vy**2 + vz**2)','v') + '\n')
outputfile.write('ddVz_dVz = ' + str(sympy.diff(dV[2], vz)).replace('(vx**2 + vy**2 + vz**2)**(3/2)','v^3').replace('sqrt(vx**2 + vy**2 + vz**2)','v') + '\n')
outputfile.write('\n')


outputfile.close()