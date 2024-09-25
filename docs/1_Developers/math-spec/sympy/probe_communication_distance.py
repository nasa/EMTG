import sympy

outputfile = open('probe_communication_distance.out', 'w')

xsc, ysc, zsc, xp, yp, zp = sympy.symbols('xsc ysc zsc xp yp zp', real=True)
Rsc = sympy.Matrix([xsc, ysc, zsc])
Rp = sympy.Matrix([xp, yp, zp])
R = Rsc - Rp
r = (R[0]*R[0] + R[1]*R[1] + R[2]*R[2])**0.5


outputfile.write('R[0] = ' + str(R[0]) + '\n')
outputfile.write('R[1] = ' + str(R[1]) + '\n')
outputfile.write('R[2] = ' + str(R[2]) + '\n')
outputfile.write('r = ' + str(r) + '\n')
outputfile.write('\n')
outputfile.write('dr_dxsc = ' + str(sympy.diff(r, xsc)).replace('*((-xp + xsc)**2 + (-yp + ysc)**2 + (-zp + zsc)**2)**(-0.5)','/r').replace('1.0*','') + '\n')
outputfile.write('dr_dysc = ' + str(sympy.diff(r, ysc)).replace('*((-xp + xsc)**2 + (-yp + ysc)**2 + (-zp + zsc)**2)**(-0.5)','/r').replace('1.0*','') + '\n')
outputfile.write('dr_dzsc = ' + str(sympy.diff(r, zsc)).replace('*((-xp + xsc)**2 + (-yp + ysc)**2 + (-zp + zsc)**2)**(-0.5)','/r').replace('1.0*','') + '\n')
outputfile.write('dr_dxp = ' + str(sympy.diff(r, xp)).replace('*((-xp + xsc)**2 + (-yp + ysc)**2 + (-zp + zsc)**2)**(-0.5)','/r').replace('1.0*','') + '\n')
outputfile.write('dr_dyp = ' + str(sympy.diff(r, yp)).replace('*((-xp + xsc)**2 + (-yp + ysc)**2 + (-zp + zsc)**2)**(-0.5)','/r').replace('1.0*','') + '\n')
outputfile.write('dr_dzp = ' + str(sympy.diff(r, zp)).replace('*((-xp + xsc)**2 + (-yp + ysc)**2 + (-zp + zsc)**2)**(-0.5)','/r').replace('1.0*','') + '\n')
outputfile.write('\n')


outputfile.close()