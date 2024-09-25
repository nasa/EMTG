import sympy
from sympy import cos, sin, sqrt, asin, atan2, acos

outputfile = open('SphericalRADEC_to_Cartesian.out', 'w')

#actual angle
x, y, z, vx, vy, vz = sympy.symbols('x y z vx vy vz')

r = sqrt(x * x + y * y + z * z);
RA = atan2(y, x);
DEC = asin(z / r);
v = sqrt(vx * vx + vy * vy + vz * vz);
vRA = atan2(vy, vx);
vDEC = asin(vz / r);

outputfile.write('double dr_dx = ('  + str(sympy.diff(r, x)).replace('**2','2') + ');\n')
outputfile.write('double dr_dx = ('  + str(sympy.diff(r, y)).replace('**2','2') + ');\n')
outputfile.write('double dr_dx = ('  + str(sympy.diff(r, z)).replace('**2','2') + ');\n')
outputfile.write('double dr_dx = ('  + str(sympy.diff(r, vx)).replace('**2','2') + ');\n')
outputfile.write('double dr_dx = ('  + str(sympy.diff(r, vy)).replace('**2','2') + ');\n')
outputfile.write('double dr_dx = ('  + str(sympy.diff(r, vz)).replace('**2','2') + ');\n')
outputfile.write('\n')
outputfile.write('double dRA_dx = ('  + str(sympy.diff(RA, x)).replace('**2','2') + ');\n')
outputfile.write('double dRA_dx = ('  + str(sympy.diff(RA, y)).replace('**2','2') + ');\n')
outputfile.write('double dRA_dx = ('  + str(sympy.diff(RA, z)).replace('**2','2') + ');\n')
outputfile.write('double dRA_dx = ('  + str(sympy.diff(RA, vx)).replace('**2','2') + ');\n')
outputfile.write('double dRA_dx = ('  + str(sympy.diff(RA, vy)).replace('**2','2') + ');\n')
outputfile.write('double dRA_dx = ('  + str(sympy.diff(RA, vz)).replace('**2','2') + ');\n')
outputfile.write('\n')
outputfile.write('double dDEC_dx = ('  + str(sympy.diff(DEC, x)).replace('**2','2') + ');\n')
outputfile.write('double dDEC_dx = ('  + str(sympy.diff(DEC, y)).replace('**2','2') + ');\n')
outputfile.write('double dDEC_dx = ('  + str(sympy.diff(DEC, z)).replace('**2','2') + ');\n')
outputfile.write('double dDEC_dx = ('  + str(sympy.diff(DEC, vx)).replace('**2','2') + ');\n')
outputfile.write('double dDEC_dx = ('  + str(sympy.diff(DEC, vy)).replace('**2','2') + ');\n')
outputfile.write('double dDEC_dx = ('  + str(sympy.diff(DEC, vz)).replace('**2','2') + ');\n')
outputfile.write('\n')
outputfile.write('double dv_dx = ('  + str(sympy.diff(v, x)).replace('**2','2') + ');\n')
outputfile.write('double dv_dx = ('  + str(sympy.diff(v, y)).replace('**2','2') + ');\n')
outputfile.write('double dv_dx = ('  + str(sympy.diff(v, z)).replace('**2','2') + ');\n')
outputfile.write('double dv_dx = ('  + str(sympy.diff(v, vx)).replace('**2','2') + ');\n')
outputfile.write('double dv_dx = ('  + str(sympy.diff(v, vy)).replace('**2','2') + ');\n')
outputfile.write('double dv_dx = ('  + str(sympy.diff(v, vz)).replace('**2','2') + ');\n')
outputfile.write('\n')
outputfile.write('double dvRA_dx = ('  + str(sympy.diff(vRA, x)).replace('**2','2') + ');\n')
outputfile.write('double dvRA_dx = ('  + str(sympy.diff(vRA, y)).replace('**2','2') + ');\n')
outputfile.write('double dvRA_dx = ('  + str(sympy.diff(vRA, z)).replace('**2','2') + ');\n')
outputfile.write('double dvRA_dx = ('  + str(sympy.diff(vRA, vx)).replace('**2','2') + ');\n')
outputfile.write('double dvRA_dx = ('  + str(sympy.diff(vRA, vy)).replace('**2','2') + ');\n')
outputfile.write('double dvRA_dx = ('  + str(sympy.diff(vRA, vz)).replace('**2','2') + ');\n')
outputfile.write('\n')
outputfile.write('double dvDEC_dx = ('  + str(sympy.diff(vDEC, x)).replace('**2','2') + ');\n')
outputfile.write('double dvDEC_dx = ('  + str(sympy.diff(vDEC, y)).replace('**2','2') + ');\n')
outputfile.write('double dvDEC_dx = ('  + str(sympy.diff(vDEC, z)).replace('**2','2') + ');\n')
outputfile.write('double dvDEC_dx = ('  + str(sympy.diff(vDEC, vx)).replace('**2','2') + ');\n')
outputfile.write('double dvDEC_dx = ('  + str(sympy.diff(vDEC, vy)).replace('**2','2') + ');\n')
outputfile.write('double dvDEC_dx = ('  + str(sympy.diff(vDEC, vz)).replace('**2','2') + ');\n')


outputfile.close()