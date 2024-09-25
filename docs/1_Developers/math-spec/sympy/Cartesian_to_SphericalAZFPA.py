import sympy
from sympy import cos, sin, sqrt, asin, atan2, acos
from sympy import pi

outputfile = open('Cartesian_to_SphericalAZFPA.out', 'w')

#actual angle
x, y, z, vx, vy, vz = sympy.symbols('x y z vx vy vz')

r = (x * x + y * y + z * z)**(1/2);
RA = atan2(y, x);
DEC = asin(z / r);
v = sqrt(vx * vx + vy * vy + vz * vz);
FPA = acos( (x*vx + y*vy + z*vz) / r / v )
        
#azimuth is complicated
xhat = sympy.Matrix([cos(RA)*cos(DEC), sin(RA)*cos(DEC), sin(DEC)])
yhat = sympy.Matrix([cos(RA + pi / 2.0), sin(RA + pi / 2.0), 0.0])
zhat = sympy.Matrix([-cos(RA)*sin(DEC), -sin(RA)*sin(DEC), cos(DEC)])
R = (xhat.row_join(yhat).row_join(zhat))
V = sympy.Matrix([vx, vy, vz])
Vprime = R * V
AZ = atan2(Vprime[1], Vprime[2])

u = (x*vx + y*vy + z*vz) / r / v
outputfile.write('du_dx = ' + str(sympy.simplify(sympy.diff(u, x))) + '\n')
outputfile.write('\n')

outputfile.write('double dr_dx = ('  + str(sympy.diff(r, x)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dr_dx = ('  + str(sympy.diff(r, y)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dr_dx = ('  + str(sympy.diff(r, z)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dr_dx = ('  + str(sympy.diff(r, vx)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dr_dx = ('  + str(sympy.diff(r, vy)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dr_dx = ('  + str(sympy.diff(r, vz)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('\n')
outputfile.write('double dRA_dx = ('  + str(sympy.diff(RA, x)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dRA_dx = ('  + str(sympy.diff(RA, y)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dRA_dx = ('  + str(sympy.diff(RA, z)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dRA_dx = ('  + str(sympy.diff(RA, vx)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dRA_dx = ('  + str(sympy.diff(RA, vy)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dRA_dx = ('  + str(sympy.diff(RA, vz)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('\n')
outputfile.write('double dDEC_dx = ('  + str(sympy.diff(DEC, x)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dDEC_dx = ('  + str(sympy.diff(DEC, y)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dDEC_dx = ('  + str(sympy.diff(DEC, z)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dDEC_dx = ('  + str(sympy.diff(DEC, vx)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dDEC_dx = ('  + str(sympy.diff(DEC, vy)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('double dDEC_dx = ('  + str(sympy.diff(DEC, vz)).replace('**2','2').replace('(x2 + y2 + z2)','(r2)') + ');\n')
outputfile.write('\n')
outputfile.write('double dv_dx = ('  + str(sympy.diff(v, x)).replace('**2','2').replace('(vx2 + vy2 + vz2)','v') + ');\n')
outputfile.write('double dv_dx = ('  + str(sympy.diff(v, y)).replace('**2','2').replace('(vx2 + vy2 + vz2)','v') + ');\n')
outputfile.write('double dv_dx = ('  + str(sympy.diff(v, z)).replace('**2','2').replace('(vx2 + vy2 + vz2)','v') + ');\n')
outputfile.write('double dv_dx = ('  + str(sympy.diff(v, vx)).replace('**2','2').replace('(vx2 + vy2 + vz2)','v') + ');\n')
outputfile.write('double dv_dx = ('  + str(sympy.diff(v, vy)).replace('**2','2').replace('(vx2 + vy2 + vz2)','v') + ');\n')
outputfile.write('double dv_dx = ('  + str(sympy.diff(v, vz)).replace('**2','2').replace('(vx2 + vy2 + vz2)','v') + ');\n')
outputfile.write('\n')
outputfile.write('double dAZ_dx  = ('  + str(sympy.diff(AZ, x )) + ');\n')#.replace('**2','2').replace('(vx2 + vy2 + vz2)','v').replace('(x2 + y2 + z2)','(r2)').replace('(vx*x + vy*y + vz*z)','rdotv') + ');\n')
outputfile.write('double dAZ_dy  = ('  + str(sympy.diff(AZ, y )) + ');\n')#.replace('**2','2').replace('(vx2 + vy2 + vz2)','v').replace('(x2 + y2 + z2)','(r2)').replace('(vx*x + vy*y + vz*z)','rdotv') + ');\n')
outputfile.write('double dAZ_dz  = ('  + str(sympy.diff(AZ, z )) + ');\n')#.replace('**2','2').replace('(vx2 + vy2 + vz2)','v').replace('(x2 + y2 + z2)','(r2)').replace('(vx*x + vy*y + vz*z)','rdotv') + ');\n')
outputfile.write('double dAZ_dvx = ('  + str(sympy.diff(AZ, vx)).replace('**2','2').replace('(vx2 + vy2 + vz2)','v').replace('(x2 + y2 + z2)','(r2)').replace('(vx*x + vy*y + vz*z)','rdotv') + ');\n')
outputfile.write('double dAZ_dvy = ('  + str(sympy.diff(AZ, vy)).replace('**2','2').replace('(vx2 + vy2 + vz2)','v').replace('(x2 + y2 + z2)','(r2)').replace('(vx*x + vy*y + vz*z)','rdotv') + ');\n')
outputfile.write('double dAZ_dvz = ('  + str(sympy.diff(AZ, vz)).replace('**2','2').replace('(vx2 + vy2 + vz2)','v').replace('(x2 + y2 + z2)','(r2)').replace('(vx*x + vy*y + vz*z)','rdotv') + ');\n')
outputfile.write('\n')
outputfile.write('double dFPA_dx  = ('  + str(sympy.diff(FPA, x )).replace('(x**2 + y**2 + z**2)','r2').replace('sqrt(vx**2 + vy**2 + vz**2)','v').replace('(vx*x + vy*y + vz*z)','rdotv').replace('(vx**2 + vy**2 + vz**2)','v2').replace('rdotv**2','rdotv2').replace('*r2**(-0.5)','/sqrt(r)').replace('*r2**(-1.5)','/sqrt(r3)') + ');\n')
outputfile.write('double dFPA_dy  = ('  + str(sympy.diff(FPA, y )).replace('(x**2 + y**2 + z**2)','r2').replace('sqrt(vx**2 + vy**2 + vz**2)','v').replace('(vx*x + vy*y + vz*z)','rdotv').replace('(vx**2 + vy**2 + vz**2)','v2').replace('rdotv**2','rdotv2').replace('*r2**(-0.5)','/sqrt(r)').replace('*r2**(-1.5)','/sqrt(r3)') + ');\n')
outputfile.write('double dFPA_dz  = ('  + str(sympy.diff(FPA, z )).replace('(x**2 + y**2 + z**2)','r2').replace('sqrt(vx**2 + vy**2 + vz**2)','v').replace('(vx*x + vy*y + vz*z)','rdotv').replace('(vx**2 + vy**2 + vz**2)','v2').replace('rdotv**2','rdotv2').replace('*r2**(-0.5)','/sqrt(r)').replace('*r2**(-1.5)','/sqrt(r3)') + ');\n')
outputfile.write('double dFPA_dvx = ('  + str(sympy.diff(FPA, vx)).replace('(x**2 + y**2 + z**2)','r2').replace('sqrt(vx**2 + vy**2 + vz**2)','v').replace('(vx*x + vy*y + vz*z)','rdotv').replace('(vx**2 + vy**2 + vz**2)','v2').replace('rdotv**2','rdotv2').replace('*r2**(-0.5)','/sqrt(r)').replace('*r2**(-1.5)','/sqrt(r3)') + ');\n')
outputfile.write('double dFPA_dvy = ('  + str(sympy.diff(FPA, vy)).replace('(x**2 + y**2 + z**2)','r2').replace('sqrt(vx**2 + vy**2 + vz**2)','v').replace('(vx*x + vy*y + vz*z)','rdotv').replace('(vx**2 + vy**2 + vz**2)','v2').replace('rdotv**2','rdotv2').replace('*r2**(-0.5)','/sqrt(r)').replace('*r2**(-1.5)','/sqrt(r3)') + ');\n')
outputfile.write('double dFPA_dvz = ('  + str(sympy.diff(FPA, vz)).replace('(x**2 + y**2 + z**2)','r2').replace('sqrt(vx**2 + vy**2 + vz**2)','v').replace('(vx*x + vy*y + vz*z)','rdotv').replace('(vx**2 + vy**2 + vz**2)','v2').replace('rdotv**2','rdotv2').replace('*r2**(-0.5)','/sqrt(r)').replace('*r2**(-1.5)','/sqrt(r3)') + ');\n')

#outputfile.write('double dFPA_dx  = ('  + str(sympy.diff(FPA, x )).replace('**2','2').replace('(vx2 + vy2 + vz2)','v').replace('(x2 + y2 + z2)','(r2)').replace('(vx*x + vy*y + vz*z)','rdotv') + ');\n')
#outputfile.write('double dFPA_dy  = ('  + str(sympy.diff(FPA, y )).replace('**2','2').replace('(vx2 + vy2 + vz2)','v').replace('(x2 + y2 + z2)','(r2)').replace('(vx*x + vy*y + vz*z)','rdotv') + ');\n')
#outputfile.write('double dFPA_dz  = ('  + str(sympy.diff(FPA, z )).replace('**2','2').replace('(vx2 + vy2 + vz2)','v').replace('(x2 + y2 + z2)','(r2)').replace('(vx*x + vy*y + vz*z)','rdotv') + ');\n')
#outputfile.write('double dFPA_dvx = ('  + str(sympy.diff(FPA, vx)).replace('**2','2').replace('(vx2 + vy2 + vz2)','v').replace('(x2 + y2 + z2)','(r2)').replace('(vx*x + vy*y + vz*z)','rdotv') + ');\n')
#outputfile.write('double dFPA_dvy = ('  + str(sympy.diff(FPA, vy)).replace('**2','2').replace('(vx2 + vy2 + vz2)','v').replace('(x2 + y2 + z2)','(r2)').replace('(vx*x + vy*y + vz*z)','rdotv') + ');\n')
#outputfile.write('double dFPA_dvz = ('  + str(sympy.diff(FPA, vz)).replace('**2','2').replace('(vx2 + vy2 + vz2)','v').replace('(x2 + y2 + z2)','(r2)').replace('(vx*x + vy*y + vz*z)','rdotv') + ');\n')


outputfile.close()