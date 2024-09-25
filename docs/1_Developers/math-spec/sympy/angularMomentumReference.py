import sympy

outputfile = open('angularMomentumReference.out', 'w')

rx, ry, rz, vx, vy, vz = sympy.symbols('rx ry rz vx vy vz')
rrx, rry, rrz = sympy.symbols('rrx rry rrz')
A, B = sympy.symbols('A B')

R = sympy.Matrix([rx, ry, rz])
V = sympy.Matrix([vx, vy, vz])

# r = sympy.Function('r')(rx, ry, rz)
# v = sympy.Function('v')(vx, vy, vz)

# Rref = sympy.Matrix([rrx, rry, rrz])

# R_probe_reference = Rref - R

# H = sympy.Cross(R, V)

# zhat = sympy.Matrix([0, 0, 1])

# angle = sympy.atan2( sympy.Dot(sympy.Cross(H, R_probe_reference), zhat), sympy.Dot(H, R_probe_reference) )

# outputfile.write('angle = ' + str(angle).replace('r(rx, ry, rz)','r').replace('v(vx, vy, vz)','v') + '\n')
outputfile.write('\n')
outputfile.write('diff(atan2(A,B),A) = ' + str(sympy.diff(sympy.atan2(A,B),A)) + '\n')
outputfile.write('diff(atan2(A,B),B) = ' + str(sympy.diff(sympy.atan2(A,B),B)) + '\n')
outputfile.write('\n')


outputfile.close()