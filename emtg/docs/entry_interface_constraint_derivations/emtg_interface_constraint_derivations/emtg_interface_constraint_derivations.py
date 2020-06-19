from sympy import *
r_bcf, r_pa = symbols('r_bcf r_pa')
th1, th2, th3 = symbols('th1 th2 th3')
x_pa, y_pa, z_pa = symbols('x_pa y_pa z_pa')
x_bcf, y_bcf, z_bcf = symbols('x_bcf y_bcf z_bcf')
R_bcf2fprime = Array([[cos(th1), sin(th1), 0], [-sin(th1), cos(th1), 0], [0, 0, 1]]).tomatrix() # 3 axis rotation
R_fprime2fdprime = Array([[1, 0, 0], [0, cos(th2), sin(th2)], [0, -sin(th2), cos(th2)]]).tomatrix() # 1 axis rotation
R_fdprime2pa = Array([[cos(th3), sin(th3), 0], [-sin(th3), cos(th3), 0], [0, 0, 1]]).tomatrix() # 3 axis rotation
R_bcf2pa = R_fdprime2pa*R_fprime2fdprime*R_bcf2fprime # composite rotation
R_pa2bcf = R_bcf2pa.T # PA->BCF
r_pa = Matrix([[x_pa], [y_pa], [z_pa]])
r_bcf = Matrix([[x_bcf], [y_bcf], [z_bcf]])
