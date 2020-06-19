#Kepler solver for PyEMTG
#Ryne Beeson 7-9-2014

#  Orbital Elements -> r, v Arrays
def coe2rv(oe, mu):
    #  import statements
    from math  import cos, sin, sqrt
    from numpy import matrix, array, zeros
    
    #  if 'a' is set to zero, then body is at rest in the
    #+ current frame of reference.
    #+ return a zero vector for r and v
    if oe[0] == 0.0: return zeros(3,1), zeros(3,1)
    
    a  = oe[0]
    e  = oe[1]
    i  = oe[2]
    Om = oe[3]
    om = oe[4]
    f  = oe[5]
    
    p  = a*(1 - e*e)
    r  = p/(1 + e*cos(f))
    rv = matrix([r*cos(f), r*sin(f),   0])
    vv = matrix([-sin(f),  e + cos(f), 0])
    vv = sqrt(mu/p)*vv
    
    c0 = cos(Om); s0 = sin(Om)
    co = cos(om); so = sin(om)
    ci = cos(i);  si = sin(i)
    
    R  = matrix([[c0*co - s0*so*ci, -c0*so - s0*co*ci,  s0*si],
                 [s0*co + c0*so*ci, -s0*so + c0*co*ci, -c0*si],
                 [so*si,             co*si,             ci]])
                 
    ri = array(R*rv.T); ri = ri.reshape(3)
    vi = array(R*vv.T); vi = vi.reshape(3)
    
    return ri, vi


#  laguerre_conway solver
def laguerre_conway(e, M):
    #  import statements
    from math  import sqrt, sin, cos
    from numpy import abs
    #  generate an initial guess for eccentric anomaly (E)
    En = (M*(1 - sin(M + e)) + (M + e)*sin(M)) / (1 + sin(M) - sin(M + e))
    n  = 4
    #  for-loop
    for i in range(1,20):
        f      =  M - En + e*sin(En)
        fdash  = -1 + e*cos(En)
        fddash = -e*sin(En)
        g = sqrt(((n - 1)**2) * (fdash**2) - n*(n - 1)*f*fddash)
        if fdash > 0:
            En1 = En - (n*f/(fdash + g))
        else:
            En1 = En - (n*f/(fdash - g))
        #  calculate error
        error = abs(En1 - En)
        En = En1
        if error <= 1E-4:
            break
    
    #  return the eccentric anomaly
    return En

def safe_acos(num):
    from math import acos, pi
    if num > 1.0:
        return 0.0
    elif num < -1.0:
        return pi
    else:
        return acos(num)


#  solve kepler's equation using laguerre_conway
def kepler(sma, ecc, inc, RAAN, AOP, MA, reference_time, epoch, mu):
    #  input delta_time in (secs)
    #  import statements
    from math import pi, sqrt, tan, atan
    #  delta_time
    delta_time = (epoch - reference_time)*86400.0
    n = 1.0 / sqrt((sma**3)/mu)
    MA = (MA + n*delta_time) % (2.0*pi)
    #  update E
    E   = laguerre_conway(ecc, MA) % (2.0*pi)
    #  calculate f
    f = 2.0*atan(tan(E/2.0) * sqrt((1 + ecc) / (1 - ecc))) % (2.0*pi)
    #  coe2rv
    r, v = coe2rv([sma, ecc, inc, RAAN, AOP, f], mu)
    #  return r and v

    return r, v

def cart2kep(r,v,mu):
    from math  import pi, acos
    from numpy import linalg, cross, dot
    
    r_mag = linalg.norm(r)
    v_mag = linalg.norm(v)
    h = cross(r,v)    
    h_mag = linalg.norm(h)    
    h_hat = [x/h_mag for x in h]
    n = cross([0,0,1],h_hat)    
    n_mag = linalg.norm(n)
    crossVH = cross(v,h)
    e_vec = [crossVH[0]/mu-r[0]/r_mag,crossVH[1]/mu-r[1]/r_mag,crossVH[2]/mu-r[2]/r_mag]    
    E = v_mag*v_mag/2 - mu/r_mag
    
    a = -mu/(2*E)  
    e = linalg.norm(e_vec) 
    i = 0;
    OM = 0;
    omega = 0;
    nu = 0;

    delta = 1e-7
    if (e > delta and n_mag>delta):
        i = safe_acos(h[2]/h_mag)
        OM = safe_acos(n[0]/n_mag)
        if (n[1] < 0):
            OM = 2*pi - OM
            
        omega = safe_acos(dot(n,e_vec)/(n_mag*e))
        if (e_vec[2] < 0):
            omega = 2*pi - omega
            
        nu = safe_acos(dot(e_vec,r)/(e*r_mag))
        if (dot(r,v) < 0):
            nu = 2*pi - nu
            
    elif (e < delta and n_mag > delta):
        i = safe_acos(h[2]/h_mag)
        omega = 0
    
        OM = safe_acos(n[0]/n_mag)
        if (n[1] < 0):
            OM = 2*pi - OM   
    
        nu = safe_acos(dot(n,r)/(n_mag*r_mag))
        if (r[2] < 0):
            nu = 2*pi - nu

    elif (e < delta and n_mag < delta):
        i = 0
        omega = 0
        OM = 0
        nu = safe_acos(r[2]/r_mag)
        if (r[2] < 0):
            nu = 2*pi - nu
            
    elif (1 > e and e > delta and n_mag < delta):
        omega = safe_acos(e_vec[0]/e)
        if (e_vec[1] < 0):
            omega = 2*pi - omega
            
        i = 0
        OM = 0
    
        nu = safe_acos(dot(e_vec,r,3)/(e*r_mag))
        if (dot(r,v,3) < 0):
            nu = 2*pi - nu

    elif (e > 1 and n_mag < delta):
        omega = safe_acos(e_vec[0]/e)
        if (e_vec[1] < 0):    
            omega = 2*pi - omega
        
        i = 0
        OM = 0
        nu = safe_acos(dot(e_vec,r)/(e*r_mag))
        if (dot(r,v,3) < 0):
            nu = 2*pi - nu
    
    oe = [a,e,i,OM,omega,nu]
    
    return oe