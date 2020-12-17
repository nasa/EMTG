#functions to convert between state representations
#takes a vector, writes a vector
#Jacob Englander, borrowing lots of things from Noble Hatten and Ryne Beeson
#6-11-2020

#All state representations convert to and from Cartesian. Cross-conversions go through conversion
#at the bottom we have a utility to convert a decision vector
#
#need the following:
#
#CartesiantoCOE                             XXXXXXXXXX
#CartesiantoSphericalAZFPA                  XXXXXXXXXX
#CartesiantoSphericalRADEC                  XXXXXXXXXX
#CartesiantoIncomingBplane                  XXXXXXXXXX
#CartesiantoOutgoingBplane                  XXXXXXXXXX
#CartesiantoMEE                             XXXXXXXXXX
#COEtoCartesian                             XXXXXXXXXX
#COEtoSphericalAZFPA                        XXXXXXXXXX
#COEtoSphericalRADEC                        XXXXXXXXXX
#COEtoIncomingBplane                        XXXXXXXXXX
#COEtoOutgoingBplane                        XXXXXXXXXX
#COEtoMEE                                   XXXXXXXXXX
#SphericalAZFPAtoCOE                        XXXXXXXXXX
#SphericalAZFPAtoCartesian                  XXXXXXXXXX
#SphericalAZFPAtoSphericalRADEC             XXXXXXXXXX
#SphericalAZFPAtoIncomingBplane             XXXXXXXXXX
#SphericalAZFPAtoOutgoingBplane             XXXXXXXXXX
#SphericalAZFPAtoMEE                        XXXXXXXXXX
#SphericalRADECtoCOE                        XXXXXXXXXX
#SphericalRADECtoSphericalAZFPA             XXXXXXXXXX
#SphericalRADECtoCartesian                  XXXXXXXXXX
#SphericalRADECtoIncomingBplane             XXXXXXXXXX
#SphericalRADECtoOutgoingBplane             XXXXXXXXXX
#SphericalRADECtoMEE                        XXXXXXXXXX
#IncomingBplanetoCOE                        XXXXXXXXXX
#IncomingBplanetoSphericalAZFPA             XXXXXXXXXX
#IncomingBplanetoSphericalRADEC             XXXXXXXXXX
#IncomingBplanetoCartesian                  XXXXXXXXXX
#IncomingBplanetoOutgoingBplane             XXXXXXXXXX
#IncomingBplanetoMEE                        XXXXXXXXXX
#OutgoingBplanetoCOE                        XXXXXXXXXX
#OutgoingBplanetoSphericalAZFPA             XXXXXXXXXX
#OutgoingBplanetoSphericalRADEC             XXXXXXXXXX
#OutgoingBplanetoIncomingBplane             XXXXXXXXXX
#OutgoingBplanetoCartesian                  XXXXXXXXXX
#OutgoingBplanetoMEE                        XXXXXXXXXX
#MEEtoCOE                                   XXXXXXXXXX
#MEEtoSphericalAZFPA                        XXXXXXXXXX
#MEEtoSphericalRADEC                        XXXXXXXXXX
#MEEtoIncomingBplane                        XXXXXXXXXX
#MEEtoOutgoingBplane                        XXXXXXXXXX
#MEEtoCartesian                             XXXXXXXXXX

class StateConverter(object):
    def COEtoCartesian(self, stateCOE, mu=1.0):
        #  import statements
        from math  import cos, sin, sqrt
        from numpy import matrix, array, zeros
    
        #  if 'a' is set to zero, then body is at rest in the
        #+ current frame of reference.
        #+ return a zero vector for r and v
        if stateCOE[0] == 0.0: return zeros(6,1)
    
        a  = stateCOE[0]
        e  = stateCOE[1]
        i  = stateCOE[2]
        Om = stateCOE[3]
        om = stateCOE[4]
        f  = stateCOE[5]
    
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

        stateCartesian = [ri[0], ri[1], ri[2], vi[0], vi[1], vi[2]]
    
        return stateCartesian

    def CartesiantoCOE(self, stateCartesian,mu=1.0):
        from kepler import safe_acos
        from math  import pi, acos
        from numpy import linalg, cross, dot, array

        r = array([stateCartesian[0], stateCartesian[1], stateCartesian[2]])
        v = array([stateCartesian[3], stateCartesian[4], stateCartesian[5]])
    
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

        #deal with the "true anomaly is zero but truncation error says not" case
        from math import fabs
        if fabs(nu - 2*pi) < 1.0e-9:
            nu = 0.0
    
        stateCOE = [a,e,i,OM,omega,nu]
    
        return stateCOE

    def SphericalRADECtoCartesian(self, stateSphericalRADEC): #doesn't need mu
        from math import sin, cos
        r = float(stateSphericalRADEC[0])                                                                               
        RA = float(stateSphericalRADEC[1])                                                                              
        DEC = float(stateSphericalRADEC[2])                                                                             
        v = float(stateSphericalRADEC[3])                                                                               
        vRA = float(stateSphericalRADEC[4])                                                                                 
        vDEC = float(stateSphericalRADEC[5])                                                                            
        #convert to cartesian                                                                                               
        cosRA = cos(RA)                                                                                                     
        sinRA = sin(RA)                                                                                                     
        cosDEC = cos(DEC)                                                                                                   
        sinDEC = sin(DEC)                                                                                                   
        cosvRA = cos(vRA)                                                                                                   
        sinvRA = sin(vRA)                                                                                                   
        cosvDEC = cos(vDEC)                                                                                                 
        sinvDEC = sin(vDEC)                                                                                                 
                                                                                                                                        
        x = r * cosRA * cosDEC                                                                                              
        y = r * sinRA * cosDEC                                                                                              
        z = r * sinDEC;                                                                                                     
        vx = v * cosvRA * cosvDEC                                                                                         
        vy = v * sinvRA * cosvDEC                                                                                         
        vz = v * sinvDEC;    

        return [x, y, z, vx, vy, vz]

    def CartesiantoSphericalRADEC(self, stateCartesian): #doesn't need mu
        from math import asin, atan2
        x =  float(stateCartesian[0])                                                                                                
        y =  float(stateCartesian[1])                                                                                                
        z =  float(stateCartesian[2])                                                                                                
        vx = float(stateCartesian[3])                                                                                               
        vy = float(stateCartesian[4])                                                                                               
        vz = float(stateCartesian[5])                                                                                               
                                                                                                                                     
        #compute the correct state                                                                                       
        r = (x**2 + y**2 + z**2)**0.5                                                                                    
        v = (vx**2 + vy**2 + vz**2)**0.5                                                                                 
        RA = atan2(y, x)                                                                                                 
        DEC = asin(z / r)                                                                                                
        vRA = atan2(vy, vx)                                                                                            
        vDEC = asin(vz / v)    

        return [r, RA, DEC, v, vRA, vDEC]

    def SphericalAZFPAtoCartesian(self, stateSphericalAZFPA): #doesn't need mu
        from math import sin, cos
        r = float(stateSphericalAZFPA[0])                                                                               
        RA = float(stateSphericalAZFPA[1])                                                                              
        DEC = float(stateSphericalAZFPA[2])                                                                             
        v = float(stateSphericalAZFPA[3])                                                                               
        AZ = float(stateSphericalAZFPA[4])                                                                                 
        FPA = float(stateSphericalAZFPA[5])
                                                                                                   
        cosRA = cos(RA)                                                                                                     
        sinRA = sin(RA)                                                                                                     
        cosDEC = cos(DEC)                                                                                                   
        sinDEC = sin(DEC)                                                                                                   
        cosAZ = cos(AZ)                                                                                                     
        sinAZ = sin(AZ)                                                                                                     
        cosFPA = cos(FPA)                                                                                                   
        sinFPA = sin(FPA)                                                                                                   
                                                                                                                                        
        x = r * cosRA * cosDEC                                                                                              
        y = r * sinRA * cosDEC                                                                                              
        z = r * sinDEC;                                                                                                     
        vx = -v * (sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA);                                  
        vy =  v * (sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA);                                  
        vz =  v * (cosFPA*sinDEC + cosDEC * cosAZ*sinFPA);   

        return [x, y, z, vx, vy, vz]

    def CartesiantoSphericalAZFPA(self, stateCartesian): #doesn't need mu
        from math import asin, atan2, acos, cos, sin, pi
        import numpy

        x =  float(stateCartesian[0])                                                                                                
        y =  float(stateCartesian[1])                                                                                                
        z =  float(stateCartesian[2])                                                                                                
        vx = float(stateCartesian[3])                                                                                               
        vy = float(stateCartesian[4])                                                                                               
        vz = float(stateCartesian[5]) 
    
        #compute the correct state                                                                                  
        r = (x**2 + y**2 + z**2)**0.5                                                                               
        v = (vx**2 + vy**2 + vz**2)**0.5                                                                            
        RA = atan2(y, x)                                                                                            
        DEC = asin(z / r)                                                                                           
                                                                                                                                
        #compute AZ and FPA                                                                                         
        FPA = acos( (x*vx + y*vy + z*vz) / r / v )                                                            
                                                                                                                                
        #azimuth is complicated                                                                                     
        xhat = numpy.matrix([cos(RA)*cos(DEC), sin(RA)*cos(DEC), sin(DEC)]).T                                       
        yhat = numpy.matrix([cos(RA + pi / 2.0), sin(RA + pi / 2.0), 0.0]).T                                        
        zhat = numpy.matrix([-cos(RA)*sin(DEC), -sin(RA)*sin(DEC), cos(DEC)]).T                                     
        R = numpy.hstack([xhat, yhat, zhat]).T                                                                      
        V = numpy.matrix([vx, vy, vz]).T / v                                                                      
        Vprime = R * V                                                                                              
        AZ = atan2(Vprime[1], Vprime[2])    

        return [r, RA, DEC, v, AZ, FPA]

    def CartesiantoIncomingBplane(self, stateCartesian, mu=1.0):
        import numpy as np
        r = np.array([stateCartesian[0], stateCartesian[1], stateCartesian[2]])
        v = np.array([stateCartesian[3], stateCartesian[4], stateCartesian[5]])

        from posVel2BPlane import posVel2BPlane
        bplane = posVel2BPlane()

        e = bplane.eVector(r, v, mu)
        h = bplane.hVector(r, v)
        S = bplane.sVector(r, v, mu)
        B = bplane.bVector(r, v, mu)
        R = bplane.rVector(r, v, mu)
        T = bplane.tVector(r, v, mu)
        BR = bplane.bDotR(r, v, mu)
        BT = bplane.bDotT(r, v, mu)
        theta = bplane.bTheta(r, v, mu)
        Bmag = bplane.bScalar(r, v, mu)
        rp = bplane.rPeri(r, v, mu)
        vinf = bplane.vInfMag(r, v, mu)
        RA = bplane.vInfRA(r, v, mu)
        Dec = bplane.vInfDec(r, v, mu)
        trueAnomaly = bplane.trueAnomaly(r, v, mu)

        return [vinf, RA, Dec, Bmag, theta, trueAnomaly]

    def CartesiantoOutgoingBplane(self, stateCartesian, mu=1.0):
        import numpy as np
        r = np.array([stateCartesian[0], stateCartesian[1], stateCartesian[2]])
        v = np.array([stateCartesian[3], stateCartesian[4], stateCartesian[5]])

        from posVel2BPlaneOut import posVel2BPlaneOut
        bplane = posVel2BPlaneOut()

        e = bplane.eVector(r, v, mu)
        h = bplane.hVector(r, v)
        S = bplane.sVector(r, v, mu)
        B = bplane.bVector(r, v, mu)
        R = bplane.rVector(r, v, mu)
        T = bplane.tVector(r, v, mu)
        BR = bplane.bDotR(r, v, mu)
        BT = bplane.bDotT(r, v, mu)
        theta = bplane.bTheta(r, v, mu)
        Bmag = bplane.bScalar(r, v, mu)
        rp = bplane.rPeri(r, v, mu)
        vinf = bplane.vInfMag(r, v, mu)
        RA = bplane.vInfRA(r, v, mu)
        Dec = bplane.vInfDec(r, v, mu)
        trueAnomaly = bplane.trueAnomaly(r, v, mu)

        return [vinf, RA, Dec, Bmag, theta, trueAnomaly]

    def CartesiantoMEE(self, stateCartesian, mu=1.0):
        stateCOE = self.CartesiantoCOE(stateCartesian, mu)

        return self.COEtoMEE(stateCOE, mu)

    def IncomingBplanetoCartesian(self, stateIncomingBplane, mu=1.0):
        import numpy as np
        from BPlane2PosVel import BPlane2PosVel
        bplane = BPlane2PosVel()
    
        eMag = bplane.eMag(stateIncomingBplane, mu)
        s = bplane.sVector(stateIncomingBplane)
        hMag = bplane.hMag(stateIncomingBplane)
        t = bplane.tVector(stateIncomingBplane)
        r = bplane.rVector(stateIncomingBplane)
        BR = bplane.bDotR(stateIncomingBplane)
        BT = bplane.bDotT(stateIncomingBplane)
        B = bplane.bVector(stateIncomingBplane)
        hUnit = bplane.hUnit(stateIncomingBplane)
        h = bplane.hVector(stateIncomingBplane)
        e = bplane.eVector(stateIncomingBplane, mu)
        rp = bplane.positionVector(stateIncomingBplane, mu)
        vp = bplane.velocityVector(stateIncomingBplane, mu)
        r = bplane.positionVector(stateIncomingBplane, mu)
        v = bplane.velocityVector(stateIncomingBplane, mu)

        return [r[0], r[1], r[2], v[0], v[1], v[2]]

    def OutgoingBplanetoCartesian(self, stateOutgoingBplane, mu=1.0):
        import numpy as np
        from BPlane2PosVelOut import BPlane2PosVelOut
        bplane = BPlane2PosVelOut()
    
        eMagBack = bplane.eMag(stateOutgoingBplane, mu)
        sBack = bplane.sVector(stateOutgoingBplane)
        hMagBack = bplane.hMag(stateOutgoingBplane)
        tBack = bplane.tVector(stateOutgoingBplane)
        rBack = bplane.rVector(stateOutgoingBplane)
        BRBack = bplane.bDotR(stateOutgoingBplane)
        BTBack = bplane.bDotT(stateOutgoingBplane)
        BBack = bplane.bVector(stateOutgoingBplane)
        hUnitBack = bplane.hUnit(stateOutgoingBplane)
        hBack = bplane.hVector(stateOutgoingBplane)
        eBack = bplane.eVector(stateOutgoingBplane, mu)
        rpBack = bplane.positionVector(stateOutgoingBplane, mu)
        vpBack = bplane.velocityVector(stateOutgoingBplane, mu)
        r = bplane.positionVector(stateOutgoingBplane, mu)
        v = bplane.velocityVector(stateOutgoingBplane, mu)

        return [r[0], r[1], r[2], v[0], v[1], v[2]]

    def MEEtoCOE(self, stateMEE, mu=1.0, prograde=True):
        from math import sqrt, atan2, atan, pi
        P = stateMEE[0]
        F = stateMEE[1]
        G = stateMEE[2]
        H = stateMEE[3]
        K = stateMEE[4]
        L = stateMEE[5]

        SMA = P / (1.0 - F * F - G * G)
        ECC = sqrt(F * F + G * G)

        A = sqrt(H * H + K * K)
        B = atan2(G / ECC, F / ECC)

        RAAN = atan2(K / A, H / A)

        #how do you check an MEE state for whether it is prograde or retrograde?

        INC = []
        AOP = []
        if prograde:
            INC = 2.0 * atan(A);
            AOP = B - RAAN;
        else:
            INC = pi - 2.0 * atan(A);
            AOP = B + RAAN

        TA = L - B

        return [SMA, ECC, INC, RAAN, AOP, TA]

    def COEtoMEE(self, stateCOE, mu=1.0):
        from math import sin, cos, tan
        SMA =  stateCOE[0]
        ECC =  stateCOE[1]
        INC =  stateCOE[2]
        RAAN = stateCOE[3]
        AOP =  stateCOE[4]
        TA =   stateCOE[5]

        tanINCover2 = tan(INC / 2.0)
        cosINCover2 = cos(INC / 2.0)
        cosRAANplusAOP = cos(RAAN + AOP)
        sinRAANplusAOP = sin(RAAN + AOP)
        cosRAAN = cos(RAAN)
        sinRAAN = sin(RAAN)
            
        P = SMA * (1 - ECC * ECC)
        F = ECC * cosRAANplusAOP
        G = ECC * sinRAANplusAOP
        H = tanINCover2 * cosRAAN
        K = tanINCover2 * sinRAAN
        L = RAAN + AOP + TA

        return [P, F, G, H, K, L]

    def COEtoSphericalAZFPA(self, stateCOE, mu=1.0):
        stateCartesian = self.COEtoCartesian(stateCOE, mu)

        stateSphericalAZFPA = self.CartesiantoSphericalAZFPA(stateCartesian)

        return stateSphericalAZFPA

    def COEtoSphericalRADEC(self, stateCOE, mu=1.0):
        stateCartesian = self.COEtoCartesian(stateCOE, mu)

        stateSphericalRADEC = self.CartesiantoSphericalRADEC(stateCartesian)

        return stateSphericalRADEC

    def COEtoIncomingBplane(self, stateCOE, mu=1.0):
        stateCartesian = self.COEtoCartesian(stateCOE, mu)

        stateIncomingBplane = self.CartesiantoIncomingBplane(stateCartesian, mu)
    
        return stateIncomingBplane

    def COEtoOutgoingBplane(self, stateCOE, mu=1.0):
        stateCartesian = self.COEtoCartesian(stateCOE, mu)

        stateOutgoingBplane = self.CartesiantoOutgoingBplane(stateCartesian, mu)
    
        return stateOutgoingBplane

    def SphericalAZFPAtoCOE(self, stateSphericalAZFPA, mu=1.0):
        stateCartesian = self.SphericalAZFPAtoCartesian(stateSphericalAZFPA)

        stateCOE = self.CartesiantoCOE(stateCartesian, mu)

        return stateCOE

    def SphericalAZFPAtoSphericalRADEC(self, stateSphericalAZFPA):
        stateCartesian = self.SphericalAZFPAtoCartesian(stateSphericalAZFPA)

        stateSphericalRADEC = self.CartesiantoSphericalRADEC(stateCartesian)

        return stateSphericalRADEC

    def SphericalAZFPAtoIncomingBplane(self, stateSphericalAZFPA, mu=1.0):
        stateCartesian = self.SphericalAZFPAtoCartesian(stateSphericalAZFPA, mu)

        stateIncomingBplane = self.CartesiantoIncomingBplane(stateCartesian, mu)
    
        return stateIncomingBplane

    def SphericalAZFPAtoOutgoingBplane(self, stateSphericalAZFPA, mu=1.0):
        stateCartesian = self.SphericalAZFPAtoCartesian(stateSphericalAZFPA)

        stateOutgoingBplane = self.CartesiantoOutgoingBplane(stateCartesian, mu)
    
        return stateOutgoingBplane
    
    def SphericalAZFPAtoMEE(self, stateSphericalAZFPA, mu=1.0):
        stateCOE = self.SphericalAZFPAtoCOE(stateSphericalAZFPA, mu)

        return self.COEtoMEE(stateCOE, mu)

    def SphericalRADECtoCOE(self, stateSphericalRADEC, mu=1.0):
        stateCartesian = self.SphericalRADECtoCartesian(stateSphericalRADEC)

        stateCOE = self.CartesiantoCOE(stateCartesian, mu)

        return stateCOE

    def SphericalRADECtoSphericalAZFPA(self, stateSphericalRADEC):
        stateCartesian = self.SphericalRADECtoCartesian(stateSphericalRADEC)

        stateSphericalAZFPA = self.CartesiantoSphericalAZFPA(stateCartesian)

        return stateSphericalAZFPA

    def SphericalRADECtoIncomingBplane(self, stateSphericalRADEC, mu=1.0):
        stateCartesian = self.SphericalRADECtoCartesian(stateSphericalRADEC)

        stateIncomingBplane = self.CartesiantoIncomingBplane(stateCartesian, mu)
    
        return stateIncomingBplane

    def SphericalRADECtoOutgoingBplane(self, stateSphericalRADEC, mu=1.0):
        stateCartesian = self.SphericalRADECtoCartesian(stateSphericalRADEC)

        stateOutgoingBplane = self.CartesiantoOutgoingBplane(stateCartesian, mu)
    
        return stateOutgoingBplane
    
    def SphericalRADECtoMEE(self, stateSphericalRADEC, mu=1.0):
        stateCOE = self.SphericalRADECtoCOE(stateSphericalRADEC, mu)

        return self.COEtoMEE(stateCOE, mu)

    def IncomingBplanetoCOE(self, stateIncomingBplane, mu=1.0):
        stateCartesian = self.IncomingBplanetoCartesian(stateIncomingBplane, mu)

        stateCOE = self.CartesiantoCOE(stateCartesian, mu)
    
        return stateCOE

    def IncomingBplanetoSphericalAZFPA(self, stateIncomingBplane, mu=1.0):
        stateCartesian = self.IncomingBplanetoCartesian(stateIncomingBplane, mu)

        stateSphericalAZFPA = self.CartesiantoSphericalAZFPA(stateCartesian)
    
        return stateSphericalAZFPA

    def IncomingBplanetoSphericalRADEC(self, stateIncomingBplane, mu=1.0):
        stateCartesian = self.IncomingBplanetoCartesian(stateIncomingBplane, mu)

        stateSphericalRADEC = self.CartesiantoSphericalRADEC(stateCartesian)
    
        return stateSphericalRADEC

    def IncomingBplanetoOutgoingBplane(self, stateIncomingBplane, mu=1.0):
        stateCartesian = self.IncomingBplanetoCartesian(stateIncomingBplane, mu)

        stateOutgoingBplane = self.CartesiantoOutgoingBplane(stateCartesian)
    
        return stateOutgoingBplane    

    def IncomingBplanetoMEE(self, stateIncomingBplane, mu=1.0):
        stateCOE = self.IncomingBplanetoCOE(stateIncomingBplane, mu)

        return self.COEtoMEE(stateCOE, mu)

    def OutgoingBplanetoCOE(self, stateOutgoingBplane, mu=1.0):
        stateCartesian = self.OutgoingBplanetoCartesian(stateOutgoingBplane, mu)

        stateCOE = self.CartesiantoCOE(stateCartesian, mu)
    
        return stateCOE

    def OutgoingBplanetoSphericalAZFPA(self, stateOutgoingBplane, mu=1.0):
        stateCartesian = self.OutgoingBplanetoCartesian(stateOutgoingBplane, mu)

        stateSphericalAZFPA = self.CartesiantoSphericalAZFPA(stateCartesian)
    
        return stateSphericalAZFPA

    def OutgoingBplanetoSphericalRADEC(self, stateOutgoingBplane, mu=1.0):
        stateCartesian = self.OutgoingBplanetoCartesian(stateOutgoingBplane, mu)

        stateSphericalRADEC = self.CartesiantoSphericalRADEC(stateCartesian)
    
        return stateSphericalRADEC

    def OutgoingBplanetoIncomingBplane(self, stateOutgoingBplane, mu=1.0):
        stateCartesian = self.OutgoingBplanetoCartesian(stateOutgoingBplane, mu)

        stateIncomingBplane = self.CartesiantoIncomingBplane(stateCartesian)
    
        return stateIncomingBplane
    
    def OutgoingBplanetoMEE(self, stateOutgoingBplane, mu=1.0):
        stateCOE = self.OutgoingBplanetoCOE(stateOutgoingBplane, mu)

        return self.COEtoMEE(stateCOE, mu)
    
    def MEEtoCartesian(self, stateMEE, mu=1.0, prograde=True):
        stateCOE = self.MEEtoCOE(stateMEE, mu, prograde)

        return self.COEtoCartesian(stateCOE, mu)

    def MEEtoSphericalRADEC(self, stateMEE, mu=1.0, prograde=True):
        stateCOE = self.MEEtoCOE(stateMEE, mu, prograde)

        return self.COEtoSphericalRADEC(stateCOE, mu)
    
    def MEEtoSphericalAZFPA(self, stateMEE, mu=1.0, prograde=True):
        stateCOE = self.MEEtoCOE(stateMEE, mu, prograde)

        return self.COEtoSphericalAZFPA(stateCOE, mu)
    
    def MEEtoIncomingAsymptote(self, stateMEE, mu=1.0, prograde=True):
        stateCOE = self.MEEtoCOE(stateMEE, mu, prograde)

        return self.COEtoIncomingAsymptote(stateCOE, mu)
    
    def MEEtoOutgoingAsymptote(self, stateMEE, mu=1.0, prograde=True):
        stateCOE = self.MEEtoCOE(stateMEE, mu, prograde)

        return self.COEtoOutgoingAsymptote(stateCOE, mu)
    ################################################################################
    def convertDecisionVector(self, X, desiredStateRep, keywords, mu=1.0): #X is the decision vector, keywords is a list of strings that may include 'and' or 'or', desiredStateRep is a string, mu is mu
        CartesianNames = ['x','y','z','vx','vy','vz']
        COENames = ['SMA','ECC','INC','RAAN','AOP','TA']
        SphericalRADECNames = ['r','RA','DEC','v','vRA','vDEC']
        SphericalAZFPANames = ['r','RA','DEC','v','AZ','FPA']
        IncomingBplaneNames= ['VINFin','RHAin','DHAin','BRADIUSin','BTHETAin','TAin']
        OutgoingBplaneNames= ['VINFout','RHAout','DHAout','BRADIUSout','BTHETAout','TAout']
        MEENames = ['P', 'F', 'G', 'H', 'K', 'L']

        #set the counter to zero
        Xindex = 0

        while Xindex < len(X):                                                                                   
            description = X[Xindex][0]                                                                                    
            prefix = description.split(':')[0]

            if any(keyword in description for keyword in keywords) and not 'PeriapseLaunch' in description: #don't convert PeriapseLaunch
                #something to do here
                if 'left state x' in description: #it's cartesian
                    stateCartesian = [float(X[Xindex][1]), float(X[Xindex+1][1]), float(X[Xindex+2][1]), float(X[Xindex+3][1]), float(X[Xindex+4][1]), float(X[Xindex+5][1])]

                    if desiredStateRep == 'Cartesian':
                        Xindex += 6
                        continue
                    elif desiredStateRep == 'COE':
                        stateCOE = self.CartesiantoCOE(stateCartesian, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + CartesianNames[i], 'left state ' + COENames[i])
                            X[Xindex + i][1] = stateCOE[i]
                        Xindex += 6
                    elif desiredStateRep == 'SphericalRADEC':
                        stateSphericalRADEC = self.CartesiantoSphericalRADEC(stateCartesian)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + CartesianNames[i], 'left state ' + SphericalRADECNames[i])
                            X[Xindex + i][1] = stateSphericalRADEC[i]
                        Xindex += 6
                    elif desiredStateRep == 'SphericalAZFPA':
                        stateSphericalAZFPA = self.CartesiantoSphericalAZFPA(stateCartesian)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + CartesianNames[i], 'left state ' + SphericalAZFPANames[i])
                            X[Xindex + i][1] = stateSphericalAZFPA[i]
                        Xindex += 6
                    elif desiredStateRep == 'MEE':
                        stateMEE = self.CartesiantoMEE(stateCartesian, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + CartesianNames[i], 'left state ' + MEENames[i])
                            X[Xindex + i][1] = stateMEE[i]
                        Xindex += 6
                    elif desiredStateRep == 'IncomingBplane':
                        stateIncomingBplane = self.CartesiantoIncomingBplane(stateCartesian, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + CartesianNames[i], 'left state ' + IncomingBplaneNames[i])
                            X[Xindex + i][1] = stateIncomingBplane[i]
                        Xindex += 6
                    elif desiredStateRep == 'OutgoingBplane':
                        stateOutgoingBplane = self.CartesiantoOutgoingBplane(stateCartesian, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + CartesianNames[i], 'left state ' + OutgoingBplaneNames[i])
                            X[Xindex + i][1] = stateOutgoingBplane[i]
                        Xindex += 6
                        
                elif 'left state AZ' in description: #it's SphericalAZFPA, but we actually have to step back four indices
                    stateSphericalAZFPA = [float(X[Xindex-4][1]), float(X[Xindex-3][1]), float(X[Xindex-2][1]), float(X[Xindex-1][1]), float(X[Xindex][1]), float(X[Xindex+1][1])]

                    if desiredStateRep == 'SphericalAZFPA':
                        Xindex += 2
                        continue
                    elif desiredStateRep == 'COE':
                        stateCOE = self.SphericalAZFPAtoCOE(stateSphericalAZFPA, mu)
                        for i in range(-4, 2):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + SphericalAZFPANames[i + 4], 'left state ' + COENames[i + 4])
                            X[Xindex + i][1] = stateCOE[i + 4]
                        Xindex += 2
                    elif desiredStateRep == 'SphericalRADEC':
                        stateSphericalRADEC = self.SphericalAZFPAtoSphericalRADEC(stateSphericalAZFPA)
                        for i in range(-4, 2):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + SphericalAZFPANames[i + 4], 'left state ' + SphericalRADECNames[i + 4])
                            X[Xindex + i][1] = stateSphericalRADEC[i + 4]
                        Xindex += 2
                    elif desiredStateRep == 'IncomingBplane':
                        stateIncomingBplane = self.SphericalAZFPAtoIncomingBplane(stateSphericalAZFPA, mu)
                        for i in range(-4, 2):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + SphericalAZFPANames[i + 4], 'left state ' + IncomingBplaneNames[i + 4])
                            X[Xindex + i][1] = stateIncomingBplane[i + 4]
                        Xindex += 2
                    elif desiredStateRep == 'OutgoingBplane':
                        stateOutgoingBplane = self.SphericalAZFPAtoOutgoingBplane(stateSphericalAZFPA, mu)
                        for i in range(-4, 2):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + SphericalAZFPANames[i + 4], 'left state ' + OutgoingBplaneNames[i + 4])
                            X[Xindex + i][1] = stateOutgoingBplane[i + 4]
                        Xindex += 2
                    elif desiredStateRep == 'Cartesian':
                        stateCartesian = self.SphericalAZFPAtoCartesian(stateSphericalAZFPA)
                        for i in range(-4, 2):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + SphericalAZFPANames[i + 4], 'left state ' + CartesianNames[i + 4])
                            X[Xindex + i][1] = stateCartesian[i + 4]
                        Xindex += 2
                    elif desiredStateRep == 'MEE':
                        stateMEE = self.SphericalAZFPAtoMEE(stateSphericalAZFPA)
                        for i in range(-4, 2):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + SphericalAZFPANames[i + 4], 'left state ' + MEENames[i + 4])
                            X[Xindex + i][1] = stateCartesian[i + 4]
                        Xindex += 2
                    
                elif 'left state vRA' in description: #it's SphericalRADEC, but we actually have to step back four indices
                    stateSphericalRADEC = [float(X[Xindex-4][1]), float(X[Xindex-3][1]), float(X[Xindex-2][1]), float(X[Xindex-1][1]), float(X[Xindex][1]), float(X[Xindex+1][1])]
                
                    if desiredStateRep == 'SphericalRADEC':
                        Xindex += 2
                        continue                
                    elif desiredStateRep == 'COE':
                        stateCOE = self.SphericalRADECtoCOE(stateSphericalRADEC, mu)
                        for i in range(-4, 2):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + SphericalRADECNames[i + 4], 'left state ' + COENames[i + 4])
                            X[Xindex + i][1] = stateCOE[i + 4]
                        Xindex += 2
                    elif desiredStateRep == 'SphericalAZFPA':
                        stateSphericalAZFPA = self.SphericalRADECtoSphericalAZFPA(stateSphericalRADEC)
                        for i in range(-4, 2):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + SphericalRADECNames[i + 4], 'left state ' + SphericalAZFPANames[i + 4])
                            X[Xindex + i][1] = stateSphericalAZFPA[i + 4]
                        Xindex += 2
                    elif desiredStateRep == 'IncomingBplane':
                        stateIncomingBplane = self.SphericalRADECtoIncomingBplane(stateSphericalRADEC, mu)
                        for i in range(-4, 2):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + SphericalRADECNames[i + 4], 'left state ' + IncomingBplaneNames[i + 4])
                            X[Xindex + i][1] = stateIncomingBplane[i + 4]
                        Xindex += 2
                    elif desiredStateRep == 'OutgoingBplane':
                        stateOutgoingBplane = self.SphericalRADECtoOutgoingBplane(stateSphericalRADEC, mu)
                        for i in range(-4, 2):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + SphericalRADECNames[i + 4], 'left state ' + OutgoingBplaneNames[i + 4])
                            X[Xindex + i][1] = stateOutgoingBplane[i + 4]
                        Xindex += 2
                    elif desiredStateRep == 'Cartesian':
                        stateCartesian = self.SphericalRADECtoCartesian(stateSphericalRADEC)
                        for i in range(-4, 2):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + SphericalRADECNames[i + 4], 'left state ' + CartesianNames[i + 4])
                            X[Xindex + i][1] = stateCartesian[i + 4]
                        Xindex += 2
                    elif desiredStateRep == 'MEE':
                        stateMEE = self.SphericalRADECtoMEE(stateSphericalRADEC)
                        for i in range(-4, 2):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + SphericalRADECNames[i + 4], 'left state ' + MEENames[i + 4])
                            X[Xindex + i][1] = stateCartesian[i + 4]
                        Xindex += 2
            
                elif 'left state SMA' in description: #it's COE
                    stateCOE = [float(X[Xindex][1]), float(X[Xindex+1][1]), float(X[Xindex+2][1]), float(X[Xindex+3][1]), float(X[Xindex+4][1]), float(X[Xindex+5][1])]
                
                    if desiredStateRep == 'COE':
                        Xindex += 6
                        continue
                    elif desiredStateRep == 'Cartesian':
                        stateCartesian = self.COEtoCartesian(stateCOE, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + COENames[i], 'left state ' + CartesianNames[i])
                            X[Xindex + i][1] = stateCartesian[i]
                        Xindex += 6
                    elif desiredStateRep == 'SphericalRADEC':
                        stateSphericalRADEC = self.COEtoSphericalRADEC(stateCOE, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + COENames[i], 'left state ' + SphericalRADECNames[i])
                            X[Xindex + i][1] = stateSphericalRADEC[i]
                        Xindex += 6
                    elif desiredStateRep == 'SphericalAZFPA':
                        stateSphericalAZFPA = self.COEtoSphericalAZFPA(stateCOE, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + COENames[i], 'left state ' + SphericalAZFPANames[i])
                            X[Xindex + i][1] = stateSphericalAZFPA[i]
                        Xindex += 6
                    elif desiredStateRep == 'MEE':
                        stateMEE = self.COEtoMEE(stateCOE, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + COENames[i], 'left state ' + MEENames[i])
                            X[Xindex + i][1] = stateMEE[i]
                        Xindex += 6
                    elif desiredStateRep == 'IncomingBplane':
                        stateIncomingBplane = self.COEtoIncomingBplane(stateCOE, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + COENames[i], 'left state ' + IncomingBplaneNames[i])
                            X[Xindex + i][1] = stateIncomingBplane[i]
                        Xindex += 6
                    elif desiredStateRep == 'OutgoingBplane':
                        stateOutgoingBplane = self.COEtoOutgoingBplane(stateCOE, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + COENames[i], 'left state ' + OutgoingBplaneNames[i])
                            X[Xindex + i][1] = stateOutgoingBplane[i]
                        Xindex += 6

                elif 'left state VINFin' in description: #it's IncomingBplane
                    stateIncomingBplane = [float(X[Xindex][1]), float(X[Xindex+1][1]), float(X[Xindex+2][1]), float(X[Xindex+3][1]), float(X[Xindex+4][1]), float(X[Xindex+5][1])]
                    if desiredStateRep == 'IncomingBplane':
                        Xindex += 6
                        continue
                    elif desiredStateRep == 'COE':
                        stateCOE = self.IncomingBplanetoCOE(stateIncomingBplane, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + IncomingBplaneNames[i], 'left state ' + COENames[i])
                            X[Xindex + i][1] = stateCOE[i]
                        Xindex += 6
                    elif desiredStateRep == 'SphericalRADEC':
                        stateSphericalRADEC = self.IncomingBplanetoSphericalRADEC(stateIncomingBplane, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + IncomingBplaneNames[i], 'left state ' + SphericalRADECNames[i])
                            X[Xindex + i][1] = stateSphericalRADEC[i]
                        Xindex += 6
                    elif desiredStateRep == 'SphericalAZFPA':
                        stateSphericalAZFPA = self.IncomingBplanetoSphericalAZFPA(stateIncomingBplane, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + IncomingBplaneNames[i], 'left state ' + SphericalAZFPANames[i])
                            X[Xindex + i][1] = stateSphericalAZFPA[i]
                        Xindex += 6
                    elif desiredStateRep == 'MEE':
                        stateMEE = self.IncomingBplanetoMEEtoMEE(stateIncomingBplane, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + IncomingBplaneNames[i], 'left state ' + MEENames[i])
                            X[Xindex + i][1] = stateMEE[i]
                        Xindex += 6
                    elif desiredStateRep == 'Cartesian':
                        stateCartesian = self.IncomingBplanetoCartesian(stateIncomingBplane, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + IncomingBplaneNames[i], 'left state ' + CartesianNames[i])
                            X[Xindex + i][1] = stateCartesian[i]
                        Xindex += 6
                    elif desiredStateRep == 'OutgoingBplane':
                        stateOutgoingBplane = self.IncomingBplanetoOutgoingBplane(stateIncomingBplane, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + IncomingBplaneNames[i], 'left state ' + OutgoingBplaneNames[i])
                            X[Xindex + i][1] = stateOutgoingBplane[i]
                        Xindex += 6
            
                elif 'left state VINFout' in description: #it's OutgoingBplane
                    stateOutgoingBplane = [float(X[Xindex][1]), float(X[Xindex+1][1]), float(X[Xindex+2][1]), float(X[Xindex+3][1]), float(X[Xindex+4][1]), float(X[Xindex+5][1])]
                    if desiredStateRep == 'OutgoingBplane':
                        Xindex += 6
                        continue
                    elif desiredStateRep == 'COE':
                        stateCOE = self.OutgoingBplanetoCOE(stateOutgoingBplane, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + OutgoingBplaneNames[i], 'left state ' + COENames[i])
                            X[Xindex + i][1] = stateCOE[i]
                        Xindex += 6
                    elif desiredStateRep == 'SphericalRADEC':
                        stateSphericalRADEC = self.OutgoingBplanetoSphericalRADEC(stateOutgoingBplane, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + OutgoingBplaneNames[i], 'left state ' + SphericalRADECNames[i])
                            X[Xindex + i][1] = stateSphericalRADEC[i]
                        Xindex += 6
                    elif desiredStateRep == 'SphericalAZFPA':
                        stateSphericalAZFPA = self.OutgoingBplanetoSphericalAZFPA(stateOutgoingBplane, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + OutgoingBplaneNames[i], 'left state ' + SphericalAZFPANames[i])
                            X[Xindex + i][1] = stateSphericalAZFPA[i]
                        Xindex += 6
                    elif desiredStateRep == 'MEE':
                        stateMEE = self.OutgoingBplanetoMEEtoMEE(stateOutgoingBplane, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + OutgoingBplaneNames[i], 'left state ' + MEENames[i])
                            X[Xindex + i][1] = stateMEE[i]
                        Xindex += 6
                    elif desiredStateRep == 'Cartesian':
                        stateCartesian = self.OutgoingBplanetoCartesian(stateOutgoingBplane, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + OutgoingBplaneNames[i], 'left state ' + CartesianNames[i])
                            X[Xindex + i][1] = stateCartesian[i]
                        Xindex += 6
                    elif desiredStateRep == 'IncomingBplane':
                        stateIncomingBplane = self.OutgoingBplanetoIncomingBplane(stateOutgoingBplane, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + OutgoingBplaneNames[i], 'left state ' + IncomingBplaneNames[i])
                            X[Xindex + i][1] = stateIncomingBplane[i]
                        Xindex += 6

                elif 'left state P' in description: #it's MEE
                    stateMEE = [float(X[Xindex][1]), float(X[Xindex+1][1]), float(X[Xindex+2][1]), float(X[Xindex+3][1]), float(X[Xindex+4][1]), float(X[Xindex+5][1])]
                
                    if desiredStateRep == 'MEE':
                        Xindex += 6
                        continue
                    elif desiredStateRep == 'Cartesian':
                        stateCartesian = self.MEEtoCartesian(stateMEE, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + MEENames[i], 'left state ' + CartesianNames[i])
                            X[Xindex + i][1] = stateCartesian[i]
                        Xindex += 6
                    elif desiredStateRep == 'SphericalRADEC':
                        stateSphericalRADEC = self.MEEtoSphericalRADEC(stateMEE, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + MEENames[i], 'left state ' + SphericalRADECNames[i])
                            X[Xindex + i][1] = stateSphericalRADEC[i]
                        Xindex += 6
                    elif desiredStateRep == 'SphericalAZFPA':
                        stateSphericalAZFPA = self.MEEtoSphericalAZFPA(stateMEE, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + MEENames[i], 'left state ' + SphericalAZFPANames[i])
                            X[Xindex + i][1] = stateSphericalAZFPA[i]
                        Xindex += 6
                    elif desiredStateRep == 'COE':
                        stateCOE = self.MEEtoCOE(stateMEE, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + MEENames[i], 'left state ' + COENames[i])
                            X[Xindex + i][1] = stateCOE[i]
                        Xindex += 6
                    elif desiredStateRep == 'IncomingBplane':
                        stateIncomingBplane = self.MEEtoIncomingBplane(stateMEE, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + MEENames[i], 'left state ' + IncomingBplaneNames[i])
                            X[Xindex + i][1] = stateIncomingBplane[i]
                        Xindex += 6
                    elif desiredStateRep == 'OutgoingBplane':
                        stateOutgoingBplane = self.MEEtoOutgoingBplane(stateMEE, mu)
                        for i in range(0, 6):
                            X[Xindex + i][0] = X[Xindex + i][0].replace('left state ' + MEENames[i], 'left state ' + OutgoingBplaneNames[i])
                            X[Xindex + i][1] = stateOutgoingBplane[i]
                        Xindex += 6

                else:#Not a 6-state. Nothing to see here, so just advance
                    Xindex += 1

            else:#Doesn't match the string we're searching for. Nothing to see here, so just advance
                Xindex += 1

        return X