import spiceypy
try:
    import spiceypy as spice
    import spiceypy.utils.support_types as stypes
except:
    print("spiceypy not available")
    
import numpy as np

class MissionEventFunc(object):
    import spiceypy as spice
    import numpy as np
    def evaluate(self, epoch):
        print("please implement an evaluate method for your MissionEventFunc")
        
class speFunc(MissionEventFunc):
    """
    Used to calculate SPE angle
    """
    def __init__(self, spiceId, units = "deg"):
        self.spiceId = spiceId
        self.units = units
        return
    
    def evaluate(self, epoch):
        # get position of object w.r.t. Earth
        stateWrtEarth, lt = spice.spkez(self.spiceId, epoch, "J2000", 'NONE', 399)

        # get position of Sun w.r.t. Earth
        sunStateWrtEarth, lt = spice.spkez(10, epoch, "J2000", 'NONE', 399)

        earthStateWrtSpacecraft = [ -x for x in stateWrtEarth]
        sunStateWrtSpacecraft = sunStateWrtEarth - stateWrtEarth

        cosAngle_SPE = np.dot(sunStateWrtSpacecraft[0:3], earthStateWrtSpacecraft[0:3])/np.linalg.norm(sunStateWrtSpacecraft[0:3])/np.linalg.norm(earthStateWrtSpacecraft[0:3])
        
        spe = np.arccos(cosAngle_SPE)
        
        if self.units == "deg":
            spe *= (180. / np.pi)
        return spe
        
class sepFunc(MissionEventFunc):
    """
    Used to calculate SEP angle
    """
    def __init__(self, spiceId, units = "deg"):
        self.spiceId = spiceId
        self.units = units
        return
    
    def evaluate(self, epoch):
        # get position of object w.r.t. Earth
        stateWrtEarth, lt = spice.spkez(self.spiceId, epoch, "J2000", 'NONE', 399)

        # get position of Sun w.r.t. Earth
        sunStateWrtEarth, lt = spice.spkez(10, epoch, "J2000", 'NONE', 399)

        earthStateWrtSpacecraft = [ -x for x in stateWrtEarth]
        sunStateWrtSpacecraft = sunStateWrtEarth - stateWrtEarth

        cosAngle_SEP = np.dot(sunStateWrtEarth[0:3], stateWrtEarth[0:3])/np.linalg.norm(sunStateWrtEarth[0:3])/np.linalg.norm(stateWrtEarth[0:3])
        
        sep = np.arccos(cosAngle_SEP)
        
        if self.units == "deg":
            sep *= (180. / np.pi)
        return sep
        
class rmagFunc(MissionEventFunc):
    """
    Used to calculate how far a body is from another body
    """
    def __init__(self, spiceId, body, units = "km"):
        self.spiceId = spiceId
        self.body = body
        self.units = units
        self.unitsDict = {
                            "km": 1.0,
                            "AU": 149597870.691,
                            "Rj": 71492.0
                        }
        return
        
    def addUnit(self, unitName, unitValueInKm):
        self.unitsDict[unitName] = unitValueInKm
        return
    
    def evaluate(self, epoch):
        state, lt = self.spice.spkez(self.spiceId, epoch, 'J2000', 'NONE', self.spice.bodn2c(self.body.name))
        radius, lon, lat = spice.reclat(state[0:3])
        radius /= self.unitsDict[self.units]
        return radius
        
class radialVelocityFunc(MissionEventFunc):
    """
    Calculate radial velocity relative to some body in some frame
    """
    
    def __init__(self, spiceId, body, frame):
        self.spiceId = spiceId
        self.body = body
        self.frame = frame
        return
        
    def evaluate(self, epoch):
        import numpy as np
        state, lt = self.spice.spkez(self.spiceId, epoch, self.frame, 'NONE', self.spice.bodn2c(self.body.name))
        r = state[0:3]
        v = state[3:6]
        rMag = np.linalg.norm(r)
        rUnit = r / rMag
        vRadial = np.dot(v, rUnit)
        return vRadial
    
class latitudeFunc(MissionEventFunc):
    """
    Used to calculate CENTRIC latitude
    """
    def __init__(self, spiceId, centralBody, frame, angleUnits = 'deg'):
        self.spiceId = spiceId
        self.centralBody = centralBody
        self.frame = frame
        self.angleUnits = angleUnits
        return
    
    def evaluate(self, epoch):
        state, lt = self.spice.spkez(self.spiceId, epoch, self.frame, 'NONE', self.spice.bodn2c(self.centralBody.name))
        radius, lon, lat = spice.reclat(state[0:3])
        if self.angleUnits == 'deg':
            lat *= 180. / np.pi
        return lat
        
class deticLatitudeFunc(MissionEventFunc):
    """
    Used to calculate DETIC latitude
    """
    def __init__(self, spiceId, centralBody, frame, angleUnits = 'deg'):
        self.spiceId = spiceId
        self.centralBody = centralBody
        self.frame = frame
        self.angleUnits = angleUnits
        return
    
    def evaluate(self, epoch):
        state, lt = self.spice.spkez(self.spiceId, epoch, self.frame, 'NONE', self.spice.bodn2c(self.centralBody.name))
        lon, lat, alt = spice.recgeo(state[0:3], self.centralBody.equatorialRadius, self.centralBody.flatteningCoef)
        if self.angleUnits == 'deg':
            lat *= 180. / np.pi
        return lat
        
class longitudeFunc(MissionEventFunc):
    """
    Used to longitude
    """
    def __init__(self, spiceId, centralBody, frame, angleUnits = 'deg'):
        self.spiceId = spiceId
        self.centralBody = centralBody
        self.frame = frame
        self.angleUnits = angleUnits
        return
    
    def evaluate(self, epoch):
        state, lt = self.spice.spkez(self.spiceId, epoch, self.frame, 'NONE', self.spice.bodn2c(self.centralBody.name))
        radius, lon, lat = spice.reclat(state[0:3])
        if self.angleUnits == 'deg':
            lon *= 180. / np.pi
        return lon
        
class dLatitudeDTFunc(MissionEventFunc):
    """
    Used to calculate d/dt [centric latitude]. With this, we can use an event to locate when centric latitude itself hits an extremum
    """
    def __init__(self, spiceId, centralBody, frame):
        self.spiceId = spiceId
        self.centralBody = centralBody
        self.frame = frame
        self.rmagFunction = rmagFunc(self.spiceId, self.centralBody)
        return
        
    def evaluate(self, epoch):
        #print('hello')
        state, lt = self.spice.spkez(self.spiceId, epoch, self.frame, 'NONE', self.spice.bodn2c(self.centralBody.name))
        
        # get state in j2000, then convert to IAU_SUN but DO NOT take into account rotation of frame.
        # so, just use pxform on both position and velocity
        mat = self.spice.pxform(self.frame, "IAU_SUN", epoch)
        
        r = np.array(state[0:3])
        v = np.array(state[3:6])
        r = np.dot(mat, r)
        v = np.dot(mat, v)
        rmag = self.rmagFunction.evaluate(epoch)
        vmag = np.linalg.norm(v)
        
        c1 = 1. / np.sqrt(1. - (r[2] / rmag)**2)
        c2 = (1. / rmag) * v[2] - (r[2] / rmag**3) * np.dot(r, v)
        val = c2
        #print((1. / rmag) * v[2], (r[2] / rmag**3) * np.dot(r, v))
        #print(val)
        return val
    
class altitudeFunc(MissionEventFunc):
    """
    Calculates altitude assuming a spherical central body
    """
    def __init__(self, spiceId, centralBody):
        self.centralBody = centralBody
        self.spiceId = spiceId

    def evaluate(self, epoch):
        state, lt = self.spice.spkez(self.spiceId, epoch, "J2000", 'NONE', self.spice.bodn2c(self.centralBody.name))
        r = self.np.array([state[0], state[1], state[2]])
        return self.np.linalg.norm(r) - self.centralBody.radius
        
class rmagEventFunc(MissionEventFunc):
    """
    Used to determine when one body passes through a certain distance from another body
    """
    def __init__(self, spiceId, rmag, body):
        self.spiceId = spiceId
        self.targetRmag = rmag
        self.body = body
        self.evalRmag = rmagFunc(spiceId, body)
        return
    
    def evaluate(self, epoch):
        f = self.evalRmag.evaluate(epoch) - self.targetRmag
        return f
        
class latitudeEventFunc(MissionEventFunc):
    def __init__(self, spiceId, latitude, centralBody, frame, angleUnits = 'deg'):
        self.spiceId = spiceId
        self.targetLatitude = latitude
        self.centralBody = centralBody
        self.frame = frame
        self.angleUnits = angleUnits
        self.evalLatitude = latitudeFunc(spiceId, centralBody, frame, angleUnits)
        return
    
    def evaluate(self, epoch):
        f = self.evalLatitude.evaluate(epoch) - self.targetLatitude
        return f
        
class latitudeExtremumEventFunc(MissionEventFunc):
    def __init__(self, spiceId, centralBody, frame):
        self.spiceId = spiceId
        self.centralBody = centralBody
        self.frame = frame
        self.evalDLatitudeDTFunc = dLatitudeDTFunc(self.spiceId, self.centralBody, self.frame)
        return
        
    def evaluate(self, epoch):
        f = self.evalDLatitudeDTFunc.evaluate(epoch) # extremum is when this evaluation hits zero, so no subtraction
        return f
        

class apseEventFunc(MissionEventFunc):
    def __init__(self, SPICE_ID, central_body):
        self.central_body = central_body
        self.SPICE_ID = SPICE_ID

    def evaluate(self, epoch):
        spacecraft_state, light_time = self.spice.spkez(self.SPICE_ID, epoch, "J2000", 'NONE', self.spice.bodn2c(self.central_body.name))
        r = self.np.array([spacecraft_state[0], spacecraft_state[1], spacecraft_state[2]])
        v = self.np.array([spacecraft_state[3], spacecraft_state[4], spacecraft_state[5]])
        return self.np.dot(r, v)

class altitudeEventFunc(MissionEventFunc):
    def __init__(self, SPICE_ID, altitude, central_body):
        self.central_body = central_body
        self.SPICE_ID = SPICE_ID
        self.target_altitude = altitude
        self.evalAltitude = altitudeFunc(SPICE_ID, central_body)

    def evaluate(self, epoch):
        #spacecraft_state, light_time = self.spice.spkez(self.SPICE_ID, epoch, "J2000", 'NONE', self.spice.bodn2c(self.central_body.name))
        #r = self.np.array([spacecraft_state[0], spacecraft_state[1], spacecraft_state[2]])
        altitude = self.evalAltitude.evaluate(epoch)
        return altitude - self.target_altitude

class elevationFromGroundStationEventFunc(MissionEventFunc):
    def __init__(self, central_body, target_SPICE_ID, origin_SPICE_ID=None, latitude=None, longitude=None, altitude=None, max_or_min='minimize'):
        self.central_body = central_body
        self.target_SPICE_ID = target_SPICE_ID
        self.origin_SPICE_ID = origin_SPICE_ID
        self.altitude = altitude
        if latitude != None:
            self.latitude = latitude * self.np.pi / 180.0
        if longitude != None:
            self.longitude = longitude * self.np.pi / 180.0
        self.mult = 1.0
        if 'max' in max_or_min.lower():
            self.mult = -1.0
     
    def evaluate(self, epoch):
        
        # compute the origin Cartesian state
        gs_position_bcf = self.np.empty(3)
        gs_position_ICRF = self.np.empty(3)
        # if we provided an SPK for the origin
        if self.origin_SPICE_ID != None:
            gs_state_ICRF, light_time = spice.spkez(self.origin_SPICE_ID, epoch, "J2000", 'NONE', self.central_body.SPICE_ID)
            rot_M = self.spice.tipbod("J2000", self.central_body.SPICE_ID, epoch)
            #rot_M = self.np.eye(3)
            gs_position_ICRF = gs_state_ICRF[0:3]
            gs_position_bcf = self.np.dot(rot_M, gs_position_ICRF)

            #gs_position_bcf_norm = self.np.linalg.norm(gs_position_bcf)
            #gs_position_bcf[0] /= gs_position_bcf_norm
            #gs_position_bcf[1] /= gs_position_bcf_norm
            #gs_position_bcf[2] /= gs_position_bcf_norm
            #print(gs_position_bcf)

        else: # we need to compute the body fixed state from the provided lat/lon/alt
            sLat = self.np.sin(self.latitude);
            cLat = self.np.cos(self.latitude);
            sLon = self.np.sin(self.longitude);
            Rp = self.central_body.radius * (1.0 - self.central_body.flattening);
            ex2 = (self.central_body.radius**2.0 - Rp**2.0) / self.central_body.radius**2;
            ee2 = 0.0; # true for oblate spheroid, not true for triaxial ellipsoid
            v = self.central_body.radius / self.np.sqrt(1.0 - ex2 * sLat**2.0 - ee2 * cLat**2.0 * sLon**2.0);

            gs_position_bcf[0] = (v + self.altitude) * cLat * self.np.cos(self.longitude);
            gs_position_bcf[1] = (v * (1.0 - ee2) + self.altitude) * cLat * sLon;
            gs_position_bcf[2] = (v * (1.0 - ex2) + self.altitude) * sLat;
            rot_M = self.spice.tipbod("J2000", self.central_body.SPICE_ID, epoch)
            gs_position_ICRF = self.np.dot(self.np.transpose(rot_M), gs_position_bcf)


        # compute the surface normal at the origin
        f2 = (1.0 - self.central_body.flattening) * (1.0 - self.central_body.flattening);
        rx = gs_position_bcf[0];
        ry = gs_position_bcf[1];
        rz = gs_position_bcf[2];
        spheroid_major_axis = self.np.sqrt(rx*rx + ry*ry + rz*rz / f2);
        spheroid_major_axis2 = spheroid_major_axis * spheroid_major_axis;

        nx = 2.0 * rx / (spheroid_major_axis2);
        ny = 2.0 * ry / (spheroid_major_axis2);
        nz = 2.0 * rz / (spheroid_major_axis2 * f2);
        nnorm = self.np.sqrt(nx * nx + ny * ny + nz * nz);
        surface_normal_vec = self.np.empty(3)
        surface_normal_vec[0] = nx
        surface_normal_vec[1] = ny
        surface_normal_vec[2] = nz
        surface_normal_vec_norm = self.np.linalg.norm(surface_normal_vec)
        surface_normal_vec[0] /= surface_normal_vec_norm
        surface_normal_vec[1] /= surface_normal_vec_norm
        surface_normal_vec[2] /= surface_normal_vec_norm
        #print(surface_normal_vec)

        # get the target state in ICRF
        spacecraft_state_ICRF, light_time = self.spice.spkez(self.target_SPICE_ID, epoch, "J2000", 'NONE', self.central_body.SPICE_ID)
        spacecraft_position_ICRF = spacecraft_state_ICRF[0:3]

        # get the relative state vector in ICRF
        relative_position_ICRF = self.np.empty(3)
        relative_position_ICRF[0] = spacecraft_position_ICRF[0] - gs_position_ICRF[0]
        relative_position_ICRF[1] = spacecraft_position_ICRF[1] - gs_position_ICRF[1]
        relative_position_ICRF[2] = spacecraft_position_ICRF[2] - gs_position_ICRF[2]

        rot_M = self.spice.tipbod("J2000", self.central_body.SPICE_ID, epoch)
        #rot_M = self.np.eye(3)
        rel_pos_body_fixed = self.np.dot(rot_M, relative_position_ICRF)
        rel_pos_norm = self.np.linalg.norm(rel_pos_body_fixed)
        rel_pos_body_fixed[0] /= rel_pos_norm
        rel_pos_body_fixed[1] /= rel_pos_norm
        rel_pos_body_fixed[2] /= rel_pos_norm

        #print(rel_pos_body_fixed)

        #print(spice.timout(epoch, "DD-MON-YYYY HR:MN:SC.### UTC ::UTC"), self.np.dot(rot_M, spacecraft_position_ICRF), gs_position_bcf)

        # compute the elevation angle
        #cos_zenith_angle = self.np.dot(surface_normal_vec, rel_pos_body_fixed) / (surface_normal_vec_norm * rel_pos_norm)
        cos_zenith_angle = self.np.dot(surface_normal_vec, rel_pos_body_fixed)
        sin_detic_elevation = cos_zenith_angle
        #return self.mult * self.np.arccos(cos_zenith_angle)
        #print(self.np.arcsin(sin_detic_elevation) * 180.0 / self.np.pi)
        return self.mult * self.np.arcsin(sin_detic_elevation)
