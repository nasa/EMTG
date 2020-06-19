# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:48:29 2020

@author: nhatten
"""

"""
Revision history
2020-06-11; Noble Hatten; created
"""

from mathUtilities import mathUtilities
import autograd.numpy as np # use the autograd version of np

class posVel2BPlaneOut(object):
    def __init__(self):
        self.mathUtil = mathUtilities()
        return

    # the stuff
    def eVector(self, r, v, mu):
        """
        eccentricity vector
        """
        rmag = np.linalg.norm(r)
        vmag = np.linalg.norm(v)
        e = (1.0 / mu) * ((vmag**2 - (mu / rmag)) * r - np.dot(r, v) * v)
        return e

    def nVector(self, r, v, mu):
        """
        not the node vector: h cross e; i think battin calls it p
        """
        h = self.hVector(r, v)
        e = self.eVector(r, v, mu)
        n = np.cross(h, e)
        return n

    def hVector(self, r, v):
        """
        angular momentum vector
        """

        h = np.cross(r, v)
        return h

    def sVector(self, r, v, mu):
        """
        s vector
        """
        e = self.eVector(r, v, mu)
        emag = np.linalg.norm(e)
        eunit = e / emag
        n = self.nVector(r, v, mu)
        nmag = np.linalg.norm(n)
        nunit = n / nmag
        s = (-1.0 / emag) * eunit + (1.0 - (1.0 / emag**2))**0.5 * nunit
        return s

    def tVector(self, r, v, mu):
        """
        t vector as a unit vector
        """
        t = self.tVectorNotUnit(r, v, mu)
        tmag = np.linalg.norm(t)
        t = t / tmag
        return t

    def tVectorNotUnit(self, r, v, mu):
        """
        t vector
        """

        # reference vector is chosen here
        k = self.mathUtil.k_unit()
        referenceVector = k

        h = self.hVector(r, v)
        s = self.sVector(r, v, mu)
        t = np.cross(s, referenceVector)
        return t

    def rVector(self, r, v, mu):
        """
        bplane R vector as a unit vector
        """
        R = self.rVectorNotUnit(r, v, mu)
        rmag = np.linalg.norm(R)
        R = R / rmag
        return R

    def rVectorNotUnit(self, r, v, mu):
        """
        bplane R vector
        """
        s = self.sVector(r, v, mu)
        t = self.tVector(r, v, mu)
        R = np.cross(s, t)
        return R

    def bVector(self, r, v, mu):
        """
        b vector
        """
        bmag = self.bScalar(r, v, mu)
        e = self.eVector(r, v, mu)
        emag = np.linalg.norm(e)
        eunit = e / emag
        n = self.nVector(r, v, mu)
        nmag = np.linalg.norm(n)
        nunit = n / nmag
        b = bmag * ((1.0 - (1.0 / emag**2))**0.5 * eunit + (1.0 / emag) * nunit)
        return b

    def bScalar(self, r, v, mu):
        """
        magnitude of b vector
        also called Delta or miss distance
        """
        h = self.hVector(r, v)
        hmag = np.linalg.norm(h)
        e = self.eVector(r, v, mu)
        emag = np.linalg.norm(e)
        bmag = hmag**2 / (mu * (emag**2 - 1.0)**0.5)
        return bmag

    def bDotR(self, r, v, mu):
        """
        B dot R
        """
        b = self.bVector(r, v, mu)
        R = self. rVector(r, v, mu)
        BR = np.dot(b, R)
        return BR

    def bDotT(self, r, v, mu):
        """
        B dot T
        """
        b = self.bVector(r, v, mu)
        t = self. tVector(r, v, mu)
        BT = np.dot(b, t)
        return BT

    def bTheta(self, r, v, mu):
        """
        b plane theta (clock angle)
        """
        BR = self.bDotR(r, v, mu)
        BT = self.bDotT(r, v, mu)
        theta = np.arctan2(BR, BT)
        return theta

    def rPeri(self, r, v, mu):
        """
        periapsis radius
        """

        e = self.eVector(r, v, mu)
        emag = np.linalg.norm(e)
        vInfMag = self.vInfMag(r, v, mu)
        rp = (mu * (emag - 1.0)) / vInfMag**2
        return rp

    def vInfMag(self, r, v, mu):
        """
        magnitude of v infinity
        """
        vmag = np.linalg.norm(v)
        rmag = np.linalg.norm(r)

        vInf = (vmag**2 - (2.0 * mu) / rmag)**0.5
        return vInf

    def vInfRA(self, r, v, mu):
        """
        right ascension of v infinity
        """
        s = self.sVector(r, v, mu)
        RA = np.arctan2(s[1], s[0])
        return RA

    def vInfDec(self, r, v, mu):
        """
        declination of v infinity
        """
        s = self.sVector(r, v, mu)
        Dec = np.arcsin(s[2])
        return Dec

    # the derivatives of the stuff
    def eVector_derivs(self, r, v, mu):
        """
        checked
        """
        drdx = self.drdx(r)
        dvdx = self.dvdx(v)
        drvdx = np.dot(r, dvdx) + np.dot(v, drdx)
        drvvdx = np.outer(v, drvdx) + np.dot(np.dot(r, v), dvdx)

        rmag = np.linalg.norm(r)
        drmagdr = self.mathUtil.column_vector_norm2_deriv(r)
        dOneByRmagdx = (-1.0 / rmag**2) * np.dot(drmagdr, drdx)

        vmag = np.linalg.norm(v)
        dvmagdv = self.mathUtil.column_vector_norm2_deriv(v)
        dvmag2dx = 2.0 * vmag * np.dot(dvmagdv, dvdx)

        term1 = np.outer(r, (dvmag2dx - mu * dOneByRmagdx)) + (vmag**2 - mu / rmag) * drdx

        term2 = drvvdx
        deVectordxUnscaled = term1 - term2
        # final scaling
        deVectordx = (1.0 / mu) * deVectordxUnscaled
        return deVectordx

    def nVector_derivs(self, r, v, mu):
        """
        checked
        """
        e = self.eVector(r, v, mu)
        h = self.hVector(r, v)
        hCross = self.mathUtil.crossmat(h)
        deVectordx = self.eVector_derivs(r, v, mu)
        dhVectordx = self.hVector_derivs(r, v)
        dnVectordx = np.zeros((3,6))
        dnVectordx[0,:] = dhVectordx[1,:] * e[2] + h[1] * deVectordx[2,:] - dhVectordx[2,:] * e[1] - h[2] * deVectordx[1,:]
        dnVectordx[1,:] = dhVectordx[2,:] * e[0] + h[2] * deVectordx[0,:] - dhVectordx[0,:] * e[2] - h[0] * deVectordx[2,:]
        dnVectordx[2,:] = dhVectordx[0,:] * e[1] + h[0] * deVectordx[1,:] - dhVectordx[1,:] * e[0] - h[1] * deVectordx[0,:]

        # this vector/matrix/tensor formulation isn't working for some reason, so i'm just going with the scalar version
        #dhCrossdx = np.zeros((3,3,6))
        #dhCrossdx[0,1,:] = -dhVectordx[2,:]
        #dhCrossdx[0,2,:] = dhVectordx[1,:]
        #dhCrossdx[1,0,:] = dhVectordx[2,:]
        #dhCrossdx[1,2,:] = -dhVectordx[0,:]
        #dhCrossdx[2,0,:] = -dhVectordx[1,:]
        #dhCrossdx[2,1,:] = dhVectordx[0,:]
        #dnVectordx = self.mathUtil.tensor_bullet2_vector(dhCrossdx, e) + np.dot(hCross, deVectordx)
        return dnVectordx

    def hCross(self, r, v):
        h = self.hVector(r, v)
        hCrossMat = self.mathUtil.crossmat(h)
        return hCrossMat

    def hVector_derivs(self, r, v):
        """
        checked
        typed
        """
        dvdx = self.dvdx(v)
        rCross = self.mathUtil.crossmat(r)
        drCrossdx = np.zeros((3,3,6))
        drCrossdx[2,1,0] = 1.0
        drCrossdx[1,2,0] = -1.0
        drCrossdx[0,2,1] = 1.0
        drCrossdx[2,0,1] = -1.0
        drCrossdx[0,1,2] = -1.0
        drCrossdx[1,0,2] = 1.0
        dhVectordx = self.mathUtil.tensor_bullet2_vector(drCrossdx, v) + np.dot(rCross, dvdx)
        return dhVectordx
    
    def hVector_derivs_v2(self, r, v):
        """
        """
        rCross = self.mathUtil.crossmat(r)
        vCross = self.mathUtil.crossmat(v)
        dhVectordx = np.zeros((3,6))
        dhVectordr = -vCross
        dhVectordv = rCross
        dhVectordx[:,0:3] = dhVectordr
        dhVectordx[:,3:6] = dhVectordv
        return dhVectordx
    
    def trueAnomaly_derivs(self, r, v, mu):
        e = self.eVector(r, v, mu)
        argx = np.dot(e, r)
        eCrossR = self.eCrossR(r, v, mu)
        argy = np.linalg.norm(eCrossR)
        denom = argx**2 + argy**2
        dTadArgx = -argy / denom
        dTadArgy = argx / denom
        drdx = self.drdx(r)
        dedx = self.eVector_derivs(r, v, mu)
        dArgxdx = np.dot(e, drdx) + np.dot(r, dedx) 
        dECrossRdx = self.eCrossR_derivs(r, v, mu)
        dArgydx = (1./argy) * np.dot(eCrossR, dECrossRdx)
        d_eCrossR_dx = self.eCrossR_derivs(r, v, mu)
        dTadx = np.dot(dTadArgx, dArgxdx) + np.dot(dTadArgy, dArgydx)
        if (np.dot(r, v) < 0):
            dTadx = -dTadx
        return dTadx
    
    def eCrossR(self, r, v, mu):
        e = self.eVector(r, v, mu)
        eCross = self.mathUtil.crossmat(e)
        result = np.dot(eCross, r)
        return result
    
    def eCrossR_derivs(self, r, v, mu):
        """
        """
        e = self.eVector(r, v, mu)
        deCrossRde = -self.mathUtil.crossmat(r)
        deCrossRdr = self.mathUtil.crossmat(e)
        dedx = self.eVector_derivs(r, v, mu)
        drdx = self.drdx(r)
        deCrossRdx = np.dot(deCrossRde, dedx) + np.dot(deCrossRdr, drdx)
        return deCrossRdx

    def sVector_derivs(self, r, v, mu):
        """
        checked
        typed
        """
        e = self.eVector(r, v, mu)
        emag = np.linalg.norm(e)
        eUnit = e / emag
        n = self.nVector(r, v, mu)
        nmag = np.linalg.norm(n)
        nUnit = n / nmag
        rootTerm = np.sqrt(1.0 - (1.0 / emag**2))

        deVectordx = self.eVector_derivs(r, v, mu)
        deUnitdx = np.dot(self.mathUtil.unit_vector_deriv(e), deVectordx)
        dOneByemagdx = (-1.0 / emag**2) * np.dot(eUnit, deVectordx)

        dnVectordx = self.nVector_derivs(r, v, mu)
        dnUnitdx = np.dot(self.mathUtil.unit_vector_deriv(n), dnVectordx)

        drootTermdx = 0.5 * (1.0 - emag**(-2))**(-0.5) * 2.0 * emag**(-3) * np.dot(eUnit, deVectordx)

        term1 = np.outer(eUnit, dOneByemagdx)
        term2 = (1.0 / emag) * deUnitdx
        term3 = np.outer(nUnit, drootTermdx)
        term4 = rootTerm * dnUnitdx
        dsVectordx = -term1 - term2 + term3 + term4
        return dsVectordx

    def tVector_derivs(self, r, v, mu):
        """
        checked
        typed
        """

        # reference vector is chosen here
        k = self.mathUtil.k_unit()
        dkdx = np.zeros((3,6))
        referenceVector = k
        dreferenceVectordx = dkdx

        h = self.hVector(r, v)
        dhdx = self.hVector_derivs(r, v)
        s = self.sVector(r, v, mu)
        dsVectordx = self.sVector_derivs(r, v, mu)

        dtVectordx = np.zeros((3,6))
        dtVectordx[0,:] = dsVectordx[1,:] * referenceVector[2] + s[1] * dreferenceVectordx[2,:] - dsVectordx[2,:] * referenceVector[1] - s[2] * dreferenceVectordx[1,:]
        dtVectordx[1,:] = dsVectordx[2,:] * referenceVector[0] + s[2] * dreferenceVectordx[0,:] - dsVectordx[0,:] * referenceVector[2] - s[0] * dreferenceVectordx[2,:]
        dtVectordx[2,:] = dsVectordx[0,:] * referenceVector[1] + s[0] * dreferenceVectordx[1,:] - dsVectordx[1,:] * referenceVector[0] - s[1] * dreferenceVectordx[0,:]

        # simplify because k is simple
        #dtVectordx[0,:] = dsVectordx[1,:]
        #dtVectordx[1,:] = -dsVectordx[0,:]
        
        t = self.tVectorNotUnit(r, v, mu)
        dtVectordx = np.dot(self.mathUtil.unit_vector_deriv(t), dtVectordx)
        return dtVectordx

    def rVector_derivs(self, r, v, mu):
        """
        checked
        typed
        """
        s = self.sVector(r, v, mu)
        dsVectordx = self.sVector_derivs(r, v, mu)
        t = self.tVector(r, v, mu)
        dtVectordx = self.tVector_derivs(r, v, mu)

        drVectordx = np.zeros((3,6))
        drVectordx[0,:] = dsVectordx[1,:] * t[2] + s[1] * dtVectordx[2,:] - dsVectordx[2,:] * t[1] - s[2] * dtVectordx[1,:]
        drVectordx[1,:] = dsVectordx[2,:] * t[0] + s[2] * dtVectordx[0,:] - dsVectordx[0,:] * t[2] - s[0] * dtVectordx[2,:]
        drVectordx[2,:] = dsVectordx[0,:] * t[1] + s[0] * dtVectordx[1,:] - dsVectordx[1,:] * t[0] - s[1] * dtVectordx[0,:]

        R = self.rVectorNotUnit(r, v, mu)
        drVectordx = np.dot(self.mathUtil.unit_vector_deriv(R), drVectordx)

        return drVectordx

    def bVector_derivs(self, r, v, mu):
        """
        checked
        typed
        """
        B = self.bVector(r, v, mu)
        bScalar = self.bScalar(r, v, mu)
        dbScalardx = self.bScalar_derivs(r, v, mu)

        e = self.eVector(r, v, mu)
        emag = np.linalg.norm(e)
        eUnit = e / emag
        n = self.nVector(r, v, mu)
        nmag = np.linalg.norm(n)
        nUnit = n / nmag
        rootTerm = np.sqrt(1.0 - (1.0 / emag**2))

        deVectordx = self.eVector_derivs(r, v, mu)
        deUnitdx = np.dot(self.mathUtil.unit_vector_deriv(e), deVectordx)
        dOneByemagdx = (-1.0 / emag**2) * np.dot(eUnit, deVectordx)

        dnVectordx = self.nVector_derivs(r, v, mu)
        dnUnitdx = np.dot(self.mathUtil.unit_vector_deriv(n), dnVectordx)

        drootTermdx = 0.5 * (1.0 - emag**(-2))**(-0.5) * 2.0 * emag**(-3) * np.dot(eUnit, deVectordx)

        # the unit vector part
        term1 = rootTerm * deUnitdx
        term2 = np.outer(eUnit, drootTermdx)
        term3 = -(1.0 / emag) * dnUnitdx
        term4 = -np.outer(nUnit, dOneByemagdx)
        dbUnitVectordx = term1 + term2 - term3 - term4

        # adding the magnitude
        dbVectordx = np.outer(B / bScalar, dbScalardx) + bScalar * dbUnitVectordx

        return dbVectordx

    def bScalar_derivs(self, r, v, mu):
        """
        checked
        typed
        """
        h = self.hVector(r, v)
        hmag = np.linalg.norm(h)
        dhVectordx = self.hVector_derivs(r, v)
        e = self.eVector(r, v, mu)
        emag = np.linalg.norm(e)
        deVectordx = self.eVector_derivs(r, v, mu)
        term1 = -(emag**2 - 1.0)**(-1.5) * np.dot(e, deVectordx)
        term2 = 2.0 * np.dot(h, dhVectordx)
        dbScalardx = (1.0 / mu) * (hmag**2 * term1 + (emag**2 - 1.0)**(-0.5) * term2)

        return dbScalardx

    def bDotR_derivs(self, r, v, mu):
        """
        checked
        typed
        """
        B = self.bVector(r, v, mu)
        dbVectordx = self.bVector_derivs(r, v, mu)
        R = self. rVector(r, v, mu)
        drVectordx = self.rVector_derivs(r, v, mu)

        dbDotRdx = np.dot(R, dbVectordx) + np.dot(B, drVectordx)
        return dbDotRdx

    def bDotT_derivs(self, r, v, mu):
        """
        checked
        typed
        """
        B = self.bVector(r, v, mu)
        dbVectordx = self.bVector_derivs(r, v, mu)
        t = self. tVector(r, v, mu)
        dtVectordx = self.tVector_derivs(r, v, mu)

        dbDotTdx = np.dot(t, dbVectordx) + np.dot(B, dtVectordx)
        return dbDotTdx

    def bTheta_derivs(self, r, v, mu):
        """
        checked
        typed
        """
        BR = self.bDotR(r, v, mu)
        BT = self.bDotT(r, v, mu)
        dbDotRdx = self.bDotR_derivs(r, v, mu)
        dbDotTdx = self.bDotT_derivs(r, v, mu)
        datan2 = self.mathUtil.d_atan2(BR, BT)
        dbThetadx = datan2[0] * dbDotRdx + datan2[1] * dbDotTdx
        return dbThetadx

    def rPeri_derivs(self, r, v, mu):
        """
        checked, but some are off by ~1e-13
        typed
        """
        e = self.eVector(r, v, mu)
        emag = np.linalg.norm(e)
        deVectordx = self.eVector_derivs(r, v, mu)

        vInfMag = self.vInfMag(r, v, mu)
        dvInfMagdx = self.vInfMag_derivs(r, v, mu)

        drPeridx = mu * ((1.0 / (emag * vInfMag**2)) * np.dot(e, deVectordx) - (2.0 / vInfMag**3) * (emag - 1.0) *dvInfMagdx)
        return drPeridx

    def vInfMag_derivs(self, r, v, mu):
        """
        checked
        typed
        """
        drdx = self.drdx(r)
        dvdx = self.dvdx(v)
        rmag = np.linalg.norm(r)
        vmag = np.linalg.norm(v)

        dvInfMagdx = 0.5 * (vmag**2 - 2.0 * mu / rmag)**(-0.5) * (2.0 * np.dot(v, dvdx) + (2.0 * mu / rmag**3) * np.dot(r, drdx))
        return dvInfMagdx

    def vInfRA_derivs(self, r, v, mu):
        """
        checked
        typed
        """
        s = self.sVector(r, v, mu)
        temp = self.mathUtil.d_atan2(s[1], s[0])
        dRAds = np.array([temp[1], temp[0], 0.0])
        dsVectordx = self.sVector_derivs(r, v, mu)
        dRAdx = np.dot(dRAds, dsVectordx)
        return dRAdx

    def vInfDec_derivs(self, r, v, mu):
        """
        checked
        typed
        """
        s = self.sVector(r, v, mu)
        dDecds = np.array([0.0, 0.0, (1.0 - s[2]**2)**(-0.5)])
        dsVectordx = self.sVector_derivs(r, v, mu)
        dDecdx = np.dot(dDecds, dsVectordx)
        return dDecdx

    def drdx(self, r):
        """
        typed
        """
        drdx = np.zeros((3,6))
        drdx[0,0] = 1.0
        drdx[1,1] = 1.0
        drdx[2,2] = 1.0
        return drdx

    def dvdx(self, v):
        """
        typed
        """
        dvdx = np.zeros((3,6))
        dvdx[0,3] = 1.0
        dvdx[1,4] = 1.0
        dvdx[2,5] = 1.0
        return dvdx

    def trueAnomaly(self, r, v, mu):
        e = self.eVector(r, v, mu)
        TA = self.mathUtil.angle_between_2_vectors(e, r)
        if (np.dot(r, v) < 0):
            TA = 2. * np.pi - TA
        return TA

    def periapsisPositionVector(self, r, v, mu):
        # in direction of eccentricity vector
        e = self.eVector(r, v, mu)
        eMag = np.linalg.norm(self.eVector(r, v, mu))
        rpDirection = e / eMag

        # magnitude
        h = self.hVector(r, v)
        hMag = np.linalg.norm(h)
        p = hMag**2 / mu
        rpMag = p / (1.0 + eMag)
        r = rpMag * rpDirection
        return r
    def periapsisVelocityVector(self, r, v, mu):
        e = self.eVector(r, v, mu)
        eMag = np.linalg.norm(self.eVector(r, v, mu))
        h = self.hVector(r, v)
        hMag = np.linalg.norm(h)
        p = hMag**2 / mu
        vpDir = np.dot(self.mathUtil.crossmat(h), e)
        vpDir = vpDir / np.linalg.norm(vpDir)
        vpMag = np.sqrt(mu / p) * (eMag + 1.0)
        v = vpMag * vpDir
        return v