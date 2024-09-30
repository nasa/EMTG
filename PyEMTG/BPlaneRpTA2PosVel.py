
from mathUtilities import mathUtilities
import autograd.numpy as np # use the autograd version of np

class BPlaneRpTA2PosVel(object):
    # x = [vInfMag, vInfRA, vInfDec, rp, bTheta, TA]
    def __init__(self):
        self.mathUtil = mathUtilities()
        return

    # the stuff
    def eMag(self, x, mu):
        """
        magnitude of eccentricity
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        rp = x[3]
        bTheta = x[4]
        TA = x[5]
        
        eMag = (rp * vInfMag**2) / mu + 1.
        return eMag

    def eMag_derivs(self, x, mu):
        """
        derivatives of magnitude of eccentricity
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        rp = x[3]
        bTheta = x[4]
        TA = x[5]

        term = vInfMag / mu
        deMagdx = term * np.array([2. * rp, 
                                   0.0, 
                                   0.0,
                                   vInfMag,
                                   0.0,
                                   0.0])
        return deMagdx

    def bScalar_derivs(self, x, mu):
        """
        derivatives of b scalar parameter
        when rp is in the state rep instead of b,
        b depends on rp and vinf
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        rp = x[3]
        bTheta = x[4]
        TA = x[5]
        vp = self.vMagPeri(x, mu)
        
        dvPeridx = self.vMagPeri_derivs(x, mu)
        dbScalardx = np.zeros(6)
        dbScalardx[0] = (rp / vInfMag) * (dvPeridx[0] - vp / vInfMag)
        dbScalardx[3] = (vp / vInfMag) + (rp / vInfMag) * dvPeridx[3]
        return dbScalardx

    def bTheta_derivs(self, x):
        """
        derivatives of b theta parameter
        typed
        """
        dbThetadx = np.zeros(6)
        dbThetadx[4] = 1.0
        return dbThetadx

    def TA_derivs(self, x):
        """
        derivatives of true anomaly
        typed
        """
        dTAdx = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0])
        return dTAdx

    def rpMag(self, x, mu):
        """
        periapsis magnitude
        """
        rpMag = x[3]
        return rpMag

    def vInfVector(self, x):
        """
        v infinity vector
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        rp = x[3]
        bTheta = x[4]
        TA = x[5]

        vInfVector = vInfMag * np.array([np.cos(vInfDec) * np.cos(vInfRA), np.cos(vInfDec) * np.sin(vInfRA), np.sin(vInfDec)])
        return vInfVector

    def sVector(self, x):
        """
        s vector
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        rp = x[3]
        bTheta = x[4]
        TA = x[5]

        s = self.vInfVector(x) / vInfMag
        return s

    def sVector_derivs(self, x):
        """
        derivatives of s vector
        typed
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        rp = x[3]
        bTheta = x[4]
        TA = x[5]

        ca = np.cos(vInfRA)
        sa = np.sin(vInfRA)
        cd = np.cos(vInfDec)
        sd = np.sin(vInfDec)

        dsdx = np.zeros((3,6))
        dsdx[0,1] = -cd * sa
        dsdx[0,2] = -sd * ca
        dsdx[1,1] = cd * ca
        dsdx[1,2] = -sd * sa
        dsdx[2,2] = cd

        return dsdx

    def hMag(self, x, mu):
        """
        angular momentum magnitude
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        rp = x[3]
        bTheta = x[4]
        TA = x[5]
        
        vMagPeri = self.vMagPeri(x, mu)

        hMag = rp * vMagPeri
        return hMag
    
    def vMagPeri(self, x, mu):
        """
        velocity magnitude at periapse
        """
        vInfMag = x[0]
        rp = x[3]
        vMagPeri = (vInfMag**2 + 2. * mu / rp)**0.5
        return vMagPeri
    
    def vMagPeri_derivs(self, x, mu):
        """
        derivatives of velocity magnitude at periapse w/r/t/ b plane state
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        rp = x[3]
        bTheta = x[4]
        TA = x[5]
        
        vp = self.vMagPeri(x, mu)
        dvMagPeridx = np.zeros((6))
        dvMagPeridx[0] = vInfMag / vp
        dvMagPeridx[3] = -mu / (vp * rp**2)
        return dvMagPeridx

    def hMag_derivs(self, x, mu):
        """
        derivatives of angular momentum magnitude
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        rp = x[3]
        bTheta = x[4]
        TA = x[5]
        
        vMagPeri = self.vMagPeri(x, mu)
        
        dhMagdrp = vMagPeri - mu / (rp * vMagPeri)
        dhMagdvInfMag = (rp * vInfMag) / vMagPeri
        
        dhMagdx = np.array([dhMagdvInfMag, 0.0, 0.0, dhMagdrp, 0.0, 0.0])

        return dhMagdx
    
    def bMag(self, x, mu):
        """
        magnitude of b vector
        """
        vInfMag = x[0]
        rp = x[3]
        vMagPeri = self.vMagPeri(x, mu)
        bMag = (rp * vMagPeri) / vInfMag
        return bMag

    def tVector(self, x):
        """
        T vector as a unit vector
        """
        t = self.tVectorNotUnit(x)
        t = t / np.linalg.norm(t)
        return t

    def tVector_derivs(self, x):
        """
        derivatives of the T unit vector
        typed
        """
        # choose the reference vector here
        k = self.mathUtil.k_unit()
        referenceVector = k
        dreferenceVectordx = np.zeros((3,6)) # for k reference vector (or any constant reference vector)

        s = self.sVector(x)
        crossProduct = np.cross(s, referenceVector)
        crossProductMag = np.linalg.norm(crossProduct)
        dcrossProductds, dcrossProductdreferenceVector = self.mathUtil.d_crossproduct(s, referenceVector)
        
        dsdx = self.sVector_derivs(x)
        dtdxNotUnit = np.dot(dcrossProductds, dsdx) + np.dot(dcrossProductdreferenceVector, dreferenceVectordx)
        dOneByCrossProductMagdx = (-1.0 / crossProductMag**3) * np.dot(crossProduct, dtdxNotUnit)
        dtdx = (1.0 / crossProductMag) * dtdxNotUnit + np.outer(crossProduct, dOneByCrossProductMagdx)

        return dtdx

    def tVectorNotUnit(self, x):
        """
        T vector, not necessarily as a unit vector
        """
        s = self.sVector(x)
        k = self.mathUtil.k_unit()
        #h = self.hVector(x)
        referenceVector = k
        t = np.cross(s, referenceVector)
        return t

    def rVector(self, x):
        """
        bplane R vector as a unit vector
        """
        R = self.rVectorNotUnit(x)
        rmag = np.linalg.norm(R)
        R = R / rmag
        return R

    def rVector_derivs(self, x):
        """
        bplane R unit vector derivatives
        typed
        """

        s = self.sVector(x)
        t = self.tVector(x)
        crossProduct = np.cross(s, t)
        crossProductMag = np.linalg.norm(crossProduct)
        dcrossProductds, dcrossProductdt = self.mathUtil.d_crossproduct(s, t)
        dsdx = self.sVector_derivs(x)
        dtdx = self.tVector_derivs(x)
        dRdxNotUnit = np.dot(dcrossProductds, dsdx) + np.dot(dcrossProductdt, dtdx)
        dOneByCrossProductMagdx = (-1.0 / crossProductMag**3) * np.dot(crossProduct, dRdxNotUnit)
        dRdx = (1.0 / crossProductMag) * dRdxNotUnit + np.outer(crossProduct, dOneByCrossProductMagdx)

        #term1 = (1.0 / crossProductMag) * dRdxNotUnit
        #term2 = np.outer(crossProduct, dOneByCrossProductMagdx)
        #print('***')
        #print('dRdx')
        #print('crossProduct: ', crossProduct)
        #print('dRdxNotUnit column1: ', dRdxNotUnit[0:3,1])
        #print(crossProduct[0]*dRdxNotUnit[0,1], crossProduct[1]*dRdxNotUnit[1,1])
        #print(crossProduct[2], dOneByCrossProductMagdx[1])
        #print(dRdxNotUnit[2,1], term2[2,1])
        #print('***')
        return dRdx

    def rVectorNotUnit(self, x):
        """
        bplane R vector
        """
        s = self.sVector(x)
        t = self.tVector(x)
        R = np.cross(s, t)
        return R

    def bDotR(self, x, mu):
        """
        B dot R
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        bScalar = self.bMag(x, mu)
        bTheta = x[4]
        TA = x[5]

        BR = bScalar * np.sin(bTheta)
        return BR

    def bDotT(self, x, mu):
        """
        B dot T
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        bScalar = self.bMag(x, mu)
        bTheta = x[4]
        TA = x[5]

        BT = bScalar * np.cos(bTheta)
        return BT

    def bVector(self, x, mu):
        """
        b plane vector B
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        rp = x[3]
        bTheta = x[4]
        TA = x[5]

        BR = self.bDotR(x, mu)
        BT = self.bDotT(x, mu)
        R = self.rVector(x)
        T = self.tVector(x)
        b = BR * R + BT * T
        return b

    def bVector_derivs(self, x, mu):
        """
        derivatives of b plane vector B
        typed
        """
        vInfMag = x[0]
        vInfRA = x[1]
        vInfDec = x[2]
        rp = x[3]
        bTheta = x[4]
        TA = x[5]
        
        bScalar = self.bMag(x, mu)

        sTheta = np.sin(bTheta)
        cTheta = np.cos(bTheta)

        dbScalardx = self.bScalar_derivs(x, mu)
        dbThetadx = self.bTheta_derivs(x)
        dsThetadx = cTheta * dbThetadx
        dcThetadx = -sTheta * dbThetadx
        R = self.rVector(x)
        dRdx = self.rVector_derivs(x)
        T = self.tVector(x)
        dTdx = self.tVector_derivs(x)
        dBdx = sTheta * np.outer(R, dbScalardx) + bScalar * np.outer(R, dsThetadx) + bScalar * sTheta * dRdx + cTheta * np.outer(T, dbScalardx) + bScalar * np.outer(T, dcThetadx) + bScalar * cTheta * dTdx
        
        #term1 = sTheta * np.outer(R, dbScalardx)
        #term2 = bScalar * np.outer(R, dsThetadx)
        #term3 = bScalar * sTheta * dRdx
        #term4 = cTheta * np.outer(T, dbScalardx)
        #term5 = bScalar * np.outer(T, dcThetadx)
        #term6 = bScalar * cTheta * dTdx
        #print('***')
        #print('term1: ', term1[2,1])
        #print('term2: ', term2[2,1]) 
        #print('term3: ', term3[2,1]) 
        #print('term4: ', term4[2,1]) 
        #print('term5: ', term5[2,1]) 
        #print('term6: ', term6[2,1])
        #print('***')
        return dBdx

    def hUnit(self, x, mu):
        """
        angular momentum unit vector
        """
        b = self.bVector(x, mu)
        s = self.sVector(x)
        hUnit = np.cross(b, s)
        hUnit = hUnit / np.linalg.norm(hUnit)
        return hUnit

    def hUnit_derivs(self, x, mu):
        """
        derivatives of angular momentum unit vector
        typed
        """
        b = self.bVector(x, mu)
        s = self.sVector(x)
        crossProduct = np.cross(b, s)
        crossProductMag = np.linalg.norm(crossProduct)
        dcrossProductdb, dcrossProductds = self.mathUtil.d_crossproduct(b, s)
        dbdx = self.bVector_derivs(x, mu)
        dsdx = self.sVector_derivs(x)
        dhdxNotUnit = np.dot(dcrossProductdb, dbdx) + np.dot(dcrossProductds, dsdx)
        dOneByCrossProductMagdx = (-1.0 / crossProductMag**3) * np.dot(crossProduct, dhdxNotUnit)        
        dhUnitdx = (1.0 / crossProductMag) * dhdxNotUnit + np.outer(crossProduct, dOneByCrossProductMagdx)
        return dhUnitdx

    def hVector(self, x, mu):
        """
        angular momentum vector
        """
        hUnit = self.hUnit(x, mu)
        hMag = self.hMag(x, mu)
        h = hMag * hUnit
        return h

    def hVector_derivs(self, x, mu):
        """
        angular momentum vector derivatives
        typed
        """
        hMag = self.hMag(x, mu)
        dhMagdx = self.hMag_derivs(x, mu)
        hUnit = self.hUnit(x, mu)
        dhUnitdx = self.hUnit_derivs(x, mu)
        dhdx = hMag * dhUnitdx + np.outer(hUnit, dhMagdx)
        return dhdx

    def TAinf(self, x, mu):
        """
        true anomaly of the incoming asymptote
        """
        eMag = self.eMag(x, mu)
        TAinf = -np.arccos(-1.0/eMag)
        return TAinf

    def TAinf_derivs(self, x, mu):
        """
        derivatives of true anomaly of the incoming asymptote
        typed
        """
        eMag = self.eMag(x, mu)
        deMagdx = self.eMag_derivs(x, mu)
        dTAinfdx = (1.0 / (eMag * np.sqrt(eMag**2 - 1.0))) * deMagdx
        return dTAinfdx

    def eVector(self, x, mu):
        """
        eccentricity vector
        """
        eMag = self.eMag(x, mu)
        TAinf = self.TAinf(x, mu)
        s = self.sVector(x)
        bUnit = self.bVector(x, mu)
        bUnit = bUnit / np.linalg.norm(bUnit)
        e = np.cos(np.pi - TAinf) * s - np.sin(np.pi - TAinf) * bUnit
        e = e / np.linalg.norm(e)
        e = eMag * e
        return e

    def eVector_derivs(self, x, mu):
        """
        eccentricity vector derivatives
        typed
        """
        eMag = self.eMag(x, mu)
        deMagdx = self.eMag_derivs(x, mu)
        eUnit = self.eUnitVector(x, mu)
        deUnitdx = self.eUnitVector_derivs(x, mu)
        dedx = eMag * deUnitdx + np.outer(eUnit, deMagdx)

        return dedx

    def eUnitVector(self, x, mu):
        """
        unit vector in direction of eccentricity vector
        """
        TAinf = self.TAinf(x, mu)
        s = self.sVector(x)
        bUnit = self.bVector(x, mu)
        bUnit = bUnit / np.linalg.norm(bUnit)
        e = self.eVector(x, mu)
        e = e / np.linalg.norm(e)
        return e

    def eUnitVector_derivs(self, x, mu):
        """
        derivatives of unit vector in direction of eccentricity vector
        typed
        """
        TAinf = self.TAinf(x, mu)
        dTAinfdx = self.TAinf_derivs(x, mu)
        s = self.sVector(x)
        dsdx = self.sVector_derivs(x)
        b = self.bVector(x, mu)
        bUnit = b / np.linalg.norm(b)
        dbdx = self.bVector_derivs(x, mu)
        dbUnitdb = self.mathUtil.unit_vector_deriv(b)
        dbUnitdx = np.dot(dbUnitdb, dbdx)

        arg = np.pi - TAinf
        cArg = np.cos(arg)
        sArg = np.sin(arg)
        dcArgdx = sArg * dTAinfdx
        dsArgdx = -cArg * dTAinfdx

        term = cArg * s - sArg * bUnit
        
        # term is a unit vector because s and bUnit are unit vectors
        # so we don't need to add the magnitude terms
        #termMag = np.linalg.norm(term)

        dTermdx = dsdx * cArg + np.outer(s, dcArgdx) - dbUnitdx * sArg - np.outer(bUnit, dsArgdx)
        
        #dOneByTermMagdx = (-1.0 / termMag**3) * np.dot(term, dTermdx)
        #deUnitdx = (1.0 / termMag) * dTermdx + np.outer(term, dOneByTermMagdx)
        
        deUnitdx = dTermdx

        return deUnitdx

    def positionVector(self, x, mu):
        """
        position vector
        """
        TA = x[5]
        h = self.hVector(x, mu)
        e = self.eVector(x, mu)
        hMag = self.hMag(x, mu)
        eMag = self.eMag(x, mu)
        crossp = np.dot(self.mathUtil.crossmat(h), e)
        crosspMag = np.linalg.norm(crossp)

        rMag = hMag**2 / (mu * (1.0 + eMag * np.cos(TA)))
        term1 = (np.cos(TA) / eMag) * e
        term2 = (np.sin(TA) / crosspMag) * crossp
        r = rMag * (term1 + term2)
        return r

    def positionVectorFromheTA(self, h, e, TA, mu):
        """
        position from from h, e, TA for testing purposes
        """
        hMag = np.linalg.norm(h)
        eMag = np.linalg.norm(e)
        crossp = np.dot(self.mathUtil.crossmat(h), e)
        crosspMag = np.linalg.norm(crossp)

        rMag = hMag**2 / (mu * (1.0 + eMag * np.cos(TA)))
        term1 = (np.cos(TA) / eMag) * e
        term2 = (np.sin(TA) / crosspMag) * crossp
        r = rMag * (term1 + term2)
        return r

    def dPositionVectordh(self, x, mu):
        """
        derivatives of position vector wrt angular momentum vector
        typed
        """
        TA = x[5]
        cTA = np.cos(TA)
        sTA = np.sin(TA)
        h = self.hVector(x, mu)
        e = self.eVector(x, mu)
        hMag = self.hMag(x, mu)
        eMag = self.eMag(x, mu)
        eUnit = e / eMag
        crossProduct = np.cross(h, e)
        crossProductMag = np.linalg.norm(crossProduct)

        factor = 1.0 / (mu * (1.0 + eMag * cTA))

        term1 = cTA * np.outer(eUnit, 2.0 * h)
        #term2 = sTA * np.dot((1.0 / crossProductMag) * (np.identity(3) - (1.0 / crossProductMag**2) * np.outer(crossProduct, crossProduct)), -self.mathUtil.crossmat(e))
        
        term2 = (sTA / crossProductMag) * (
            2. * np.outer(crossProduct, h) + 
            hMag**2 * np.dot((-np.identity(3) + (1. / crossProductMag**2) * 
                       np.outer(crossProduct, crossProduct)), 
                             self.mathUtil.crossmat(e)))
        drdh = factor * (term1 + term2)
        return drdh

    def dPositionVectorde(self, x, mu):
        """
        derivatives of position vector wrt eccentricity vector
        typed
        """
        TA = x[5]
        cTA = np.cos(TA)
        sTA = np.sin(TA)
        h = self.hVector(x, mu)
        e = self.eVector(x, mu)
        hMag = self.hMag(x, mu)
        eMag = self.eMag(x, mu)
        eUnit = e / eMag
        crossProduct = np.cross(h, e)
        crossProductMag = np.linalg.norm(crossProduct)
        crossProductUnit = crossProduct / crossProductMag
        factor = hMag**2 / mu
        #term1_old = np.outer((-cTA / ((1.0 + eMag * cTA)**2)) * eUnit, eUnit * cTA + crossProduct / crossProductMag * sTA)
        
        term1 = (-cTA / (1. + eMag * cTA)**2) * np.outer(
            cTA * eUnit + sTA * crossProductUnit, eUnit)
        
        factor2 = 1.0 / (1.0 + eMag * cTA)
        #term2_old = factor2 * ((cTA / eMag) * (np.identity(3) - (1.0 / eMag**2) * np.outer(e, e)) + (sTA / crossProductMag) * np.dot((np.identity(3) - (1.0 / crossProductMag**2) * np.outer(crossProduct, crossProduct)), self.mathUtil.crossmat(h)))
                
        term2 = factor2 * (
            (cTA / eMag) * (np.identity(3) - (1. / eMag**2) * np.outer(e, e)) +
            (sTA / crossProductMag) * np.dot((np.identity(3) - (1. / crossProductMag**2) * np.outer(crossProduct, crossProduct)), self.mathUtil.crossmat(h)))
        
        #term2 = term2_old
        # print(term1_old)
        # print("\n")
        # print(term1)
        # print("\n")
        # print("\n")
        # print("\n")
        # print("\n")
        # print(term2_old)
        # print("\n")
        # print(term2)
        # print("\n")
        # print(term3)
        # print("\n")
        
        drde = factor * (term1 + term2)
        return drde

    def dPositionVectordTA(self, x, mu):
        """
        derivatives of position vector wrt true anomaly
        typed
        """
        TA = x[5]
        cTA = np.cos(TA)
        sTA = np.sin(TA)
        h = self.hVector(x, mu)
        e = self.eVector(x, mu)
        hMag = self.hMag(x, mu)
        eMag = self.eMag(x, mu)
        eUnit = e / eMag
        crossProduct = np.cross(h, e)
        crossProductMag = np.linalg.norm(crossProduct)
        factor = hMag**2 / mu
        factor1 = (eMag * sTA) / ((1.0 + eMag * cTA)**2)
        term1 = factor1 * (eUnit * cTA + (sTA / crossProductMag) * crossProduct)
        factor2 = 1.0 / (1.0 + eMag * cTA)
        term2 = factor2 * (-sTA * eUnit + (cTA / crossProductMag) * crossProduct)
        drdTA = factor * (term1 + term2)
        return drdTA

    def positionVector_derivs(self, x, mu):
        """
        derivatives of position vector
        typed
        """
        drdh = self.dPositionVectordh(x, mu)
        drde = self.dPositionVectorde(x, mu)
        drdTA = self.dPositionVectordTA(x, mu)
        dhdx = self.hVector_derivs(x, mu)
        dedx = self.eVector_derivs(x, mu)
        dTAdx = self.TA_derivs(x)
        drdx = np.dot(drdh, dhdx) + np.dot(drde, dedx) + np.outer(drdTA, dTAdx)
        return drdx

    def velocityVector(self, x, mu):
        """
        velocity vector
        """

        TA = x[5]
        h = self.hVector(x, mu)
        e = self.eVector(x, mu)
        hMag = np.linalg.norm(h)
        eMag = np.linalg.norm(e)
        crossp = np.dot(self.mathUtil.crossmat(h), e)
        crosspMag = np.linalg.norm(crossp)

        factor = -mu / hMag
        term1 = (np.sin(TA) / eMag) * e
        term2 = -((eMag + np.cos(TA)) / crosspMag) * crossp
        v = factor * (term1 + term2)

        return v

    def velocityVectorFromheTA(self, h, e, TA, mu):
        """
        velocity vector as function of h, e, TA
        """

        hMag = np.linalg.norm(h)
        eMag = np.linalg.norm(e)
        crossp = np.dot(self.mathUtil.crossmat(h), e)
        crosspMag = np.linalg.norm(crossp)

        factor = -mu / hMag
        term1 = (np.sin(TA) / eMag) * e
        term2 = -((eMag + np.cos(TA)) / crosspMag) * crossp
        v = factor * (term1 + term2)

        return v

    def dVelocityVectordh(self, x, mu):
        """
        derivatives of Velocity vector wrt angular momentum vector
        typed
        """
        TA = x[5]
        sTA = np.sin(TA)
        cTA = np.cos(TA)
        h = self.hVector(x, mu)
        e = self.eVector(x, mu)
        eCross = self.mathUtil.crossmat(e)
        hMag = np.linalg.norm(h)
        eMag = np.linalg.norm(e)
        eUnit = e / eMag
        crossProduct = np.cross(h, e)
        crossProductMag = np.linalg.norm(crossProduct)

        factor1 = 1.0 / hMag**3
        term1 = factor1 * (-np.outer(sTA * eUnit - ((eMag + cTA) / crossProductMag) * crossProduct, h))
        factor2 = (-(eMag + cTA)) / (hMag * crossProductMag)
        term2 = factor2 * (-eCross + (1.0 / crossProductMag**2) * np.dot(np.outer(crossProduct, crossProduct), eCross))

        dvdh = -mu * (term1 + term2)
        return dvdh

    def dVelocityVectorde(self, x, mu):
        """
        derivatives of Velocity vector wrt eccentricity vector
        typed
        """
        TA = x[5]
        sTA = np.sin(TA)
        cTA = np.cos(TA)
        h = self.hVector(x, mu)
        e = self.eVector(x, mu)
        hCross = self.mathUtil.crossmat(h)
        hMag = np.linalg.norm(h)
        eMag = np.linalg.norm(e)
        crossProduct = np.cross(h, e)
        crossProductMag = np.linalg.norm(crossProduct)

        term1 = sTA * (1.0 / eMag) * (np.identity(3) - (1.0 / eMag**2) * np.outer(e, e))
        term21 = np.outer((1.0 / crossProductMag) * crossProduct, (1.0 / eMag) * e)
        term22 = (eMag + cTA) * (1.0 / crossProductMag) * (hCross - (1.0 / crossProductMag**2) * np.dot(np.outer(crossProduct, crossProduct), hCross))
        term2 = -(term21 + term22)

        factor = -mu / hMag

        dvde = factor * (term1 + term2)
        return dvde

    def dVelocityVectordTA(self, x, mu):
        """
        derivatives of Velocity vector wrt true anomaly
        typed
        """
        TA = x[5]
        h = self.hVector(x, mu)
        e = self.eVector(x, mu)
        hMag = np.linalg.norm(h)
        eMag = np.linalg.norm(e)
        crossp = np.dot(self.mathUtil.crossmat(h), e)
        crosspMag = np.linalg.norm(crossp)

        factor = -mu / hMag
        term1 = (np.cos(TA) / eMag) * e
        term2 = ((np.sin(TA)) / crosspMag) * crossp
        dvdTA = factor * (term1 + term2)
        return dvdTA

    def velocityVector_derivs(self, x, mu):
        """
        derivatives of Velocity vector
        typed
        """
        dvdh = self.dVelocityVectordh(x, mu)
        dvde = self.dVelocityVectorde(x, mu)
        dvdTA = self.dVelocityVectordTA(x, mu)
        dhdx = self.hVector_derivs(x, mu)
        dedx = self.eVector_derivs(x, mu)
        dTAdx = self.TA_derivs(x)
        dvdx = np.dot(dvdh, dhdx) + np.dot(dvde, dedx) + np.outer(dvdTA, dTAdx)
        return dvdx
        