# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 16:17:04 2020

@author: nhatten
"""

"""

Revision history
2018-12-18; Noble Hatten; begun

"""

from mathUtilities import mathUtilities
from posVel2BPlaneOut import posVel2BPlaneOut
from BPlane2PosVelOut import BPlane2PosVelOut
import autograd.numpy as np # use the autograd version of np
from autograd import grad
from autograd import jacobian

if __name__ == "__main__":
    np.set_printoptions(precision=20)
    mathUtil = mathUtilities()

    # initial conditions
    #r = np.array([7000.0, 1000.0, -5000.0])
    #v = np.array([-1.0, 22.0, 7.0])
    
    # outgoing asypmtote (TA = ~10 deg)
    r = np.array([-1389.182993860459, -6830.153393624416, -1183.630020141228])
    v = np.array([10.080073981938, -3.269814704898709, 1.400774599839738])
    
    
    # incoming asymptote (TA = ~350 deg)
    #r = np.array([-3636.648315948511, -5882.749052462134, -1466.239205353441])
    #v = np.array([9.416328692291867, -4.948766622733928, 1.050814939799521])
    
    
    mu = 3.986e5

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

    # autograd derivatives
    
    # true anomaly
    d_trueAnomaly_d_r_func = jacobian(bplane.trueAnomaly, 0)
    d_trueAnomaly_d_r_ad = d_trueAnomaly_d_r_func(r, v, mu)
    d_trueAnomaly_d_v_func = jacobian(bplane.trueAnomaly, 1)
    d_trueAnomaly_d_v_ad = d_trueAnomaly_d_v_func(r, v, mu)
    
    d_trueAnomaly_d_x = bplane.trueAnomaly_derivs(r, v, mu)
    
    # e cross r
    d_eCrossR_d_r_func = jacobian(bplane.eCrossR, 0)
    d_eCrossR_d_r_ad = d_eCrossR_d_r_func(r, v, mu)
    d_eCrossR_d_v_func = jacobian(bplane.eCrossR, 1)
    d_eCrossR_d_v_ad = d_eCrossR_d_v_func(r, v, mu)
    
    d_eCrossR_d_x = bplane.eCrossR_derivs(r, v, mu)

    # eVector
    d_eVector_d_r_func = jacobian(bplane.eVector, 0)
    d_eVector_d_r_ad = d_eVector_d_r_func(r, v, mu)
    d_eVector_d_v_func = jacobian(bplane.eVector, 1)
    d_eVector_d_v_ad = d_eVector_d_v_func(r, v, mu)

    d_eVector_d_x = bplane.eVector_derivs(r, v, mu)

    # nVector
    d_nVector_d_r_func = jacobian(bplane.nVector, 0)
    d_nVector_d_r_ad = d_nVector_d_r_func(r, v, mu)
    d_nVector_d_v_func = jacobian(bplane.nVector, 1)
    d_nVector_d_v_ad = d_nVector_d_v_func(r, v, mu)

    d_nVector_d_x = bplane.nVector_derivs(r, v, mu)

    # hVector
    d_hVector_d_r_func = jacobian(bplane.hVector, 0)
    d_hVector_d_r_ad = d_hVector_d_r_func(r, v)
    d_hVector_d_v_func = jacobian(bplane.hVector, 1)
    d_hVector_d_v_ad = d_hVector_d_v_func(r, v)

    d_hVector_d_x = bplane.hVector_derivs(r, v)
    
    # hVector v2
    d_hVector_d_x_v2 = bplane.hVector_derivs_v2(r, v)
    

    # sVector
    d_sVector_d_r_func = jacobian(bplane.sVector, 0)
    d_sVector_d_r_ad = d_sVector_d_r_func(r, v, mu)
    d_sVector_d_v_func = jacobian(bplane.sVector, 1)
    d_sVector_d_v_ad = d_sVector_d_v_func(r, v, mu)

    d_sVector_d_x = bplane.sVector_derivs(r, v, mu)

    # tVector
    d_tVector_d_r_func = jacobian(bplane.tVector, 0)
    d_tVector_d_r_ad = d_tVector_d_r_func(r, v, mu)
    d_tVector_d_v_func = jacobian(bplane.tVector, 1)
    d_tVector_d_v_ad = d_tVector_d_v_func(r, v, mu)

    d_tVector_d_x = bplane.tVector_derivs(r, v, mu)

    # rVector
    d_rVector_d_r_func = jacobian(bplane.rVector, 0)
    d_rVector_d_r_ad = d_rVector_d_r_func(r, v, mu)
    d_rVector_d_v_func = jacobian(bplane.rVector, 1)
    d_rVector_d_v_ad = d_rVector_d_v_func(r, v, mu)

    d_rVector_d_x = bplane.rVector_derivs(r, v, mu)

    # bVector
    d_bVector_d_r_func = jacobian(bplane.bVector, 0)
    d_bVector_d_r_ad = d_bVector_d_r_func(r, v, mu)
    d_bVector_d_v_func = jacobian(bplane.bVector, 1)
    d_bVector_d_v_ad = d_bVector_d_v_func(r, v, mu)

    d_bVector_d_x = bplane.bVector_derivs(r, v, mu)

    # bScalar
    d_bScalar_d_r_func = jacobian(bplane.bScalar, 0)
    d_bScalar_d_r_ad = d_bScalar_d_r_func(r, v, mu)
    d_bScalar_d_v_func = jacobian(bplane.bScalar, 1)
    d_bScalar_d_v_ad = d_bScalar_d_v_func(r, v, mu)

    d_bScalar_d_x = bplane.bScalar_derivs(r, v, mu)

    # bDotR
    d_bDotR_d_r_func = jacobian(bplane.bDotR, 0)
    d_bDotR_d_r_ad = d_bDotR_d_r_func(r, v, mu)
    d_bDotR_d_v_func = jacobian(bplane.bDotR, 1)
    d_bDotR_d_v_ad = d_bDotR_d_v_func(r, v, mu)

    d_BdotR_d_x = bplane.bDotR_derivs(r, v, mu)

    # bDotT
    d_bDotT_d_r_func = jacobian(bplane.bDotT, 0)
    d_bDotT_d_r_ad = d_bDotT_d_r_func(r, v, mu)
    d_bDotT_d_v_func = jacobian(bplane.bDotT, 1)
    d_bDotT_d_v_ad = d_bDotT_d_v_func(r, v, mu)

    d_BdotT_d_x = bplane.bDotT_derivs(r, v, mu)

    # bTheta
    d_bTheta_d_r_func = jacobian(bplane.bTheta, 0)
    d_bTheta_d_r_ad = d_bTheta_d_r_func(r, v, mu)
    d_bTheta_d_v_func = jacobian(bplane.bTheta, 1)
    d_bTheta_d_v_ad = d_bTheta_d_v_func(r, v, mu)

    d_bTheta_dx = bplane.bTheta_derivs(r, v, mu)

    # rPeri
    d_rPeri_d_r_func = jacobian(bplane.rPeri, 0)
    d_rPeri_d_r_ad = d_rPeri_d_r_func(r, v, mu)
    d_rPeri_d_v_func = jacobian(bplane.rPeri, 1)
    d_rPeri_d_v_ad = d_rPeri_d_v_func(r, v, mu)

    d_rPeri_d_x = bplane.rPeri_derivs(r, v, mu)

    # vInfMag
    d_vInfMag_d_r_func = jacobian(bplane.vInfMag, 0)
    d_vInfMag_d_r_ad = d_vInfMag_d_r_func(r, v, mu)
    d_vInfMag_d_v_func = jacobian(bplane.vInfMag, 1)
    d_vInfMag_d_v_ad = d_vInfMag_d_v_func(r, v, mu)

    d_vInfMag_d_x = bplane.vInfMag_derivs(r, v, mu)

    # vInfRA
    d_vInfRA_d_r_func = jacobian(bplane.vInfRA, 0)
    d_vInfRA_d_r_ad = d_vInfRA_d_r_func(r, v, mu)
    d_vInfRA_d_v_func = jacobian(bplane.vInfRA, 1)
    d_vInfRA_d_v_ad = d_vInfRA_d_v_func(r, v, mu)

    d_vInfRA_d_x = bplane.vInfRA_derivs(r, v, mu)

    # vInfDec
    d_vInfDec_d_r_func = jacobian(bplane.vInfDec, 0)
    d_vInfDec_d_r_ad = d_vInfDec_d_r_func(r, v, mu)
    d_vInfDec_d_v_func = jacobian(bplane.vInfDec, 1)
    d_vInfDec_d_v_ad = d_vInfDec_d_v_func(r, v, mu)

    d_vInfDec_d_x = bplane.vInfDec_derivs(r, v, mu)

    # hCross
    d_hCross_d_r_func = jacobian(bplane.hCross, 0)
    d_hCross_d_r_ad = d_hCross_d_r_func(r, v)
    d_hCross_d_v_func = jacobian(bplane.hCross, 1)
    d_hCross_d_v_ad = d_hCross_d_v_func(r, v)

    # periapsis position
    rp = bplane.periapsisPositionVector(r, v, mu)

    # periapsis velocity
    vp =  bplane.periapsisVelocityVector(r, v, mu)

    ## print values
    #print('S = ', S)
    #print('B = ', B)
    #print('R = ', R)
    #print('T = ', T)
    #print('BR = ', BR)
    #print('BT = ', BT)
    #print('theta (deg) = ', np.rad2deg(theta))
    #print('b = ', Bmag)
    #print('rp = ', rp)
    #print('vinf = ', vinf)
    #print('vinf RA = ', np.rad2deg(RA))
    #print('vinf Dec = ', np.rad2deg(Dec))
    #print('')

    ## print derivatives
    # print('')
    # print('true anomaly = ', trueAnomaly)
    # print('d_trueAnomaly_d_r difference = \n', d_trueAnomaly_d_r_ad - d_trueAnomaly_d_x[0:3]) # correct
    # print('d_trueAnomaly_d_v difference = \n', d_trueAnomaly_d_v_ad - d_trueAnomaly_d_x[3:6]) # correct
    # print('d_trueAnomaly_d_r_ad = \n', d_trueAnomaly_d_r_ad)
    # print('d_trueAnomaly_d_v_ad = \n', d_trueAnomaly_d_v_ad)
    # print('d_trueAnomaly_d_r = \n', d_trueAnomaly_d_x[0:3])
    # print('d_trueAnomaly_d_v = \n', d_trueAnomaly_d_x[3:6])
    # print('')
    # print('d_eCrossR_d_r_ad = \n', d_eCrossR_d_r_ad - d_eCrossR_d_x[0:3,0:3]) # correct
    # print('d_eCrossR_d_v_ad = \n', d_eCrossR_d_v_ad - d_eCrossR_d_x[0:3,3:6]) # correct
    # print('d_eCrossR_d_r_ad = \n', d_eCrossR_d_r_ad)
    # print('d_eCrossR_d_v_ad = \n', d_eCrossR_d_v_ad)
    # print('d_eCrossR_d_r = \n', d_eCrossR_d_x[0:3,0:3])
    # print('d_eCrossR_d_v = \n', d_eCrossR_d_x[0:3,3:6])
    # print('')
    #print('d_eVector_d_r_ad = \n', d_eVector_d_r_ad - d_eVector_d_x[0:3,0:3])
    #print('d_eVector_d_v_ad = \n', d_eVector_d_v_ad - d_eVector_d_x[0:3,3:6])
    #print('d_eVector_d_r = \n', d_eVector_d_x[0:3,0:3])
    #print('d_eVector_d_v = \n', d_eVector_d_x[0:3,3:6])
    #print('')
    #print('d_nVector_d_r_ad = \n', d_nVector_d_r_ad - d_nVector_d_x[0:3,0:3])
    #print('d_nVector_d_v_ad = \n', d_nVector_d_v_ad - d_nVector_d_x[0:3,3:6])
    #print('d_nVector_d_r = \n', d_nVector_d_x[0:3,0:3])
    #print('d_nVector_d_v = \n', d_nVector_d_x[0:3,3:6])
    #print('')
    #print('d_hVector_d_r_ad = \n', d_hVector_d_r_ad - d_hVector_d_x[0:3,0:3])
    #print('d_hVector_d_v_ad = \n', d_hVector_d_v_ad - d_hVector_d_x[0:3,3:6])
    #print('d_hVector_d_r = \n', d_hVector_d_x[0:3,0:3])
    #print('d_hVector_d_v = \n', d_hVector_d_x[0:3,3:6])
    #print('d_hVector_d_r v2 = \n', d_hVector_d_x_v2[0:3,0:3])
    #print('d_hVector_d_v v2 = \n', d_hVector_d_x_v2[0:3,3:6])
    # print('')
    # print('d_sVector_d_r diff = \n', d_sVector_d_r_ad - d_sVector_d_x[0:3,0:3])
    # print('d_sVector_d_v diff = \n', d_sVector_d_v_ad - d_sVector_d_x[0:3,3:6])
    # print('d_sVector_d_r_ad = \n', d_sVector_d_r_ad)
    # print('d_sVector_d_v_ad = \n', d_sVector_d_v_ad)
    # print('d_sVector_d_r = \n', d_sVector_d_x[0:3,0:3])
    # print('d_sVector_d_v = \n', d_sVector_d_x[0:3,3:6])
    # print('')
    #print('d_tVector_d_r_ad = \n', d_tVector_d_r_ad - d_tVector_d_x[0:3,0:3])
    #print('d_tVector_d_v_ad = \n', d_tVector_d_v_ad - d_tVector_d_x[0:3,3:6])
    #print('d_tVector_d_r = \n', d_tVector_d_x[0:3,0:3])
    #print('d_tVector_d_v = \n', d_tVector_d_x[0:3,3:6])
    #print('')
    #print('d_rVector_d_r_ad = \n', d_rVector_d_r_ad - d_rVector_d_x[0:3,0:3])
    #print('d_rVector_d_v_ad = \n', d_rVector_d_v_ad - d_rVector_d_x[0:3,3:6])
    #print('d_rVector_d_r = \n', d_rVector_d_x[0:3,0:3])
    #print('d_rVector_d_v = \n', d_rVector_d_x[0:3,3:6])
    # print('')
    # print('d_bVector_d_r diff = \n', d_bVector_d_r_ad - d_bVector_d_x[0:3,0:3])
    # print('d_bVector_d_v diff = \n', d_bVector_d_v_ad - d_bVector_d_x[0:3,3:6])
    # print('d_bVector_d_r_ad = \n', d_bVector_d_r_ad)
    # print('d_bVector_d_v_ad = \n', d_bVector_d_v_ad)
    # print('d_bVector_d_r = \n', d_bVector_d_x[0:3,0:3])
    # print('d_bVector_d_v = \n', d_bVector_d_x[0:3,3:6])
    # print('')
    #print('d_bScalar_d_r_ad = \n', d_bScalar_d_r_ad - d_bScalar_d_x[0:3])
    #print('d_bScalar_d_v_ad = \n', d_bScalar_d_v_ad - d_bScalar_d_x[3:6])
    #print('d_bScalar_d_r = \n', d_bScalar_d_x[0:3])
    #print('d_bScalar_d_v = \n', d_bScalar_d_x[3:6]) 
    #print('')
    #print('d_bDotR_d_r_ad = \n', d_bDotR_d_r_ad - d_BdotR_d_x[0:3])
    #print('d_bDotR_d_v_ad = \n', d_bDotR_d_v_ad - d_BdotR_d_x[3:6])
    #print('d_bDotR_d_r = \n', d_BdotR_d_x[0:3])
    #print('d_bDotR_d_v = \n', d_BdotR_d_x[3:6])
    #print('')
    #print('d_bDotT_d_r_ad = \n', d_bDotT_d_r_ad - d_BdotT_d_x[0:3])
    #print('d_bDotT_d_v_ad = \n', d_bDotT_d_v_ad - d_BdotT_d_x[3:6])
    #print('d_bDotT_d_r = \n', d_BdotT_d_x[0:3])
    #print('d_bDotT_d_v = \n', d_BdotT_d_x[3:6])
    # print('')
    # print('d_bTheta_d_r_ad = \n', d_bTheta_d_r_ad - d_bTheta_dx[0:3]) # correct
    # print('d_bTheta_d_v_ad = \n', d_bTheta_d_v_ad - d_bTheta_dx[3:6]) # correct
    # print('d_bTheta_d_r = \n', d_bTheta_dx[0:3])
    # print('d_bTheta_d_v = \n', d_bTheta_dx[3:6])
    # print('')
    # print('d_rPeri_d_r_ad = \n', d_rPeri_d_r_ad - d_rPeri_d_x[0:3]) # correct
    # print('d_rPeri_d_v_ad = \n', d_rPeri_d_v_ad - d_rPeri_d_x[3:6]) # correct
    # print('d_rPeri_d_r = \n', d_rPeri_d_x[0:3])
    # print('d_rPeri_d_v = \n', d_rPeri_d_x[3:6])
    #print('')
    # print('d_vInfMag_d_r_ad = \n', d_vInfMag_d_r_ad - d_vInfMag_d_x[0:3]) # correct
    # print('d_vInfMag_d_v_ad = \n', d_vInfMag_d_v_ad - d_vInfMag_d_x[3:6]) # correct
    # print('d_vInfMag_d_r = \n', d_vInfMag_d_x[0:3])
    # print('d_vInfMag_d_v = \n', d_vInfMag_d_x[3:6])
    #print('')
    # print('d_vInfRA_d_r_ad = \n', d_vInfRA_d_r_ad - d_vInfRA_d_x[0:3]) # correct
    # print('d_vInfRA_d_v_ad = \n', d_vInfRA_d_v_ad - d_vInfRA_d_x[3:6]) # correct
    # print('d_vInfRA_d_r = \n', d_vInfRA_d_x[0:3])
    # print('d_vInfRA_d_v = \n', d_vInfRA_d_x[3:6])
    #print('')
    # print('d_vInfDec_d_r_ad = \n', d_vInfDec_d_r_ad - d_vInfDec_d_x[0:3]) # correct
    # print('d_vInfDec_d_v_ad = \n', d_vInfDec_d_v_ad - d_vInfDec_d_x[3:6]) # correct
    # print('d_vInfDec_d_r = \n', d_vInfDec_d_x[0:3])
    # print('d_vInfDec_d_v = \n', d_vInfDec_d_x[3:6])
    

    # now, go back the other way
    TA0 = 0.0 # fix at periapsis
    
    # true anomaly not fixed at periapsis, but taken from initial Cartesian state
    TA = trueAnomaly
    
    x = np.array([vinf, RA, Dec, Bmag, theta, TA])
    print(x)
    xp = np.array([vinf, RA, Dec, Bmag, theta, TA0]) # periapsis
    #x = xp
    
    xDeg = np.zeros((6))
    xDeg[0] = x[0]
    xDeg[1] = x[1] * 180./np.pi
    xDeg[2] = x[2] * 180./np.pi
    xDeg[3] = x[3]
    xDeg[4] = x[4] * 180./np.pi
    xDeg[5] = x[5] * 180./np.pi
    
    # print('B plane state = \n', x)
    # print("\n")
    # print('B plane state (deg) = \n', xDeg)
    # print("\n")
    
    back = BPlane2PosVelOut()

    # calculate the stuff
    eMagBack = back.eMag(x, mu)
    sBack = back.sVector(x)
    hMagBack = back.hMag(x)
    tBack = back.tVector(x)
    rBack = back.rVector(x)
    BRBack = back.bDotR(x)
    BTBack = back.bDotT(x)
    BBack = back.bVector(x)
    hUnitBack = back.hUnit(x)
    hBack = back.hVector(x)
    eBack = back.eVector(x, mu)
    rpBack = back.positionVector(xp, mu)
    vpBack = back.velocityVector(xp, mu)
    rBack = back.positionVector(x, mu)
    vBack = back.velocityVector(x, mu)

    # autograd derivatives

    # eMag
    d_eMag_d_x_func = jacobian(back.eMag, 0)
    d_eMag_d_x_ad = d_eMag_d_x_func(x, mu)
    d_eMag_d_x = back.eMag_derivs(x, mu)

    # hMag
    d_hMag_d_x_func = jacobian(back.hMag, 0)
    d_hMag_d_x_ad = d_hMag_d_x_func(x)
    d_hMag_d_x = back.hMag_derivs(x)

    # sVector
    d_sVector_d_x_func = jacobian(back.sVector, 0)
    d_sVector_d_x_ad = d_sVector_d_x_func(x)
    d_sVector_d_x = back.sVector_derivs(x)

    # tVector
    d_tVector_d_x_func = jacobian(back.tVector, 0)
    d_tVector_d_x_ad = d_tVector_d_x_func(x)
    d_tVector_d_x = back.tVector_derivs(x)

    # rVector
    d_rVector_d_x_func = jacobian(back.rVector, 0)
    d_rVector_d_x_ad = d_rVector_d_x_func(x)
    d_rVector_d_x = back.rVector_derivs(x)

    # bVector
    d_bVector_d_x_func = jacobian(back.bVector, 0)
    d_bVector_d_x_ad = d_bVector_d_x_func(x)
    d_bVector_d_x = back.bVector_derivs(x)

    # hUnit
    d_hUnit_d_x_func = jacobian(back.hUnit, 0)
    d_hUnit_d_x_ad = d_hUnit_d_x_func(x)
    d_hUnit_d_x = back.hUnit_derivs(x)

    # hVector
    d_hVector_d_x_func = jacobian(back.hVector, 0)
    d_hVector_d_x_ad = d_hVector_d_x_func(x)
    d_hVector_d_x = back.hVector_derivs(x)

    # eUnitVector
    d_eUnitVector_d_x_func = jacobian(back.eUnitVector, 0)
    d_eUnitVector_d_x_ad = d_eUnitVector_d_x_func(x, mu)
    d_eUnitVector_d_x = back.eUnitVector_derivs(x, mu)

    # eVector
    d_eVector_d_x_func = jacobian(back.eVector, 0)
    d_eVector_d_x_ad = d_eVector_d_x_func(x, mu)
    d_eVector_d_x = back.eVector_derivs(x, mu)

    # TAinf
    d_TAinf_d_x_func = jacobian(back.TAinf, 0)
    d_TAinf_d_x_ad = d_TAinf_d_x_func(x, mu)
    d_TAinf_d_x = back.TAinf_derivs(x, mu)

    # position wrt ang mo
    d_position_d_h_func = jacobian(back.positionVectorFromheTA, 0)
    d_position_d_h_ad = d_position_d_h_func(hBack, eBack, x[5], mu)
    d_position_d_h = back.dPositionVectordh(x, mu)

    # position wrt ecc vec
    d_position_d_e_func = jacobian(back.positionVectorFromheTA, 1)
    d_position_d_e_ad = d_position_d_e_func(hBack, eBack, x[5], mu)
    d_position_d_e = back.dPositionVectorde(x, mu)

    # position wrt ecc vec
    d_position_d_TA_func = jacobian(back.positionVectorFromheTA, 2)
    d_position_d_TA_ad = d_position_d_TA_func(hBack, eBack, x[5], mu)
    d_position_d_TA = back.dPositionVectordTA(x, mu)

    # velocity wrt ang mo
    d_velocity_d_h_func = jacobian(back.velocityVectorFromheTA, 0)
    d_velocity_d_h_ad = d_velocity_d_h_func(hBack, eBack, x[5], mu)
    d_velocity_d_h = back.dVelocityVectordh(x, mu)

    # velocity wrt ecc vec
    d_velocity_d_e_func = jacobian(back.velocityVectorFromheTA, 1)
    d_velocity_d_e_ad = d_velocity_d_e_func(hBack, eBack, x[5], mu)
    d_velocity_d_e = back.dVelocityVectorde(x, mu)

    # velocity wrt TA
    d_velocity_d_TA_func = jacobian(back.velocityVectorFromheTA, 2)
    d_velocity_d_TA_ad = d_velocity_d_TA_func(hBack, eBack, x[5], mu)
    d_velocity_d_TA = back.dVelocityVectordTA(x, mu)

    # positionVector
    d_positionVector_d_x_func = jacobian(back.positionVector, 0)
    d_positionVector_d_x_ad = d_positionVector_d_x_func(x, mu)
    d_positionVector_d_x = back.positionVector_derivs(x, mu)

    # velocityVector
    d_velocityVector_d_x_func = jacobian(back.velocityVector, 0)
    d_velocityVector_d_x_ad = d_velocityVector_d_x_func(x, mu)
    d_velocityVector_d_x = back.velocityVector_derivs(x, mu)

    ## print derivatives
    print('BPlane2PosVel derivatives')
    print('')
    #print('d_eMag_d_x_ad = \n', d_eMag_d_x_ad - d_eMag_d_x) # correct
    #print('d_eMag_d_x = \n', d_eMag_d_x) # correct
    #print('')
    #print('d_hMag_d_x_ad = \n', d_hMag_d_x_ad - d_hMag_d_x) # correct
    #print('d_hMag_d_x = \n', d_hMag_d_x) # correct
    #print('')
    #print('d_sVector_d_x_ad = \n', d_sVector_d_x_ad - d_sVector_d_x) # correct
    #print('d_sVector_d_x = \n', d_sVector_d_x) # correct
    #print('')
    #print('d_tVector_d_x_ad = \n', d_tVector_d_x_ad - d_tVector_d_x) # correct
    #print('d_tVector_d_x = \n', d_tVector_d_x) # correct
    #print('')
    #print('d_rVector_d_x_ad = \n', d_rVector_d_x_ad - d_rVector_d_x) # correct
    #print('d_rVector_d_x = \n', d_rVector_d_x) # correct
    # print('')
    # print('d_bVector_d_x_ad = \n', d_bVector_d_x_ad - d_bVector_d_x) # doesn't look exactly correct, but that is because of cross product terms that should cancel not EXACTLY canceling
    # ## this occurs in rVector_derivs, actually: print(crossProduct[0]*dRdxNotUnit[0,1], crossProduct[1]*dRdxNotUnit[1,1])
    # print('d_bVector_d_x = \n', d_bVector_d_x) # correct
    # print('')
    # print('d_TAinf_d_x diff = \n', d_TAinf_d_x_ad - d_TAinf_d_x) # correct
    # print('d_TAinf_d_x_ad = \n', d_TAinf_d_x_ad) # correct
    # print('d_TAinf_d_x = \n', d_TAinf_d_x) # correct
    # print('')
    #print('d_hUnit_d_x_ad = \n', d_hUnit_d_x_ad - d_hUnit_d_x) # correct
    #print('d_hUnit_d_x = \n', d_hUnit_d_x) # correct
    # print('')
    # print('d_hVector_d_x_ad = \n', d_hVector_d_x_ad - d_hVector_d_x) # should be correct, except there are scaling issues that result in small errors in cross product cancellations in hUnitVector_derivs resulting in slightly larger errors later
    # print('d_hVector_d_x = \n', d_hVector_d_x) # correct
    # print('')
    # print('d_eUnitVector_d_x diff = \n', d_eUnitVector_d_x_ad - d_eUnitVector_d_x) # correct
    # print('d_eUnitVector_d_x_ad = \n', d_eUnitVector_d_x_ad) # correct
    # print('d_eUnitVector_d_x = \n', d_eUnitVector_d_x) # correct
    # print('')
    # print('d_eVector_d_x_ad = \n', d_eVector_d_x_ad - d_eVector_d_x) # correct
    # print('d_eVector_d_x = \n', d_eVector_d_x) # correct
    # print('')
    # print('d_positionVector_d_h difference = \n', d_position_d_h_ad - d_position_d_h) # correct
    # print('d_positionVector_d_h_ad = \n', d_position_d_h_ad) # correct
    # print('d_positionVector_d_h = \n', d_position_d_h) # correct
    # print('')
    # print('d_positionVector_d_e difference = \n', d_position_d_e_ad - d_position_d_e) # correct
    # print('d_positionVector_d_e_ad = \n', d_position_d_e_ad) # correct
    # print('d_positionVector_d_e = \n', d_position_d_e) # correct
    # print('')
    # print('d_positionVector_d_TA_ad = \n', d_position_d_TA_ad - d_position_d_TA) # correct
    # print('d_positionVector_d_TA = \n', d_position_d_TA) # correct
    print('')
    print('d_positionVector_d_x difference = \n', d_positionVector_d_x_ad - d_positionVector_d_x) # wrong
    print('d_positionVector_d_x_ad = \n', d_positionVector_d_x_ad) # wrong
    print('d_positionVector_d_x = \n', d_positionVector_d_x) # wrong
    print('')
    # print('d_velocityVector_d_h_ad = \n', d_velocity_d_h_ad - d_velocity_d_h) # correct
    # print('d_velocityVector_d_h = \n', d_velocity_d_h) # correct
    #print('')
    #print('d_velocityVector_d_e_ad = \n', d_velocity_d_e_ad - d_velocity_d_e) # correct
    #print('d_velocityVector_d_e = \n', d_velocity_d_e) # correct
    #print('')
    #print('d_velocityVector_d_TA_ad = \n', d_velocity_d_TA_ad - d_velocity_d_TA) # correct
    #print('d_velocityVector_d_TA = \n', d_velocity_d_TA) # correct
    print('')
    print('dvelocityVector_d_x difference = \n', d_velocityVector_d_x_ad - d_velocityVector_d_x) # correct
    print('dvelocityVector_d_x_ad = \n', d_velocityVector_d_x_ad) # correct
    print('dvelocityVector_d_x = \n', d_velocityVector_d_x) # correct


    # print the stuff
    #print('')
    #print('*****************************')
    #print('*****************************')
    #print('*****************************')
    #print('')
    #print('eMag forward = ', np.linalg.norm(e) - eMagBack)
    #print('eMag back = ', eMagBack)
    #print('')
    #print('s vector forward = ', S - sBack)
    #print('s vector backward = ', sBack)
    #print('')
    #print('hMag forward = ', np.linalg.norm(h) - hMagBack)
    #print('hMag backward = ', hMagBack)
    #print('')
    #print('tVector forward = ', T - tBack)
    #print('tVector backward = ', tBack)
    #print('')
    #print('rVector forward = ', R - rBack)
    #print('rVector backward = ', rBack)
    #print('')
    #print('bDotT forward = ', BT - BTBack)
    #print('bDotT backward = ', BTBack)
    #print('')
    #print('bDotR forward = ', BR - BRBack)
    #print('bDotR backward = ', BRBack)
    #print('')
    #print('BVector forward = ', B - BBack)
    #print('BVector backward = ', BBack)
    #print('')
    #print('hUnit forward = ', h/np.linalg.norm(h) - hUnitBack)
    #print('hUnit backward = ', hUnitBack)
    #print('')
    #print('h forward = ', h - hBack)
    #print('h backward = ', hBack)
    #print('')
    #print('e forward = ', e - eBack)
    #print('e backward = ', eBack)
    # print('')
    # print('periapsis position fwd - bwd = ', rp - rpBack)
    # print('periapsis position backward = ', rpBack)
    # print('')
    # print('periapsis velocity fwd - bwd = ', vp - vpBack)
    # print('periapsis velocity backward = ', vpBack)
    # print('')
    # print('position fwd - bwd = ', r - rBack)
    # print('position forward = ', r)
    # print('position backward = ', rBack)
    # print('')
    # print('velocity fwd - bwd = ', v - vBack)
    # print('velocity forward = ', v)
    # print('velocity backward = ', vBack)
    # print('')

    print("Done!")
