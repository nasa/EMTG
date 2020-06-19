import numpy as np

def EOM_inertial_2bodyconstant_thrust(t, X, Thrust, Mdot, mu):
    r = np.linalg.norm(X[0:3])
    r3 = r**3
    dX = np.zeros(7)

    dX[0] = X[3]
    dX[1] = X[4]
    dX[2] = X[5]
    dX[3] = -mu*X[0]/r3 + Thrust[0]/X[6]
    dX[4] = -mu*X[1]/r3 + Thrust[1]/X[6]
    dX[5] = -mu*X[2]/r3 + Thrust[2]/X[6]
    dX[6] = -Mdot

    return dX

def EOM_jacobian_inertial_2bodyconstant_thrust(t, X, Thrust, Mdot, mu):
    r = np.linalg.norm(X[0:3])

    r3 = r**3
    r5 = r3*r*r
    ddX = np.zeros((7,7))

    ddX[0,3] = 1.0
    ddX[1,4] = 1.0
    ddX[2,5] = 1.0
    ddX[3,0] = -mu * (1/r3 - 3 * X[0] * X[0] / r5)
    ddX[3,6] = -Thrust[0] / (X[6]**2)
    ddX[4,1] = -mu * (1/r3 - 3 * X[1] * X[1] / r5)
    ddX[4,6] = -Thrust[1] / (X[6]**2)
    ddX[5,2] = -mu * (1/r3 - 3 * X[2] * X[2] / r5)
    ddX[5,6] = -Thrust[2] / (X[6]**2)

    return ddX

def EOM_inertial_2body(t, X, mu):
    r = np.linalg.norm(X[0:3])
    r3 = r**3
    dX = np.zeros(6)

    dX[0] = X[3]
    dX[1] = X[4]
    dX[2] = X[5]
    dX[3] = -mu*X[0]/r3
    dX[4] = -mu*X[1]/r3
    dX[5] = -mu*X[2]/r3

    return dX

def EOM_jacobian_intertial_2body(t, X, mu):
    r = np.linalg.norm(X[0:3])

    r3 = r**3
    r5 = r3*r*r
    ddX = np.zeros((6,6))

    ddX[0,3] = 1.0
    ddX[1,4] = 1.0
    ddX[2,5] = 1.0
    ddX[3,0] = -mu * (1/r3 - 3 * X[0] * X[0] / r5)
    ddX[4,1] = -mu * (1/r3 - 3 * X[1] * X[1] / r5)
    ddX[5,2] = -mu * (1/r3 - 3 * X[2] * X[2] / r5)

    return ddX