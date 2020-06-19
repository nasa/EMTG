import numpy as np
import math

def rotate_from_ecliptic_to_equatorial6(state):
    R3x3 = np.array([1.0,            0.0,            0.0,
                    0.0,     0.91748206,   -0.397777156,
                    0.0,     0.397777156,  0.9174820621]).reshape((3,3))

    Zero3x3 = np.zeros(9).reshape((3,3))

    R6x6 = np.concatenate((np.concatenate((R3x3, Zero3x3), 0), np.concatenate((Zero3x3, R3x3), 0)), 1)

    rotatedstate = np.dot(R6x6, state)

    return rotatedstate

def rotate_from_ecliptic_to_equatorial3(state):
    R3x3 = np.array([1.0,            0.0,            0.0,
                    0.0,     0.91748206,   -0.397777156,
                    0.0,     0.397777156,  0.9174820621]).reshape((3,3))

    rotatedstate = np.dot(R3x3, state)

    return rotatedstate