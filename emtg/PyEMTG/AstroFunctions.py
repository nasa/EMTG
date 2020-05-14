# EMTG: Evolutionary Mission Trajectory Generator
# An open-source global optimization tool for preliminary mission design
# Provided by NASA Goddard Space Flight Center
#
# Copyright (c) 2013 - 2020 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Other Rights Reserved.
#
# Licensed under the NASA Open Source License (the "License"); 
# You may not use this file except in compliance with the License. 
# You may obtain a copy of the License at:
# https://opensource.org/licenses/NASA-1.3
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
# express or implied.   See the License for the specific language
# governing permissions and limitations under the License.

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