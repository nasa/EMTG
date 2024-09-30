#EMTG: Evolutionary Mission Trajectory Generator
#An open-source global optimization tool for preliminary mission design
#Provided by NASA Goddard Space Flight Center
#
#Copyright (c) 2014 - 2024 United States Government as represented by the
#Administrator of the National Aeronautics and Space Administration.
#All Other Rights Reserved.
#
#Licensed under the NASA Open Source License (the "License"); 
#You may not use this file except in compliance with the License. 
#You may obtain a copy of the License at:
#https://opensource.org/license/nasa1-3-php
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
#express or implied.   See the License for the specific language
#governing permissions and limitations under the License.

"""
PyEMTG.py
========================
PyEMTG GUI is a graphical user interface (GUI) written in Python to process
EMTG input and output files. It provides users easier access to most of the 
Core EMTG application functionality and contains EMTG post-processing 
functionality.

"""

#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('wxagg') 
import wx

import PyEMTG_interface
import numpy as np
import MissionOptions as MO


if __name__ == '__main__':
    application = wx.App(redirect=False)
    PyEMTG_interface.PyEMTG_interface(None)
    application.MainLoop()