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