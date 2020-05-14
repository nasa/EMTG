#import matplotlib.pyplot as plt
import wx
import PyEMTG_interface
import numpy as np
import MissionOptions as MO


if __name__ == '__main__':
    application = wx.App()
    PyEMTG_interface.PyEMTG_interface(None)
    application.MainLoop()