#EMTG: Evolutionary Mission Trajectory Generator
#An open-source global optimization tool for preliminary mission design
#Provided by NASA Goddard Space Flight Center
#
#Copyright (c) 2014 - 2018 United States Government as represented by the
#Administrator of the National Aeronautics and Space Administration.
#All Other Rights Reserved.
#
#Licensed under the NASA Open Source License (the "License"); 
#You may not use this file except in compliance with the License. 
#You may obtain a copy of the License at:
#https://opensource.org/licenses/NASA-1.3
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
#express or implied.   See the License for the specific language
#governing permissions and limitations under the License.


import sys
import os
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(currentdir + "/" + 'SpacecraftOptionsPanel')
sys.path.append(currentdir + "/" +'JourneyOptionsPanel')
from SpacecraftOptionsPanel import SpacecraftOptionsPanel
from GlobalOptionsPanel import GlobalOptionsPanel
from JourneyOptionsPanel import JourneyOptionsPanel
from SolverOptionsPanel import SolverOptionsPanel
from PhysicsOptionsPanel import PhysicsOptionsPanel
from OutputOptionsPanel import OutputOptionsPanel

import wx

class OptionsBook(wx.Notebook):
    #class for Options notebook
    def __init__(self, parent, options):
        wx.Notebook.__init__(self, parent=parent, id=wx.ID_ANY, style=
                             wx.BK_DEFAULT
                             #wx.BK_TOP 
                             #wx.BK_BOTTOM
                             #wx.BK_LEFT
                             #wx.BK_RIGHT
                             )
                             
        font = self.GetFont()
        font.SetPointSize(10)
        self.SetFont(font)


        #create tabs
        self.tabGlobal = GlobalOptionsPanel(self, options)
        self.AddPage(self.tabGlobal, "Global Mission Options")
        self.tabSpacecraft = SpacecraftOptionsPanel(self, options)
        self.AddPage(self.tabSpacecraft, "Spacecraft Options")
        self.tabJourney = JourneyOptionsPanel(self, options)
        self.AddPage(self.tabJourney, "Journey Options")
        self.tabSolver = SolverOptionsPanel(self, options)
        self.AddPage(self.tabSolver, "Solver Options")
        self.tabPhysics = PhysicsOptionsPanel(self, options)
        self.AddPage(self.tabPhysics, "Physics Options")
        self.tabOutput = OutputOptionsPanel(self, options)
        self.AddPage(self.tabOutput, "Output Options")

    def update(self):
        if hasattr(self, 'tabGlobal'):
            self.tabGlobal.update()
        
        if hasattr(self, 'tabSpacecraft'):
            self.tabSpacecraft.update()
        
        if hasattr(self, 'tabJourney'):
            self.tabJourney.update()
        
        if hasattr(self, 'tabSolver'):
            self.tabSolver.update()
        
        if hasattr(self, 'tabPhysics'):
            self.tabPhysics.update()
        
        if hasattr(self, 'tabOutput'):
            self.tabOutput.update()