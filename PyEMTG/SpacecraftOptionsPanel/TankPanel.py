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

import wx
import wx.adv
import platform

class TankPanel(wx.Panel):
    def __init__(self, parent, missionoptions):
        self.missionoptions = missionoptions
        self.parent = parent
        
        wx.Panel.__init__(self, parent)

        #tank-related fields
        self.tanksgrid = wx.FlexGridSizer(20,2,5,5)
        self.tanksgridtitle = wx.StaticText(self, -1, "Tanks")

        self.lblenable_electric_propellant_tank_constraint = wx.StaticText(self, -1, "Enable electric propulsion propellant tank constraint?")
        self.chkenable_electric_propellant_tank_constraint = wx.CheckBox(self, -1)
        self.lblmaximum_electric_propellant = wx.StaticText(self, -1, "Maximum electric propulsion propellant (kg)")
        self.txtmaximum_electric_propellant = wx.TextCtrl(self, -1, "txtmaximum_electric_propellant")
        self.lblenable_chemical_propellant_tank_constraint = wx.StaticText(self, -1, "Enable chemical propulsion tank constraints?")
        self.chkenable_chemical_propellant_tank_constraint = wx.CheckBox(self, -1)
        self.lblmaximum_chemical_fuel = wx.StaticText(self, -1, "Maximum chemical fuel (kg)")
        self.txtmaximum_chemical_fuel = wx.TextCtrl(self, -1, "txtmaximum_chemical_fuel")
        self.lblmaximum_chemical_oxidizer = wx.StaticText(self, -1, "Maximum chemical oxidizer (kg)")
        self.txtmaximum_chemical_oxidizer = wx.TextCtrl(self, -1, "txtmaximum_chemical_oxidizer")
        self.lblbipropellant_mixture_ratio = wx.StaticText(self, -1, "Bipropellant mixture ratio")
        self.txtbipropellant_mixture_ratio = wx.TextCtrl(self, -1, "txtbipropellant_mixture_ratio")

        self.tanksgrid.AddMany([ self.lblenable_electric_propellant_tank_constraint, self.chkenable_electric_propellant_tank_constraint,
                                 self.lblmaximum_electric_propellant, self.txtmaximum_electric_propellant,
                                 self.lblenable_chemical_propellant_tank_constraint, self.chkenable_chemical_propellant_tank_constraint,
                                 self.lblmaximum_chemical_fuel, self.txtmaximum_chemical_fuel,
                                 self.lblmaximum_chemical_oxidizer, self.txtmaximum_chemical_oxidizer,
                                 self.lblbipropellant_mixture_ratio, self.txtbipropellant_mixture_ratio])

        tanksbox = wx.BoxSizer(wx.VERTICAL)
        tanksbox.AddMany([self.tanksgridtitle, self.tanksgrid])

        #now tie everything together
        self.mainbox = wx.BoxSizer(wx.VERTICAL)
        self.mainbox.Add(tanksbox)

        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.tanksgridtitle.SetFont(font)

        self.SetSizer(self.mainbox)

        #bindings

        #propellant tank options
        self.chkenable_chemical_propellant_tank_constraint.Bind(wx.EVT_CHECKBOX, self.Changeenable_chemical_propellant_tank_constraint)
        self.txtmaximum_chemical_fuel.Bind(wx.EVT_KILL_FOCUS, self.Changemaximum_chemical_fuel)
        self.txtmaximum_chemical_oxidizer.Bind(wx.EVT_KILL_FOCUS, self.Changemaximum_chemical_oxidizer)
        self.txtbipropellant_mixture_ratio.Bind(wx.EVT_KILL_FOCUS, self.Changebipropellant_mixture_ratio)
        self.chkenable_electric_propellant_tank_constraint.Bind(wx.EVT_CHECKBOX, self.Changeenable_electric_propellant_tank_constraint)
        self.txtmaximum_electric_propellant.Bind(wx.EVT_KILL_FOCUS, self.Changemaximum_electric_propellant)

    def update(self):
        self.chkenable_electric_propellant_tank_constraint.SetValue(self.missionoptions.enable_electric_propellant_tank_constraint)
        self.txtmaximum_electric_propellant.SetValue(str(self.missionoptions.maximum_electric_propellant))
        self.chkenable_chemical_propellant_tank_constraint.SetValue(self.missionoptions.enable_chemical_propellant_tank_constraint)
        self.txtmaximum_chemical_fuel.SetValue(str(self.missionoptions.maximum_chemical_fuel))
        self.txtmaximum_chemical_oxidizer.SetValue(str(self.missionoptions.maximum_chemical_oxidizer))
        self.txtbipropellant_mixture_ratio.SetValue(str(self.missionoptions.bipropellant_mixture_ratio))

        if self.missionoptions.mission_type in [0, 1, 2, 3, 5, 9]:
            self.lblenable_electric_propellant_tank_constraint.Show(True)
            self.chkenable_electric_propellant_tank_constraint.Show(True)
        else:
            self.lblenable_electric_propellant_tank_constraint.Show(False)
            self.chkenable_electric_propellant_tank_constraint.Show(False)

        if self.missionoptions.enable_chemical_propellant_tank_constraint:
            self.lblmaximum_chemical_fuel.Show(True)
            self.lblmaximum_chemical_oxidizer.Show(True)
            self.txtmaximum_chemical_fuel.Show(True)
            self.txtmaximum_chemical_oxidizer.Show(True)
        else:
            self.lblmaximum_chemical_fuel.Show(False)
            self.lblmaximum_chemical_oxidizer.Show(False)
            self.txtmaximum_chemical_fuel.Show(False)
            self.txtmaximum_chemical_oxidizer.Show(False)
            

        if self.missionoptions.enable_electric_propellant_tank_constraint:
            self.lblmaximum_electric_propellant.Show(True)
            self.txtmaximum_electric_propellant.Show(True)
        else:
            self.lblmaximum_electric_propellant.Show(False)
            self.txtmaximum_electric_propellant.Show(False)

        #re-size the panel
        self.Layout()
        self.parent.update(updateTanks=False)

    #event handlers
    def Changeenable_chemical_propellant_tank_constraint(self, e):
        e.Skip()
        self.missionoptions.enable_chemical_propellant_tank_constraint = int(self.chkenable_chemical_propellant_tank_constraint.GetValue())
        self.update()

    def Changemaximum_chemical_fuel(self, e):
        e.Skip()
        self.missionoptions.maximum_chemical_fuel = eval(self.txtmaximum_chemical_fuel.GetValue())
        self.update()
        
    def Changemaximum_chemical_oxidizer(self, e):
        e.Skip()
        self.missionoptions.maximum_chemical_oxidizer = eval(self.txtmaximum_chemical_oxidizer.GetValue())
        self.update()
    
    def Changebipropellant_mixture_ratio(self, e):
        e.Skip()
        self.missionoptions.bipropellant_mixture_ratio = eval(self.txtbipropellant_mixture_ratio.GetValue())
        self.update()

    def Changeenable_electric_propellant_tank_constraint(self, e):
        e.Skip()
        self.missionoptions.enable_electric_propellant_tank_constraint = int(self.chkenable_electric_propellant_tank_constraint.GetValue())
        self.update()

    def Changemaximum_electric_propellant(self, e):
        e.Skip()
        self.missionoptions.maximum_electric_propellant = eval(self.txtmaximum_electric_propellant.GetValue())
        self.update()