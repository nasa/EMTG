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

import TankPanel

class LibraryHardwarePanel(wx.Panel):
    def __init__(self, parent, missionoptions):
        self.missionoptions = missionoptions
        self.parent = parent
        
        wx.Panel.__init__(self, parent)
        
        self.librarygridtitle = wx.StaticText(self, -1, "Hardware library options")
        self.btnSetDefaultLibraries = wx.Button(self, -1, "Reset to default library files")

        #path selector

        self.lblHardwarePath = wx.StaticText(self, -1, "Hardware library path")
        self.txtHardwarePath = wx.TextCtrl(self, -1, "HardwarePath", size=(400,-1))
        self.btnGetHardwarePath = wx.Button(self, -1, "...")
        HardwarePathSizerSizer = wx.BoxSizer(wx.HORIZONTAL)
        HardwarePathSizerSizer.AddMany([self.txtHardwarePath, self.btnGetHardwarePath])

        #file selectors

        self.lblLaunchVehicleLibraryFile = wx.StaticText(self, -1, "Launch vehicle library file")
        self.txtLaunchVehicleLibraryFile = wx.TextCtrl(self, -1, "LaunchVehicleLibraryFile", size=(400,-1))
        self.btnGetLaunchVehicleLibraryFile = wx.Button(self, -1, "...")
        LaunchVehicleLibraryFileSizer = wx.BoxSizer(wx.HORIZONTAL)
        LaunchVehicleLibraryFileSizer.AddMany([self.txtLaunchVehicleLibraryFile, self.btnGetLaunchVehicleLibraryFile])

        self.lblSpacecraftOptionsFile = wx.StaticText(self, -1, "Spacecraft file")
        self.txtSpacecraftOptionsFile = wx.TextCtrl(self, -1, "SpacecraftOptionsFile", size=(400,-1))
        self.btnGetSpacecraftOptionsFile = wx.Button(self, -1, "...")
        SpacecraftOptionsFileSizer = wx.BoxSizer(wx.HORIZONTAL)
        SpacecraftOptionsFileSizer.AddMany([self.txtSpacecraftOptionsFile, self.btnGetSpacecraftOptionsFile])

        self.lblPowerSystemsLibraryFile = wx.StaticText(self, -1, "Power systems library file")
        self.txtPowerSystemsLibraryFile = wx.TextCtrl(self, -1, "PowerSystemsLibraryFile", size=(400,-1))
        self.btnGetPowerSystemsLibraryFile = wx.Button(self, -1, "...")
        PowerSystemsLibraryFileSizer = wx.BoxSizer(wx.HORIZONTAL)
        PowerSystemsLibraryFileSizer.AddMany([self.txtPowerSystemsLibraryFile, self.btnGetPowerSystemsLibraryFile])

        self.lblPropulsionSystemsLibraryFile = wx.StaticText(self, -1, "Propulsion systems library file")
        self.txtPropulsionSystemsLibraryFile = wx.TextCtrl(self, -1, "PropulsionSystemsLibraryFile", size=(400,-1))
        self.btnGetPropulsionSystemsLibraryFile = wx.Button(self, -1, "...")
        PropulsionSystemsLibraryFileSizer = wx.BoxSizer(wx.HORIZONTAL)
        PropulsionSystemsLibraryFileSizer.AddMany([self.txtPropulsionSystemsLibraryFile, self.btnGetPropulsionSystemsLibraryFile])

        #item selectors
        self.lblLaunchVehicleKey = wx.StaticText(self, -1, "Launch vehicle")
        self.lblPowerSystemKey = wx.StaticText(self, -1, "Power system")
        self.lblElectricPropulsionSystemKey = wx.StaticText(self, -1, "Electric propulsion system")
        self.lblChemicalPropulsionSystemKey = wx.StaticText(self, -1, "Chemical propulsion system")
        
        self.LaunchVehicleChoices = [self.missionoptions.LaunchVehicleKey]
        self.cmbLaunchVehicleKey = wx.ComboBox(self, -1, choices=self.LaunchVehicleChoices, style=wx.CB_READONLY)
        self.PowerSystemChoices = [self.missionoptions.PowerSystemKey]
        self.cmbPowerSystemKey = wx.ComboBox(self, -1, choices=self.PowerSystemChoices, style=wx.CB_READONLY)
        self.ElectricPropulsionSystemChoices = [self.missionoptions.ElectricPropulsionSystemKey]
        self.cmbElectricPropulsionSystemKey = wx.ComboBox(self, -1, choices=self.ElectricPropulsionSystemChoices, style=wx.CB_READONLY)
        self.ChemicalPropulsionSystemChoices = [self.missionoptions.ChemicalPropulsionSystemKey]
        self.cmbChemicalPropulsionSystemKey = wx.ComboBox(self, -1, choices=self.ChemicalPropulsionSystemChoices, style=wx.CB_READONLY)
              

        self.lblnumber_of_electric_propulsion_systems = wx.StaticText(self, -1, "Number of thrusters")
        self.txtnumber_of_electric_propulsion_systems = wx.TextCtrl(self, -1, "number_of_electric_propulsion_systems")

        librarygrid = wx.FlexGridSizer(20, 2, 5, 5)
        librarygrid.AddMany([self.lblHardwarePath, HardwarePathSizerSizer,
                             self.lblLaunchVehicleLibraryFile, LaunchVehicleLibraryFileSizer,
                             self.lblSpacecraftOptionsFile, SpacecraftOptionsFileSizer,
                             self.lblPowerSystemsLibraryFile, PowerSystemsLibraryFileSizer,
                             self.lblPropulsionSystemsLibraryFile, PropulsionSystemsLibraryFileSizer,
                             self.lblLaunchVehicleKey, self.cmbLaunchVehicleKey,
                             self.lblPowerSystemKey, self.cmbPowerSystemKey,
                             self.lblElectricPropulsionSystemKey, self.cmbElectricPropulsionSystemKey,
                             self.lblChemicalPropulsionSystemKey, self.cmbChemicalPropulsionSystemKey,
                             self.lblnumber_of_electric_propulsion_systems, self.txtnumber_of_electric_propulsion_systems])

        librarybox = wx.BoxSizer(wx.VERTICAL)
        librarybox.AddMany([self.librarygridtitle, librarygrid, self.btnSetDefaultLibraries])

        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.librarygridtitle.SetFont(font)

        

        self.TankPanel = TankPanel.TankPanel(self, self.missionoptions)
        tanksbox = wx.BoxSizer(wx.VERTICAL)
        tanksbox.Add(self.TankPanel)

        #now tie everything together
        self.mainbox = wx.BoxSizer(wx.VERTICAL)
        self.mainbox.AddMany([librarybox, tanksbox])

        self.SetSizer(self.mainbox)

        #bindings
        #library file options
        self.btnSetDefaultLibraries.Bind(wx.EVT_BUTTON, self.ClickSetDefaultLibraries)    
        
        self.txtHardwarePath.Bind(wx.EVT_KILL_FOCUS,self.ChangeHardwarePath)
        self.btnGetHardwarePath.Bind(wx.EVT_BUTTON,self.GetHardwarePath)
        self.txtLaunchVehicleLibraryFile.Bind(wx.EVT_KILL_FOCUS,self.ChangeLaunchVehicleLibraryFile)
        self.btnGetLaunchVehicleLibraryFile.Bind(wx.EVT_BUTTON,self.GetLaunchVehicleLibraryFile)
        self.txtSpacecraftOptionsFile.Bind(wx.EVT_KILL_FOCUS,self.ChangeSpacecraftOptionsFile)
        self.btnGetSpacecraftOptionsFile.Bind(wx.EVT_BUTTON,self.GetSpacecraftOptionsFile)
        self.txtPowerSystemsLibraryFile.Bind(wx.EVT_KILL_FOCUS,self.ChangePowerSystemsLibraryFile)
        self.btnGetPowerSystemsLibraryFile.Bind(wx.EVT_BUTTON,self.GetPowerSystemsLibraryFile)
        self.txtPropulsionSystemsLibraryFile.Bind(wx.EVT_KILL_FOCUS,self.ChangePropulsionSystemsLibraryFile)
        self.btnGetPropulsionSystemsLibraryFile.Bind(wx.EVT_BUTTON,self.GetPropulsionSystemsLibraryFile)

        self.cmbLaunchVehicleKey.Bind(wx.EVT_COMBOBOX, self.ChangeLaunchVehicleKey)
        self.cmbPowerSystemKey.Bind(wx.EVT_COMBOBOX, self.ChangePowerSystemKey)
        self.cmbElectricPropulsionSystemKey.Bind(wx.EVT_COMBOBOX, self.ChangeElectricPropulsionSystemKey)
        self.cmbChemicalPropulsionSystemKey.Bind(wx.EVT_COMBOBOX, self.ChangeChemicalPropulsionSystemKey)

        self.txtnumber_of_electric_propulsion_systems.Bind(wx.EVT_KILL_FOCUS, self.Changenumber_of_electric_propulsion_systems)

    def update(self, updateTanks=True):
        import os
        import re

        self.txtHardwarePath.SetValue(self.missionoptions.HardwarePath)
        self.txtLaunchVehicleLibraryFile.SetValue(self.missionoptions.LaunchVehicleLibraryFile)
        self.txtSpacecraftOptionsFile.SetValue(self.missionoptions.SpacecraftOptionsFile)
        self.txtPowerSystemsLibraryFile.SetValue(self.missionoptions.PowerSystemsLibraryFile)
        self.txtPropulsionSystemsLibraryFile.SetValue(self.missionoptions.PropulsionSystemsLibraryFile)

        self.txtnumber_of_electric_propulsion_systems.SetValue(str(self.missionoptions.number_of_electric_propulsion_systems))

        #now, ingest the hardware files and update the choices for each library
        #first launch vehicle
        if os.path.isdir(self.missionoptions.HardwarePath):
            if os.path.exists(os.path.join(self.missionoptions.HardwarePath, self.missionoptions.LaunchVehicleLibraryFile)):
                self.LaunchVehicleChoices = []
                with open(os.path.join(self.missionoptions.HardwarePath, self.missionoptions.LaunchVehicleLibraryFile), 'r') as library:
                    for line in library:
                        if not line.startswith('#'):
                            self.LaunchVehicleChoices.append(re.split('[, ]', line)[0])
                self.cmbLaunchVehicleKey.Set(self.LaunchVehicleChoices)
                self.cmbLaunchVehicleKey.SetSelection(self.LaunchVehicleChoices.index(self.missionoptions.LaunchVehicleKey))
            else:
                self.append = ['FILE NOT FOUND']
                self.cmbLaunchVehicleKey.Set(self.append)
                self.cmbLaunchVehicleKey.SetSelection(0)
        else:
            self.append = ['PATH NOT FOUND']
            self.cmbLaunchVehicleKey.Set(self.append)
            self.cmbLaunchVehicleKey.SetSelection(0)

        #then the others, if appropriate
        if self.missionoptions.SpacecraftModelInput == 0:
            self.lblSpacecraftOptionsFile.Show(False)
            self.txtSpacecraftOptionsFile.Show(False)
            self.btnGetSpacecraftOptionsFile.Show(False)
            self.lblPowerSystemsLibraryFile.Show(True)
            self.txtPowerSystemsLibraryFile.Show(True)
            self.btnGetPowerSystemsLibraryFile.Show(True)
            self.lblPropulsionSystemsLibraryFile.Show(True)
            self.txtPropulsionSystemsLibraryFile.Show(True)
            self.btnGetPropulsionSystemsLibraryFile.Show(True)
            self.lblPowerSystemKey.Show(True)
            self.lblElectricPropulsionSystemKey.Show(True)
            self.lblChemicalPropulsionSystemKey.Show(True)
            self.cmbPowerSystemKey.Show(True)
            self.cmbElectricPropulsionSystemKey.Show(True)
            self.cmbChemicalPropulsionSystemKey.Show(True)
            self.lblnumber_of_electric_propulsion_systems.Show(True)
            self.txtnumber_of_electric_propulsion_systems.Show(True)
            self.TankPanel.Show(True)
            
        
            if updateTanks:
                self.TankPanel.update()
            
            
            if os.path.isdir(self.missionoptions.HardwarePath):
                if os.path.exists(os.path.join(self.missionoptions.HardwarePath, self.missionoptions.PowerSystemsLibraryFile)):
                    self.PowerSystemChoices = []
                    with open(os.path.join(self.missionoptions.HardwarePath, self.missionoptions.PowerSystemsLibraryFile), 'r') as library:
                        for line in library:
                            if not line.startswith('#'):
                                self.PowerSystemChoices.append(re.split('[, ]', line)[0])
                    self.cmbPowerSystemKey.Set(self.PowerSystemChoices)
                    self.cmbPowerSystemKey.SetSelection(self.PowerSystemChoices.index(self.missionoptions.PowerSystemKey))
                else:
                    self.PowerSystemChoices = ['FILE NOT FOUND']
                    self.cmbPowerSystemKey.Set(self.PowerSystemChoices)
                    self.cmbPowerSystemKey.SetSelection(0)
                                    
                if os.path.exists(os.path.join(self.missionoptions.HardwarePath, self.missionoptions.PropulsionSystemsLibraryFile)):            
                    self.ElectricPropulsionSystemChoices = []
                    with open(os.path.join(self.missionoptions.HardwarePath, self.missionoptions.PropulsionSystemsLibraryFile), 'r') as library:
                        for line in library:
                            if not line.startswith('#'):
                                self.ElectricPropulsionSystemChoices.append(re.split('[, ]', line)[0])
                    self.cmbElectricPropulsionSystemKey.Set(self.ElectricPropulsionSystemChoices)
                    self.cmbElectricPropulsionSystemKey.SetSelection(self.ElectricPropulsionSystemChoices.index(self.missionoptions.ElectricPropulsionSystemKey))
            
                    self.ChemicalPropulsionSystemChoices = []
                    with open(os.path.join(self.missionoptions.HardwarePath, self.missionoptions.PropulsionSystemsLibraryFile), 'r') as library:
                        for line in library:
                            if not line.startswith('#'):
                                self.ChemicalPropulsionSystemChoices.append(re.split('[, ]', line)[0])
                    self.cmbChemicalPropulsionSystemKey.Set(self.ChemicalPropulsionSystemChoices)
                    self.cmbChemicalPropulsionSystemKey.SetSelection(self.ChemicalPropulsionSystemChoices.index(self.missionoptions.ChemicalPropulsionSystemKey))
                else:
                    self.ElectricPropulsionSystemChoices = ['FILE NOT FOUND']
                    self.ChemicalPropulsionSystemChoices = ['FILE NOT FOUND']
                    self.cmbElectricPropulsionSystemKey.Set(self.ElectricPropulsionSystemChoices)
                    self.cmbElectricPropulsionSystemKey.SetSelection(0)
                    self.cmbChemicalPropulsionSystemKey.Set(self.ChemicalPropulsionSystemChoices)
                    self.cmbChemicalPropulsionSystemKey.SetSelection(0)
            else: #bad path
                    self.PowerSystemChoices = ['PATH NOT FOUND']
                    self.cmbPowerSystemKey.Set(self.PowerSystemChoices)
                    self.cmbPowerSystemKey.SetSelection(0)
                    self.ElectricPropulsionSystemChoices = ['PATH NOT FOUND']
                    self.ChemicalPropulsionSystemChoices = ['PATH NOT FOUND']
                    self.cmbElectricPropulsionSystemKey.Set(self.ElectricPropulsionSystemChoices)
                    self.cmbElectricPropulsionSystemKey.SetSelection(0)
                    self.cmbChemicalPropulsionSystemKey.Set(self.ChemicalPropulsionSystemChoices)
                    self.cmbChemicalPropulsionSystemKey.SetSelection(0)
        
        else: #read spacecraft file
            self.lblSpacecraftOptionsFile.Show(True)
            self.txtSpacecraftOptionsFile.Show(True)
            self.btnGetSpacecraftOptionsFile.Show(True)
            self.lblPowerSystemsLibraryFile.Show(False)
            self.txtPowerSystemsLibraryFile.Show(False)
            self.btnGetPowerSystemsLibraryFile.Show(False)
            self.lblPropulsionSystemsLibraryFile.Show(False)
            self.txtPropulsionSystemsLibraryFile.Show(False)
            self.btnGetPropulsionSystemsLibraryFile.Show(False)
            self.lblPowerSystemKey.Show(False)
            self.lblElectricPropulsionSystemKey.Show(False)
            self.lblChemicalPropulsionSystemKey.Show(False)
            self.cmbPowerSystemKey.Show(False)
            self.cmbElectricPropulsionSystemKey.Show(False)
            self.cmbChemicalPropulsionSystemKey.Show(False)
            self.lblnumber_of_electric_propulsion_systems.Show(False)
            self.txtnumber_of_electric_propulsion_systems.Show(False)
            self.TankPanel.Show(False)


        #re-size the panel
        self.Layout()
        if platform.system() == 'Windows':
            self.parent.SetupScrolling(scrollToTop=False)

    #event handlers for spacecraft options

    def ClickSetDefaultLibraries(self, e):
        e.Skip()
        self.missionoptions.HardwarePath = self.parent.parent.Parent.default_HardwarePath
        self.missionoptions.ThrottleTableFile = self.parent.parent.Parent.default_ThrottleTableFile
        self.missionoptions.LaunchVehicleLibraryFile = self.parent.parent.Parent.default_LaunchVehicleLibraryFile
        self.missionoptions.PowerSystemsLibraryFile = self.parent.parent.Parent.default_PowerSystemsLibraryFile
        self.missionoptions.PropulsionSystemsLibraryFile = self.parent.parent.Parent.default_PropulsionSystemsLibraryFile
        self.missionoptions.SpacecraftOptionsFile = self.parent.parent.Parent.default_SpacecraftOptionsFile

        self.update()

        
    def ChangeHardwarePath(self, e):
        e.Skip()
        self.missionoptions.HardwarePath = self.txtHardwarePath.GetValue()
        if self.missionoptions.HardwarePath[-1] not in ['/','\\']:
            self.missionoptions.HardwarePath = self.missionoptions.HardwarePath + '/'
        self.update()

    def ChangeLaunchVehicleLibraryFile(self, e):
        e.Skip()
        self.missionoptions.LaunchVehicleLibraryFile = self.txtLaunchVehicleLibraryFile.GetValue()
        self.update()

    def ChangeSpacecraftOptionsFile(self, e):
        e.Skip()
        self.missionoptions.SpacecraftOptionsFile = self.txtSpacecraftOptionsFile.GetValue()
        self.update()

    def ChangePowerSystemsLibraryFile(self, e):
        e.Skip()
        self.missionoptions.PowerSystemsLibraryFile = self.txtPowerSystemsLibraryFile.GetValue()
        self.update()

    def ChangePropulsionSystemsLibraryFile(self, e):
        e.Skip()
        self.missionoptions.PropulsionSystemsLibraryFile = self.txtPropulsionSystemsLibraryFile.GetValue()
        self.update()

    def GetHardwarePath(self, e):
        #file load dialog to get name of hardware folder
        dlg = wx.DirDialog(self, "Choose a hardware folder", self.missionoptions.HardwarePath)
        if dlg.ShowModal() == wx.ID_OK:
            self.missionoptions.HardwarePath = dlg.GetPath()
            self.txtHardwarePath.SetValue(self.missionoptions.HardwarePath)
        dlg.Destroy()

    def GetLaunchVehicleLibraryFile(self, e):
        dlg = wx.FileDialog(self, "Choose a launch vehicle library file", self.missionoptions.HardwarePath, self.missionoptions.LaunchVehicleLibraryFile, "*.emtg_launchvehicleopt", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
        else:
            filename = self.missionoptions.LaunchVehicleLibraryFile

        dlg.Destroy()

        self.missionoptions.LaunchVehicleLibraryFile = filename

        self.update()

    def GetSpacecraftOptionsFile(self, e):
        dlg = wx.FileDialog(self, "Choose a spacecraft file", self.missionoptions.HardwarePath, self.missionoptions.SpacecraftOptionsFile, "*.emtg_spacecraftopt", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
        else:
            filename = self.missionoptions.SpacecraftOptionsFile

        dlg.Destroy()

        self.missionoptions.SpacecraftOptionsFile = filename

        self.update()

    def GetPowerSystemsLibraryFile(self, e):
        dlg = wx.FileDialog(self, "Choose a power systems library file", self.missionoptions.HardwarePath, self.missionoptions.PowerSystemsLibraryFile, "*.emtg_powersystemsopt", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
        else:
            filename = self.missionoptions.PowerSystemsLibraryFile

        dlg.Destroy()

        self.missionoptions.PowerSystemsLibraryFile = filename

        self.update()

    def GetPropulsionSystemsLibraryFile(self, e):
        dlg = wx.FileDialog(self, "Choose a propulsion systems library file", self.missionoptions.HardwarePath, self.missionoptions.PropulsionSystemsLibraryFile, "*.emtg_propulsionsystemopt", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
        else:
            filename = self.missionoptions.PropulsionSystemsLibraryFile

        dlg.Destroy()

        self.missionoptions.PropulsionSystemsLibraryFile = filename

        self.update()

    def ChangeLaunchVehicleKey(self, e):
        e.Skip()
        self.missionoptions.LaunchVehicleKey = self.cmbLaunchVehicleKey.GetStringSelection()

    def ChangePowerSystemKey(self, e):
        e.Skip()
        self.missionoptions.PowerSystemKey = self.cmbPowerSystemKey.GetStringSelection()

    def ChangeElectricPropulsionSystemKey(self, e):
        e.Skip()
        self.missionoptions.ElectricPropulsionSystemKey = self.cmbElectricPropulsionSystemKey.GetStringSelection()

    def ChangeChemicalPropulsionSystemKey(self, e):
        e.Skip()
        self.missionoptions.ChemicalPropulsionSystemKey = self.cmbChemicalPropulsionSystemKey.GetStringSelection()

    def Changenumber_of_electric_propulsion_systems(self, e):
        e.Skip()
        self.missionoptions.number_of_electric_propulsion_systems = eval(self.txtnumber_of_electric_propulsion_systems.GetValue())
        self.update()