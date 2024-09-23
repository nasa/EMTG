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

class ConventionalHardwarePanel(wx.Panel):
    def __init__(self, parent, missionoptions):
        self.missionoptions = missionoptions
        self.parent = parent
        
        wx.Panel.__init__(self, parent)
        self.btnSetDefaultLibraries = wx.Button(self, -1, "Reset to default library files")

        #path selector
        
        #spacecraft and launch vehicle fields
        LaunchVehicleGrid = wx.FlexGridSizer(20,2,5,5)
        self.LaunchVehicleGridTitle = wx.StaticText(self, -1, "Launch Vehicle options")

        self.lblHardwarePath = wx.StaticText(self, -1, "Hardware library path")
        self.txtHardwarePath = wx.TextCtrl(self, -1, "HardwarePath", size=(400,-1))
        self.btnGetHardwarePath = wx.Button(self, -1, "...")
        HardwarePathSizerSizer = wx.BoxSizer(wx.HORIZONTAL)
        HardwarePathSizerSizer.AddMany([self.txtHardwarePath, self.btnGetHardwarePath])

        self.lblLaunchVehicleLibraryFile = wx.StaticText(self, -1, "Launch vehicle library file")
        self.txtLaunchVehicleLibraryFile = wx.TextCtrl(self, -1, "LaunchVehicleLibraryFile", size=(400,-1))
        self.btnGetLaunchVehicleLibraryFile = wx.Button(self, -1, "...")
        LaunchVehicleLibraryFileSizer = wx.BoxSizer(wx.HORIZONTAL)
        LaunchVehicleLibraryFileSizer.AddMany([self.txtLaunchVehicleLibraryFile, self.btnGetLaunchVehicleLibraryFile])

        
        #item selectors
        self.lblLaunchVehicleKey = wx.StaticText(self, -1, "Launch vehicle")
        self.LaunchVehicleChoices = [self.missionoptions.LaunchVehicleKey]
        self.cmbLaunchVehicleKey = wx.ComboBox(self, -1, choices=self.LaunchVehicleChoices, style=wx.CB_READONLY)

        LaunchVehicleGrid.AddMany([self.lblHardwarePath, HardwarePathSizerSizer,
                                   self.lblLaunchVehicleLibraryFile, LaunchVehicleLibraryFileSizer,
                                   self.lblLaunchVehicleKey, self.cmbLaunchVehicleKey,
                                   self.btnSetDefaultLibraries])

        LaunchVehicleBox = wx.BoxSizer(wx.VERTICAL)
        LaunchVehicleBox.AddMany([self.LaunchVehicleGridTitle, LaunchVehicleGrid])

        #propulsion
        propulsiongrid = wx.FlexGridSizer(30,2,5,5)
        self.propulsiongridtitle = wx.StaticText(self, -1, "Propulsion options")

        self.lblIspChem = wx.StaticText(self, -1, "Chemical Isp (s)")
        self.txtIspChem = wx.TextCtrl(self, -1, "IspChem")
        self.lblTCM_Isp = wx.StaticText(self, -1, "TCM Isp (s)")
        self.txtTCM_Isp = wx.TextCtrl(self, -1, "TCM_Isp")

        self.lblengine_type = wx.StaticText(self, -1, "Engine type")
        enginetypes = ['0: fixed thrust/Isp','1: constant Isp, efficiency, EMTG computes input power','2: choice of power model, constant efficiency, EMTG chooses Isp',
                       '3: choice of power model, constant efficiency and Isp','4: continuously-varying specific impulse','5: custom thrust and mass flow rate polynomial']

        self.cmbengine_type = wx.ComboBox(self, -1, choices = enginetypes, style=wx.CB_READONLY)

        self.lblnumber_of_electric_propulsion_systems = wx.StaticText(self, -1, "Number of thrusters")
        self.txtnumber_of_electric_propulsion_systems = wx.TextCtrl(self, -1, "number_of_electric_propulsion_systems")

        self.lblthrust_scale_factor = wx.StaticText(self, -1, "Thrust scale factor")
        self.txtthrust_scale_factor = wx.TextCtrl(self, -1, "thrust_scale_factor")

        self.lblThrust = wx.StaticText(self, -1, "Electric thruster thrust (N)")
        self.txtThrust = wx.TextCtrl(self, -1, "Thrust")

        self.lblIspLT = wx.StaticText(self, -1, "Electric thruster Isp (s)")
        self.txtIspLT = wx.TextCtrl(self, -1, "IspLT")

        self.lblIspLT_minimum = wx.StaticText(self, -1, "Minimum Isp for VSI systems (s)")
        self.txtIspLT_minimum = wx.TextCtrl(self, -1, "IspLT_minimum")

        self.lbluser_defined_engine_efficiency = wx.StaticText(self, -1, "Thruster efficiency")
        self.txtuser_defined_engine_efficiency = wx.TextCtrl(self, -1, "user_defined_engine_efficiency")

        self.lblengine_coefficient_spacer = wx.StaticText(self, -1, "")
        self.lblengine_coefficient0 = wx.StaticText(self, -1, "1.0")
        self.lblengine_coefficient1 = wx.StaticText(self, -1, "P")
        self.lblengine_coefficient2 = wx.StaticText(self, -1, "P^2")
        self.lblengine_coefficient3 = wx.StaticText(self, -1, "P^3")
        self.lblengine_coefficient4 = wx.StaticText(self, -1, "P^4")
        self.lblengine_coefficient5 = wx.StaticText(self, -1, "P^5")
        self.lblengine_coefficient6 = wx.StaticText(self, -1, "P^6")

        self.lblengine_input_thrust_coefficients = wx.StaticText(self, -1, "Custom thrust coefficients (mN)")
        self.txtengine_input_thrust_coefficients0 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[0]")
        self.txtengine_input_thrust_coefficients1 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[1]")
        self.txtengine_input_thrust_coefficients2 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[2]")
        self.txtengine_input_thrust_coefficients3 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[3]")
        self.txtengine_input_thrust_coefficients4 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[4]")
        self.txtengine_input_thrust_coefficients5 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[5]")
        self.txtengine_input_thrust_coefficients6 = wx.TextCtrl(self, -1, "engine_input_thrust_coefficients[6]")

        self.lblengine_input_mass_flow_rate_coefficients = wx.StaticText(self, -1, "Custom mass flow rate coefficients (mg/s)")
        self.txtengine_input_mass_flow_rate_coefficients0 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[0]")
        self.txtengine_input_mass_flow_rate_coefficients1 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[1]")
        self.txtengine_input_mass_flow_rate_coefficients2 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[2]")
        self.txtengine_input_mass_flow_rate_coefficients3 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[3]")
        self.txtengine_input_mass_flow_rate_coefficients4 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[4]")
        self.txtengine_input_mass_flow_rate_coefficients5 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[5]")
        self.txtengine_input_mass_flow_rate_coefficients6 = wx.TextCtrl(self, -1, "engine_input_mass_flow_rate_coefficients[6]")

        self.lblengine_input_power_bounds = wx.StaticText(self, -1, "Thruster input power bounds")
        self.txtengine_input_power_bounds_lower = wx.TextCtrl(self, -1, "engine_input_power_bounds[0]")
        self.txtengine_input_power_bounds_upper = wx.TextCtrl(self, -1, "engine_input_power_bounds[1]")
        enginepowerbox = wx.BoxSizer(wx.HORIZONTAL)
        enginepowerbox.AddMany([self.txtengine_input_power_bounds_lower, self.txtengine_input_power_bounds_upper])

        self.lblthrottle_table = wx.StaticText(self, -1, "Throttle table file")
        self.txtthrottle_table = wx.TextCtrl(self, -1, "throttle_table", size=(500,-1))
        self.btnthrottle_table_default = wx.Button(self, -1, "Default")
        throttle_table_box = wx.BoxSizer(wx.HORIZONTAL)
        throttle_table_box.AddMany([self.txtthrottle_table, self.btnthrottle_table_default])

        customthrustergrid = wx.FlexGridSizer(3, 8, 5, 5)
        customthrustergrid.AddMany([self.lblengine_coefficient_spacer, self.lblengine_coefficient0, self.lblengine_coefficient1, self.lblengine_coefficient2, self.lblengine_coefficient3, self.lblengine_coefficient4, self.lblengine_coefficient5, self.lblengine_coefficient6,
                                    self.lblengine_input_thrust_coefficients, self.txtengine_input_thrust_coefficients0, self.txtengine_input_thrust_coefficients1, self.txtengine_input_thrust_coefficients2, self.txtengine_input_thrust_coefficients3, self.txtengine_input_thrust_coefficients4, self.txtengine_input_thrust_coefficients5, self.txtengine_input_thrust_coefficients6,
                                    self.lblengine_input_mass_flow_rate_coefficients, self.txtengine_input_mass_flow_rate_coefficients0, self.txtengine_input_mass_flow_rate_coefficients1, self.txtengine_input_mass_flow_rate_coefficients2, self.txtengine_input_mass_flow_rate_coefficients3, self.txtengine_input_mass_flow_rate_coefficients4, self.txtengine_input_mass_flow_rate_coefficients5, self.txtengine_input_mass_flow_rate_coefficients6])
        


        propulsiongrid.AddMany([self.lblIspChem, self.txtIspChem,
                                self.lblTCM_Isp, self.txtTCM_Isp,
                                self.lblengine_type, self.cmbengine_type,
                                self.lblnumber_of_electric_propulsion_systems, self.txtnumber_of_electric_propulsion_systems,
                                self.lblthrust_scale_factor, self.txtthrust_scale_factor,
                                self.lblThrust, self.txtThrust,
                                self.lblIspLT, self.txtIspLT,
                                self.lblIspLT_minimum, self.txtIspLT_minimum,
                                self.lbluser_defined_engine_efficiency, self.txtuser_defined_engine_efficiency,
                                self.lblengine_input_power_bounds, enginepowerbox,
                                self.lblthrottle_table, throttle_table_box])



        propulsionbox = wx.BoxSizer(wx.VERTICAL)
        propulsionbox.AddMany([self.propulsiongridtitle, propulsiongrid, customthrustergrid])

        #power
        self.powergrid = wx.FlexGridSizer(20,2,5,5)
        self.powergridtitle = wx.StaticText(self, -1, "Power options")

        self.lblpower_at_1_AU = wx.StaticText(self, -1, "Power at 1 AU (kW)")
        self.txtpower_at_1_AU = wx.TextCtrl(self, -1, "power_at_1_AU")

        self.lblpower_source_type = wx.StaticText(self, -1, "Power source type")
        power_source_choices = ['0: fixed','1: solar','2: proprietary (if compiled)']
        self.cmbpower_source_type = wx.ComboBox(self, -1, choices=power_source_choices, style=wx.CB_READONLY)

        self.lblsolar_power_model_type = wx.StaticText(self, -1, "Solar power model type")
        solar_power_model_type_choices = ['0: classic Sauer','1: polynomial']
        self.cmbsolar_power_model_type = wx.ComboBox(self, -1, choices=solar_power_model_type_choices, style=wx.CB_READONLY)

        self.lblsolar_power_gamma = wx.StaticText(self, -1, "Solar power coefficients")
        self.txtsolar_power_gamma0 = wx.TextCtrl(self, -1, "solar_power_gamma[0]")
        self.txtsolar_power_gamma1 = wx.TextCtrl(self, -1, "solar_power_gamma[1]")
        self.txtsolar_power_gamma2 = wx.TextCtrl(self, -1, "solar_power_gamma[2]")
        self.txtsolar_power_gamma3 = wx.TextCtrl(self, -1, "solar_power_gamma[3]")
        self.txtsolar_power_gamma4 = wx.TextCtrl(self, -1, "solar_power_gamma[4]")
        self.txtsolar_power_gamma5 = wx.TextCtrl(self, -1, "solar_power_gamma[5]")
        self.txtsolar_power_gamma6 = wx.TextCtrl(self, -1, "solar_power_gamma[6]")
        solarpowerbox = wx.BoxSizer(wx.HORIZONTAL)
        solarpowerbox.AddMany([self.txtsolar_power_gamma0, self.txtsolar_power_gamma1, self.txtsolar_power_gamma2, self.txtsolar_power_gamma3, self.txtsolar_power_gamma4, self.txtsolar_power_gamma5, self.txtsolar_power_gamma6])

        self.lblspacecraft_power_model_type = wx.StaticText(self, -1, "Spacecraft power model type")
        power_model_choices = ['0: P_sc = A + B/r + C/r^2','1: P_sc = A if P > A, A + B(C - P) otherwise']
        self.cmbspacecraft_power_model_type = wx.ComboBox(self, -1, choices=power_model_choices, style = wx.CB_READONLY)
        
        self.lblspacecraft_power_coefficients = wx.StaticText(self, -1, "Spacecraft power coefficients (kW)")
        self.txtspacecraft_power_coefficients0 = wx.TextCtrl(self, -1, "spacecraft_power_coefficients[0]")
        self.txtspacecraft_power_coefficients1 = wx.TextCtrl(self, -1, "spacecraft_power_coefficients[1]")
        self.txtspacecraft_power_coefficients2 = wx.TextCtrl(self, -1, "spacecraft_power_coefficients[2]")
        spacecraftpowerbox = wx.BoxSizer(wx.HORIZONTAL)
        spacecraftpowerbox.AddMany([self.txtspacecraft_power_coefficients0, self.txtspacecraft_power_coefficients1, self.txtspacecraft_power_coefficients2])

        self.lblpower_decay_rate = wx.StaticText(self, -1, "Power decay rate (fraction per year)")
        self.txtpower_decay_rate = wx.TextCtrl(self, -1, "power_decay_rate")

        self.powergrid.AddMany([self.lblpower_at_1_AU, self.txtpower_at_1_AU,
                           self.lblpower_source_type, self.cmbpower_source_type,
                           self.lblsolar_power_model_type, self.cmbsolar_power_model_type,
                           self.lblsolar_power_gamma, solarpowerbox,
                           self.lblspacecraft_power_model_type, self.cmbspacecraft_power_model_type,
                           self.lblspacecraft_power_coefficients, spacecraftpowerbox,
                           self.lblpower_decay_rate, self.txtpower_decay_rate])

        powerbox = wx.BoxSizer(wx.VERTICAL)
        powerbox.AddMany([self.powergridtitle, self.powergrid])

        

        self.TankPanel = TankPanel.TankPanel(self, self.missionoptions)
        tanksbox = wx.BoxSizer(wx.VERTICAL)
        tanksbox.Add(self.TankPanel)

        #now tie everything together
        self.mainbox = wx.BoxSizer(wx.VERTICAL)
        self.mainbox.AddMany([LaunchVehicleBox, propulsionbox, powerbox, tanksbox])

        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.LaunchVehicleGridTitle.SetFont(font)
        self.powergridtitle.SetFont(font)
        self.propulsiongridtitle.SetFont(font)

        self.SetSizer(self.mainbox)

        #bindings
        #regular spacecraft tab
        
        self.txtHardwarePath.Bind(wx.EVT_KILL_FOCUS,self.ChangeHardwarePath)
        self.btnGetHardwarePath.Bind(wx.EVT_BUTTON,self.GetHardwarePath)
        self.txtLaunchVehicleLibraryFile.Bind(wx.EVT_KILL_FOCUS,self.ChangeLaunchVehicleLibraryFile)
        self.btnGetLaunchVehicleLibraryFile.Bind(wx.EVT_BUTTON,self.GetLaunchVehicleLibraryFile)
        self.cmbLaunchVehicleKey.Bind(wx.EVT_COMBOBOX, self.ChangeLaunchVehicleKey)

        self.txtIspChem.Bind(wx.EVT_KILL_FOCUS, self.ChangeIspChem)
        self.txtTCM_Isp.Bind(wx.EVT_KILL_FOCUS, self.ChangeTCM_Isp)
        self.cmbengine_type.Bind(wx.EVT_COMBOBOX, self.Changeengine_type)
        self.txtnumber_of_electric_propulsion_systems.Bind(wx.EVT_KILL_FOCUS, self.Changenumber_of_electric_propulsion_systems)
        self.txtthrust_scale_factor.Bind(wx.EVT_KILL_FOCUS, self.Changethrust_scale_factor)
        self.txtThrust.Bind(wx.EVT_KILL_FOCUS, self.ChangeThrust)
        self.txtIspLT.Bind(wx.EVT_KILL_FOCUS, self.ChangeIspLT)
        self.txtIspLT_minimum.Bind(wx.EVT_KILL_FOCUS, self.ChangeIspLT_minimum)
        self.txtuser_defined_engine_efficiency.Bind(wx.EVT_KILL_FOCUS, self.Changeuser_defined_engine_efficiency)
        self.txtengine_input_thrust_coefficients0.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients0)
        self.txtengine_input_thrust_coefficients1.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients1)
        self.txtengine_input_thrust_coefficients2.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients2)
        self.txtengine_input_thrust_coefficients3.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients3)
        self.txtengine_input_thrust_coefficients4.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients4)
        self.txtengine_input_thrust_coefficients5.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients5)
        self.txtengine_input_thrust_coefficients6.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_thrust_coefficients6)
        self.txtengine_input_mass_flow_rate_coefficients0.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients0)
        self.txtengine_input_mass_flow_rate_coefficients1.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients1)
        self.txtengine_input_mass_flow_rate_coefficients2.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients2)
        self.txtengine_input_mass_flow_rate_coefficients3.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients3)
        self.txtengine_input_mass_flow_rate_coefficients4.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients4)
        self.txtengine_input_mass_flow_rate_coefficients5.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients5)
        self.txtengine_input_mass_flow_rate_coefficients6.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_mass_flow_rate_coefficients6)
        self.txtengine_input_power_bounds_lower.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_power_bounds_lower)
        self.txtengine_input_power_bounds_upper.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_input_power_bounds_upper)
        self.txtthrottle_table.Bind(wx.EVT_KILL_FOCUS, self.Changethrottle_table)
        self.btnthrottle_table_default.Bind(wx.EVT_BUTTON, self.ClickThrottleTableDefaultButton)
        self.txtpower_at_1_AU.Bind(wx.EVT_KILL_FOCUS, self.Changepower_at_1_AU)
        self.cmbpower_source_type.Bind(wx.EVT_COMBOBOX, self.Changepower_source_type)
        self.cmbsolar_power_model_type.Bind(wx.EVT_COMBOBOX, self.Changesolar_power_model_type)
        self.txtsolar_power_gamma0.Bind(wx.EVT_KILL_FOCUS, self.Changesolar_power_gamma0)
        self.txtsolar_power_gamma1.Bind(wx.EVT_KILL_FOCUS, self.Changesolar_power_gamma1)
        self.txtsolar_power_gamma2.Bind(wx.EVT_KILL_FOCUS, self.Changesolar_power_gamma2)
        self.txtsolar_power_gamma3.Bind(wx.EVT_KILL_FOCUS, self.Changesolar_power_gamma3)
        self.txtsolar_power_gamma4.Bind(wx.EVT_KILL_FOCUS, self.Changesolar_power_gamma4)
        self.txtsolar_power_gamma5.Bind(wx.EVT_KILL_FOCUS, self.Changesolar_power_gamma5)
        self.txtsolar_power_gamma6.Bind(wx.EVT_KILL_FOCUS, self.Changesolar_power_gamma6)
        self.cmbspacecraft_power_model_type.Bind(wx.EVT_COMBOBOX, self.Changespacecraft_power_model_type)
        self.txtspacecraft_power_coefficients0.Bind(wx.EVT_KILL_FOCUS, self.Changespacecraft_power_coefficients0)
        self.txtspacecraft_power_coefficients1.Bind(wx.EVT_KILL_FOCUS, self.Changespacecraft_power_coefficients1)
        self.txtspacecraft_power_coefficients2.Bind(wx.EVT_KILL_FOCUS, self.Changespacecraft_power_coefficients2)
        self.txtpower_decay_rate.Bind(wx.EVT_KILL_FOCUS, self.Changepower_decay_rate)

    def update(self, updateTanks=True):
        import os
        import re

        if updateTanks:
            self.TankPanel.update()
        
        self.txtHardwarePath.SetValue(self.missionoptions.HardwarePath)
        self.txtLaunchVehicleLibraryFile.SetValue(self.missionoptions.LaunchVehicleLibraryFile)

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

        self.txtIspChem.SetValue(str(self.missionoptions.IspChem))
        self.txtTCM_Isp.SetValue(str(self.missionoptions.TCM_Isp))
        self.cmbengine_type.SetSelection(self.missionoptions.engine_type)
        self.txtnumber_of_electric_propulsion_systems.SetValue(str(self.missionoptions.number_of_electric_propulsion_systems))
        self.txtthrust_scale_factor.SetValue(str(self.missionoptions.thrust_scale_factor))
        self.txtThrust.SetValue(str(self.missionoptions.Thrust))
        self.txtIspLT.SetValue(str(self.missionoptions.IspLT))
        self.txtIspLT_minimum.SetValue(str(self.missionoptions.IspLT_minimum))
        self.txtuser_defined_engine_efficiency.SetValue(str(self.missionoptions.user_defined_engine_efficiency))
        self.txtengine_input_thrust_coefficients0.SetValue(str(self.missionoptions.engine_input_thrust_coefficients[0]))
        self.txtengine_input_thrust_coefficients1.SetValue(str(self.missionoptions.engine_input_thrust_coefficients[1]))
        self.txtengine_input_thrust_coefficients2.SetValue(str(self.missionoptions.engine_input_thrust_coefficients[2]))
        self.txtengine_input_thrust_coefficients3.SetValue(str(self.missionoptions.engine_input_thrust_coefficients[3]))
        self.txtengine_input_thrust_coefficients4.SetValue(str(self.missionoptions.engine_input_thrust_coefficients[4]))
        self.txtengine_input_thrust_coefficients5.SetValue(str(self.missionoptions.engine_input_thrust_coefficients[5]))
        self.txtengine_input_thrust_coefficients6.SetValue(str(self.missionoptions.engine_input_thrust_coefficients[6]))
        self.txtengine_input_mass_flow_rate_coefficients0.SetValue(str(self.missionoptions.engine_input_mass_flow_rate_coefficients[0]))
        self.txtengine_input_mass_flow_rate_coefficients1.SetValue(str(self.missionoptions.engine_input_mass_flow_rate_coefficients[1]))
        self.txtengine_input_mass_flow_rate_coefficients2.SetValue(str(self.missionoptions.engine_input_mass_flow_rate_coefficients[2]))
        self.txtengine_input_mass_flow_rate_coefficients3.SetValue(str(self.missionoptions.engine_input_mass_flow_rate_coefficients[3]))
        self.txtengine_input_mass_flow_rate_coefficients4.SetValue(str(self.missionoptions.engine_input_mass_flow_rate_coefficients[4]))
        self.txtengine_input_mass_flow_rate_coefficients5.SetValue(str(self.missionoptions.engine_input_mass_flow_rate_coefficients[5]))
        self.txtengine_input_mass_flow_rate_coefficients6.SetValue(str(self.missionoptions.engine_input_mass_flow_rate_coefficients[6]))
        self.txtengine_input_power_bounds_lower.SetValue(str(self.missionoptions.engine_input_power_bounds[0]))
        self.txtengine_input_power_bounds_upper.SetValue(str(self.missionoptions.engine_input_power_bounds[1]))
        self.txtthrottle_table.SetValue(self.missionoptions.ThrottleTableFile)
        self.txtpower_at_1_AU.SetValue(str(self.missionoptions.power_at_1_AU))
        self.cmbpower_source_type.SetSelection(self.missionoptions.power_source_type)
        self.cmbsolar_power_model_type.SetSelection(self.missionoptions.solar_power_model_type)
        self.txtsolar_power_gamma0.SetValue(str(self.missionoptions.solar_power_gamma[0]))
        self.txtsolar_power_gamma1.SetValue(str(self.missionoptions.solar_power_gamma[1]))
        self.txtsolar_power_gamma2.SetValue(str(self.missionoptions.solar_power_gamma[2]))
        self.txtsolar_power_gamma3.SetValue(str(self.missionoptions.solar_power_gamma[3]))
        self.txtsolar_power_gamma4.SetValue(str(self.missionoptions.solar_power_gamma[4]))
        self.txtsolar_power_gamma5.SetValue(str(self.missionoptions.solar_power_gamma[5]))
        self.txtsolar_power_gamma6.SetValue(str(self.missionoptions.solar_power_gamma[6]))
        self.cmbspacecraft_power_model_type.SetSelection(self.missionoptions.spacecraft_power_model_type)
        self.txtspacecraft_power_coefficients0.SetValue(str(self.missionoptions.spacecraft_power_coefficients[0]))
        self.txtspacecraft_power_coefficients1.SetValue(str(self.missionoptions.spacecraft_power_coefficients[1]))
        self.txtspacecraft_power_coefficients2.SetValue(str(self.missionoptions.spacecraft_power_coefficients[2]))
        self.txtpower_decay_rate.SetValue(str(self.missionoptions.power_decay_rate))

        #switch various spacecraft fields visible and invisible depending on the mission type

        #impulsive vs low-thrust missions
        if self.missionoptions.mission_type in [6, 9]:
            #impulsive mission
            self.powergridtitle.Show(False)
            self.lblIspChem.Show(True)
            self.lblengine_type.Show(False)
            self.lblnumber_of_electric_propulsion_systems.Show(False)
            self.lblthrust_scale_factor.Show(False)
            self.lblThrust.Show(False)
            self.lblIspLT.Show(False)
            self.lblIspLT_minimum.Show(False)
            self.lbluser_defined_engine_efficiency.Show(False)
            self.lblengine_input_thrust_coefficients.Show(False)
            self.lblengine_input_mass_flow_rate_coefficients.Show(False)
            self.lblengine_input_power_bounds.Show(False)
            self.lblpower_at_1_AU.Show(False)
            self.lblpower_source_type.Show(False)
            self.lblsolar_power_gamma.Show(False)
            self.lblspacecraft_power_model_type.Show(False)
            self.lblspacecraft_power_coefficients.Show(False)
            self.lblpower_decay_rate.Show(False)
            self.lblengine_coefficient0.Show(False)
            self.lblengine_coefficient1.Show(False)
            self.lblengine_coefficient2.Show(False)
            self.lblengine_coefficient3.Show(False)
            self.lblengine_coefficient4.Show(False)
            self.lblengine_coefficient5.Show(False)
            self.lblengine_coefficient6.Show(False)
            self.txtIspChem.Show(True)
            self.cmbengine_type.Show(False)
            self.txtnumber_of_electric_propulsion_systems.Show(False)
            self.txtthrust_scale_factor.Show(False)
            self.txtThrust.Show(False)
            self.txtIspLT.Show(False)
            self.txtIspLT_minimum.Show(False)
            self.txtuser_defined_engine_efficiency.Show(False)
            self.txtengine_input_thrust_coefficients0.Show(False)
            self.txtengine_input_thrust_coefficients1.Show(False)
            self.txtengine_input_thrust_coefficients2.Show(False)
            self.txtengine_input_thrust_coefficients3.Show(False)
            self.txtengine_input_thrust_coefficients4.Show(False)
            self.txtengine_input_thrust_coefficients5.Show(False)
            self.txtengine_input_thrust_coefficients6.Show(False)
            self.txtengine_input_mass_flow_rate_coefficients0.Show(False)
            self.txtengine_input_mass_flow_rate_coefficients1.Show(False)
            self.txtengine_input_mass_flow_rate_coefficients2.Show(False)
            self.txtengine_input_mass_flow_rate_coefficients3.Show(False)
            self.txtengine_input_mass_flow_rate_coefficients4.Show(False)
            self.txtengine_input_mass_flow_rate_coefficients5.Show(False)
            self.txtengine_input_mass_flow_rate_coefficients6.Show(False)
            self.txtengine_input_power_bounds_lower.Show(False)
            self.txtengine_input_power_bounds_upper.Show(False)
            self.txtpower_at_1_AU.Show(False)
            self.cmbpower_source_type.Show(False)
            self.lblsolar_power_model_type.Show(False)
            self.cmbsolar_power_model_type.Show(False)
            self.txtsolar_power_gamma0.Show(False)
            self.txtsolar_power_gamma1.Show(False)
            self.txtsolar_power_gamma2.Show(False)
            self.txtsolar_power_gamma3.Show(False)
            self.txtsolar_power_gamma4.Show(False)
            self.txtsolar_power_gamma5.Show(False)
            self.txtsolar_power_gamma6.Show(False)
            self.cmbspacecraft_power_model_type.Show(False)
            self.txtspacecraft_power_coefficients0.Show(False)
            self.txtspacecraft_power_coefficients1.Show(False)
            self.txtspacecraft_power_coefficients2.Show(False)
            self.txtpower_decay_rate.Show(False)
            self.lblthrottle_table.Show(False)
            self.txtthrottle_table.Show(False)
            self.btnthrottle_table_default.Show(False)
        else:
            #low-thrust missions
            self.lblengine_type.Show(True)
            self.cmbengine_type.Show(True)

            if self.missionoptions.engine_type == 0:
                #fixed thrust/Isp, no power information required
                self.powergridtitle.Show(False)
                self.lblIspChem.Show(False)
                self.lblnumber_of_electric_propulsion_systems.Show(False)
                self.lblthrust_scale_factor.Show(True)
                self.lblThrust.Show(True)
                self.lblIspLT.Show(True)
                self.lblIspLT_minimum.Show(False)
                self.lbluser_defined_engine_efficiency.Show(False)
                self.lblengine_input_thrust_coefficients.Show(False)
                self.lblengine_input_mass_flow_rate_coefficients.Show(False)
                self.lblengine_input_power_bounds.Show(False)
                self.lblpower_at_1_AU.Show(False)
                self.lblpower_source_type.Show(False)
                self.lblsolar_power_gamma.Show(False)
                self.lblspacecraft_power_model_type.Show(False)
                self.lblspacecraft_power_coefficients.Show(False)
                self.lblpower_decay_rate.Show(False)
                self.lblengine_coefficient0.Show(False)
                self.lblengine_coefficient1.Show(False)
                self.lblengine_coefficient2.Show(False)
                self.lblengine_coefficient3.Show(False)
                self.lblengine_coefficient4.Show(False)
                self.lblengine_coefficient5.Show(False)
                self.lblengine_coefficient6.Show(False)
                self.txtIspChem.Show(False)
                self.txtnumber_of_electric_propulsion_systems.Show(False)
                self.txtthrust_scale_factor.Show(True)
                self.txtThrust.Show(True)
                self.txtIspLT.Show(True)
                self.txtIspLT_minimum.Show(False)
                self.txtuser_defined_engine_efficiency.Show(False)
                self.txtengine_input_thrust_coefficients0.Show(False)
                self.txtengine_input_thrust_coefficients1.Show(False)
                self.txtengine_input_thrust_coefficients2.Show(False)
                self.txtengine_input_thrust_coefficients3.Show(False)
                self.txtengine_input_thrust_coefficients4.Show(False)
                self.txtengine_input_thrust_coefficients5.Show(False)
                self.txtengine_input_thrust_coefficients6.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients0.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients1.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients2.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients3.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients4.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients5.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients6.Show(False)
                self.txtengine_input_power_bounds_lower.Show(False)
                self.txtengine_input_power_bounds_upper.Show(False)
                self.txtpower_at_1_AU.Show(False)
                self.cmbpower_source_type.Show(False)
                self.txtsolar_power_gamma0.Show(False)
                self.txtsolar_power_gamma1.Show(False)
                self.txtsolar_power_gamma2.Show(False)
                self.txtsolar_power_gamma3.Show(False)
                self.txtsolar_power_gamma4.Show(False)
                self.txtsolar_power_gamma5.Show(False)
                self.txtsolar_power_gamma6.Show(False)
                self.cmbspacecraft_power_model_type.Show(False)
                self.txtspacecraft_power_coefficients0.Show(False)
                self.txtspacecraft_power_coefficients1.Show(False)
                self.txtspacecraft_power_coefficients2.Show(False)
                self.txtpower_decay_rate.Show(False)
                self.lblthrottle_table.Show(False)
                self.txtthrottle_table.Show(False)
                self.btnthrottle_table_default.Show(False)
            elif self.missionoptions.engine_type == 1:
                #constant Isp, efficiency, EMTG computes input power
                #do not need anything except Isp, efficiency
                self.powergridtitle.Show(False)
                self.lblIspChem.Show(False)
                self.lblnumber_of_electric_propulsion_systems.Show(False)
                self.lblthrust_scale_factor.Show(True)
                self.lblThrust.Show(False)
                self.lblIspLT.Show(True)
                self.lblIspLT_minimum.Show(False)
                self.lbluser_defined_engine_efficiency.Show(True)
                self.lblengine_input_thrust_coefficients.Show(False)
                self.lblengine_input_mass_flow_rate_coefficients.Show(False)
                self.lblengine_input_power_bounds.Show(False)
                self.lblpower_at_1_AU.Show(False)
                self.lblpower_source_type.Show(False)
                self.lblsolar_power_gamma.Show(False)
                self.lblspacecraft_power_model_type.Show(False)
                self.lblspacecraft_power_coefficients.Show(False)
                self.lblpower_decay_rate.Show(False)
                self.lblengine_coefficient0.Show(False)
                self.lblengine_coefficient1.Show(False)
                self.lblengine_coefficient2.Show(False)
                self.lblengine_coefficient3.Show(False)
                self.lblengine_coefficient4.Show(False)
                self.lblengine_coefficient5.Show(False)
                self.lblengine_coefficient6.Show(False)
                self.txtIspChem.Show(False)
                self.txtnumber_of_electric_propulsion_systems.Show(False)
                self.txtthrust_scale_factor.Show(True)
                self.txtThrust.Show(False)
                self.txtIspLT.Show(True)
                self.txtIspLT_minimum.Show(False)
                self.txtuser_defined_engine_efficiency.Show(True)
                self.txtengine_input_thrust_coefficients0.Show(False)
                self.txtengine_input_thrust_coefficients1.Show(False)
                self.txtengine_input_thrust_coefficients2.Show(False)
                self.txtengine_input_thrust_coefficients3.Show(False)
                self.txtengine_input_thrust_coefficients4.Show(False)
                self.txtengine_input_thrust_coefficients5.Show(False)
                self.txtengine_input_thrust_coefficients6.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients0.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients1.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients2.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients3.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients4.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients5.Show(False)
                self.txtengine_input_mass_flow_rate_coefficients6.Show(False)
                self.txtengine_input_power_bounds_lower.Show(False)
                self.txtengine_input_power_bounds_upper.Show(False)
                self.txtpower_at_1_AU.Show(False)
                self.cmbpower_source_type.Show(False)
                self.txtsolar_power_gamma0.Show(False)
                self.txtsolar_power_gamma1.Show(False)
                self.txtsolar_power_gamma2.Show(False)
                self.txtsolar_power_gamma3.Show(False)
                self.txtsolar_power_gamma4.Show(False)
                self.txtsolar_power_gamma5.Show(False)
                self.txtsolar_power_gamma6.Show(False)
                self.cmbspacecraft_power_model_type.Show(False)
                self.txtspacecraft_power_coefficients0.Show(False)
                self.txtspacecraft_power_coefficients1.Show(False)
                self.txtspacecraft_power_coefficients2.Show(False)
                self.txtpower_decay_rate.Show(False)
                self.lblthrottle_table.Show(False)
                self.txtthrottle_table.Show(False)
                self.btnthrottle_table_default.Show(False)
            else:
                #engine types 2 and greater all require power information but NOT chemical Isp information
                self.powergridtitle.Show(True)
                self.lblpower_at_1_AU.Show(True)
                self.lblpower_source_type.Show(True)
                self.lblspacecraft_power_model_type.Show(True)
                self.lblspacecraft_power_coefficients.Show(True)
                self.lblpower_decay_rate.Show(True)
                self.txtpower_at_1_AU.Show(True)
                self.cmbpower_source_type.Show(True)
                if self.missionoptions.power_source_type == 1:
                    self.lblpower_at_1_AU.SetLabel("Power at BOL, 1 AU (kW)")
                    self.lblsolar_power_model_type.Show(True)
                    self.cmbsolar_power_model_type.Show(True)
                    self.lblsolar_power_gamma.Show(True)
                    self.txtsolar_power_gamma0.Show(True)
                    self.txtsolar_power_gamma1.Show(True)
                    self.txtsolar_power_gamma2.Show(True)
                    self.txtsolar_power_gamma3.Show(True)
                    self.txtsolar_power_gamma4.Show(True)
                    self.txtsolar_power_gamma5.Show(True)
                    self.txtsolar_power_gamma6.Show(True)
                else:
                    self.lblpower_at_1_AU.SetLabel("Power at BOL (kW)")
                    self.lblsolar_power_model_type.Show(False)
                    self.cmbsolar_power_model_type.Show(False)
                    self.lblsolar_power_gamma.Show(False)
                    self.txtsolar_power_gamma0.Show(False)
                    self.txtsolar_power_gamma1.Show(False)
                    self.txtsolar_power_gamma2.Show(False)
                    self.txtsolar_power_gamma3.Show(False)
                    self.txtsolar_power_gamma4.Show(False)
                    self.txtsolar_power_gamma5.Show(False)
                    self.txtsolar_power_gamma6.Show(False)

                self.cmbspacecraft_power_model_type.Show(True)
                self.txtspacecraft_power_coefficients0.Show(True)
                self.txtspacecraft_power_coefficients1.Show(True)
                self.txtspacecraft_power_coefficients2.Show(True)
                self.txtpower_decay_rate.Show(True)

                if self.missionoptions.engine_type == 2:
                    #choice of power model, constant efficiency, EMTG chooses Isp
                    #all other options off
                    self.lblnumber_of_electric_propulsion_systems.Show(False)
                    self.lblthrust_scale_factor.Show(True)
                    self.lblThrust.Show(False)
                    self.lblIspLT.Show(True)
                    self.lblIspLT_minimum.Show(True)
                    self.lbluser_defined_engine_efficiency.Show(True)
                    self.lblengine_input_thrust_coefficients.Show(False)
                    self.lblengine_input_mass_flow_rate_coefficients.Show(False)
                    self.lblengine_input_power_bounds.Show(True)
                    self.txtnumber_of_electric_propulsion_systems.Show(False)
                    self.txtthrust_scale_factor.Show(True)
                    self.txtThrust.Show(False)
                    self.txtIspLT.Show(True)
                    self.txtIspLT_minimum.Show(True)
                    self.txtuser_defined_engine_efficiency.Show(True)
                    self.lblengine_coefficient_spacer.Show(False)
                    self.lblengine_coefficient0.Show(False)
                    self.lblengine_coefficient1.Show(False)
                    self.lblengine_coefficient2.Show(False)
                    self.lblengine_coefficient3.Show(False)
                    self.lblengine_coefficient4.Show(False)
                    self.lblengine_coefficient5.Show(False)
                    self.lblengine_coefficient6.Show(False)
                    self.txtengine_input_thrust_coefficients0.Show(False)
                    self.txtengine_input_thrust_coefficients1.Show(False)
                    self.txtengine_input_thrust_coefficients2.Show(False)
                    self.txtengine_input_thrust_coefficients3.Show(False)
                    self.txtengine_input_thrust_coefficients4.Show(False)
                    self.txtengine_input_thrust_coefficients5.Show(False)
                    self.txtengine_input_thrust_coefficients6.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients0.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients1.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients2.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients3.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients4.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients5.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients6.Show(False)
                    self.txtengine_input_power_bounds_lower.Show(True)
                    self.txtengine_input_power_bounds_upper.Show(True)
                    self.lblthrottle_table.Show(False)
                    self.txtthrottle_table.Show(False)
                    self.btnthrottle_table_default.Show(False)

                elif self.missionoptions.engine_type == 3:
                    #choice of power model, constant efficiency and Isp
                    #all other options off
                    self.lblnumber_of_electric_propulsion_systems.Show(False)
                    self.lblthrust_scale_factor.Show(True)
                    self.lblThrust.Show(False)
                    self.lblIspLT.Show(True)
                    self.lblIspLT_minimum.Show(False)
                    self.lbluser_defined_engine_efficiency.Show(True)
                    self.lblengine_input_thrust_coefficients.Show(False)
                    self.lblengine_input_mass_flow_rate_coefficients.Show(False)
                    self.lblengine_input_power_bounds.Show(True)
                    self.txtnumber_of_electric_propulsion_systems.Show(False)
                    self.txtthrust_scale_factor.Show(True)
                    self.txtThrust.Show(False)
                    self.txtIspLT.Show(True)
                    self.txtIspLT_minimum.Show(False)
                    self.txtuser_defined_engine_efficiency.Show(True)
                    self.lblengine_coefficient_spacer.Show(False)
                    self.lblengine_coefficient0.Show(False)
                    self.lblengine_coefficient1.Show(False)
                    self.lblengine_coefficient2.Show(False)
                    self.lblengine_coefficient3.Show(False)
                    self.lblengine_coefficient4.Show(False)
                    self.lblengine_coefficient5.Show(False)
                    self.lblengine_coefficient6.Show(False)
                    self.txtengine_input_thrust_coefficients0.Show(False)
                    self.txtengine_input_thrust_coefficients1.Show(False)
                    self.txtengine_input_thrust_coefficients2.Show(False)
                    self.txtengine_input_thrust_coefficients3.Show(False)
                    self.txtengine_input_thrust_coefficients4.Show(False)
                    self.txtengine_input_thrust_coefficients5.Show(False)
                    self.txtengine_input_thrust_coefficients6.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients0.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients1.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients2.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients3.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients4.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients5.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients6.Show(False)
                    self.txtengine_input_power_bounds_lower.Show(True)
                    self.txtengine_input_power_bounds_upper.Show(True)
                    self.lblthrottle_table.Show(False)
                    self.txtthrottle_table.Show(False)
                    self.btnthrottle_table_default.Show(False)
                
                elif self.missionoptions.engine_type == 4:
                    #continuously-varying specific impulse
                    #requires efficiency, Isp min and max
                    self.lblnumber_of_electric_propulsion_systems.Show(False)
                    self.lblthrust_scale_factor.Show(True)
                    self.lblThrust.Show(False)
                    self.lblIspLT.Show(True)
                    self.lblIspLT_minimum.Show(True)
                    self.lbluser_defined_engine_efficiency.Show(True)
                    self.lblengine_input_thrust_coefficients.Show(False)
                    self.lblengine_input_mass_flow_rate_coefficients.Show(False)
                    self.lblengine_input_power_bounds.Show(True)
                    self.txtnumber_of_electric_propulsion_systems.Show(False)
                    self.txtthrust_scale_factor.Show(True)
                    self.txtThrust.Show(False)
                    self.txtIspLT.Show(True)
                    self.txtIspLT_minimum.Show(True)
                    self.txtuser_defined_engine_efficiency.Show(True)
                    self.lblengine_coefficient_spacer.Show(False)
                    self.lblengine_coefficient0.Show(False)
                    self.lblengine_coefficient1.Show(False)
                    self.lblengine_coefficient2.Show(False)
                    self.lblengine_coefficient3.Show(False)
                    self.lblengine_coefficient4.Show(False)
                    self.lblengine_coefficient5.Show(False)
                    self.lblengine_coefficient6.Show(False)
                    self.txtengine_input_thrust_coefficients0.Show(False)
                    self.txtengine_input_thrust_coefficients1.Show(False)
                    self.txtengine_input_thrust_coefficients2.Show(False)
                    self.txtengine_input_thrust_coefficients3.Show(False)
                    self.txtengine_input_thrust_coefficients4.Show(False)
                    self.txtengine_input_thrust_coefficients5.Show(False)
                    self.txtengine_input_thrust_coefficients6.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients0.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients1.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients2.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients3.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients4.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients5.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients6.Show(False)
                    self.txtengine_input_power_bounds_lower.Show(True)
                    self.txtengine_input_power_bounds_upper.Show(True)
                    self.lblthrottle_table.Show(False)
                    self.txtthrottle_table.Show(False)
                    self.btnthrottle_table_default.Show(False)

                elif self.missionoptions.engine_type == 5:
                    #custom thrust and mass flow rate polynomial
                    self.lblnumber_of_electric_propulsion_systems.Show(True)
                    self.lblthrust_scale_factor.Show(True)
                    self.lblThrust.Show(False)
                    self.lblIspLT.Show(False)
                    self.lblIspLT_minimum.Show(False)
                    self.lbluser_defined_engine_efficiency.Show(False)
                    self.lblengine_input_thrust_coefficients.Show(True)
                    self.lblengine_input_mass_flow_rate_coefficients.Show(True)
                    self.lblengine_input_power_bounds.Show(True)
                    self.txtnumber_of_electric_propulsion_systems.Show(True)
                    self.txtthrust_scale_factor.Show(True)
                    self.txtThrust.Show(False)
                    self.txtIspLT.Show(False)
                    self.txtIspLT_minimum.Show(False)
                    self.txtuser_defined_engine_efficiency.Show(False)
                    self.lblengine_coefficient_spacer.Show(True)
                    self.lblengine_coefficient0.Show(True)
                    self.lblengine_coefficient1.Show(True)
                    self.lblengine_coefficient2.Show(True)
                    self.lblengine_coefficient3.Show(True)
                    self.lblengine_coefficient4.Show(True)
                    self.lblengine_coefficient5.Show(True)
                    self.lblengine_coefficient6.Show(True)
                    self.txtengine_input_thrust_coefficients0.Show(True)
                    self.txtengine_input_thrust_coefficients1.Show(True)
                    self.txtengine_input_thrust_coefficients2.Show(True)
                    self.txtengine_input_thrust_coefficients3.Show(True)
                    self.txtengine_input_thrust_coefficients4.Show(True)
                    self.txtengine_input_thrust_coefficients5.Show(True)
                    self.txtengine_input_thrust_coefficients6.Show(True)
                    self.txtengine_input_mass_flow_rate_coefficients0.Show(True)
                    self.txtengine_input_mass_flow_rate_coefficients1.Show(True)
                    self.txtengine_input_mass_flow_rate_coefficients2.Show(True)
                    self.txtengine_input_mass_flow_rate_coefficients3.Show(True)
                    self.txtengine_input_mass_flow_rate_coefficients4.Show(True)
                    self.txtengine_input_mass_flow_rate_coefficients5.Show(True)
                    self.txtengine_input_mass_flow_rate_coefficients6.Show(True)
                    self.txtengine_input_power_bounds_lower.Show(True)
                    self.txtengine_input_power_bounds_upper.Show(True)
                    self.lblthrottle_table.Show(False)
                    self.txtthrottle_table.Show(False)
                    self.btnthrottle_table_default.Show(False)
                else:
                    #hard-coded thrust model
                    self.lblnumber_of_electric_propulsion_systems.Show(True)
                    self.lblthrust_scale_factor.Show(True)
                    self.lblThrust.Show(False)
                    self.lblIspLT.Show(False)
                    self.lblIspLT_minimum.Show(False)
                    self.lbluser_defined_engine_efficiency.Show(False)
                    self.lblengine_input_thrust_coefficients.Show(False)
                    self.lblengine_input_mass_flow_rate_coefficients.Show(False)
                    self.lblengine_input_power_bounds.Show(False)
                    self.txtnumber_of_electric_propulsion_systems.Show(True)
                    self.txtthrust_scale_factor.Show(True)
                    self.txtThrust.Show(False)
                    self.txtIspLT.Show(False)
                    self.txtIspLT_minimum.Show(False)
                    self.txtuser_defined_engine_efficiency.Show(False)
                    self.lblengine_coefficient_spacer.Show(False)
                    self.lblengine_coefficient0.Show(False)
                    self.lblengine_coefficient1.Show(False)
                    self.lblengine_coefficient2.Show(False)
                    self.lblengine_coefficient3.Show(False)
                    self.lblengine_coefficient4.Show(False)
                    self.lblengine_coefficient5.Show(False)
                    self.lblengine_coefficient6.Show(False)
                    self.txtengine_input_thrust_coefficients0.Show(False)
                    self.txtengine_input_thrust_coefficients1.Show(False)
                    self.txtengine_input_thrust_coefficients2.Show(False)
                    self.txtengine_input_thrust_coefficients3.Show(False)
                    self.txtengine_input_thrust_coefficients4.Show(False)
                    self.txtengine_input_thrust_coefficients5.Show(False)
                    self.txtengine_input_thrust_coefficients6.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients0.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients1.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients2.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients3.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients4.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients5.Show(False)
                    self.txtengine_input_mass_flow_rate_coefficients6.Show(False)
                    self.txtengine_input_power_bounds_lower.Show(False)
                    self.txtengine_input_power_bounds_upper.Show(False)

                    #throttle table options
                    if self.missionoptions.engine_type in [29, 30, 31, 32]:
                        self.lblthrottle_table.Show(True)
                        self.txtthrottle_table.Show(True)
                        self.btnthrottle_table_default.Show(True)
                    else:
                        self.lblthrottle_table.Show(False)
                        self.txtthrottle_table.Show(False)
                        self.btnthrottle_table_default.Show(False)
            
            needChemIsp = False
            for Journey in self.missionoptions.Journeys:
                if Journey.arrival_type in [1, 4, 0]:
                    needChemIsp = True
                elif Journey.enable_periapse_burns:
                            needChemIsp = True
            self.txtIspChem.Show(needChemIsp)
            self.lblIspChem.Show(needChemIsp)


        #re-size the panel
        self.Layout()
        if platform.system() == 'Windows':
            self.parent.SetupScrolling(scrollToTop=False)

    #event handlers for spacecraft options
    

    def ClickSetDefaultLibraries(self, e):
        e.Skip()
        self.missionoptions.HardwarePath = self.parent.parent.Parent.default_HardwarePath
        self.missionoptions.LaunchVehicleLibraryFile = self.parent.parent.Parent.default_LaunchVehicleLibraryFile

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

    def ChangeLaunchVehicleKey(self, e):
        e.Skip()
        self.missionoptions.LaunchVehicleKey = self.cmbLaunchVehicleKey.GetStringSelection()

        
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


    def ChangeIspChem(self, e):
        e.Skip()
        self.missionoptions.IspChem = eval(self.txtIspChem.GetValue())
        self.update()

    def ChangeTCM_Isp(self, e):
        e.Skip()
        self.missionoptions.TCM_Isp = eval(self.txtTCM_Isp.GetValue())
        self.update()

    def Changeengine_type(self, e):
        self.missionoptions.engine_type = self.cmbengine_type.GetSelection()
        self.update()

    def Changenumber_of_electric_propulsion_systems(self, e):
        e.Skip()
        self.missionoptions.number_of_electric_propulsion_systems = eval(self.txtnumber_of_electric_propulsion_systems.GetValue())
        self.update()

    def Changethrust_scale_factor(self, e):
        e.Skip()
        self.missionoptions.thrust_scale_factor = eval(self.txtthrust_scale_factor.GetValue())
        self.update()
        

    def ChangeThrust(self, e):
        e.Skip()
        self.missionoptions.Thrust = eval(self.txtThrust.GetValue())
        self.update()

    def ChangeIspLT(self, e):
        e.Skip()
        self.missionoptions.IspLT = eval(self.txtIspLT.GetValue())
        self.update()

    def ChangeIspLT_minimum(self, e):
        e.Skip()
        self.missionoptions.IspLT_minimum = eval(self.txtIspLT_minimum.GetValue())
        self.update()

    def Changeuser_defined_engine_efficiency(self, e):
        e.Skip()
        self.missionoptions.user_defined_engine_efficiency = eval(self.txtuser_defined_engine_efficiency.GetValue())
        self.update()

    def Changeengine_input_thrust_coefficients0(self, e):
        e.Skip()
        self.missionoptions.engine_input_thrust_coefficients[0] = eval(self.txtengine_input_thrust_coefficients0.GetValue())
        self.update()

    def Changeengine_input_thrust_coefficients1(self, e):
        e.Skip()
        self.missionoptions.engine_input_thrust_coefficients[1] = eval(self.txtengine_input_thrust_coefficients1.GetValue())
        self.update()

    def Changeengine_input_thrust_coefficients2(self, e):
        e.Skip()
        self.missionoptions.engine_input_thrust_coefficients[2] = eval(self.txtengine_input_thrust_coefficients2.GetValue())
        self.update()

    def Changeengine_input_thrust_coefficients3(self, e):
        e.Skip()
        self.missionoptions.engine_input_thrust_coefficients[3] = eval(self.txtengine_input_thrust_coefficients3.GetValue())
        self.update()

    def Changeengine_input_thrust_coefficients4(self, e):
        e.Skip()
        self.missionoptions.engine_input_thrust_coefficients[4] = eval(self.txtengine_input_thrust_coefficients4.GetValue())
        self.update()

    def Changeengine_input_thrust_coefficients5(self, e):
        e.Skip()
        self.missionoptions.engine_input_thrust_coefficients[5] = eval(self.txtengine_input_thrust_coefficients5.GetValue())
        self.update()

    def Changeengine_input_thrust_coefficients6(self, e):
        e.Skip()
        self.missionoptions.engine_input_thrust_coefficients[6] = eval(self.txtengine_input_thrust_coefficients6.GetValue())
        self.update()

    def Changeengine_input_mass_flow_rate_coefficients0(self, e):
        e.Skip()
        self.missionoptions.engine_input_mass_flow_rate_coefficients[0] = eval(self.txtengine_input_mass_flow_rate_coefficients0.GetValue())
        self.update()

    def Changeengine_input_mass_flow_rate_coefficients1(self, e):
        e.Skip()
        self.missionoptions.engine_input_mass_flow_rate_coefficients[1] = eval(self.txtengine_input_mass_flow_rate_coefficients1.GetValue())
        self.update()

    def Changeengine_input_mass_flow_rate_coefficients2(self, e):
        e.Skip()
        self.missionoptions.engine_input_mass_flow_rate_coefficients[2] = eval(self.txtengine_input_mass_flow_rate_coefficients2.GetValue())
        self.update()

    def Changeengine_input_mass_flow_rate_coefficients3(self, e):
        e.Skip()
        self.missionoptions.engine_input_mass_flow_rate_coefficients[3] = eval(self.txtengine_input_mass_flow_rate_coefficients3.GetValue())
        self.update()

    def Changeengine_input_mass_flow_rate_coefficients4(self, e):
        e.Skip()
        self.missionoptions.engine_input_mass_flow_rate_coefficients[4] = eval(self.txtengine_input_mass_flow_rate_coefficients4.GetValue())
        self.update()

    def Changeengine_input_mass_flow_rate_coefficients5(self, e):
        e.Skip()
        self.missionoptions.engine_input_mass_flow_rate_coefficients[5] = eval(self.txtengine_input_mass_flow_rate_coefficients5.GetValue())
        self.update()

    def Changeengine_input_mass_flow_rate_coefficients6(self, e):
        e.Skip()
        self.missionoptions.engine_input_mass_flow_rate_coefficients[6] = eval(self.txtengine_input_mass_flow_rate_coefficients6.GetValue())
        self.update()

    def Changeengine_input_power_bounds_lower(self, e):
        e.Skip()
        self.missionoptions.engine_input_power_bounds[0] = eval(self.txtengine_input_power_bounds_lower.GetValue())
        self.update()

    def Changeengine_input_power_bounds_upper(self, e):
        e.Skip()
        self.missionoptions.engine_input_power_bounds[1] = eval(self.txtengine_input_power_bounds_upper.GetValue())
        self.update()

    def Changethrottle_table(self, e):
        e.Skip()
        self.missionoptions.throttle_table = self.txtthrottle_table.GetValue()
        self.update()

    def ClickThrottleTableDefaultButton(self, e):
        e.Skip()
        self.missionoptions.throttle_table = self.parent.parent.Parent.default_thruster_file
        self.update()

    def Changepower_at_1_AU(self, e):
        e.Skip()
        self.missionoptions.power_at_1_AU = eval(self.txtpower_at_1_AU.GetValue())
        self.update()

    def Changepower_source_type(self, e):
        e.Skip()
        self.missionoptions.power_source_type = self.cmbpower_source_type.GetSelection()
        self.update()

    def Changesolar_power_model_type(self, e):
        e.Skip()
        self.missionoptions.solar_power_model_type = self.cmbsolar_power_model_type.GetSelection()
        self.update()

    def Changesolar_power_gamma0(self, e):
        e.Skip()
        self.missionoptions.solar_power_gamma[0] = eval(self.txtsolar_power_gamma0.GetValue())
        self.update()

    def Changesolar_power_gamma1(self, e):
        e.Skip()
        self.missionoptions.solar_power_gamma[1] = eval(self.txtsolar_power_gamma1.GetValue())
        self.update()

    def Changesolar_power_gamma2(self, e):
        e.Skip()
        self.missionoptions.solar_power_gamma[2] = eval(self.txtsolar_power_gamma2.GetValue())
        self.update()

    def Changesolar_power_gamma3(self, e):
        e.Skip()
        self.missionoptions.solar_power_gamma[3] = eval(self.txtsolar_power_gamma3.GetValue())
        self.update()

    def Changesolar_power_gamma4(self, e):
        e.Skip()
        self.missionoptions.solar_power_gamma[4] = eval(self.txtsolar_power_gamma4.GetValue())
        self.update()

    def Changesolar_power_gamma5(self, e):
        e.Skip()
        self.missionoptions.solar_power_gamma[5] = eval(self.txtsolar_power_gamma5.GetValue())
        self.update()
        
    def Changesolar_power_gamma6(self, e):
        e.Skip()
        self.missionoptions.solar_power_gamma[6] = eval(self.txtsolar_power_gamma6.GetValue())
        self.update()

    def Changespacecraft_power_model_type(self, e):
        e.Skip()
        self.missionoptions.spacecraft_power_model_type = self.cmbspacecraft_power_model_type.GetSelection()
        self.update()

    def Changespacecraft_power_coefficients0(self, e):
        e.Skip()
        self.missionoptions.spacecraft_power_coefficients[0] = eval(self.txtspacecraft_power_coefficients0.GetValue())
        self.update()

    def Changespacecraft_power_coefficients1(self, e):
        e.Skip()
        self.missionoptions.spacecraft_power_coefficients[1] = eval(self.txtspacecraft_power_coefficients1.GetValue())
        self.update()

    def Changespacecraft_power_coefficients2(self, e):
        e.Skip()
        self.missionoptions.spacecraft_power_coefficients[2] = eval(self.txtspacecraft_power_coefficients2.GetValue())
        self.update()

    def Changepower_decay_rate(self, e):
        e.Skip()
        self.missionoptions.power_decay_rate = eval(self.txtpower_decay_rate.GetValue())
        self.update()