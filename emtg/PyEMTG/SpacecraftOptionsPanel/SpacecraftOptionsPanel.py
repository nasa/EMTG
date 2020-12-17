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
import wx.lib.scrolledpanel
import platform

import ConventionalHardwarePanel
import LibraryHardwarePanel

class SpacecraftOptionsPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, missionoptions):
        self.missionoptions = missionoptions
        self.parent = parent
        
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)
        
        #common fields
        topgrid = wx.FlexGridSizer(20,2,5,5)
        topgridtitle = wx.StaticText(self, -1, "Common hardware options")

        self.lblmaximum_mass = wx.StaticText(self, -1, "Maximum mass (kg)")
        self.txtmaximum_mass = wx.TextCtrl(self, -1, "maximum_mass")
        
        self.lblallow_initial_mass_to_vary = wx.StaticText(self, -1, "Allow initial mass to vary")
        self.chkallow_initial_mass_to_vary = wx.CheckBox(self, -1)

        self.lblSpacecraftModelInput = wx.StaticText(self, -1, "Spacecraft model type")
        SpacecraftModelInputChoices = ['0: Assemble from libraries', '1: Read .emtg_spacecraftoptions file', '2: Assemble from missionoptions object']
        self.cmbSpacecraftModelInput = wx.ComboBox(self, -1, choices=SpacecraftModelInputChoices, style=wx.CB_READONLY)

        topgrid.AddMany([self.lblmaximum_mass, self.txtmaximum_mass,
                         self.lblallow_initial_mass_to_vary, self.chkallow_initial_mass_to_vary,
                         self.lblSpacecraftModelInput, self.cmbSpacecraftModelInput])

        topbox = wx.BoxSizer(wx.VERTICAL)
        topbox.AddMany([topgridtitle, topgrid])

        self.ConventionalHardwarePanel = ConventionalHardwarePanel.ConventionalHardwarePanel(self, self.missionoptions)
        ConventionalHardwareBox = wx.BoxSizer(wx.VERTICAL)
        ConventionalHardwareBox.Add(self.ConventionalHardwarePanel)

        self.LibraryHardwarePanel = LibraryHardwarePanel.LibraryHardwarePanel(self, self.missionoptions)
        LibraryHardwareBox = wx.BoxSizer(wx.VERTICAL)
        LibraryHardwareBox.Add(self.LibraryHardwarePanel)

        #margining fields
        margingrid = wx.FlexGridSizer(20,2,5,5)
        margingridtitle = wx.StaticText(self, -1, "Margins")

        self.lblpower_margin = wx.StaticText(self, -1, "Power margin (fraction)")
        self.txtpower_margin = wx.TextCtrl(self, -1, "power_margin")

        self.lblLV_margin = wx.StaticText(self, -1, "Launch vehicle margin (fraction)")
        self.txtLV_margin = wx.TextCtrl(self, -1, "LV_margin")

        self.lblelectric_propellant_margin = wx.StaticText(self, -1, "Electric propulsion propellant margin")
        self.txtelectric_propellant_margin = wx.TextCtrl(self, -1, "txtelectric_propellant_margin")
        self.lblchemical_propellant_margin = wx.StaticText(self, -1, "Chemical propulsion propellant margin")
        self.txtchemical_propellant_margin = wx.TextCtrl(self, -1, "txtchemical_propellant_margin")

        self.lblengine_duty_cycle = wx.StaticText(self, -1, "Thruster duty cycle")
        self.txtengine_duty_cycle = wx.TextCtrl(self, -1, "engine_duty_cycle")

        self.lblduty_cycle_type = wx.StaticText(self, -1, "Duty cycle type")
        self.cmbduty_cycle_type = wx.ComboBox(self, -1, choices=['0: Averaged','1: Realistic'], style=wx.CB_READONLY)

        margingrid.AddMany([self.lblLV_margin, self.txtLV_margin,
                            self.lblpower_margin, self.txtpower_margin,
                            self.lblengine_duty_cycle, self.txtengine_duty_cycle,
                            self.lblduty_cycle_type, self.cmbduty_cycle_type,
                            self.lblelectric_propellant_margin, self.txtelectric_propellant_margin,
                            self.lblchemical_propellant_margin, self.txtchemical_propellant_margin])

        marginbox = wx.BoxSizer(wx.VERTICAL)
        marginbox.AddMany([margingridtitle, margingrid])

        
        #ACS
        ACSgrid = wx.FlexGridSizer(20,2,5,5)
        ACSgridtitle = wx.StaticText(self, -1, "ACS options")

        self.lbltrackACS = wx.StaticText(self, -1, "Track ACS propellant?")
        self.chktrackACS = wx.CheckBox(self, -1)
        self.lblACS_kg_per_day = wx.StaticText(self, -1, "ACS propellant use per day (kg)")
        self.txtACS_kg_per_day = wx.TextCtrl(self, -1, "ACS_kg_per_day")

        ACSgrid.AddMany([self.lbltrackACS, self.chktrackACS,
                         self.lblACS_kg_per_day, self.txtACS_kg_per_day])

        ACSbox = wx.BoxSizer(wx.VERTICAL)
        ACSbox.AddMany([ACSgridtitle, ACSgrid])

        #throttle grid
        ThrottleLogicGrid = wx.FlexGridSizer(20,2,5,5)
        ThrottleLogictitle = wx.StaticText(self, -1, "Throttle grid options")

        self.lblthrottle_logic_mode = wx.StaticText(self, -1, "Throttle logic mode")
        throttle_logic_types = ['maximum number of thrusters','minimum number of thrusters']
        self.cmbthrottle_logic_mode = wx.ComboBox(self, -1, choices = throttle_logic_types, style = wx.CB_READONLY)

        self.lblthrottle_sharpness = wx.StaticText(self, -1, "Throttle sharpness")
        self.txtthrottle_sharpness = wx.TextCtrl(self, -1, "throttle_sharpness")

        ThrottleLogicGrid.AddMany([self.lblthrottle_logic_mode, self.cmbthrottle_logic_mode,
                                   self.lblthrottle_sharpness, self.txtthrottle_sharpness])

        ThrottleLogicBox = wx.BoxSizer(wx.VERTICAL)
        ThrottleLogicBox.AddMany([ThrottleLogictitle, ThrottleLogicGrid])

        # Power Source Decay Reference Epoch
        Decaygrid = wx.FlexGridSizer(20,2,5,5)
        self.lblpower_system_decay_reference_epoch = wx.StaticText(self, -1, "Power Source Decay Reference Epoch")
        self.txtpower_system_decay_reference_epoch = wx.TextCtrl(self, -1, "power_system_decay_reference_epoch")
        self.powerdecayCalendar = wx.adv.CalendarCtrl(self, -1)

        Decaygrid.AddMany([self.lblpower_system_decay_reference_epoch, self.txtpower_system_decay_reference_epoch, self.powerdecayCalendar])

        Decaybox = wx.BoxSizer(wx.VERTICAL)
        Decaybox.AddMany([Decaygrid])

        #now tie everything together
        leftvertsizer = wx.BoxSizer(wx.VERTICAL)
        leftvertsizer.AddMany([topbox, ConventionalHardwareBox, LibraryHardwareBox])
        rightvertsizer = wx.BoxSizer(wx.VERTICAL)
        rightvertsizer.AddMany([marginbox, ACSbox, ThrottleLogicBox, Decaybox])
        
        self.mainbox = wx.BoxSizer(wx.HORIZONTAL)
        self.mainbox.AddMany([leftvertsizer, rightvertsizer]) 

        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        topgridtitle.SetFont(font)
        margingridtitle.SetFont(font)
        ACSgridtitle.SetFont(font)
        ThrottleLogictitle.SetFont(font)

        self.SetSizer(self.mainbox)
        self.SetupScrolling()

        #bindings
        
        #globals
        self.txtmaximum_mass.Bind(wx.EVT_KILL_FOCUS, self.Changemaximum_mass)
        self.chkallow_initial_mass_to_vary.Bind(wx.EVT_CHECKBOX, self.Changeallow_initial_mass_to_vary)
        self.cmbSpacecraftModelInput.Bind(wx.EVT_COMBOBOX, self.ChangeSpacecraftModelInput)

        #margins
        self.txtpower_margin.Bind(wx.EVT_KILL_FOCUS, self.Changepower_margin)
        self.txtLV_margin.Bind(wx.EVT_KILL_FOCUS, self.ChangeLV_margin)
        self.txtengine_duty_cycle.Bind(wx.EVT_KILL_FOCUS, self.Changeengine_duty_cycle)
        self.cmbduty_cycle_type.Bind(wx.EVT_COMBOBOX, self.Changeduty_cycle_type)
        self.txtchemical_propellant_margin.Bind(wx.EVT_KILL_FOCUS, self.Changechemical_propellant_margin)
        self.txtelectric_propellant_margin.Bind(wx.EVT_KILL_FOCUS, self.Changeelectric_propellant_margin)

        #ACS
        self.chktrackACS.Bind(wx.EVT_CHECKBOX, self.ChangetrackACS)
        self.txtACS_kg_per_day.Bind(wx.EVT_KILL_FOCUS, self.ChangeACS_kg_per_day)

        #throttle grid
        self.cmbthrottle_logic_mode.Bind(wx.EVT_COMBOBOX, self.Changethrottle_logic_mode)
        self.txtthrottle_sharpness.Bind(wx.EVT_KILL_FOCUS, self.Changethrottle_sharpness)

        # decay grid
        self.powerdecayCalendar.Bind(wx.adv.EVT_CALENDAR_SEL_CHANGED, self.ChangepowerdecayCalendar)
        self.txtpower_system_decay_reference_epoch.Bind(wx.EVT_KILL_FOCUS, self.Changepower_system_reference_epoch)

    def update(self):
        self.txtmaximum_mass.SetValue(str(self.missionoptions.maximum_mass))
        self.chkallow_initial_mass_to_vary.SetValue(self.missionoptions.allow_initial_mass_to_vary)
        self.cmbSpacecraftModelInput.SetSelection(self.missionoptions.SpacecraftModelInput)
        
        self.txtpower_margin.SetValue(str(self.missionoptions.power_margin))
        self.txtLV_margin.SetValue(str(self.missionoptions.LV_margin))        
        self.txtengine_duty_cycle.SetValue(str(self.missionoptions.engine_duty_cycle))       
        self.cmbduty_cycle_type.SetSelection(self.missionoptions.duty_cycle_type)
        self.txtelectric_propellant_margin.SetValue(str(self.missionoptions.electric_propellant_margin))
        self.txtchemical_propellant_margin.SetValue(str(self.missionoptions.chemical_propellant_margin))

        if self.missionoptions.SpacecraftModelInput == 2: #assemble from GUI
            self.ConventionalHardwarePanel.Show(True)
            self.LibraryHardwarePanel.Show(False)
            self.ConventionalHardwarePanel.update()

        else:
            self.ConventionalHardwarePanel.Show(False)
            self.LibraryHardwarePanel.Show(True)
            self.LibraryHardwarePanel.update()

        self.chktrackACS.SetValue(self.missionoptions.trackACS)
        self.txtACS_kg_per_day.SetValue(str(self.missionoptions.ACS_kg_per_day))

        if self.missionoptions.trackACS:
            self.lblACS_kg_per_day.Show(True)
            self.txtACS_kg_per_day.Show(True)
        else:
            self.lblACS_kg_per_day.Show(False)
            self.txtACS_kg_per_day.Show(False)
 
        self.cmbthrottle_logic_mode.SetSelection(self.missionoptions.throttle_logic_mode)
        self.txtthrottle_sharpness.SetValue(str(self.missionoptions.throttle_sharpness))

        self.txtpower_system_decay_reference_epoch.SetValue(str(self.missionoptions.power_system_decay_reference_epoch))
        DecayreferenceEpoch = wx.DateTime.FromJDN(self.missionoptions.power_system_decay_reference_epoch + 2400000.5)
        self.powerdecayCalendar.SetDate(DecayreferenceEpoch.MakeUTC())


        #re-size the panel
        self.Layout()
        if platform.system() == 'Windows':
            self.SetupScrolling(scrollToTop=False)

    #event handlers for spacecraft options
    def Changemaximum_mass(self, e):
        e.Skip()
        self.missionoptions.maximum_mass = eval(self.txtmaximum_mass.GetValue())
        self.update()

    def Changeallow_initial_mass_to_vary(self, e):
        e.Skip()
        self.missionoptions.allow_initial_mass_to_vary = int(self.chkallow_initial_mass_to_vary.GetValue())
        self.update()

    def ChangeSpacecraftModelInput(self, e):
        e.Skip()
        self.missionoptions.SpacecraftModelInput = self.cmbSpacecraftModelInput.GetSelection()
        self.update()

    def Changechemical_propellant_margin(self, e):
        e.Skip()
        self.missionoptions.chemical_propellant_margin = eval(self.txtchemical_propellant_margin.GetValue())
        self.update()

    def Changeelectric_propellant_margin(self, e):
        e.Skip()
        self.missionoptions.electric_propellant_margin = eval(self.txtelectric_propellant_margin.GetValue())
        self.update()

    def Changepower_margin(self, e):
        e.Skip()
        self.missionoptions.power_margin = eval(self.txtpower_margin.GetValue())
        self.update()

    def ChangeLV_margin(self, e):
        e.Skip()
        self.missionoptions.LV_margin = eval(self.txtLV_margin.GetValue())
        self.update()

    def Changeengine_duty_cycle(self, e):
        e.Skip()
        self.missionoptions.engine_duty_cycle = eval(self.txtengine_duty_cycle.GetValue())
        self.update()

    def Changeduty_cycle_type(self, e):
        e.Skip()
        self.missionoptions.duty_cycle_type = self.cmbduty_cycle_type.GetSelection()
        self.update()

    def ChangetrackACS(self, e):
        e.Skip()
        self.missionoptions.trackACS = self.chktrackACS.GetValue()
        self.update()

    def ChangeACS_kg_per_day(self, e):
        e.Skip()
        self.missionoptions.ACS_kg_per_day = eval(self.txtACS_kg_per_day.GetValue())
        self.update()

    def Changethrottle_logic_mode(self, e):
        self.missionoptions.throttle_logic_mode = self.cmbthrottle_logic_mode.GetSelection()
        self.update()

    def Changethrottle_sharpness(self, e):
        e.Skip()
        self.missionoptions.throttle_sharpness = eval(self.txtthrottle_sharpness.GetValue())
        self.update()

    def Changepower_system_reference_epoch(self, e):
        e.Skip()

        dateString = self.txtpower_system_decay_reference_epoch.GetValue()

        from timeUtilities import stringToJD

        self.missionoptions.power_system_decay_reference_epoch = stringToJD(dateString, self.missionoptions.universe_folder)

        self.update()

    def ChangepowerdecayCalendar(self, e):
        powerdecayReferenceEpoch = self.powerdecayCalendar.GetDate()
        powerdecayReferenceEpoch = powerdecayReferenceEpoch.FromTimezone(wx.DateTime.TimeZone(offset=0))
        self.missionoptions.power_system_decay_reference_epoch = powerdecayReferenceEpoch.GetMJD()
        self.update()