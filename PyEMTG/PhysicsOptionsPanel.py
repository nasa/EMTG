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

class PhysicsOptionsPanel(wx.lib.scrolledpanel.ScrolledPanel):    
    def __init__(self, parent, missionoptions):
        self.missionoptions = missionoptions
        self.parent = parent

        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)

        ephemerisgrid = wx.FlexGridSizer(9,2,5,5)
        perturbgrid = wx.GridSizer(8,2,8,8)
        integratorgrid = wx.GridSizer(10,2,5,5)
        
        self.lblephemeris_source = wx.StaticText(self, -1, "Ephemeris Source")
        ephemeris_source_typestypes = ['Static','SPICE','SplineEphem']
        self.cmbephemeris_source = wx.ComboBox(self, -1, choices = ephemeris_source_typestypes, style=wx.CB_READONLY)

        self.lblSPICE_leap_seconds_kernel = wx.StaticText(self, -1, "Leap seconds kernel")
        self.txtSPICE_leap_seconds_kernel = wx.TextCtrl(self, -1, "SPICE_leap_seconds_kernel", size=(200,-1))

        self.lblSPICE_reference_frame_kernel = wx.StaticText(self, -1, "Frame kernel")
        self.txtSPICE_reference_frame_kernel = wx.TextCtrl(self, -1, "SPICE_reference_frame_kernel", size=(200,-1))

        self.lbluniverse_folder = wx.StaticText(self, -1, "Universe folder")
        self.txtuniverse_folder = wx.TextCtrl(self, -1, "universe_folder", size=(400,-1))
        self.btnGetNewUniverseFolder = wx.Button(self, -1, "...")
        self.btnSetDefaultUniverse = wx.Button(self, -1, "Default")
        UniverseButtonSizer = wx.BoxSizer(wx.HORIZONTAL)
        UniverseButtonSizer.AddMany([self.txtuniverse_folder, self.btnGetNewUniverseFolder, self.btnSetDefaultUniverse])

        self.lblSplineEphem_points_per_period = wx.StaticText(self, -1, "SplineEphem sample points per orbit period")
        self.txtSplineEphem_points_per_period = wx.TextCtrl(self, -1, "SplineEphem_points_per_period", size=(100,-1))

        self.lblSplineEphem_non_central_body_sun_points_per_period = wx.StaticText(self, -1, "SplineEphem sample points of the sun relative to the central body")
        self.txtSplineEphem_non_central_body_sun_points_per_period = wx.TextCtrl(self, -1, "SplineEphem_non_central_body_sun_points_per_period", size=(100,-1))

        self.lblSplineEphem_truncate_ephemeris_at_maximum_mission_epoch = wx.StaticText(self, -1, "Shorten SplineEphem to maximum mission epoch? (less memory but impedes MBH)")
        self.chkSplineEphem_truncate_ephemeris_at_maximum_mission_epoch = wx.CheckBox(self, -1)

        
        self.lblearliestPossibleEpoch = wx.StaticText(self, -1, "Earliest possible SplineEphem epoch")
        self.txtearliestPossibleEpoch = wx.TextCtrl(self, -1, "earliestPossibleEpoch")
        self.earliestPossibleEpochCalendar = wx.adv.CalendarCtrl(self, -1)
        earliestcalendarbox = wx.BoxSizer(wx.HORIZONTAL)
        earliestcalendarbox.AddMany([self.txtearliestPossibleEpoch, self.earliestPossibleEpochCalendar])

        self.lbllatestPossibleEpoch = wx.StaticText(self, -1, "latest possible SplineEphem epoch")
        self.txtlatestPossibleEpoch = wx.TextCtrl(self, -1, "latestPossibleEpoch")
        self.latestPossibleEpochCalendar = wx.adv.CalendarCtrl(self, -1)
        latestcalendarbox = wx.BoxSizer(wx.HORIZONTAL)
        latestcalendarbox.AddMany([self.txtlatestPossibleEpoch, self.latestPossibleEpochCalendar])

        self.lblperturb_SRP = wx.StaticText(self, -1, "Enable SRP")
        self.chkperturb_SRP = wx.CheckBox(self, -1)

        self.lblperturb_thirdbody = wx.StaticText(self, -1, "Enable third body")
        self.chkperturb_thirdbody = wx.CheckBox(self, -1)

        self.lblperturb_J2 = wx.StaticText(self, -1, "Enable central-body J2")
        self.chkperturb_J2 = wx.CheckBox(self, -1)

        self.lblspacecraft_area = wx.StaticText(self, -1, "Spacecraft area (in m^2)")
        self.txtspacecraft_area = wx.TextCtrl(self, -1, "spacecraft_area")

        self.lblcoefficient_of_reflectivity = wx.StaticText(self, -1, "Coefficient of reflectivity")
        self.txtcoefficient_of_reflectivity = wx.TextCtrl(self, -1, "coefficient_of_reflectivity")

        self.lblsolar_percentage = wx.StaticText(self, -1, "Solar percentage [0, 1]")
        self.txtsolar_percentage = wx.TextCtrl(self, -1, "solar_percentage")

        self.lblsolar_flux = wx.StaticText(self, -1, "Solar constant (flux at 1 AU)")
        self.txtsolar_flux = wx.TextCtrl(self, -1, "solar_flux")

        self.lblspeed_of_light_vac = wx.StaticText(self, -1, "Speed of light in a vacuum (m/s)")
        self.txtspeed_of_light_vac = wx.TextCtrl(self, -1, "speed_of_light_vac")

        self.lblintegrator_tolerance = wx.StaticText(self, -1, "Integrator tolerance")
        self.txtintegrator_tolerance = wx.TextCtrl(self, -1, "integrator_tolerance")

        self.lblpropagatorType = wx.StaticText(self, -1, "Propagator type")
        propagatorType_choices = ["Keplerian", "Integrator"]
        self.cmbpropagatorType = wx.ComboBox(self, -1, choices=propagatorType_choices, style=wx.CB_READONLY)

        self.lblintegratorType = wx.StaticText(self, -1, "Integrator type")
        integratorType_choices = ["rk7813M adaptive step", "rk8 fixed step"]
        self.cmbintegratorType = wx.ComboBox(self, -1, choices=integratorType_choices, style=wx.CB_READONLY)

        self.lblintegration_time_step_size = wx.StaticText(self, -1, "Integrator time step size (seconds)")
        self.txtintegration_time_step_size = wx.TextCtrl(self, -1, "integration_time_step_size")

        ephemerisgrid.AddMany([self.lblephemeris_source, self.cmbephemeris_source,
                              self.lblSPICE_leap_seconds_kernel, self.txtSPICE_leap_seconds_kernel,
                              self.lblSPICE_reference_frame_kernel, self.txtSPICE_reference_frame_kernel,
                              self.lbluniverse_folder, UniverseButtonSizer,
                              self.lblSplineEphem_points_per_period, self.txtSplineEphem_points_per_period,
                              self.lblSplineEphem_non_central_body_sun_points_per_period, self.txtSplineEphem_non_central_body_sun_points_per_period,
                              self.lblSplineEphem_truncate_ephemeris_at_maximum_mission_epoch, self.chkSplineEphem_truncate_ephemeris_at_maximum_mission_epoch,
                              self.lblearliestPossibleEpoch, earliestcalendarbox,
                              self.lbllatestPossibleEpoch, latestcalendarbox])
        perturbgrid.AddMany([ self.lblperturb_SRP, self.chkperturb_SRP,
                              self.lblperturb_thirdbody, self.chkperturb_thirdbody,
                              self.lblperturb_J2, self.chkperturb_J2,
                              self.lblspacecraft_area, self.txtspacecraft_area,
                              self.lblcoefficient_of_reflectivity, self.txtcoefficient_of_reflectivity,
                              self.lblsolar_percentage, self.txtsolar_percentage,
                              self.lblsolar_flux, self.txtsolar_flux,
                              self.lblspeed_of_light_vac, self.txtspeed_of_light_vac
                              ])

        integratorgrid.AddMany([self.lblintegrator_tolerance, self.txtintegrator_tolerance,
                                self.lblpropagatorType, self.cmbpropagatorType,
                                self.lblintegratorType, self.cmbintegratorType,
                                self.lblintegration_time_step_size, self.txtintegration_time_step_size])


        lblLeftTitle = wx.StaticText(self, -1, "Ephemeris settings")
        vboxleft = wx.BoxSizer(wx.VERTICAL)
        vboxleft.AddMany([lblLeftTitle, ephemerisgrid])

        lblRightTitle = wx.StaticText(self, -1, "Perturbation settings")
        vboxright = wx.BoxSizer(wx.VERTICAL)
        vboxright.AddMany([lblRightTitle, perturbgrid])

        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        lblLeftTitle.SetFont(font)
        lblRightTitle.SetFont(font)

        self.mainbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.mainbox.Add(vboxleft)
        self.mainbox.AddSpacer(20)
        self.mainbox.Add(vboxright)


        spiralgrid = wx.GridSizer(2,2,5,5)
        self.lblspiral_segments = wx.StaticText(self, -1, "Number of spiral segments")
        self.txtspiral_segments = wx.TextCtrl(self, -1, "spiral_segments")
        spiralgrid.AddMany([self.lblspiral_segments, self.txtspiral_segments])
        lblBottomTitle = wx.StaticText(self, -1, "Spiral settings")
        lblBottomTitle.SetFont(font)
        vboxspiral = wx.BoxSizer(wx.VERTICAL)
        vboxspiral.AddMany([lblBottomTitle, spiralgrid])

        StateRepresentationgrid = wx.GridSizer(3,2,5,5)
        StateRepresentationChoices = ['Cartesian', 'SphericalRADEC', 'SphericalAZFPA', 'COE', 'MEE', "IncomingBplane", "OutgoingBplane"]
        self.lblPeriapseBoundaryStateRepresentation = wx.StaticText(self, -1, "PeriapseBoundary state representation")
        self.cmbPeriapseBoundaryStateRepresentation = wx.ComboBox(self, -1, choices=StateRepresentationChoices, style=wx.CB_READONLY)
        self.lblParallelShootingStateRepresentation = wx.StaticText(self, -1, "Parallel shooting decision variable state representation")
        self.cmbParallelShootingStateRepresentation = wx.ComboBox(self, -1, choices=StateRepresentationChoices[0:5], style=wx.CB_READONLY) #parallel shooting can't use the asymptotic coordinate sets
        self.lblParallelShootingConstraintStateRepresentation = wx.StaticText(self, -1, "Parallel shooting constraint state representation")
        self.cmbParallelShootingConstraintStateRepresentation = wx.ComboBox(self, -1, choices=['Cartesian','same as encoded state representation'], style=wx.CB_READONLY)
        StateRepresentationgrid.AddMany([self.lblPeriapseBoundaryStateRepresentation, self.cmbPeriapseBoundaryStateRepresentation,
                                         self.lblParallelShootingStateRepresentation, self.cmbParallelShootingStateRepresentation,
                                         self.lblParallelShootingConstraintStateRepresentation, self.cmbParallelShootingConstraintStateRepresentation])
        lblStateRepresentationBottomTitle = wx.StaticText(self, -1, "State Representation settings")
        lblStateRepresentationBottomTitle.SetFont(font)
        vboxStateRepresentation = wx.BoxSizer(wx.VERTICAL)
        vboxStateRepresentation.AddMany([lblStateRepresentationBottomTitle, StateRepresentationgrid])

        self.mainvbox = wx.BoxSizer(wx.VERTICAL)
        self.mainvbox.Add(self.mainbox)
        self.mainvbox.AddSpacer(20)
        self.mainvbox.AddMany([vboxspiral, integratorgrid, vboxStateRepresentation])

        self.SetSizer(self.mainvbox)
        self.SetupScrolling()

        #bindings
        self.cmbephemeris_source.Bind(wx.EVT_COMBOBOX,self.Changeephemeris_source)
        self.txtSPICE_leap_seconds_kernel.Bind(wx.EVT_KILL_FOCUS,self.ChangeSPICE_leap_seconds_kernel)
        self.txtSPICE_reference_frame_kernel.Bind(wx.EVT_KILL_FOCUS,self.ChangeSPICE_reference_frame_kernel)
        self.txtuniverse_folder.Bind(wx.EVT_KILL_FOCUS,self.Changeuniverse_folder)
        self.txtSplineEphem_points_per_period.Bind(wx.EVT_KILL_FOCUS, self.ChangeSplineEphem_points_per_period)
        self.txtSplineEphem_non_central_body_sun_points_per_period.Bind(wx.EVT_KILL_FOCUS, self.ChangeSplineEphem_non_central_body_sun_points_per_period)
        self.chkSplineEphem_truncate_ephemeris_at_maximum_mission_epoch.Bind(wx.EVT_CHECKBOX, self.ChangeSplineEphem_truncate_ephemeris_at_maximum_mission_epoch)
        self.txtearliestPossibleEpoch.Bind(wx.EVT_KILL_FOCUS,self.ChangeearliestPossibleEpoch)
        self.earliestPossibleEpochCalendar.Bind(wx.adv.EVT_CALENDAR_SEL_CHANGED, self.ChangeearliestPossibleEpochCalendar)
        self.txtlatestPossibleEpoch.Bind(wx.EVT_KILL_FOCUS,self.ChangelatestPossibleEpoch)
        self.latestPossibleEpochCalendar.Bind(wx.adv.EVT_CALENDAR_SEL_CHANGED, self.ChangelatestPossibleEpochCalendar)
        self.btnGetNewUniverseFolder.Bind(wx.EVT_BUTTON,self.GetNewUniverseFolder)
        self.btnSetDefaultUniverse.Bind(wx.EVT_BUTTON,self.SetDefaultUniverse)
        self.chkperturb_SRP.Bind(wx.EVT_CHECKBOX,self.Changeperturb_SRP)
        self.chkperturb_thirdbody.Bind(wx.EVT_CHECKBOX,self.Changeperturb_thirdbody)
        self.chkperturb_J2.Bind(wx.EVT_CHECKBOX,self.Changeperturb_J2)
        self.txtspacecraft_area.Bind(wx.EVT_KILL_FOCUS,self.Changespacecraft_area)
        self.txtcoefficient_of_reflectivity.Bind(wx.EVT_KILL_FOCUS,self.Changecoefficient_of_reflectivity)
        self.txtsolar_percentage.Bind(wx.EVT_KILL_FOCUS,self.Changesolar_percentage)
        self.txtsolar_flux.Bind(wx.EVT_KILL_FOCUS,self.Changesolar_flux)
        self.txtspeed_of_light_vac.Bind(wx.EVT_KILL_FOCUS,self.Changespeed_of_light_vac)
        self.txtspiral_segments.Bind(wx.EVT_KILL_FOCUS, self.Changespiral_segments)
        self.txtintegrator_tolerance.Bind(wx.EVT_KILL_FOCUS, self.ChangeIntegratorTolerance)
        self.cmbpropagatorType.Bind(wx.EVT_COMBOBOX, self.ChangepropagatorType)
        self.cmbintegratorType.Bind(wx.EVT_COMBOBOX, self.ChangeintegratorType)
        self.txtintegration_time_step_size.Bind(wx.EVT_KILL_FOCUS, self.Changeintegration_time_step_size)
        self.cmbPeriapseBoundaryStateRepresentation.Bind(wx.EVT_COMBOBOX, self.ChangePeriapseBoundaryStateRepresentation)
        self.cmbParallelShootingStateRepresentation.Bind(wx.EVT_COMBOBOX, self.ChangeParallelShootingStateRepresentation)
        self.cmbParallelShootingConstraintStateRepresentation.Bind(wx.EVT_COMBOBOX, self.ChangeParallelShootingConstraintStateRepresentation)

    def update(self):

        self.cmbephemeris_source.SetSelection(self.missionoptions.ephemeris_source)
        self.txtSPICE_leap_seconds_kernel.SetValue(str(self.missionoptions.SPICE_leap_seconds_kernel))
        self.txtSPICE_reference_frame_kernel.SetValue(str(self.missionoptions.SPICE_reference_frame_kernel))
        self.txtuniverse_folder.SetValue(self.missionoptions.universe_folder)
        self.txtSplineEphem_points_per_period.SetValue(str(self.missionoptions.SplineEphem_points_per_period))
        self.txtSplineEphem_non_central_body_sun_points_per_period.SetValue(str(self.missionoptions.SplineEphem_non_central_body_sun_points_per_period))
        self.chkSplineEphem_truncate_ephemeris_at_maximum_mission_epoch.SetValue(self.missionoptions.SplineEphem_truncate_ephemeris_at_maximum_mission_epoch)
        self.txtearliestPossibleEpoch.SetValue(str(self.missionoptions.earliestPossibleEpoch))
        earliestPossibleDate = wx.DateTime.FromJDN(self.missionoptions.earliestPossibleEpoch + 2400000.5)
        self.earliestPossibleEpochCalendar.SetDate(earliestPossibleDate.MakeUTC())        
        self.txtlatestPossibleEpoch.SetValue(str(self.missionoptions.latestPossibleEpoch))
        latestPossibleDate = wx.DateTime.FromJDN(self.missionoptions.latestPossibleEpoch + 2400000.5)
        self.latestPossibleEpochCalendar.SetDate(latestPossibleDate.MakeUTC())
        self.chkperturb_SRP.SetValue(self.missionoptions.perturb_SRP)
        self.chkperturb_thirdbody.SetValue(self.missionoptions.perturb_thirdbody)
        self.chkperturb_J2.SetValue(self.missionoptions.perturb_J2)
        self.txtspacecraft_area.SetValue(str(self.missionoptions.spacecraft_area))
        self.txtcoefficient_of_reflectivity.SetValue(str(self.missionoptions.coefficient_of_reflectivity))
        self.txtsolar_percentage.SetValue(str(self.missionoptions.solar_percentage))
        self.txtsolar_flux.SetValue(str(self.missionoptions.solar_flux))
        self.txtspeed_of_light_vac.SetValue(str(self.missionoptions.speed_of_light_vac))
        self.txtspiral_segments.SetValue(str(self.missionoptions.spiral_segments))
        self.txtintegrator_tolerance.SetValue(str(self.missionoptions.integrator_tolerance))
        self.cmbpropagatorType.SetSelection(self.missionoptions.propagatorType)
        self.cmbintegratorType.SetSelection(self.missionoptions.integratorType)
        self.txtintegration_time_step_size.SetValue(str(self.missionoptions.integration_time_step_size))
        self.cmbPeriapseBoundaryStateRepresentation.SetSelection(self.missionoptions.PeriapseBoundaryStateRepresentation)
        self.cmbParallelShootingStateRepresentation.SetSelection(self.missionoptions.ParallelShootingStateRepresentation)
        self.cmbParallelShootingConstraintStateRepresentation.SetSelection(self.missionoptions.ParallelShootingConstraintStateRepresentation)

        #if SplineEphem is active, show the SplineEphem controls
        if self.missionoptions.ephemeris_source == 2:
            self.lblSplineEphem_points_per_period.Show(True)
            self.txtSplineEphem_points_per_period.Show(True)
            self.lblSplineEphem_truncate_ephemeris_at_maximum_mission_epoch.Show(True)
            self.chkSplineEphem_truncate_ephemeris_at_maximum_mission_epoch.Show(True)
        else:            
            self.lblSplineEphem_points_per_period.Show(False)
            self.txtSplineEphem_points_per_period.Show(False)
            self.lblSplineEphem_truncate_ephemeris_at_maximum_mission_epoch.Show(False)
            self.chkSplineEphem_truncate_ephemeris_at_maximum_mission_epoch.Show(False)

        #if SRP is disabled, make the options associated with it invisible
        if self.missionoptions.perturb_SRP == 1:
            self.lblspacecraft_area.Show(True)
            self.lblcoefficient_of_reflectivity.Show(True)
            self.lblsolar_percentage.Show(True)
            self.lblsolar_flux.Show(True)
            self.lblspeed_of_light_vac.Show(True)
            self.txtspacecraft_area.Show(True)
            self.txtcoefficient_of_reflectivity.Show(True)
            self.txtsolar_percentage.Show(True)
            self.txtsolar_flux.Show(True)
            self.txtspeed_of_light_vac.Show(True)
        else:
            self.lblspacecraft_area.Show(False)
            self.lblcoefficient_of_reflectivity.Show(False)
            self.lblsolar_percentage.Show(False)
            self.lblsolar_flux.Show(False)
            self.lblspeed_of_light_vac.Show(False)
            self.txtspacecraft_area.Show(False)
            self.txtcoefficient_of_reflectivity.Show(False)
            self.txtsolar_percentage.Show(False)
            self.txtsolar_flux.Show(False)
            self.txtspeed_of_light_vac.Show(False)

        #only enable propagator switch if using a phase type that supports it
        if self.missionoptions.mission_type in [6, 7, 8, 9]:
            self.lblpropagatorType.Show(True)
            self.cmbpropagatorType.Show(True)
        else:
            self.lblpropagatorType.Show(False)
            self.cmbpropagatorType.Show(False)

        #only enable integrator tolerance for mission types that integrate
        if self.missionoptions.integratorType == 0:
            self.lblintegrator_tolerance.Show(True)
            self.txtintegrator_tolerance.Show(True)
        else:
            self.lblintegrator_tolerance.Show(False)
            self.txtintegrator_tolerance.Show(False)

        self.lblintegratorType.Show(True)
        self.cmbintegratorType.Show(True)

        #re-size the panel
        self.Layout()
        if platform.system() == 'Windows':
            self.SetupScrolling(scrollToTop=False)

    #event handlers for physics options
    def Changeephemeris_source(self, e):
        self.missionoptions.ephemeris_source = self.cmbephemeris_source.GetSelection()

    def ChangeSPICE_leap_seconds_kernel(self, e):
        e.Skip()
        self.missionoptions.SPICE_leap_seconds_kernel = self.txtSPICE_leap_seconds_kernel.GetValue()

    def ChangeSPICE_reference_frame_kernel(self, e):
        e.Skip()
        self.missionoptions.SPICE_reference_frame_kernel = self.txtSPICE_reference_frame_kernel.GetValue()
        
    def ChangeSplineEphem_points_per_period(self, e):
        e.Skip()
        self.missionoptions.SplineEphem_points_per_period = int(self.txtSplineEphem_points_per_period.GetValue())

    def ChangeSplineEphem_non_central_body_sun_points_per_period(self, e):
        e.Skip()
        self.missionoptions.SplineEphem_non_central_body_sun_points_per_period = int(self.txtSplineEphem_non_central_body_sun_points_per_period.GetValue())

    def ChangeSplineEphem_truncate_ephemeris_at_maximum_mission_epoch(self, e):
        e.Skip()
        self.missionoptions.SplineEphem_truncate_ephemeris_at_maximum_mission_epoch = int(self.chkSplineEphem_truncate_ephemeris_at_maximum_mission_epoch.GetValue())
          
    def ChangeearliestPossibleEpoch(self, e):
        e.Skip()

        dateString = self.txtearliestPossibleEpoch.GetValue()

        from timeUtilities import stringToJD

        self.missionoptions.earliestPossibleEpoch = stringToJD(dateString, self.missionoptions.universe_folder)

        self.update()

    def ChangeearliestPossibleEpochCalendar(self, e):
        epoch = self.earliestPossibleEpochCalendar.GetDate()
        epoch = epoch.FromTimezone(wx.DateTime.TimeZone(offset=0))
        self.missionoptions.earliestPossibleEpoch = epoch.GetMJD()
        self.update()
        
    def ChangelatestPossibleEpoch(self, e):
        e.Skip()

        dateString = self.txtlatestPossibleEpoch.GetValue()

        from timeUtilities import stringToJD

        self.missionoptions.latestPossibleEpoch = stringToJD(dateString, self.missionoptions.universe_folder)

        self.update()

    def ChangelatestPossibleEpochCalendar(self, e):
        epoch = self.latestPossibleEpochCalendar.GetDate()
        epoch = epoch.FromTimezone(wx.DateTime.TimeZone(offset=0))
        self.missionoptions.latestPossibleEpoch = epoch.GetMJD()
        self.update()

    def Changeuniverse_folder(self, e):
        e.Skip()
        self.missionoptions.universe_folder = self.txtuniverse_folder.GetValue()

    def GetNewUniverseFolder(self, e):
        #file load dialog to get name of universe folder
        dlg = wx.DirDialog(self, "Choose a Universe folder", self.parent.Parent.dirname)
        if dlg.ShowModal() == wx.ID_OK:
            self.missionoptions.universe_folder = dlg.GetPath()
            self.txtuniverse_folder.SetValue(self.missionoptions.universe_folder)
        dlg.Destroy()

    def SetDefaultUniverse(self, e):
        self.missionoptions.universe_folder = self.parent.Parent.default_universe_path
        self.txtuniverse_folder.SetValue(self.missionoptions.universe_folder)

    def Changeperturb_SRP(self, e):
        self.missionoptions.perturb_SRP = int(self.chkperturb_SRP.GetValue())
        self.update()

    def Changeperturb_thirdbody(self, e):
        self.missionoptions.perturb_thirdbody = int(self.chkperturb_thirdbody.GetValue())
        self.parent.update()

    def Changeperturb_J2(self, e):
        self.missionoptions.perturb_J2 = int(self.chkperturb_J2.GetValue())
        self.parent.update()

    def Changespacecraft_area(self, e):
        e.Skip()
        self.missionoptions.spacecraft_area = eval(self.txtspacecraft_area.GetValue())

    def Changecoefficient_of_reflectivity(self, e):
        e.Skip()
        self.missionoptions.coefficient_of_reflectivity = eval(self.txtcoefficient_of_reflectivity.GetValue())

    def Changesolar_percentage(self, e):
        e.Skip()
        self.missionoptions.solar_percentage = eval(self.txtsolar_percentage.GetValue())

    def Changesolar_flux(self, e):
        e.Skip()
        self.missionoptions.solar_flux = eval(self.txtsolar_flux.GetValue())

    def Changespeed_of_light_vac(self, e):
        e.Skip()
        self.missionoptions.speed_of_light_vac = eval(self.txtspeed_of_light_vac.GetValue())

    def Changespiral_segments(self, e):
        e.Skip()
        self.missionoptions.spiral_segments = int(self.txtspiral_segments.GetValue())

    def ChangeIntegratorTolerance(self, e):
        e.Skip()
        self.missionoptions.integrator_tolerance = eval(self.txtintegrator_tolerance.GetValue())

    def ChangepropagatorType(self, e):
        self.missionoptions.propagatorType = self.cmbpropagatorType.GetSelection()
        self.update()
        e.Skip()

    def ChangeintegratorType(self, e):
        self.missionoptions.integratorType = self.cmbintegratorType.GetSelection()
        self.parent.update()
        e.Skip()

    def Changeintegration_time_step_size(self, e):
        self.missionoptions.integration_time_step_size = eval(self.txtintegration_time_step_size.GetValue())
        self.parent.update()
        e.Skip()
        
    def ChangePeriapseBoundaryStateRepresentation(self, e):
        self.missionoptions.PeriapseBoundaryStateRepresentation = self.cmbPeriapseBoundaryStateRepresentation.GetSelection()
        self.parent.update()
        self.missionoptions.DisassembleMasterDecisionVector()
        self.missionoptions.ConvertDecisionVector()
        self.missionoptions.AssembleMasterDecisionVector()
        e.Skip()

    def ChangeParallelShootingStateRepresentation(self, e):
        self.missionoptions.ParallelShootingStateRepresentation = self.cmbParallelShootingStateRepresentation.GetSelection()
        self.parent.update()
        self.missionoptions.DisassembleMasterDecisionVector()
        self.missionoptions.ConvertDecisionVector()
        self.missionoptions.AssembleMasterDecisionVector()
        e.Skip()

    def ChangeParallelShootingConstraintStateRepresentation(self, e):
        self.missionoptions.ParallelShootingConstraintStateRepresentation = self.cmbParallelShootingConstraintStateRepresentation.GetSelection()
        self.parent.update()
        e.Skip()