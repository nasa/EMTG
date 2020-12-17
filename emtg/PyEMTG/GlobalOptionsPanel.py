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

class GlobalOptionsPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, missionoptions):
        self.missionoptions = missionoptions
        self.parent = parent

        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)
        
        globaloptionsgrid = wx.FlexGridSizer(20,2,5,5)
        self.lblMissionName = wx.StaticText(self, -1, "Mission Name")
        self.txtMissionName = wx.TextCtrl(self, -1, "mission_name", size=(500,-1))

        self.lblMissionType = wx.StaticText(self, -1, "Mission Type")
        phasetypes = ['MGALTS (not yet implemented)', 'FBLTS (not yet implemented)', 'MGALT', 'FBLT', 'PSBI', 'PSFB', 'MGAnDSMs', 'CoastPhase', 'SundmanCoastPhase', 'Variable phase type', 'ProbeEntryPhase', 'ControlLawThrustPhase']
        self.cmbMissionType = wx.ComboBox(self, -1, choices=phasetypes, style=wx.CB_READONLY)

        self.lblobjective_type = wx.StaticText(self, -1, "Objective function")
        objectivetypes = ['0: minimum deterministic deltaV','1: minimum time','2: maximum final mass','3: maximize initial mass',
                          '4: depart as late as possible','5: depart as early as possible',
                          '6: maximize orbit energy','7: minimize launch mass','8: arrive as early as possible',
                          '9: arrive as late as possible','10: minimum propellant (not the same as 2)','11: maximum dry/wet ratio',
                          '12: maximum arrival kinetic energy', '13: minimum BOL power', '14: maximize log_10(final mass)', '15: maximize log_e(final mass)',
                          '16: maximum dry mass margin', '17: maximum dry mass', '18: maximum log_10(dry mass)', '19: maximum log_e(dry mass)',
                          '20: minimize chemical fuel', '21: minimize chemical oxidizer', '22: minimize electric propellant', '23: minimize total propellant',
                          '24: minimize waypoint tracking error', '25: minimize initial impulse magnitude', 
                          '26: maximize distance from central body']
        self.cmbobjective_type = wx.ComboBox(self, -1, choices=objectivetypes, style = wx.CB_READONLY)

        self.lblinclude_initial_impulse_in_cost = wx.StaticText(self, -1, "Include initial impulse in total deterministic delta-v")
        self.chkinclude_initial_impulse_in_cost = wx.CheckBox(self, -1)
        
        self.lblobjective_journey = wx.StaticText(self,-1,"Which journey to optimize?")
        self.txtobjective_journey = wx.TextCtrl(self,-1,"objective_journey")

        self.lblwaypoint_file_path = wx.StaticText(self, -1, "Waypoint file")
        self.txtwaypoint_file_path = wx.TextCtrl(self, -1, "waypoint_file_path", size=(450, -1))
        self.btnwaypoint_file_path = wx.Button(self, -1, "...")
        WayPointSizer = wx.BoxSizer(wx.HORIZONTAL)
        WayPointSizer.AddMany([self.txtwaypoint_file_path, self.btnwaypoint_file_path])

        self.lblcovariance_file_path = wx.StaticText(self, -1, "Covariance file")
        self.txtcovariance_file_path = wx.TextCtrl(self, -1, "covariance_file_path", size=(450, -1))
        self.btncovariance_file_path = wx.Button(self, -1, "...")
        CovarianceSizer = wx.BoxSizer(wx.HORIZONTAL)
        CovarianceSizer.AddMany([self.txtcovariance_file_path, self.btncovariance_file_path])
        
        self.lbllaunch_window_open_date = wx.StaticText(self, -1, "Launch window open date")
        self.txtlaunch_window_open_date = wx.TextCtrl(self, -1, "launch_window_open_date", size=(450, -1))
        self.LaunchDateCalendar = wx.adv.CalendarCtrl(self, -1)
        CalendarSpacer = wx.StaticText(self, -1, "")
        
        self.lblnum_timesteps = wx.StaticText(self, -1, "Number of time-steps")
        self.txtnum_timesteps = wx.TextCtrl(self, -1, "num_timesteps")

        self.lblstop_after_journey = wx.StaticText(self, -1, "Stop after journey (indexed from 0)")
        self.txtstop_after_journey = wx.TextCtrl(self, -1, "stop_after_journey")
                
        globaloptionsgrid.AddMany(  [self.lblMissionName, self.txtMissionName,
                                    self.lblMissionType, self.cmbMissionType,
                                    self.lblobjective_type, self.cmbobjective_type,
                                    self.lblobjective_journey,self.txtobjective_journey,
                                    self.lblwaypoint_file_path, WayPointSizer,
                                    self.lblcovariance_file_path, CovarianceSizer,
                                    self.lblinclude_initial_impulse_in_cost, self.chkinclude_initial_impulse_in_cost,
                                    self.lbllaunch_window_open_date, self.txtlaunch_window_open_date,
                                    CalendarSpacer, self.LaunchDateCalendar,
                                    self.lblnum_timesteps, self.txtnum_timesteps,
                                    self.lblstop_after_journey, self.txtstop_after_journey])

        globaloptionsgrid.SetFlexibleDirection(wx.BOTH)

        #constraint fields
        constraintgrid = wx.FlexGridSizer(20, 2, 5, 5)

        self.lblRLA_bounds = wx.StaticText(self, -1, "RLA bounds (degrees)")
        self.txtRLA_bounds_lower = wx.TextCtrl(self, -1, "RLA_bounds[0]")
        self.txtRLA_bounds_upper = wx.TextCtrl(self, -1, "RLA_bounds[1]")
        RLAbox = wx.BoxSizer(wx.HORIZONTAL)
        RLAbox.AddMany([self.txtRLA_bounds_lower, self.txtRLA_bounds_upper])

        self.lblDLA_bounds = wx.StaticText(self, -1, "DLA bounds (degrees)")
        self.txtDLA_bounds_lower = wx.TextCtrl(self, -1, "DLA_bounds[0]")
        self.txtDLA_bounds_upper = wx.TextCtrl(self, -1, "DLA_bounds[1]")
        DLAbox = wx.BoxSizer(wx.HORIZONTAL)
        DLAbox.AddMany([self.txtDLA_bounds_lower, self.txtDLA_bounds_upper])

        self.lblglobal_timebounded = wx.StaticText(self, -1, "Enable mission time bounds")
        self.chkglobal_timebounded = wx.CheckBox(self, -1)

        self.lbltotal_flight_time_bounds = wx.StaticText(self, -1, "Global flight time bounds (days)")
        self.txttotal_flight_time_bounds_lower = wx.TextCtrl(self, -1, "total_flight_time_bounds[0]")
        self.txttotal_flight_time_bounds_upper = wx.TextCtrl(self, -1, "total_flight_time_bounds[1]")
        GlobalTimebox = wx.BoxSizer(wx.HORIZONTAL)
        GlobalTimebox.AddMany([self.txttotal_flight_time_bounds_lower, self.txttotal_flight_time_bounds_upper])

        self.lblforced_post_launch_coast = wx.StaticText(self, -1, "Forced post-launch coast duration (days)")
        self.txtforced_post_launch_coast = wx.TextCtrl(self, -1, "forced_post_launch_coast")

        self.lblforced_pre_flyby_coast = wx.StaticText(self, -1, "Forced pre-flyby coast duration (days)")
        self.txtforced_pre_flyby_coast = wx.TextCtrl(self, -1, "forced_pre_flyby_coast")
        self.lblforced_post_flyby_coast = wx.StaticText(self, -1, "Forced post-flyby coast duration (days)")
        self.txtforced_post_flyby_coast = wx.TextCtrl(self, -1, "forced_post_flyby_coast")

        self.lblTCM_post_launch = wx.StaticText(self, -1, "Magnitude of post-launch TCM (km/s)")
        self.txtTCM_post_launch = wx.TextCtrl(self, -1, "TCM_post_launch")
        self.lblTCM_pre_flyby = wx.StaticText(self, -1, "Magnitude of pre-flyby TCM (km/s)")
        self.txtTCM_pre_flyby = wx.TextCtrl(self, -1, "Magnitude of pre-flyby TCM (km/s)")
        self.lblTCM_maneuver_fraction = wx.StaticText(self, -1, "Magnitude of post-DSM TCMs as a fraction of DSM magnitude")
        self.txtTCM_maneuver_fraction = wx.TextCtrl(self, -1, "TCM_maneuver_fraction")

        self.lblfinal_mass_constraint= wx.StaticText(self, -1, "Final/dry mass constraint value (kg)")
        self.txtfinal_mass_constraint_LowerBound = wx.TextCtrl(self, -1, "finaL_mass_constraint_bounds[0]")
        self.txtfinal_mass_constraint_UpperBound = wx.TextCtrl(self, -1, "finaL_mass_constraint_bounds[1]")
        mass_constraint_sizer = wx.BoxSizer(wx.HORIZONTAL)
        mass_constraint_sizer.AddMany([self.txtfinal_mass_constraint_LowerBound, self.txtfinal_mass_constraint_UpperBound])

        self.lblconstrain_final_mass = wx.StaticText(self, -1, "Constrain final mass?")
        self.chkconstrain_final_mass = wx.CheckBox(self, -1)

        self.lblconstrain_dry_mass = wx.StaticText(self, -1, "Constrain dry mass?")
        self.chkconstrain_dry_mass = wx.CheckBox(self, -1)

        constraintgrid.AddMany([self.lblRLA_bounds, RLAbox,
                                self.lblDLA_bounds, DLAbox,
                                self.lblglobal_timebounded, self.chkglobal_timebounded,
                                self.lbltotal_flight_time_bounds, GlobalTimebox,
                                self.lblforced_post_launch_coast, self.txtforced_post_launch_coast,
                                self.lblforced_pre_flyby_coast, self.txtforced_pre_flyby_coast,
                                self.lblforced_post_flyby_coast, self.txtforced_post_flyby_coast,
                                self.lblTCM_post_launch, self.txtTCM_post_launch,
                                self.lblTCM_pre_flyby, self.txtTCM_pre_flyby,
                                self.lblTCM_maneuver_fraction, self.txtTCM_maneuver_fraction,
                                self.lblfinal_mass_constraint, mass_constraint_sizer,
                                self.lblconstrain_final_mass, self.chkconstrain_final_mass,
                                self.lblconstrain_dry_mass, self.chkconstrain_dry_mass])

        vboxleft = wx.BoxSizer(wx.VERTICAL)
        vboxright = wx.BoxSizer(wx.VERTICAL)
        lblLeftTitle = wx.StaticText(self, -1, "Global mission options")
        lblRightTitle = wx.StaticText(self, -1, "Global mission constraints")
        vboxleft.Add(lblLeftTitle)
        vboxleft.Add(globaloptionsgrid)
        vboxright.Add(lblRightTitle)
        vboxright.Add(constraintgrid)
        
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        lblLeftTitle.SetFont(font)
        lblRightTitle.SetFont(font)

        self.mainbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.mainbox.Add(vboxleft)
        self.mainbox.AddSpacer(20)
        self.mainbox.Add(vboxright)

        self.SetSizer(self.mainbox)
        self.SetupScrolling()

        #bindings
        self.cmbMissionType.Bind(wx.EVT_COMBOBOX, self.ChangeMissionType)
        self.cmbobjective_type.Bind(wx.EVT_COMBOBOX,self.Changeobjective_type)
        self.txtobjective_journey.Bind(wx.EVT_KILL_FOCUS, self.Changeobjective_journey)
        self.txtwaypoint_file_path.Bind(wx.EVT_KILL_FOCUS, self.Changewaypoint_file_path)
        self.txtcovariance_file_path.Bind(wx.EVT_KILL_FOCUS, self.Changecovariance_file_path)
        self.btnwaypoint_file_path.Bind(wx.EVT_BUTTON, self.Clickwaypoint_file_path_button)
        self.btncovariance_file_path.Bind(wx.EVT_BUTTON, self.Clickcovariance_file_path_button)
        self.chkinclude_initial_impulse_in_cost.Bind(wx.EVT_CHECKBOX,self.Changeinclude_initial_impulse_in_cost)
        self.txtlaunch_window_open_date.Bind(wx.EVT_KILL_FOCUS,self.Changelaunch_window_open_date)
        self.LaunchDateCalendar.Bind(wx.adv.EVT_CALENDAR_SEL_CHANGED, self.ChangeLaunchDateCalendar)
        self.txtnum_timesteps.Bind(wx.EVT_KILL_FOCUS,self.Changenum_timesteps)
        self.txtDLA_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.ChangeDLA_bounds_lower)
        self.txtDLA_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.ChangeDLA_bounds_upper)
        self.txtRLA_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.ChangeRLA_bounds_lower)
        self.txtRLA_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.ChangeRLA_bounds_upper)
        self.chkglobal_timebounded.Bind(wx.EVT_CHECKBOX,self.Changeglobal_timebounded)
        self.txttotal_flight_time_bounds_lower.Bind(wx.EVT_KILL_FOCUS,self.Changetotal_flight_time_bounds_lower)
        self.txttotal_flight_time_bounds_upper.Bind(wx.EVT_KILL_FOCUS,self.Changetotal_flight_time_bounds_upper)
        self.txtforced_post_launch_coast.Bind(wx.EVT_KILL_FOCUS,self.Changeforced_post_launch_coast)
        self.txtforced_pre_flyby_coast.Bind(wx.EVT_KILL_FOCUS,self.Changeforced_pre_flyby_coast)
        self.txtforced_post_flyby_coast.Bind(wx.EVT_KILL_FOCUS,self.Changeforced_post_flyby_coast)
        self.txtTCM_post_launch.Bind(wx.EVT_KILL_FOCUS,self.ChangeTCM_post_launch)
        self.txtTCM_pre_flyby.Bind(wx.EVT_KILL_FOCUS,self.ChangeTCM_pre_flyby)
        self.txtTCM_maneuver_fraction.Bind(wx.EVT_KILL_FOCUS,self.ChangeTCM_maneuver_fraction)
        self.txtfinal_mass_constraint_LowerBound.Bind(wx.EVT_KILL_FOCUS,self.Changefinal_mass_constraint_LowerBound)
        self.txtfinal_mass_constraint_UpperBound.Bind(wx.EVT_KILL_FOCUS,self.Changefinal_mass_constraint_UpperBound)
        self.chkconstrain_dry_mass.Bind(wx.EVT_CHECKBOX, self.Changeconstrain_dry_mass)
        self.chkconstrain_final_mass.Bind(wx.EVT_CHECKBOX, self.Changeconstrain_final_mass)
        self.txtstop_after_journey.Bind(wx.EVT_KILL_FOCUS, self.Changestop_after_journey)

        

        #put this bind statement last to prevent the annoying thing where mission_name becomes mission_name - we have to set mission_name first
        self.txtMissionName.SetValue(self.missionoptions.mission_name)
        self.txtMissionName.Bind(wx.EVT_KILL_FOCUS, self.ChangeMissionName)

    def update(self):
        self.txtMissionName.SetValue(self.missionoptions.mission_name)
        self.cmbMissionType.SetSelection(self.missionoptions.mission_type)
        self.cmbobjective_type.SetSelection(self.missionoptions.objective_type)
        self.chkinclude_initial_impulse_in_cost.SetValue(self.missionoptions.include_initial_impulse_in_cost)
        self.txtlaunch_window_open_date.SetValue(str(self.missionoptions.launch_window_open_date))
        CurrentLaunchDate = wx.DateTime.FromJDN(self.missionoptions.launch_window_open_date + 2400000.5)
        self.LaunchDateCalendar.SetDate(CurrentLaunchDate.MakeUTC())
        self.txtnum_timesteps.SetValue(str(self.missionoptions.num_timesteps))
        self.txtDLA_bounds_lower.SetValue(str(self.missionoptions.DLA_bounds[0]))
        self.txtDLA_bounds_upper.SetValue(str(self.missionoptions.DLA_bounds[1]))
        self.txtRLA_bounds_lower.SetValue(str(self.missionoptions.RLA_bounds[0]))
        self.txtRLA_bounds_upper.SetValue(str(self.missionoptions.RLA_bounds[1]))
        self.chkglobal_timebounded.SetValue(self.missionoptions.global_timebounded)
        self.txttotal_flight_time_bounds_lower.SetValue(str(self.missionoptions.total_flight_time_bounds[0]))
        self.txttotal_flight_time_bounds_upper.SetValue(str(self.missionoptions.total_flight_time_bounds[1]))
        self.txtforced_post_launch_coast.SetValue(str(self.missionoptions.forced_post_launch_coast))
        self.txtforced_pre_flyby_coast.SetValue(str(self.missionoptions.forced_pre_flyby_coast))
        self.txtforced_post_flyby_coast.SetValue(str(self.missionoptions.forced_post_flyby_coast))
        self.txtTCM_post_launch.SetValue(str(self.missionoptions.TCM_post_launch))
        self.txtTCM_pre_flyby.SetValue(str(self.missionoptions.TCM_pre_flyby))
        self.txtTCM_maneuver_fraction.SetValue(str(self.missionoptions.TCM_maneuver_fraction))
        self.txtfinal_mass_constraint_LowerBound.SetValue(str(self.missionoptions.final_mass_constraint_bounds[0]))
        self.txtfinal_mass_constraint_UpperBound.SetValue(str(self.missionoptions.final_mass_constraint_bounds[1]))
        self.chkconstrain_dry_mass.SetValue(self.missionoptions.constrain_dry_mass)
        self.chkconstrain_final_mass.SetValue(self.missionoptions.constrain_final_mass)
        self.txtstop_after_journey.SetValue(str(self.missionoptions.stop_after_journey))
        self.txtobjective_journey.SetValue(str(self.missionoptions.objective_journey))
        self.txtwaypoint_file_path.SetValue(self.missionoptions.waypoint_file_path)
        self.txtcovariance_file_path.SetValue(self.missionoptions.covariance_file_path)

        #for MGA-nDSM mission types, show the maximum number of DSMs and the TCM maneuver fraction
        if self.missionoptions.mission_type in [6, 9]:
            self.lblTCM_maneuver_fraction.Show(True)
            self.txtTCM_maneuver_fraction.Show(True)
        else:
            self.lblTCM_maneuver_fraction.Show(False)
            self.txtTCM_maneuver_fraction.Show(False)
        
        if self.missionoptions.objective_type in [8, 9, 23, 26]:
            self.txtobjective_journey.Show(True)
            self.lblobjective_journey.Show(True) 
        else:
            self.txtobjective_journey.Show(False)
            self.lblobjective_journey.Show(False) 

        if self.missionoptions.objective_type in [24]:
            self.txtwaypoint_file_path.Show(True)
            self.lblwaypoint_file_path.Show(True)
            self.btnwaypoint_file_path.Show(True)
            self.txtcovariance_file_path.Show(True)
            self.lblcovariance_file_path.Show(True)
            self.btncovariance_file_path.Show(True)
        else:            
            self.txtwaypoint_file_path.Show(False)
            self.lblwaypoint_file_path.Show(False)
            self.btnwaypoint_file_path.Show(False)
            self.txtcovariance_file_path.Show(False)
            self.lblcovariance_file_path.Show(False)
            self.btncovariance_file_path.Show(False)

        #if global time bounds are off, don't show the bounds fields
        if self.missionoptions.global_timebounded == 1:
            self.lbltotal_flight_time_bounds.Show(True)
            self.txttotal_flight_time_bounds_lower.Show(True)
            self.txttotal_flight_time_bounds_upper.Show(True)
        else:
            self.lbltotal_flight_time_bounds.Show(False)
            self.txttotal_flight_time_bounds_lower.Show(False)
            self.txttotal_flight_time_bounds_upper.Show(False)

        #if the minimum dry mass constraint or the fixed dry mass or final constraint are not active, make the mass constraint field invisible
        #likewise do not show the mass constraint field unless one of those boxes is checked
        if self.missionoptions.constrain_final_mass or self.missionoptions.constrain_dry_mass:
            self.lblfinal_mass_constraint.Show(True)
            self.txtfinal_mass_constraint_LowerBound.Show(True)
            self.txtfinal_mass_constraint_UpperBound.Show(True)
        else:
            self.lblfinal_mass_constraint.Show(False)
            self.txtfinal_mass_constraint_LowerBound.Show(False)
            self.txtfinal_mass_constraint_UpperBound.Show(False)

        self.Layout()
        if platform.system() == 'Windows':
            self.SetupScrolling(scrollToTop=False)

    #event handlers for global mission options    
    def ChangeMissionName(self, e):
        e.Skip()
        self.missionoptions.mission_name = self.txtMissionName.GetValue().replace(' ','_')
        self.parent.update()
        
        
    def ChangeMissionType(self, e):
        self.missionoptions.mission_type = self.cmbMissionType.GetSelection()

        if self.missionoptions.mission_type != 9:
            for journey in self.missionoptions.Journeys:
                journey.phase_type = self.missionoptions.mission_type       
                
        self.missionoptions.DisassembleMasterDecisionVector()
        self.missionoptions.ConvertDecisionVector()
        self.missionoptions.AssembleMasterDecisionVector()

        self.parent.update()

    def Changeobjective_type(self, e):
        self.missionoptions.objective_type = self.cmbobjective_type.GetSelection()
        self.parent.update()

    def Changeinclude_initial_impulse_in_cost(self, e):
        self.missionoptions.include_initial_impulse_in_cost = int(self.chkinclude_initial_impulse_in_cost.GetValue())
    
    def Changeobjective_journey(self,e):
        e.Skip()
        self.missionoptions.objective_journey = eval(self.txtobjective_journey.GetValue())    
        self.parent.update()
    
    def Changewaypoint_file_path(self,e):
        e.Skip()
        self.missionoptions.waypoint_file_path = self.txtwaypoint_file_path.GetValue()    
        self.parent.update()

    def Changecovariance_file_path(self,e):
        e.Skip()
        self.missionoptions.covariance_file_path = self.txtcovariance_file_path.GetValue()    
        self.parent.update()
        
    def Clickwaypoint_file_path_button(self, e):
        #call dialog to choose destination list
        dlg = wx.FileDialog(self, "Choose a waypoint file", defaultDir=self.parent.Parent.dirname, wildcard="*.ephemeris", style=wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filepath = dlg.GetPath() + '/' + dlg.GetFilename()
        else:
            filepath = self.missionoptions.waypoint_file_path

        dlg.Destroy()

        self.missionoptions.waypoint_file_path = filepath

        self.update()

    def Clickcovariance_file_path_button(self, e):
        #call dialog to choose destination list
        dlg = wx.FileDialog(self, "Choose a covariance file", defaultDir=self.parent.Parent.dirname, wildcard="*.ephemeris", style=wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filepath = dlg.GetPath() + '/' + dlg.GetFilename()
        else:
            filepath = self.missionoptions.covariance_file_path

        dlg.Destroy()

        self.missionoptions.covariance_file_path = filepath

        self.update()

    def Changelaunch_window_open_date(self, e):
        e.Skip()
        
        dateString = self.txtlaunch_window_open_date.GetValue()

        from timeUtilities import stringToJD

        self.missionoptions.launch_window_open_date = stringToJD(dateString, self.missionoptions.universe_folder)

        self.update()

    def ChangeLaunchDateCalendar(self, e):
        CurrentLaunchDate = self.LaunchDateCalendar.GetDate()
        CurrentLaunchDate = CurrentLaunchDate.FromTimezone(wx.DateTime.TimeZone(offset=0))
        self.missionoptions.launch_window_open_date = CurrentLaunchDate.GetMJD()
        self.update()

    def Changenum_timesteps(self, e):
        e.Skip()
        self.missionoptions.num_timesteps = eval(self.txtnum_timesteps.GetValue())
        self.update()

    def ChangeDLA_bounds_lower(self, e):
        e.Skip()
        self.missionoptions.DLA_bounds[0] = eval(self.txtDLA_bounds_lower.GetValue())
        self.update()
        
    def ChangeDLA_bounds_upper(self, e):
        e.Skip()
        self.missionoptions.DLA_bounds[1] = eval(self.txtDLA_bounds_upper.GetValue())
        self.update()

    def ChangeRLA_bounds_lower(self, e):
        e.Skip()
        self.missionoptions.RLA_bounds[0] = eval(self.txtRLA_bounds_lower.GetValue())
        self.update()
        
    def ChangeRLA_bounds_upper(self, e):
        e.Skip()
        self.missionoptions.RLA_bounds[1] = eval(self.txtRLA_bounds_upper.GetValue())
        self.update()

    def Changeglobal_timebounded(self, e):
        self.missionoptions.global_timebounded = int(self.chkglobal_timebounded.GetValue())
        self.parent.update()

    def Changetotal_flight_time_bounds_lower(self, e):
        e.Skip()
        self.missionoptions.total_flight_time_bounds[0] = eval(self.txttotal_flight_time_bounds_lower.GetValue())
        self.update()

    def Changetotal_flight_time_bounds_upper(self, e):
        e.Skip()
        self.missionoptions.total_flight_time_bounds[1] = eval(self.txttotal_flight_time_bounds_upper.GetValue())
        self.update()

    def Changeforced_post_launch_coast(self, e):
        e.Skip()
        self.missionoptions.forced_post_launch_coast = eval(self.txtforced_post_launch_coast.GetValue())
        self.update()

    def Changeforced_pre_flyby_coast(self, e):
        e.Skip()
        self.missionoptions.forced_pre_flyby_coast = eval(self.txtforced_pre_flyby_coast.GetValue())
        self.update()
    
    def Changeforced_post_flyby_coast(self, e):
        e.Skip()
        self.missionoptions.forced_post_flyby_coast = eval(self.txtforced_post_flyby_coast.GetValue())
        self.update()

    def ChangeTCM_post_launch(self, e):
        e.Skip()
        self.missionoptions.TCM_post_launch = eval(self.txtTCM_post_launch.GetValue())
        self.update()

    def ChangeTCM_pre_flyby(self, e):
        e.Skip()
        self.missionoptions.TCM_pre_flyby = eval(self.txtTCM_pre_flyby.GetValue())
        self.update()

    def ChangeTCM_maneuver_fraction(self, e):
        e.Skip()
        self.missionoptions.TCM_maneuver_fraction = eval(self.txtTCM_maneuver_fraction.GetValue())
        self.update()

    def Changefinal_mass_constraint_LowerBound(self, e):
        e.Skip()
        self.missionoptions.final_mass_constraint_bounds[0] = eval(self.txtfinal_mass_constraint_LowerBound.GetValue())
        self.update()

    def Changefinal_mass_constraint_UpperBound(self, e):
        e.Skip()
        self.missionoptions.final_mass_constraint_bounds[1] = eval(self.txtfinal_mass_constraint_UpperBound.GetValue())
        self.update()
        
    def Changeconstrain_final_mass(self, e):
        e.Skip()
        self.missionoptions.constrain_final_mass = self.chkconstrain_final_mass.GetValue()
        self.update()

    def Changeconstrain_dry_mass(self, e):
        e.Skip()
        self.missionoptions.constrain_dry_mass = self.chkconstrain_dry_mass.GetValue()
        self.update()

    def Changestop_after_journey(self, e):
        e.Skip()
        self.missionoptions.stop_after_journey = int(self.txtstop_after_journey.GetValue())
        self.update()