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

class OutputOptionsPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, missionoptions):
        self.missionoptions = missionoptions
        self.parent = parent

        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)

        self.mainbox = wx.FlexGridSizer(30,2,5,5)

        self.lblprint_only_non_default_options = wx.StaticText(self, -1, "Print only non-default options to .emtgopt file?")
        self.chkprint_only_non_default_options = wx.CheckBox(self, -1)
        
        self.lbloutput_file_frame = wx.StaticText(self, -1, "Output file frame")
        output_file_frame_choices = ['ICRF', 'J2000_BCI', 'J2000_BCF', 'TrueOfDate_BCI', 'TrueOfDate_BCF', 'Principle Axes', 'Topocentric', 'Polar']
        self.cmboutput_file_frame = wx.ComboBox(self, -1, style=wx.CB_READONLY, choices=output_file_frame_choices)

        self.lbloutput_dormant_journeys = wx.StaticText(self, -1, "Output journey entries for wait times at intermediate and final target?")
        self.chkoutput_dormant_journeys = wx.CheckBox(self, -1)

        self.lblpost_mission_wait_time = wx.StaticText(self, -1, "Stay time at the final target")
        self.txtpost_mission_wait_time = wx.TextCtrl(self, -1, "post_mission_wait_time")

        self.lbloverride_working_directory = wx.StaticText(self, -1, "Override working directory?")
        self.chkoverride_working_directory = wx.CheckBox(self, -1)

        self.lblforced_working_directory = wx.StaticText(self, -1, "Working directory")
        self.txtforced_working_directory = wx.TextCtrl(self, -1, "forced_working_directory", size=(600,-1))
        self.btnforced_working_directory = wx.Button(self, -1, "...")
        working_directory_sizer = wx.BoxSizer(wx.HORIZONTAL)
        working_directory_sizer.AddMany([self.txtforced_working_directory, self.btnforced_working_directory])

        self.lbloverride_mission_subfolder = wx.StaticText(self, -1, "Override mission_subfolder?")
        self.chkoverride_mission_subfolder = wx.CheckBox(self, -1)

        self.lblforced_mission_subfolder = wx.StaticText(self, -1, "Mission Subfolder")
        self.txtforced_mission_subfolder = wx.TextCtrl(self, -1, "forced_mission_subfolder", size=(600,-1))
        self.btnforced_mission_subfolder = wx.Button(self, -1, "...")
        mission_subfolder_sizer = wx.BoxSizer(wx.HORIZONTAL)
        mission_subfolder_sizer.AddMany([self.txtforced_mission_subfolder, self.btnforced_mission_subfolder])

        self.lblshort_output_file_names = wx.StaticText(self, -1, "Shorten output file names?")
        self.chkshort_output_file_names = wx.CheckBox(self, -1)

        self.lblgenerate_forward_integrated_ephemeris = wx.StaticText(self, -1, "Generate forward-integrated ephemeris")
        self.chkgenerate_forward_integrated_ephemeris = wx.CheckBox(self, -1)

        self.lblforward_integrated_ephemeris_central_body_SPICE_ID = wx.StaticText(self, -1, "SPICE ID of central body for forward integrated ephemeris")
        self.txtforward_integrated_ephemeris_central_body_SPICE_ID = wx.TextCtrl(self, -1, "forward_integrated_ephemeris_central_body_SPICE_ID")

        self.lbladd_control_switch_line_to_ephemeris = wx.StaticText(self, -1, "Add duplicate control switch lines to ephemeris?")
        self.chkadd_control_switch_line_to_ephemeris = wx.CheckBox(self, -1)

        self.lblappend_mass_to_ephemeris_output = wx.StaticText(self, -1, "Append mass to forward integrated ephemeris?")
        self.chkappend_mass_to_ephemeris_output = wx.CheckBox(self, -1)

        self.lblappend_control_to_ephemeris_output = wx.StaticText(self, -1, "Append control to forward integrated ephemeris?")
        self.chkappend_control_to_ephemeris_output = wx.CheckBox(self, -1)
        
        self.lblappend_thrust_to_ephemeris_output = wx.StaticText(self, -1, "Append thrust to forward integrated ephemeris?")
        self.chkappend_thrust_to_ephemeris_output = wx.CheckBox(self, -1)

        self.lblappend_mdot_to_ephemeris_output = wx.StaticText(self, -1, "Append mass flow rate to forward integrated ephemeris?")
        self.chkappend_mdot_to_ephemeris_output = wx.CheckBox(self, -1)

        self.lblappend_Isp_to_ephemeris_output = wx.StaticText(self, -1, "Append Isp to forward integrated ephemeris?")
        self.chkappend_Isp_to_ephemeris_output = wx.CheckBox(self, -1)

        self.lblappend_number_of_active_engines_to_ephemeris_output = wx.StaticText(self, -1, "Append number of active engines to forward integrated ephemeris?")
        self.chkappend_number_of_active_engines_to_ephemeris_output = wx.CheckBox(self, -1)
        
        self.lblappend_throttle_level_to_ephemeris_output = wx.StaticText(self, -1, "Append throttle level to forward integrated ephemeris?")
        self.chkappend_throttle_level_to_ephemeris_output = wx.CheckBox(self, -1)         
        
        self.lblappend_active_power_to_ephemeris_output = wx.StaticText(self, -1, "Append active_power to forward integrated ephemeris?")
        self.chkappend_active_power_to_ephemeris_output = wx.CheckBox(self, -1)
        
        self.lblspacecraft_SPICE_ID = wx.StaticText(self, -1, "Spacecraft SPICE ID")
        self.txtspacecraft_SPICE_ID = wx.TextCtrl(self, -1, "spacecraft_SPICE_ID", size=(300,-1))
        
        self.lblpyemtg_path = wx.StaticText(self, -1, "Path to PyEMTG")
        self.txtpyemtg_path = wx.TextCtrl(self, -1, "pyemtg_path", size=(600,-1))
        self.btnpyemtg_path = wx.Button(self, -1, "...")
        PyEMTG_path_sizer = wx.BoxSizer(wx.HORIZONTAL)
        PyEMTG_path_sizer.AddMany([self.txtpyemtg_path, self.btnpyemtg_path])
        
        self.lblspice_utility_extension = wx.StaticText(self, -1, "File extension for SPICE utilities")
        self.txtspice_utility_extension = wx.TextCtrl(self, -1, "spice_utility_extension", size=(300,-1))      
                
        self.lblspice_utilities_path = wx.StaticText(self, -1, "Path to SPICE utilities (brief, mkspk, etc.)")
        self.txtspice_utilities_path = wx.TextCtrl(self, -1, "spice_utilities_path", size=(600,-1))
        self.btnspice_utilities_path = wx.Button(self, -1, "...")
        spice_utilities_path_sizer = wx.BoxSizer(wx.HORIZONTAL)
        spice_utilities_path_sizer.AddMany([self.txtspice_utilities_path, self.btnspice_utilities_path])  

        self.lblcall_system_to_generate_bsp = wx.StaticText(self, -1, "Perform a system call to write a binary SPICE kernel?")
        self.chkcall_system_to_generate_bsp = wx.CheckBox(self, -1)
        
        self.lblbackground_mode = wx.StaticText(self, -1, "Enable background mode")
        self.chkbackground_mode = wx.CheckBox(self, -1)

        self.lbloutput_STMs = wx.StaticText(self, -1, "Write STMs?")
        self.chkoutput_STMs = wx.CheckBox(self, -1)

        self.lbloutput_maneuver_and_target_spec_files = wx.StaticText(self, -1, "Write maneuver and target spec files?")
        self.chkoutput_maneuver_and_target_spec_files = wx.CheckBox(self, -1)

        self.mainbox.AddMany([  self.lblprint_only_non_default_options, self.chkprint_only_non_default_options,
                                self.lbloutput_file_frame, self.cmboutput_file_frame,
                                self.lbloutput_dormant_journeys, self.chkoutput_dormant_journeys,
                                self.lblpost_mission_wait_time, self.txtpost_mission_wait_time,
                                self.lbloverride_working_directory, self.chkoverride_working_directory,
                                self.lblforced_working_directory, working_directory_sizer,
                                self.lbloverride_mission_subfolder, self.chkoverride_mission_subfolder,
                                self.lblforced_mission_subfolder, mission_subfolder_sizer,
                                self.lblshort_output_file_names, self.chkshort_output_file_names,
                                self.lblgenerate_forward_integrated_ephemeris, self.chkgenerate_forward_integrated_ephemeris,
                                self.lblforward_integrated_ephemeris_central_body_SPICE_ID, self.txtforward_integrated_ephemeris_central_body_SPICE_ID,
                                self.lbladd_control_switch_line_to_ephemeris, self.chkadd_control_switch_line_to_ephemeris,
                                self.lblappend_mass_to_ephemeris_output, self.chkappend_mass_to_ephemeris_output,
                                self.lblappend_control_to_ephemeris_output, self.chkappend_control_to_ephemeris_output,
                                self.lblappend_thrust_to_ephemeris_output, self.chkappend_thrust_to_ephemeris_output,
                                self.lblappend_mdot_to_ephemeris_output, self.chkappend_mdot_to_ephemeris_output,
                                self.lblappend_Isp_to_ephemeris_output, self.chkappend_Isp_to_ephemeris_output,
                                self.lblappend_number_of_active_engines_to_ephemeris_output, self.chkappend_number_of_active_engines_to_ephemeris_output,
                                self.lblappend_active_power_to_ephemeris_output, self.chkappend_active_power_to_ephemeris_output,
                                self.lblappend_throttle_level_to_ephemeris_output, self.chkappend_throttle_level_to_ephemeris_output,
                                self.lblspacecraft_SPICE_ID, self.txtspacecraft_SPICE_ID,
                                self.lblpyemtg_path, PyEMTG_path_sizer,
                                self.lblspice_utilities_path, spice_utilities_path_sizer,
                                self.lblspice_utility_extension, self.txtspice_utility_extension,
                                self.lblcall_system_to_generate_bsp, self.chkcall_system_to_generate_bsp,
                                self.lbloutput_STMs, self.chkoutput_STMs,
                                self.lbloutput_maneuver_and_target_spec_files, self.chkoutput_maneuver_and_target_spec_files,
                                self.lblbackground_mode, self.chkbackground_mode])
        
        self.presentbox = wx.FlexGridSizer(1,2,5,5)
        self.presentbox.AddMany([self.mainbox])

        self.SetSizer(self.presentbox)
        self.SetupScrolling()

        #bindings
        self.chkprint_only_non_default_options.Bind(wx.EVT_CHECKBOX, self.Changeprint_only_non_default_options)
        self.cmboutput_file_frame.Bind(wx.EVT_COMBOBOX, self.Changeoutput_file_frame)
        self.chkoutput_dormant_journeys.Bind(wx.EVT_CHECKBOX, self.Changeoutput_dormant_journeys)
        self.txtpost_mission_wait_time.Bind(wx.EVT_KILL_FOCUS, self.Changepost_mission_wait_time)
        self.chkoverride_working_directory.Bind(wx.EVT_CHECKBOX, self.Changeoverride_working_directory)
        self.txtforced_working_directory.Bind(wx.EVT_KILL_FOCUS, self.Changeforced_working_directory)
        self.btnforced_working_directory.Bind(wx.EVT_BUTTON, self.Clickforced_working_directory_button)
        self.chkoverride_mission_subfolder.Bind(wx.EVT_CHECKBOX, self.Changeoverride_mission_subfolder)
        self.txtforced_mission_subfolder.Bind(wx.EVT_KILL_FOCUS, self.Changeforced_mission_subfolder)
        self.btnforced_mission_subfolder.Bind(wx.EVT_BUTTON, self.Clickforced_mission_subfolder_button)
        self.chkshort_output_file_names.Bind(wx.EVT_CHECKBOX, self.Changeshort_output_file_names)
        self.chkgenerate_forward_integrated_ephemeris.Bind(wx.EVT_CHECKBOX, self.Changegenerate_forward_integrated_ephemeris)
        self.txtforward_integrated_ephemeris_central_body_SPICE_ID.Bind(wx.EVT_KILL_FOCUS, self.Changeforward_integrated_ephemeris_central_body_SPICE_ID)
        self.chkadd_control_switch_line_to_ephemeris.Bind(wx.EVT_CHECKBOX, self.Changeadd_control_switch_line_to_ephemeris)
        self.chkappend_mass_to_ephemeris_output.Bind(wx.EVT_CHECKBOX, self.Changeappend_mass_to_ephemeris_output)
        self.chkappend_control_to_ephemeris_output.Bind(wx.EVT_CHECKBOX, self.Changeappend_control_to_ephemeris_output)
        self.chkappend_thrust_to_ephemeris_output.Bind(wx.EVT_CHECKBOX, self.Changeappend_thrust_to_ephemeris_output)
        self.chkappend_mdot_to_ephemeris_output.Bind(wx.EVT_CHECKBOX, self.Changeappend_mdot_to_ephemeris_output)
        self.chkappend_Isp_to_ephemeris_output.Bind(wx.EVT_CHECKBOX, self.Changeappend_Isp_to_ephemeris_output)
        self.chkappend_number_of_active_engines_to_ephemeris_output.Bind(wx.EVT_CHECKBOX, self.Changeappend_number_of_active_engines_to_ephemeris_output)
        self.chkappend_active_power_to_ephemeris_output.Bind(wx.EVT_CHECKBOX, self.Changeappend_active_power_to_ephemeris_output)
        self.chkappend_throttle_level_to_ephemeris_output.Bind(wx.EVT_CHECKBOX, self.Changeappend_throttle_level_to_ephemeris_output)

        self.txtspacecraft_SPICE_ID.Bind(wx.EVT_KILL_FOCUS, self.Changespacecraft_SPICE_ID)
        self.txtpyemtg_path.Bind(wx.EVT_KILL_FOCUS, self.Changepyemtg_path)
        self.btnpyemtg_path.Bind(wx.EVT_BUTTON, self.Clickpyemtg_path_button)
        self.txtspice_utility_extension.Bind(wx.EVT_KILL_FOCUS, self.Changespice_utility_extension)
        self.txtspice_utilities_path.Bind(wx.EVT_KILL_FOCUS, self.Changespice_utilities_path)
        self.btnspice_utilities_path.Bind(wx.EVT_BUTTON, self.Clickspice_utilities_path_button)
        self.chkcall_system_to_generate_bsp.Bind(wx.EVT_CHECKBOX, self.Changecall_system_to_generate_bsp)

        self.chkoutput_STMs.Bind(wx.EVT_CHECKBOX, self.Changeoutput_STMs)
        self.chkoutput_maneuver_and_target_spec_files.Bind(wx.EVT_CHECKBOX, self.Changeoutput_maneuver_and_target_spec_files)
        self.chkbackground_mode.Bind(wx.EVT_CHECKBOX, self.Changebackground_mode)

    def update(self):
        import wx
        self.chkprint_only_non_default_options.SetValue(self.missionoptions.print_only_non_default_options)
        self.cmboutput_file_frame.SetSelection(self.missionoptions.output_file_frame)
        self.chkoutput_dormant_journeys.SetValue(self.missionoptions.output_dormant_journeys)
        self.txtpost_mission_wait_time.SetValue(str(self.missionoptions.post_mission_wait_time))
        self.chkoverride_working_directory.SetValue(self.missionoptions.override_working_directory)
        self.chkoverride_mission_subfolder.SetValue(self.missionoptions.override_mission_subfolder)
        self.txtforced_working_directory.SetValue(self.missionoptions.forced_working_directory)
        self.chkshort_output_file_names.SetValue(self.missionoptions.short_output_file_names)
        self.chkgenerate_forward_integrated_ephemeris.SetValue(self.missionoptions.generate_forward_integrated_ephemeris)
        self.txtforward_integrated_ephemeris_central_body_SPICE_ID.SetValue(str(self.missionoptions.forward_integrated_ephemeris_central_body_SPICE_ID))
        self.chkadd_control_switch_line_to_ephemeris.SetValue(self.missionoptions.add_control_switch_line_to_ephemeris)
        self.chkappend_mass_to_ephemeris_output.SetValue(self.missionoptions.append_mass_to_ephemeris_output)
        self.chkappend_control_to_ephemeris_output.SetValue(self.missionoptions.append_control_to_ephemeris_output)
        self.chkappend_thrust_to_ephemeris_output.SetValue(self.missionoptions.append_thrust_to_ephemeris_output)
        self.chkappend_mdot_to_ephemeris_output.SetValue(self.missionoptions.append_mdot_to_ephemeris_output)
        self.chkappend_Isp_to_ephemeris_output.SetValue(self.missionoptions.append_Isp_to_ephemeris_output)
        self.chkappend_number_of_active_engines_to_ephemeris_output.SetValue(self.missionoptions.append_number_of_active_engines_to_ephemeris_output)
        self.chkappend_active_power_to_ephemeris_output.SetValue(self.missionoptions.append_active_power_to_ephemeris_output)
        self.chkappend_throttle_level_to_ephemeris_output.SetValue(self.missionoptions.append_throttle_level_to_ephemeris_output)
                
        self.txtspacecraft_SPICE_ID.SetValue(str(self.missionoptions.spacecraft_SPICE_ID))
        self.txtpyemtg_path.SetValue(self.missionoptions.pyemtg_path)
        self.txtspice_utilities_path.SetValue(self.missionoptions.spice_utilities_path)
        self.txtspice_utility_extension.SetValue(self.missionoptions.spice_utility_extension)
        self.chkcall_system_to_generate_bsp.SetValue(self.missionoptions.call_system_to_generate_bsp)

        self.chkoutput_STMs.SetValue(self.missionoptions.output_STMs)
        self.chkoutput_maneuver_and_target_spec_files.SetValue(self.missionoptions.output_maneuver_and_target_spec_files)
        self.chkbackground_mode.SetValue(self.missionoptions.background_mode)

        if self.missionoptions.generate_forward_integrated_ephemeris:
            self.lblforward_integrated_ephemeris_central_body_SPICE_ID.Show(True)
            self.txtforward_integrated_ephemeris_central_body_SPICE_ID.Show(True)
            self.lbladd_control_switch_line_to_ephemeris.Show(True)
            self.chkadd_control_switch_line_to_ephemeris.Show(True)
            self.lblappend_mass_to_ephemeris_output.Show(True)
            self.chkappend_mass_to_ephemeris_output.Show(True)
            self.lblappend_control_to_ephemeris_output.Show(True)
            self.chkappend_control_to_ephemeris_output.Show(True)
            self.lblappend_thrust_to_ephemeris_output.Show(True)
            self.chkappend_thrust_to_ephemeris_output.Show(True)
            self.lblappend_mdot_to_ephemeris_output.Show(True)
            self.chkappend_mdot_to_ephemeris_output.Show(True)
            self.lblappend_Isp_to_ephemeris_output.Show(True)
            self.chkappend_Isp_to_ephemeris_output.Show(True)
            self.lblappend_number_of_active_engines_to_ephemeris_output.Show(True)
            self.chkappend_number_of_active_engines_to_ephemeris_output.Show(True)
            self.lblappend_active_power_to_ephemeris_output.Show(True)
            self.chkappend_active_power_to_ephemeris_output.Show(True)
            self.lblappend_throttle_level_to_ephemeris_output.Show(True)
            self.chkappend_throttle_level_to_ephemeris_output.Show(True)
            self.lblspacecraft_SPICE_ID.Show(True)
            self.txtspacecraft_SPICE_ID.Show(True)        
            self.lblpyemtg_path.Show(True)
            self.txtpyemtg_path.Show(True)
            self.btnpyemtg_path.Show(True)        
            self.lblspice_utility_extension.Show(True)
            self.txtspice_utility_extension.Show(True)                    
            self.lblspice_utilities_path.Show(True)
            self.txtspice_utilities_path.Show(True)
            self.btnspice_utilities_path.Show(True)
            self.lblcall_system_to_generate_bsp.Show(True)
            self.chkcall_system_to_generate_bsp.Show(True)
        else:
            self.lblforward_integrated_ephemeris_central_body_SPICE_ID.Show(False)
            self.txtforward_integrated_ephemeris_central_body_SPICE_ID.Show(False)
            self.lbladd_control_switch_line_to_ephemeris.Show(False)
            self.chkadd_control_switch_line_to_ephemeris.Show(False)
            self.lblappend_mass_to_ephemeris_output.Show(False)
            self.chkappend_mass_to_ephemeris_output.Show(False)
            self.lblappend_control_to_ephemeris_output.Show(False)
            self.chkappend_control_to_ephemeris_output.Show(False)
            self.lblappend_thrust_to_ephemeris_output.Show(False)
            self.chkappend_thrust_to_ephemeris_output.Show(False)
            self.lblappend_mdot_to_ephemeris_output.Show(False)
            self.chkappend_mdot_to_ephemeris_output.Show(False)
            self.lblappend_Isp_to_ephemeris_output.Show(False)
            self.chkappend_Isp_to_ephemeris_output.Show(False)
            self.lblappend_number_of_active_engines_to_ephemeris_output.Show(False)
            self.chkappend_number_of_active_engines_to_ephemeris_output.Show(False)
            self.lblappend_active_power_to_ephemeris_output.Show(False)
            self.chkappend_active_power_to_ephemeris_output.Show(False)
            self.lblappend_throttle_level_to_ephemeris_output.Show(False)
            self.chkappend_throttle_level_to_ephemeris_output.Show(False)
            self.lblspacecraft_SPICE_ID.Show(False)
            self.txtspacecraft_SPICE_ID.Show(False)        
            self.lblpyemtg_path.Show(False)
            self.txtpyemtg_path.Show(False)
            self.btnpyemtg_path.Show(False)        
            self.lblspice_utility_extension.Show(False)
            self.txtspice_utility_extension.Show(False)                    
            self.lblspice_utilities_path.Show(False)
            self.txtspice_utilities_path.Show(False)
            self.btnspice_utilities_path.Show(False)
            self.lblcall_system_to_generate_bsp.Show(False)
            self.chkcall_system_to_generate_bsp.Show(False)

        if self.missionoptions.output_dormant_journeys:
            self.lblpost_mission_wait_time.Show(True)
            self.txtpost_mission_wait_time.Show(True)
        else:
            self.lblpost_mission_wait_time.Show(False)
            self.txtpost_mission_wait_time.Show(False)

        if self.missionoptions.override_working_directory:
            self.lblforced_working_directory.Show(True)
            self.txtforced_working_directory.Show(True)
            self.btnforced_working_directory.Show(True)
        else:
            self.lblforced_working_directory.Show(False)
            self.txtforced_working_directory.Show(False)
            self.btnforced_working_directory.Show(False)
        
        if self.missionoptions.override_mission_subfolder:
            self.lblforced_mission_subfolder.Show(True)
            self.txtforced_mission_subfolder.Show(True)
            self.btnforced_mission_subfolder.Show(True)
        else:
            self.lblforced_mission_subfolder.Show(False)
            self.txtforced_mission_subfolder.Show(False)
            self.btnforced_mission_subfolder.Show(False)

        #re-size the panel
        self.Layout()
        if platform.system() == 'Windows':
            self.SetupScrolling(scrollToTop=False)


    #handlers for output options

    def Changeprint_only_non_default_options(self, e):
        e.Skip()
        self.missionoptions.print_only_non_default_options = int(self.chkprint_only_non_default_options.GetValue())

    def Changeoutput_file_frame(self, e):
        e.Skip()
        self.missionoptions.output_file_frame = int(self.cmboutput_file_frame.GetSelection())
        self.update()

    def Changeoutput_dormant_journeys(self, e):
        e.Skip()
        self.missionoptions.output_dormant_journeys = int(self.chkoutput_dormant_journeys.GetValue())
        self.update()

    def Changepost_mission_wait_time(self, e):
        e.Skip()
        self.missionoptions.post_mission_wait_time = eval(self.txtpost_mission_wait_time.GetValue())

    def Changeoverride_working_directory(self, e):
        e.Skip()
        self.missionoptions.override_working_directory = int(self.chkoverride_working_directory.GetValue())
        self.update()

    def Changeforced_working_directory(self, e):
        e.Skip()
        self.missionoptions.forced_working_directory = self.txtforced_working_directory.GetValue()

    def Changeoverride_mission_subfolder(self, e):
        e.Skip()
        self.missionoptions.override_mission_subfolder = int(self.chkoverride_mission_subfolder.GetValue())
        self.update()

    def Changeforced_mission_subfolder(self, e):
        e.Skip()
        self.missionoptions.forced_mission_subfolder = self.txtforced_mission_subfolder.GetValue()

    def Clickforced_working_directory_button(self, e):
        e.Skip()
        #file load dialog to get name of working directory
        dlg = wx.DirDialog(self, "Choose a working directory", self.parent.Parent.dirname)
        if dlg.ShowModal() == wx.ID_OK:
            self.missionoptions.forced_working_directory = dlg.GetPath()
            self.txtforced_working_directory.SetValue(self.missionoptions.forced_working_directory)
        dlg.Destroy()
        
    def Clickforced_mission_subfolder_button(self, e):
        e.Skip()
        #file load dialog to get name of working directory
        dlg = wx.DirDialog(self, "Choose a mission subfolder", self.parent.Parent.dirname)
        if dlg.ShowModal() == wx.ID_OK:
            self.missionoptions.forced_mission_subfolder = dlg.GetPath()
            self.txtforced_mission_subfolder.SetValue(self.missionoptions.forced_mission_subfolder)
        dlg.Destroy()

    def Changeshort_output_file_names(self, e):
        e.Skip()
        self.missionoptions.short_output_file_names = int(self.chkshort_output_file_names.GetValue())

    def Changegenerate_forward_integrated_ephemeris(self, e):
        self.missionoptions.generate_forward_integrated_ephemeris = int(self.chkgenerate_forward_integrated_ephemeris.GetValue())
        self.update()

    def Changeforward_integrated_ephemeris_central_body_SPICE_ID(self, e):
        self.missionoptions.forward_integrated_ephemeris_central_body_SPICE_ID = int(self.txtforward_integrated_ephemeris_central_body_SPICE_ID.GetValue())
        self.update()

    def Changeadd_control_switch_line_to_ephemeris(self, e):
        e.Skip()
        self.missionoptions.add_control_switch_line_to_ephemeris = int(self.chkadd_control_switch_line_to_ephemeris.GetValue())

    def Changeappend_mass_to_ephemeris_output(self, e):
        e.Skip()
        self.missionoptions.append_mass_to_ephemeris_output = int(self.chkappend_mass_to_ephemeris_output.GetValue())

    def Changeappend_control_to_ephemeris_output(self, e):
        e.Skip()
        self.missionoptions.append_control_to_ephemeris_output = int(self.chkappend_control_to_ephemeris_output.GetValue())

    def Changeappend_thrust_to_ephemeris_output(self, e):
        e.Skip()
        self.missionoptions.append_thrust_to_ephemeris_output = int(self.chkappend_thrust_to_ephemeris_output.GetValue())

    def Changeappend_mdot_to_ephemeris_output(self, e):
        e.Skip()
        self.missionoptions.append_mdot_to_ephemeris_output = int(self.chkappend_mdot_to_ephemeris_output.GetValue())

    def Changeappend_Isp_to_ephemeris_output(self, e):
        e.Skip()
        self.missionoptions.append_Isp_to_ephemeris_output = int(self.chkappend_Isp_to_ephemeris_output.GetValue())

    def Changeappend_number_of_active_engines_to_ephemeris_output(self, e):
        e.Skip()
        self.missionoptions.append_number_of_active_engines_to_ephemeris_output = int(self.chkappend_number_of_active_engines_to_ephemeris_output.GetValue())
        
    def Changeappend_active_power_to_ephemeris_output(self, e):
        e.Skip()
        self.missionoptions.append_active_power_to_ephemeris_output = int(self.chkappend_active_power_to_ephemeris_output.GetValue())
        
    def Changeappend_throttle_level_to_ephemeris_output(self, e):
        e.Skip()
        self.missionoptions.append_throttle_level_to_ephemeris_output = int(self.chkappend_throttle_level_to_ephemeris_output.GetValue())

    def Changespacecraft_SPICE_ID(self, e):
        e.Skip()
        self.missionoptions.spacecraft_SPICE_ID = int(self.txtspacecraft_SPICE_ID.GetValue())

    def Changepyemtg_path(self, e):
        e.Skip()
        self.missionoptions.pyemtg_path = self.txtpyemtg_path.GetValue()

    def Clickpyemtg_path_button(self, e):
        e.Skip()
        dlg = wx.DirDialog(self, "Choose a path to PyEMTG", self.missionoptions.pyemtg_path)
        if dlg.ShowModal() == wx.ID_OK:
            self.missionoptions.pyemtg_path = dlg.GetPath().replace('\\','/') + '/'
            self.txtpyemtg_path.SetValue(self.missionoptions.pyemtg_path)
        dlg.Destroy()
    
    def Changespice_utility_extension(self, e):
        e.Skip()
        self.missionoptions.spice_utility_extension = self.txtspice_utility_extension.GetValue()
        
    def Changespice_utilities_path(self, e):
        e.Skip()
        self.missionoptions.spice_utilities_path = self.txtspice_utilities_path.GetValue()

    def Clickspice_utilities_path_button(self, e):
        e.Skip()
        dlg = wx.DirDialog(self, "Choose a path to SPICE utilites", self.missionoptions.spice_utilities_path)
        if dlg.ShowModal() == wx.ID_OK:
            self.missionoptions.spice_utilities_path = dlg.GetPath().replace('\\','/') + '/'
            self.txtspice_utilities_path.SetValue(self.missionoptions.spice_utilities_path)
        dlg.Destroy()

    def Changecall_system_to_generate_bsp(self, e):
        e.Skip()
        self.missionoptions.call_system_to_generate_bsp = int(self.chkcall_system_to_generate_bsp.GetValue())

    def Changeoutput_STMs(self, e):
        e.Skip()
        self.missionoptions.output_STMs = int(self.chkoutput_STMs.GetValue())
        
    def Changeoutput_maneuver_and_target_spec_files(self, e):
        e.Skip()
        self.missionoptions.output_maneuver_and_target_spec_files = int(self.chkoutput_maneuver_and_target_spec_files.GetValue())
                    
    def Changebackground_mode(self, e):
        e.Skip()
        self.missionoptions.background_mode = int(self.chkbackground_mode.GetValue())