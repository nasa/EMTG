import Mission
import PlotOptions
import BubbleOptions
import wx
import wx.lib.scrolledpanel
import webbrowser
import os
import PyEMTG_interface

class MissionPanel(wx.lib.scrolledpanel.ScrolledPanel):
    #mission panel
    #contains post-processing options

    def __init__(self, parent, mission):
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)

        self.mission = mission
        self.plotoptions = PlotOptions.PlotOptions(throttletablefile=parent.default_thruster_file,
                                                  default_universe_path=parent.default_universe_path,
                                                  de_file=parent.de_file,
                                                  leapseconds_file=parent.leapseconds_file)
        self.bubbleoptions = BubbleOptions.BubbleOptions(parent.default_small_bodies_file)

        #Journey selection listbox
        #first we need to create an array of journey names
        self.journeynamelist = []
        for journey in self.mission.Journeys:
            self.journeynamelist.append(journey.journey_name)
        
        if len(self.journeynamelist) > 1:
            identicaljourneys = True
            for OtherJourney in self.mission.Journeys[1:len(self.mission.Journeys)]:
                if OtherJourney.central_body.lower() != self.mission.Journeys[0].central_body.lower():
                    identicaljourneys = False
            if identicaljourneys:
                self.journeynamelist.append("All journeys")

        self.mission.ChosenJourneys = [0]
        self.mission.ActiveJourney = 0

        self.lblJourneyList = wx.StaticText(self, -1, "Choose a journey")
        self.JourneyListBox = wx.ListBox(self, -1, choices = self.journeynamelist, size=(300, -1), style=wx.LB_EXTENDED)

        #widgets for trajectory plot
        self.TrajectoryPlotBox = wx.StaticBox(self, -1, "Trajectory plot options")
        self.chkShowPlotAxes = wx.CheckBox(self, -1, "Show plot axes", size=(300, -1))
        self.chkShowBoundaryOrbits = wx.CheckBox(self, -1, "Show boundary orbits", size=(300, -1))
        self.chkShowFreePointBoundaryOrbits = wx.CheckBox(self, -1, "Show free point boundary orbits", size=(300,-1))
        self.chkShowCentralBodyOrbits = wx.CheckBox(self, -1, "Show central body orbits", size=(300, -1))
        self.chkShowMissionEvents = wx.CheckBox(self, -1, "Show mission events", size=(300, -1))
        self.chkShowPropagatedTrajectory = wx.CheckBox(self, -1, "Show progagated trajectory", size=(300, -1))
        self.chkShowThrustVectors = wx.CheckBox(self, -1, "Show thust vectors", size=(300, -1))
        self.chkShowEventLabels = wx.CheckBox(self, -1, "Show event labels", size=(300, -1))
        self.chkShowTextDescriptions = wx.CheckBox(self, -1, "Show text descriptions", size=(300, -1))
        self.chkNumberEventLabels = wx.CheckBox(self, -1, "Number event labels", size=(300,-1))
        self.chkDisplayEventDates = wx.CheckBox(self, -1, "Display event dates", size=(300,-1))
        self.chkDisplayEventSpecs = wx.CheckBox(self, -1, "Display event specs", size=(300,-1))
        self.chkDisplayEventMass = wx.CheckBox(self, -1, "Display event mass", size=(300,-1))
        self.chkDisplayArrivalPhaseAngle = wx.CheckBox(self, -1, "Display arrival phase angle", size=(300, -1))
        self.chkAutoTableTCMcolumn = wx.CheckBox(self, -1, "Display TCMs on table", size=(300, -1))

        reference_frames = ['Journey central body','Sun','Mercury','Venus','Earth','Moon','Mars Barycenter','Jupiter Barycenter','Saturn Barycenter','Uranus Barycenter','Neptune Barycenter','Pluto Barycenter']

        self.cmbPlotReferenceFrame = wx.ComboBox(self, -1, "Plot reference frame", choices=reference_frames, style=wx.CB_READONLY)
        self.cmbPlotReferenceFrame.SetSelection(0)

        self.btnPlotJourney = wx.Button(self, -1, "Plot trajectory", size = (300, -1))
        self.btnWriteMissionSummaryTable = wx.Button(self, -1, "Write mission summary table", size=(300, -1))

        TrajectoryPlotBoxSizer = wx.StaticBoxSizer(self.TrajectoryPlotBox, wx.VERTICAL)
        TrajectoryPlotBoxSizer.AddMany([self.chkShowPlotAxes, self.chkShowBoundaryOrbits, self.chkShowFreePointBoundaryOrbits, self.chkShowCentralBodyOrbits, self.chkShowMissionEvents, self.chkShowPropagatedTrajectory, self.chkShowThrustVectors, self.chkShowEventLabels, self.chkShowTextDescriptions, self.chkNumberEventLabels, self.chkDisplayEventDates, self.chkDisplayEventSpecs, self.chkDisplayEventMass, self.chkDisplayArrivalPhaseAngle, self.chkAutoTableTCMcolumn, self.cmbPlotReferenceFrame, self.btnPlotJourney, self.btnWriteMissionSummaryTable])
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.TrajectoryPlotBox.SetFont(font)
        
        LeftHandSizer = wx.BoxSizer(wx.VERTICAL)
        LeftHandSizer.AddMany([self.lblJourneyList, self.JourneyListBox, TrajectoryPlotBoxSizer])

        #widgets for data plot
        self.DataPlotBox = wx.StaticBox(self, -1, "Data Plots", size = (300, 300))
        DataPlotSizer = wx.StaticBoxSizer(self.DataPlotBox, wx.VERTICAL)

        self.chkPlotR = wx.CheckBox(self, -1, "distance from central body")
        self.chkPlotV = wx.CheckBox(self, -1, "velocity magnitude with respect to central body")
        self.chkPlotThrust = wx.CheckBox(self, -1, "applied thrust (N)")
        self.chkPlotIsp = wx.CheckBox(self, -1, "specific impulse (s)")
        self.chkPlotMdot = wx.CheckBox(self, -1, "mass flow rate (kg/s)")
        self.chkPlotEfficiency = wx.CheckBox(self, -1, "propulsion system efficiency")
        self.chkPlotThrottle = wx.CheckBox(self, -1, "control magnitude")
        self.chkPlotPower = wx.CheckBox(self, -1, "power produced by spacecraft (kW)")
        self.chkPlotGamma = wx.CheckBox(self, -1, "in-plane control angle (degrees)")
        self.chkPlotDelta = wx.CheckBox(self, -1, "out-of-plane control angle (degrees)")
        self.chkPlotArray_Thrust_Angle = wx.CheckBox(self, -1, "angle between solar array and thrust vector")
        self.chkPlotMass = wx.CheckBox(self, -1, "mass (kg)")
        self.chkPlotNumberOfEngines = wx.CheckBox(self, -1, "number of active thrusters")
        self.chkPlotActivePower = wx.CheckBox(self, -1, "power used by propulsion system (kW)")
        self.chkPlotWasteHeat = wx.CheckBox(self, -1, "propulsion system waste heat (kW)")
        self.chkPlotCriticalEvents = wx.CheckBox(self, -1, "mark critical events")
        self.chkPlotEarthDistance = wx.CheckBox(self, -1, "distance from Earth")
        self.chkPlotSpacecraftViewingAngle = wx.CheckBox(self, -1, "Latitude of Earth-Spacecraft viewing angle")
        self.chkPlotSunSpacecraftEarthAngle = wx.CheckBox(self, -1, "Sun-Spacecraft-Earth angle")
        self.chkPlotThrottleLevel = wx.CheckBox(self, -1, "Throttle level")
        self.chkPlotSunBoresightAngle = wx.CheckBox(self, -1, "Sun-Boresight angle")
        


        DataPlotGrid = wx.GridSizer(11, 2, 5, 5)
        DataPlotGrid.AddMany([  self.chkPlotR, self.chkPlotV,
                                self.chkPlotThrust, self.chkPlotIsp,
                                self.chkPlotMdot, self.chkPlotEfficiency,
                                self.chkPlotThrottle, self.chkPlotPower,
                                self.chkPlotGamma, self.chkPlotDelta,
                                self.chkPlotArray_Thrust_Angle, self.chkPlotMass,
                                self.chkPlotNumberOfEngines, self.chkPlotActivePower,
                                self.chkPlotWasteHeat, self.chkPlotCriticalEvents,
                                self.chkPlotEarthDistance, self.chkPlotSunSpacecraftEarthAngle,
                                self.chkPlotThrottleLevel, self.chkPlotSunBoresightAngle,
                                self.chkPlotSpacecraftViewingAngle])

        self.btnGenerateDataPlot = wx.Button(self, -1, "Generate plot")
        self.btnWriteDataReport = wx.Button(self, -1, "Write report")

        self.chkIncludeStateVectorInReport = wx.CheckBox(self, -1, "Include state vector in report")

        DataButtonSizer = wx.BoxSizer(wx.HORIZONTAL)
        DataButtonSizer.AddMany([self.btnGenerateDataPlot, self.btnWriteDataReport, self.chkIncludeStateVectorInReport])

        DataPlotSizer.AddMany([DataPlotGrid, DataButtonSizer])

        self.FormatBox = wx.StaticBox(self, -1, "Common plot options")
        self.FormatBox.SetFont(font)
        FormatBoxSizer = wx.StaticBoxSizer(self.FormatBox, wx.HORIZONTAL)
        self.lblFontSize = wx.StaticText(self, -1, "Font size")
        self.spnctrlFontSizeControl = wx.SpinCtrl(self, -1, min=1, max=100, initial=10, name="Font size")
        FormatBoxSizer.AddMany([self.lblFontSize, self.spnctrlFontSizeControl])

        RightHandSizer = wx.BoxSizer(wx.VERTICAL)
        RightHandSizer.AddMany([DataPlotSizer, FormatBoxSizer])

        mainplotbox = wx.BoxSizer(wx.HORIZONTAL)
        mainplotbox.AddMany([LeftHandSizer, RightHandSizer])

        #widgets for bubble search
        BubbleSearchBox = wx.StaticBox(self, -1, "Targets of Opportunity Bubble Search")
        BubbleSearchBox.SetFont(font)
        BubbleSearchBoxSizer = wx.StaticBoxSizer(BubbleSearchBox, wx.VERTICAL)

        self.lblBubbleSearchFile = wx.StaticText(self, -1, "Bubble search file")
        self.txtBubbleSearchFile = wx.TextCtrl(self, -1, "BubbleSearchFile", size=(600,-1))
        self.btnBubbleSearchFile = wx.Button(self, -1, "...")
        BubbleSearchFileSizer = wx.BoxSizer(wx.HORIZONTAL)
        BubbleSearchFileSizer.AddMany([self.txtBubbleSearchFile, self.btnBubbleSearchFile])

        self.lblLU = wx.StaticText(self, -1, "Length unit (km)")
        self.lblmu = wx.StaticText(self, -1, "mu (km^3/s^2)")
        self.lblRelativePositionFilterMagnitude = wx.StaticText(self, -1, "Relative position filter (km)")
        self.lblRelativeVelocityFilterMagnitude = wx.StaticText(self, -1, "Relative velocity filter (km/s)")
        self.lblMaximumMagnitude = wx.StaticText(self, -1, "Maximum Absolute Magnitude")
        self.txtLU = wx.TextCtrl(self, -1, "Length unit (km)")
        self.txtmu = wx.TextCtrl(self, -1, "mu (km^3/s^2)")
        self.txtRelativePositionFilterMagnitude = wx.TextCtrl(self, -1, "Relative position filter (km)")
        self.txtRelativeVelocityFilterMagnitude = wx.TextCtrl(self, -1, "Relative velocity filter (km/s)")
        self.txtMaximumMagnitude = wx.TextCtrl(self, -1, "Maximum Absolute Magnitude")
        self.lblCheckForEncountersAfterMissionEnd = wx.StaticText(self, -1, "Check for encounters after mission end?")
        self.chkCheckForEncountersAfterMissionEnd = wx.CheckBox(self, -1)
        self.lblPostMissionCheckDuration = wx.StaticText(self, -1, "Duration to propagate after mission end (days)")
        self.txtPostMissionCheckDuration = wx.TextCtrl(self, -1, "PostMissionCheckDuration")
        self.lblPostMissionCheckSteps = wx.StaticText(self, -1, "Number of time-steps for post-mission propagation")
        self.txtPostMissionCheckSteps = wx.TextCtrl(self, -1, "PostMissionCheckSteps")

        BubbleGrid = wx.FlexGridSizer(9, 2, 10, 10)
        BubbleGrid.AddMany([self.lblBubbleSearchFile, BubbleSearchFileSizer,
                            self.lblLU, self.txtLU,
                            self.lblmu, self.txtmu,
                            self.lblRelativePositionFilterMagnitude, self.txtRelativePositionFilterMagnitude,
                            self.lblRelativeVelocityFilterMagnitude, self.txtRelativeVelocityFilterMagnitude,
                            self.lblMaximumMagnitude, self.txtMaximumMagnitude,
                            self.lblCheckForEncountersAfterMissionEnd, self.chkCheckForEncountersAfterMissionEnd,
                            self.lblPostMissionCheckDuration, self.txtPostMissionCheckDuration,
                            self.lblPostMissionCheckSteps, self.txtPostMissionCheckSteps])

        self.btnGenerateBubbleSearch = wx.Button(self, -1, "Perform bubble search")

        BubbleSearchBoxSizer.AddMany([BubbleGrid, self.btnGenerateBubbleSearch])

        #widgets for throttle table matching
        ThrottleTableMatchBox = wx.StaticBox(self, -1, "Throttle table matching")
        ThrottleTableMatchBox.SetFont(font)
        ThrottleTableMatchBoxSizer = wx.StaticBoxSizer(ThrottleTableMatchBox, wx.VERTICAL)

        self.lblThrottleTableFile = wx.StaticText(self, -1, "Throttle table file")
        self.txtThrottleTableFile = wx.TextCtrl(self, -1, self.plotoptions.throttletablefile, size = (600,-1))
        self.btnThrottleTableFile = wx.Button(self, -1, "...")
        ThrottleTableFileSizer = wx.BoxSizer(wx.HORIZONTAL)
        ThrottleTableFileSizer.AddMany([self.txtThrottleTableFile, self.btnThrottleTableFile])

        self.lblThrottleSetChoices = wx.StaticText(self, -1, 'Throttle set')
        self.cmbThrottleSetChoices = wx.ComboBox(self, -1, style=wx.CB_READONLY, choices=['high-Thrust set','high-Isp set','full throttle box'])
        self.cmbThrottleSetChoices.SetSelection(self.plotoptions.throttlesetmode)

        self.btnGenerateThrottleHistogram = wx.Button(self, -1, "Generate throttle histogram")
        self.btnGenerateThrottleReport    = wx.Button(self, -1, "Generate throttle report")
        self.btnGenerateThrottlePlot      = wx.Button(self, -1, "Generate throttle plot")

        ThrottleGrid = wx.GridSizer(4, 2, 10, 10)
        ThrottleGrid.AddMany([self.lblThrottleTableFile, ThrottleTableFileSizer,
                              self.lblThrottleSetChoices, self.cmbThrottleSetChoices,
                              self.btnGenerateThrottleReport, self.btnGenerateThrottleHistogram,
                              self.btnGenerateThrottlePlot])

        ThrottleTableMatchBoxSizer.AddMany([ThrottleGrid])

        #output window
        self.lblOutputTextWindow = wx.StaticText(self, -1, "Output window")
        self.lblOutputTextWindow.SetFont(font)
        self.txtOutputTextWindow = wx.TextCtrl(self, -1, size=(800, 200), style=wx.TE_MULTILINE)
        self.btnClearOutputTextWindow = wx.Button(self, -1, "Clear output")
        self.btnClearOutputTextWindow.Bind(wx.EVT_BUTTON, self.ClickClearOutputTextWindow)
        
        #add everything to the main box
        mainbox = wx.BoxSizer(wx.VERTICAL)
        mainbox.AddMany([mainplotbox, BubbleSearchBoxSizer, ThrottleTableMatchBoxSizer, self.lblOutputTextWindow, self.txtOutputTextWindow, self.btnClearOutputTextWindow])

        self.SetSizer(mainbox)

        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.lblJourneyList.SetFont(font)
        self.DataPlotBox.SetFont(font)
        self.TrajectoryPlotBox.SetFont(font)

        self.mission.ChosenJourneys = [self.mission.ActiveJourney]
        self.JourneyListBox.SetSelection(self.mission.ChosenJourneys[0])
        self.plotoptions.update_mission_panel(self)
        self.bubbleoptions.update_mission_panel(self)
        self.SetupScrolling()

        #trajectory plot bindings
        self.JourneyListBox.Bind(wx.EVT_LISTBOX, self.ClickJourneyListBox)
        self.btnPlotJourney.Bind(wx.EVT_BUTTON, self.ClickPlotMissionButton)
        self.btnWriteMissionSummaryTable.Bind(wx.EVT_BUTTON, self.ClickWriteMissionSummaryTableButton)
        self.cmbPlotReferenceFrame.Bind(wx.EVT_COMBOBOX, self.ChangePlotReferenceFrame)
        self.chkShowPlotAxes.Bind(wx.EVT_CHECKBOX, self.ChangeShowPlotAxes)
        self.chkShowBoundaryOrbits.Bind(wx.EVT_CHECKBOX, self.ChangeShowBoundaryOrbits)
        self.chkShowFreePointBoundaryOrbits.Bind(wx.EVT_CHECKBOX, self.ChangeShowFreePointBoundaryOrbits)
        self.chkShowCentralBodyOrbits.Bind(wx.EVT_CHECKBOX, self.ChangeShowCentralBodyOrbits)
        self.chkShowMissionEvents.Bind(wx.EVT_CHECKBOX, self.ChangeShowMissionEvents)
        self.chkShowPropagatedTrajectory.Bind(wx.EVT_CHECKBOX, self.ChangeShowPropagatedTrajectory)
        self.chkShowThrustVectors.Bind(wx.EVT_CHECKBOX, self.ChangeShowThrustVectors)
        self.chkShowEventLabels.Bind(wx.EVT_CHECKBOX, self.ChangeShowEventLabels)
        self.chkShowTextDescriptions.Bind(wx.EVT_CHECKBOX, self.ChangeShowTextDescriptions)
        self.chkNumberEventLabels.Bind(wx.EVT_CHECKBOX, self.ChangeNumberEventLabels)
        self.chkDisplayEventDates.Bind(wx.EVT_CHECKBOX, self.ChangeDisplayEventDates)
        self.chkDisplayEventSpecs.Bind(wx.EVT_CHECKBOX, self.ChangeDisplayEventSpecs)
        self.chkDisplayEventMass.Bind(wx.EVT_CHECKBOX, self.ChangeDisplayEventMass)
        self.chkDisplayArrivalPhaseAngle.Bind(wx.EVT_CHECKBOX, self.ChangeDisplayArrivalPhaseAngle)
        self.chkAutoTableTCMcolumn.Bind(wx.EVT_CHECKBOX, self.ChangeAutoTableTCMColumn)

        #data plot bindings
        self.btnGenerateDataPlot.Bind(wx.EVT_BUTTON, self.ClickGenerateDataPlotButton)
        self.btnWriteDataReport.Bind(wx.EVT_BUTTON, self.ClickWriteDataReportButton)
        self.chkPlotR.Bind(wx.EVT_CHECKBOX, self.ChangePlotR)
        self.chkPlotV.Bind(wx.EVT_CHECKBOX, self.ChangePlotV)
        self.chkPlotThrust.Bind(wx.EVT_CHECKBOX, self.ChangePlotThrust)
        self.chkPlotIsp.Bind(wx.EVT_CHECKBOX, self.ChangePlotIsp)
        self.chkPlotMdot.Bind(wx.EVT_CHECKBOX, self.ChangePlotMdot)
        self.chkPlotEfficiency.Bind(wx.EVT_CHECKBOX, self.ChangePlotEfficiency)
        self.chkPlotThrottle.Bind(wx.EVT_CHECKBOX, self.ChangePlotThrottle)
        self.chkPlotPower.Bind(wx.EVT_CHECKBOX, self.ChangePlotPower)
        self.chkPlotGamma.Bind(wx.EVT_CHECKBOX, self.ChangePlotGamma)
        self.chkPlotDelta.Bind(wx.EVT_CHECKBOX, self.ChangePlotDelta)
        self.chkPlotArray_Thrust_Angle.Bind(wx.EVT_CHECKBOX, self.ChangePlotArray_Thrust_Angle)
        self.chkPlotMass.Bind(wx.EVT_CHECKBOX, self.ChangePlotMass)
        self.chkPlotNumberOfEngines.Bind(wx.EVT_CHECKBOX, self.ChangePlotNumberOfEngines)
        self.chkPlotActivePower.Bind(wx.EVT_CHECKBOX, self.ChangePlotActivePower)
        self.chkPlotWasteHeat.Bind(wx.EVT_CHECKBOX, self.ChangePlotWasteHeat)
        self.chkPlotCriticalEvents.Bind(wx.EVT_CHECKBOX, self.ChangePlotCriticalEvents)
        self.chkPlotEarthDistance.Bind(wx.EVT_CHECKBOX, self.ChangePlotEarthDistance)
        self.chkPlotSunSpacecraftEarthAngle.Bind(wx.EVT_CHECKBOX, self.ChangePlotSunSpacecraftEarthAngle)
        self.chkPlotSpacecraftViewingAngle.Bind(wx.EVT_CHECKBOX, self.ChangePlotSpacecraftViewingAngle)
        self.chkPlotThrottleLevel.Bind(wx.EVT_CHECKBOX, self.ChangePlotThrottleLevel)
        self.chkPlotSunBoresightAngle.Bind(wx.EVT_CHECKBOX, self.ChangePlotSunBoresightAngle)
        self.chkIncludeStateVectorInReport.Bind(wx.EVT_CHECKBOX, self.ChangeIncludeStateVectorInReport)

        #format bindings
        self.spnctrlFontSizeControl.Bind(wx.EVT_SPINCTRL, self.ChangeFontSize)

        #Bubble Search bindings
        self.txtBubbleSearchFile.Bind(wx.EVT_KILL_FOCUS, self.ChangeBubbleSearchFile)
        self.btnBubbleSearchFile.Bind(wx.EVT_BUTTON, self.ClickBubbleSearchFileButton)
        self.txtLU.Bind(wx.EVT_KILL_FOCUS, self.changeBubbleSearchLU)
        self.txtmu.Bind(wx.EVT_KILL_FOCUS, self.changeBubbleSearchmu)
        self.txtRelativePositionFilterMagnitude.Bind(wx.EVT_KILL_FOCUS, self.changeBubbleSearchRelativePositionFilterMagnitude)
        self.txtRelativeVelocityFilterMagnitude.Bind(wx.EVT_KILL_FOCUS, self.changeBubbleSearchRelativeVelocityFilterMagnitude)
        self.txtMaximumMagnitude.Bind(wx.EVT_KILL_FOCUS, self.changeMaximumMagnitude)
        self.chkCheckForEncountersAfterMissionEnd.Bind(wx.EVT_CHECKBOX, self.ChangeCheckForEncountersAfterMissionEnd)
        self.txtPostMissionCheckDuration.Bind(wx.EVT_KILL_FOCUS, self.ChangePostMissionCheckDuration)
        self.txtPostMissionCheckSteps.Bind(wx.EVT_KILL_FOCUS, self.ChangePostMissionCheckSteps)
        self.btnGenerateBubbleSearch.Bind(wx.EVT_BUTTON, self.ClickGenerateBubbleSearch)
        

        #Throttle table matching bindings
        self.txtThrottleTableFile.Bind(wx.EVT_KILL_FOCUS, self.ChangeThrottleTableFile)
        self.btnThrottleTableFile.Bind(wx.EVT_BUTTON, self.ClickThrottleTableFileButton)
        self.cmbThrottleSetChoices.Bind(wx.EVT_COMBOBOX, self.ChangeThrottleSetChoice)
        self.btnGenerateThrottleReport.Bind(wx.EVT_BUTTON, self.ClickGenerateThrottleReportButton)
        self.btnGenerateThrottleHistogram.Bind(wx.EVT_BUTTON, self.ClickGenerateThrottleHistogramButton)
        self.btnGenerateThrottlePlot.Bind(wx.EVT_BUTTON, self.ClickGenerateThrottlePlotButton)

    def ClickJourneyListBox(self, e):
        self.mission.ChosenJourneys = self.JourneyListBox.GetSelections()

    def ClickPlotMissionButton(self, e):
        self.mission.PlotMission(self.plotoptions)

    def ClickWriteMissionSummaryTableButton(self, e):
        dlg = wx.FileDialog(self, "Save", self.GetParent().dirname, self.mission.mission_name, '*.tex', wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            saved = True
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
        else:
            saved = False
        
        dlg.Destroy()

        if saved:
            self.mission.AutoTableMission(os.path.join(dirname, filename), self.plotoptions)
            os.system('pdflatex ' + os.path.join(dirname, filename) + ' -output-directory ' + dirname)
            webbrowser.open(os.path.join(dirname, filename.strip('.tex') + '.pdf'))
        else:
            return

    def ChangeShowPlotAxes(self, e):
        e.Skip()
        self.plotoptions.ShowPlotAxes = self.chkShowPlotAxes.GetValue()

        
    def ChangePlotReferenceFrame(self, e):
        e.Skip()
        self.plotoptions.PlotCentralBody = self.cmbPlotReferenceFrame.GetStringSelection()

    def ChangeShowMissionEvents(self, e):
        e.Skip()
        self.plotoptions.ShowMissionEvents = self.chkShowMissionEvents.GetValue()

    def ChangeShowBoundaryOrbits(self, e):
        e.Skip()
        self.GetParent().reset_file_menu() # kludgy solution to allow OSX users to get their file menu back if osx errors and deletes it after making a plot
        self.plotoptions.ShowBoundaryOrbits = self.chkShowBoundaryOrbits.GetValue()

    def ChangeShowFreePointBoundaryOrbits(self, e):
        e.Skip()
        self.GetParent().reset_file_menu() # kludgy solution to allow OSX users to get their file menu back if osx errors and deletes it after making a plot
        self.plotoptions.ShowFreePointBoundaryOrbits = self.chkShowFreePointBoundaryOrbits.GetValue()
        
    def ChangeShowCentralBodyOrbits(self, e):
        e.Skip()
        self.GetParent().reset_file_menu() # kludgy solution to allow OSX users to get their file menu back if osx errors and deletes it after making a plot
        self.plotoptions.ShowCentralBodyOrbits = self.chkShowCentralBodyOrbits.GetValue()        

    def ChangeShowPropagatedTrajectory(self, e):
        e.Skip()
        self.plotoptions.ShowPropagatedTrajectory = self.chkShowPropagatedTrajectory.GetValue()

    def ChangeShowThrustVectors(self, e):
        e.Skip()
        self.plotoptions.ShowThrustVectors = self.chkShowThrustVectors.GetValue()

    def ChangeShowEventLabels(self, e):
        e.Skip()
        self.plotoptions.LabelEvents = self.chkShowEventLabels.GetValue()

    def ChangeShowTextDescriptions(self, e):
        e.Skip()
        self.plotoptions.ShowTextDescriptions = self.chkShowTextDescriptions.GetValue()

    def ChangeNumberEventLabels(self, e):
        e.Skip()
        self.plotoptions.NumberEventLabels = self.chkNumberEventLabels.GetValue()
        
    def ChangeDisplayEventDates(self, e):
        e.Skip()
        self.plotoptions.DisplayEventDates = self.chkDisplayEventDates.GetValue()
        
    def ChangeDisplayEventSpecs(self, e):
        e.Skip()
        self.plotoptions.DisplayEventSpecs = self.chkDisplayEventSpecs.GetValue()

    def ChangeDisplayEventMass(self, e):
        e.Skip()
        self.plotoptions.DisplayEventMass = self.chkDisplayEventMass.GetValue()

    def ChangeDisplayArrivalPhaseAngle(self, e):
        e.Skip()
        self.plotoptions.DisplayArrivalPhaseAngle = self.chkDisplayArrivalPhaseAngle.GetValue()

    def ChangeAutoTableTCMColumn(self, e):
        e.Skip()
        self.plotoptions.AutoTableTCMcolumn = self.chkAutoTableTCMcolumn.GetValue()

    def ClickGenerateDataPlotButton(self, e):
        #generate a custom plot if at least one checkbox is active
        if (self.plotoptions.PlotR or self.plotoptions.PlotV or self.plotoptions.PlotThrust or self.plotoptions.PlotIsp or self.plotoptions.PlotMdot
            or self.plotoptions.PlotEfficiency or self.plotoptions.PlotThrottle or self.plotoptions.PlotPower or self.plotoptions.PlotGamma or self.plotoptions.PlotDelta
            or self.chkPlotMass or self.PlotEarthDistance or self.chkPlotSunSpacecraftEarthAngle or self.chkPlotThrottleLevel or self.chkPlotSpacecraftViewingAngle):
            self.mission.GenerateDataPlot(self.plotoptions)

    def ClickWriteDataReportButton(self, e):
        #generate a custom report if at least one checkbox is active
        if (self.plotoptions.PlotR or self.plotoptions.PlotV or self.plotoptions.PlotThrust or self.plotoptions.PlotIsp or self.plotoptions.PlotMdot
            or self.plotoptions.PlotEfficiency or self.plotoptions.PlotThrottle or self.plotoptions.PlotPower or self.plotoptions.PlotGamma or self.plotoptions.PlotDelta
            or self.chkPlotMass or self.PlotEarthDistance or self.chkPlotSunSpacecraftEarthAngle or self.chkPlotThrottleLevel or self.chkPlotSpacecraftViewingAngle):

            dlg = wx.FileDialog(self, "Save", self.GetParent().dirname, self.mission.mission_name, '*.systemsreport', wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
            if dlg.ShowModal() == wx.ID_OK:
                saved = True
                filename = dlg.GetFilename()
                dirname = dlg.GetDirectory()
            else:
                saved = False
        
            dlg.Destroy()


            if saved:
                self.mission.WriteDataReport(self.plotoptions, os.path.join(dirname, filename))
            else:
                return
            

    def ChangePlotR(self, e):
        e.Skip()
        self.plotoptions.PlotR = self.chkPlotR.GetValue()

    def ChangePlotV(self, e):
        e.Skip()
        self.plotoptions.PlotV = self.chkPlotV.GetValue()

    def ChangePlotThrust(self, e):
        e.Skip()
        self.plotoptions.PlotThrust = self.chkPlotThrust.GetValue()

    def ChangePlotIsp(self, e):
        e.Skip()
        self.plotoptions.PlotIsp = self.chkPlotIsp.GetValue()

    def ChangePlotMdot(self, e):
        e.Skip()
        self.plotoptions.PlotMdot = self.chkPlotMdot.GetValue()

    def ChangePlotEfficiency(self, e):
        e.Skip()
        self.plotoptions.PlotEfficiency = self.chkPlotEfficiency.GetValue()

    def ChangePlotThrottle(self, e):
        e.Skip()
        self.plotoptions.PlotThrottle = self.chkPlotThrottle.GetValue()

    def ChangePlotPower(self, e):
        e.Skip()
        self.plotoptions.PlotPower = self.chkPlotPower.GetValue()

    def ChangePlotGamma(self, e):
        e.Skip()
        self.plotoptions.PlotGamma = self.chkPlotGamma.GetValue()

    def ChangePlotDelta(self, e):
        e.Skip()
        self.plotoptions.PlotDelta = self.chkPlotDelta.GetValue()

    def ChangePlotArray_Thrust_Angle(self, e):
        e.Skip()
        self.plotoptions.PlotArray_Thrust_Angle = self.chkPlotArray_Thrust_Angle.GetValue()

    def ChangePlotMass(self, e):
        e.Skip()
        self.plotoptions.PlotMass = self.chkPlotMass.GetValue()

    def ChangePlotNumberOfEngines(self, e):
        e.Skip()
        self.plotoptions.PlotNumberOfEngines = self.chkPlotNumberOfEngines.GetValue()

    def ChangePlotActivePower(self, e):
        e.Skip()
        self.plotoptions.PlotActivePower = self.chkPlotActivePower.GetValue()

    def ChangePlotWasteHeat(self, e):
        e.Skip()
        self.plotoptions.PlotWasteHeat = self.chkPlotWasteHeat.GetValue()

    def ChangePlotCriticalEvents(self, e):
        e.Skip()
        self.plotoptions.PlotCriticalEvents = self.chkPlotCriticalEvents.GetValue()

    def ChangePlotEarthDistance(self, e):
        e.Skip()
        self.plotoptions.PlotEarthDistance = self.chkPlotEarthDistance.GetValue()

    def ChangePlotSunSpacecraftEarthAngle(self, e):
        e.Skip()
        self.plotoptions.PlotSunSpacecraftEarthAngle = self.chkPlotSunSpacecraftEarthAngle.GetValue()
        
    def ChangePlotSpacecraftViewingAngle(self, e):
        e.Skip()
        self.plotoptions.PlotSpacecraftViewingAngle = self.chkPlotSpacecraftViewingAngle.GetValue()

    def ChangePlotThrottleLevel(self, e):
        e.Skip()
        self.plotoptions.PlotThrottleLevel = self.chkPlotThrottleLevel.GetValue()

    def ChangePlotSunBoresightAngle(self, e):
        e.Skip()
        self.plotoptions.PlotSunBoresightAngle = self.chkPlotSunBoresightAngle.GetValue()

    def ChangeIncludeStateVectorInReport(self, e):
        e.skip()
        self.plotoptions.IncludeStateVectorInReport = self.chkIncludeStateVectorInReport.GetValue()

    def ChangeFontSize(self, e):
        e.Skip()
        self.plotoptions.FontSize = self.spnctrlFontSizeControl.GetValue()

    def ChangeBubbleSearchFile(self, e):
        self.bubbleoptions.smallbodyfile = self.txtBubbleSearchFile.GetValue()
    
    def ClickBubbleSearchFileButton(self, e):
        #file load dialog to get name of small bodies file for bubble search
        dlg = wx.FileDialog(self, "Select a small bodies file", self.GetParent().dirname, "", '*.SmallBody', wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.bubbleoptions.smallbodyfile = dlg.GetPath()
            self.txtBubbleSearchFile.SetValue(self.bubbleoptions.smallbodyfile)
        dlg.Destroy()

    def changeBubbleSearchLU(self, e):
        e.Skip()
        self.bubbleoptions.LU = eval(self.txtLU.GetValue())
        self.bubbleoptions.update_mission_panel(self)

    def changeBubbleSearchmu(self, e):
        e.Skip()
        self.bubbleoptions.mu = eval(self.txtmu.GetValue())
        self.bubbleoptions.update_mission_panel(self)

    def changeBubbleSearchRelativePositionFilterMagnitude(self, e):
        e.Skip()
        self.bubbleoptions.RelativePositionFilterMagnitude = eval(self.txtRelativePositionFilterMagnitude.GetValue())
        self.bubbleoptions.update_mission_panel(self)

    def changeBubbleSearchRelativeVelocityFilterMagnitude(self, e):
        e.Skip()
        self.bubbleoptions.RelativeVelocityFilterMagnitude = eval(self.txtRelativeVelocityFilterMagnitude.GetValue())
        self.bubbleoptions.update_mission_panel(self)

    def changeMaximumMagnitude(self, e):
        e.Skip()
        self.bubbleoptions.MaximumMagnitude = eval(self.txtMaximumMagnitude.GetValue())

    def ChangeCheckForEncountersAfterMissionEnd(self, e):
        e.Skip()
        self.bubbleoptions.CheckForEncountersAfterMissionEnd = self.chkCheckForEncountersAfterMissionEnd.GetValue()
        self.bubbleoptions.update_mission_panel(self)

    def ChangePostMissionCheckDuration(self, e):
        e.Skip()
        PostMissionCheckDuration = eval(self.txtPostMissionCheckDuration.GetValue())
        #note values too close to zero blow up the propagator
        if PostMissionCheckDuration < 1.0e-8:
            PostMissionCheckDuration = 1.0e-8

        self.bubbleoptions.PostMissionCheckDuration = PostMissionCheckDuration

    def ChangePostMissionCheckSteps(self, e):
        e.Skip()
        self.bubbleoptions.PostMissionCheckSteps = int(float(self.txtPostMissionCheckSteps.GetValue()))

    def ClickGenerateBubbleSearch(self, e):

        dlg = wx.FileDialog(self, "Save", self.GetParent().dirname, self.mission.mission_name, '*.bubble', wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            saved = True
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
        else:
            saved = False
        
        dlg.Destroy()


        if saved:
            self.mission.BubbleSearch(self.bubbleoptions, os.path.join(dirname, filename), self.txtOutputTextWindow)
        else:
            return

       #Throttle table matching bindings

    def ChangeThrottleTableFile(self, e):
        e.Skip()
        self.plotoptions.throttletablefile = self.txtThrottleTableFile.GetValue()

    def ClickThrottleTableFileButton(self, e):
        #file load dialog to get name of 
        dlg = wx.FileDialog(self, "Select a throttle table file", self.GetParent().dirname, "", '*.ThrottleTable', wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.plotoptions.throttletablefile = dlg.GetPath()
            self.txtThrottleTableFile.SetValue(self.throttletablefile)
        dlg.Destroy()

    def ChangeThrottleSetChoice(self, e):
        e.skip()
        self.plotoptions.throttlesetmode = self.cmbThrottleSetChoices.GetSelection()

    def ClickGenerateThrottleReportButton(self, e):
        dlg = wx.FileDialog(self, "Save", self.GetParent().dirname, self.mission.mission_name, '.ThrottleReport', wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            saved = True
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
        else:
            saved = False
        
        dlg.Destroy()


        if saved:
            self.mission.GenerateThrottleReport(os.path.join(dirname, filename))
        else:
            return

    def ClickGenerateThrottleHistogramButton(self, e):
        self.mission.GenerateThrottleHistogram(self.plotoptions)

    def ClickGenerateThrottlePlotButton(self, e):
        self.mission.GenerateThrottlePlot(self.plotoptions.throttletablefile, self.plotoptions)

    def ClickClearOutputTextWindow(self, e):
        e.Skip()
        self.txtOutputTextWindow.Clear()