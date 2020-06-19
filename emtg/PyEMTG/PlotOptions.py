class PlotOptions(object):
    def __init__(self, throttletablefile='none', default_universe_path='none', de_file='none', leapseconds_file='none'):
        self.ShowPlotAxes = False
        self.ShowBoundaryOrbits = True
        self.ShowFreePointBoundaryOrbits = False
        self.ShowCentralBodyOrbits = False
        self.ShowMissionEvents = True
        self.ShowPropagatedTrajectory = True
        self.ShowThrustVectors = True
        self.LabelEvents = False
        self.ShowTextDescriptions = True
        self.NumberEventLabels = True
        self.DisplayEventSpecs = True
        self.DisplayEventDates = True
        self.DisplayEventMass = True
        self.DisplayArrivalPhaseAngle = False
        self.PlotR = True
        self.PlotV = False
        self.PlotThrust = False
        self.PlotIsp = False
        self.PlotMdot = False
        self.PlotEfficiency = False
        self.PlotThrottle = False
        self.PlotPower = False
        self.PlotGamma = False
        self.PlotDelta = False
        self.PlotArray_Thrust_Angle = False
        self.PlotMass = False
        self.PlotNumberOfEngines = False
        self.PlotActivePower = False
        self.PlotWasteHeat = False
        self.PlotCriticalEvents = False
        self.PlotEarthDistance = False
        self.PlotSpacecraftViewingAngle = False
        self.PlotSunSpacecraftEarthAngle = False
        self.PlotThrottleLevel = False
        self.PlotSunBoresightAngle = False
        self.throttlesetmode = 0
        self.throttletablefile = throttletablefile
        self.FontSize = 10
        self.PlotCB_thrust_angle = False
        self.IncludeStateVectorInReport = False
        self.PlotCentralBody = 'Journey central body'        
        self.SPICEpath = default_universe_path + '/ephemeris_files/'
        self.DEfile = default_universe_path + '/ephemeris_files/' + de_file
        self.LeapSecondsFile = default_universe_path + '/ephemeris_files/' + leapseconds_file

        self.AutoTableDeltavColumn = False
        self.AutoTableAltitudeColumn = False
        self.AutoTableTCMcolumn = False

        self.EventCounter = 1

    def update_mission_panel(self, missionpanel):
        #trajectory plot options
        missionpanel.chkShowPlotAxes.SetValue(self.ShowPlotAxes)
        missionpanel.chkShowBoundaryOrbits.SetValue(self.ShowBoundaryOrbits)
        missionpanel.chkShowMissionEvents.SetValue(self.ShowMissionEvents)
        missionpanel.chkShowPropagatedTrajectory.SetValue(self.ShowPropagatedTrajectory)
        missionpanel.chkShowThrustVectors.SetValue(self.ShowThrustVectors)
        missionpanel.chkShowEventLabels.SetValue(self.LabelEvents)
        missionpanel.chkShowTextDescriptions.SetValue(self.ShowTextDescriptions)
        missionpanel.chkNumberEventLabels.SetValue(self.NumberEventLabels)          
        missionpanel.chkDisplayEventDates.SetValue(self.DisplayEventDates)
        missionpanel.chkDisplayEventSpecs.SetValue(self.DisplayEventSpecs)
        missionpanel.chkDisplayEventMass.SetValue(self.DisplayEventMass)
        missionpanel.chkDisplayArrivalPhaseAngle.SetValue(self.DisplayArrivalPhaseAngle)
        missionpanel.chkAutoTableTCMcolumn.SetValue(self.AutoTableTCMcolumn)

        #data plot options
        missionpanel.chkPlotR.SetValue(self.PlotR)
        missionpanel.chkPlotV.SetValue(self.PlotV)
        missionpanel.chkPlotThrust.SetValue(self.PlotThrust)
        missionpanel.chkPlotIsp.SetValue(self.PlotIsp)
        missionpanel.chkPlotMdot.SetValue(self.PlotMdot)
        missionpanel.chkPlotEfficiency.SetValue(self.PlotEfficiency)
        missionpanel.chkPlotThrottle.SetValue(self.PlotThrottle)
        missionpanel.chkPlotPower.SetValue(self.PlotPower)
        missionpanel.chkPlotActivePower.SetValue(self.PlotActivePower)
        missionpanel.chkPlotNumberOfEngines.SetValue(self.PlotNumberOfEngines)
        missionpanel.chkPlotGamma.SetValue(self.PlotGamma)
        missionpanel.chkPlotDelta.SetValue(self.PlotDelta)
        missionpanel.chkPlotArray_Thrust_Angle.SetValue(self.PlotArray_Thrust_Angle)
        missionpanel.chkPlotMass.SetValue(self.PlotMass)
        missionpanel.chkPlotNumberOfEngines.SetValue(self.PlotNumberOfEngines)
        missionpanel.chkPlotActivePower.SetValue(self.PlotActivePower)
        missionpanel.chkPlotWasteHeat.SetValue(self.PlotWasteHeat)
        missionpanel.chkPlotCriticalEvents.SetValue(self.PlotCriticalEvents)
        missionpanel.chkPlotEarthDistance.SetValue(self.PlotEarthDistance)
        missionpanel.chkPlotSunSpacecraftEarthAngle.SetValue(self.PlotSunSpacecraftEarthAngle)
        missionpanel.chkPlotSpacecraftViewingAngle.SetValue(self.PlotSpacecraftViewingAngle)
        missionpanel.chkPlotThrottleLevel.SetValue(self.PlotThrottleLevel)
        missionpanel.chkIncludeStateVectorInReport.SetValue(self.IncludeStateVectorInReport)

        #throttle matching options
        missionpanel.cmbThrottleSetChoices.SetSelection(self.throttlesetmode)

        #format options
        missionpanel.spnctrlFontSizeControl.SetValue(self.FontSize)