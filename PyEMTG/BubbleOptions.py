#class file for BubbleOptions class
#this class contains settings for bubble search
#bubble search is used to find flybys of opportunity from an existing trajectory

class BubbleOptions(object):
    def __init__(self):
        self.initialize()

    def __init__(self, smallbodyfile):
        self.initialize()
        if smallbodyfile != None:
            self.smallbodyfile = smallbodyfile
        else:
            self.smallbodyfile = "blah"

    def initialize(self):
        self.smallbodyfile = './MainBeltBright.SmallBody'
        self.LU = 149597870.691
        self.mu = 1.32712440018e+11
        self.RelativePositionFilterMagnitude = 0.1 * self.LU # km
        self.RelativeVelocityFilterMagnitude = 20.0 # km/s
        self.MaximumMagnitude = 14.0
        self.ifSpice = True
        self.spiceFiles = ['D:/Projects/Lucy/Universe/ephemeris_files/AllTrojans.bsp',
                           'D:/Projects/Lucy/Universe/ephemeris_files/2436157.bsp',
                           'D:/Projects/Lucy/Universe/ephemeris_files/2233683.bsp',
                           'D:/Projects/Lucy/Universe/AllMainBelt/MainBeltBright.bsp',
                           'D:/Projects/Lucy/Universe/ephemeris_files/naif0011.tls',
                           'D:/Projects/Lucy/Universe/ephemeris_files/pck00010.tpc',
                           'D:/Projects/Lucy/Universe/ephemeris_files/de430.bsp']
        self.SPICE_IDs = []
        self.SPICEID_file = 'D:/Projects/Lucy/Universe/AllTrojans.dat', 'D:/Projects/Lucy/Universe/AllMainBelt/MainBeltBrighter10.dat'
        
        self.CheckForEncountersAfterMissionEnd = False
        self.PostMissionCheckDuration = 1.0
        self.PostMissionCheckSteps = 10

    def update_mission_panel(self, missionpanel):
        missionpanel.txtBubbleSearchFile.SetValue(self.smallbodyfile)
        missionpanel.txtLU.SetValue(str(self.LU))
        missionpanel.txtmu.SetValue(str(self.mu))
        missionpanel.txtRelativePositionFilterMagnitude.SetValue(str(self.RelativePositionFilterMagnitude))
        missionpanel.txtRelativeVelocityFilterMagnitude.SetValue(str(self.RelativeVelocityFilterMagnitude))
        missionpanel.txtMaximumMagnitude.SetValue(str(self.MaximumMagnitude))
        missionpanel.chkCheckForEncountersAfterMissionEnd.SetValue(self.CheckForEncountersAfterMissionEnd)
        missionpanel.txtPostMissionCheckDuration.SetValue(str(self.PostMissionCheckDuration))
        missionpanel.txtPostMissionCheckSteps.SetValue(str(self.PostMissionCheckSteps))

        if self.CheckForEncountersAfterMissionEnd == True:
            missionpanel.lblPostMissionCheckDuration.Show(True)
            missionpanel.txtPostMissionCheckDuration.Show(True)
            missionpanel.lblPostMissionCheckSteps.Show(True)
            missionpanel.txtPostMissionCheckSteps.Show(True)
        else:
            missionpanel.lblPostMissionCheckDuration.Show(False)
            missionpanel.txtPostMissionCheckDuration.Show(False)
            missionpanel.lblPostMissionCheckSteps.Show(False)
            missionpanel.txtPostMissionCheckSteps.Show(False)

        missionpanel.Layout()