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
import platform

class JourneyStagingPanel(wx.Panel):
    def __init__(self, parent, missionoptions):
        self.missionoptions = missionoptions
        self.parent = parent
        
        wx.Panel.__init__(self, parent)
        
        self.staginggridtitle = wx.StaticText(self, -1, "Staging options")

        
        self.chkstage_after_departure = wx.CheckBox(self, -1, label='Stage after departure')
        self.chkstage_before_arrival = wx.CheckBox(self, -1, label='Stage before arrival')
        self.chkstage_after_arrival = wx.CheckBox(self, -1, label='Stage after arrival')

        stagegrid = wx.FlexGridSizer(20, 1, 5, 5)
        stagegrid.AddMany([self.chkstage_after_departure,
                           self.chkstage_before_arrival,
                           self.chkstage_after_arrival])

        stagebox = wx.BoxSizer(wx.VERTICAL)
        stagebox.AddMany([self.staginggridtitle, stagegrid])

        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.staginggridtitle.SetFont(font)

        #now tie everything together
        self.SetSizer(stagebox)

        #bindings
        
        self.chkstage_after_departure.Bind(wx.EVT_CHECKBOX, self.Changestage_after_departure)
        self.chkstage_before_arrival.Bind(wx.EVT_CHECKBOX, self.Changestage_before_arrival)
        self.chkstage_after_arrival.Bind(wx.EVT_CHECKBOX, self.Changestage_after_arrival)

    def update(self):
        self.chkstage_after_departure.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].stage_after_departure)
        self.chkstage_before_arrival.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].stage_before_arrival)
        self.chkstage_after_arrival.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].stage_after_arrival)

        #re-size the panel
        self.Layout()
        if platform.system() == 'Windows':
            self.parent.SetupScrolling(scrollToTop=False)

    #event handlers

    def Changestage_after_departure(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].stage_after_departure = self.chkstage_after_departure.GetValue()

    def Changestage_before_arrival(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].stage_before_arrival = self.chkstage_before_arrival.GetValue()

    def Changestage_after_arrival(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].stage_after_arrival = self.chkstage_after_arrival.GetValue()
