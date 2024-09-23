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

class OrbitElementsPanel(wx.Panel):
    def __init__(self, parent, missionoptions):
        self.missionoptions = missionoptions
        self.parent = parent
        
        wx.Panel.__init__(self, parent)

        #custom  elements
        self.lblelements_frame = wx.StaticText(self, -1, "Journey elements frame")
        elements_frame_choices = ['0: ICRF', '1: J2000_BCI', '2: J2000_BCF', '3: TrueOfDate_BCI', '4: TrueOfDate_BCF', '5: Principle Axes', '6: Topocentric','7: Polar', '8: SAM', '9: ObjectReferenced']
        self.cmbelements_frame = wx.ComboBox(self, -1, choices = elements_frame_choices, style=wx.CB_READONLY)
        self.elements_frame_box = wx.BoxSizer(wx.HORIZONTAL)
        self.elements_frame_box.Add(self.lblelements_frame)
        self.elements_frame_box.AddSpacer(5)
        self.elements_frame_box.Add(self.cmbelements_frame)

        self.lblelements_representation = wx.StaticText(self, -1, "Journey elements state representation")
        elements_representation_choices = ['0: Cartesian', '1: SphericalRADEC', '2: SphericalAZFPA', '3: COE', '4: MEE', '5: IncomingBplane', '6: OutgoingPlane']
        self.cmbelements_representation = wx.ComboBox(self, -1, choices = elements_representation_choices, style=wx.CB_READONLY)
        self.elements_representation_box = wx.BoxSizer(wx.HORIZONTAL)
        self.elements_representation_box.Add(self.lblelements_representation)
        self.elements_representation_box.AddSpacer(5)
        self.elements_representation_box.Add(self.cmbelements_representation)

        self.elements_epoch_box = wx.BoxSizer(wx.HORIZONTAL)
        self.lblelements_reference_epoch = wx.StaticText(self, -1, "Reference epoch")
        self.txtelements_reference_epoch = wx.TextCtrl(self, -1, "Reference epoch")
        self.chkAllowJourneyFreePointToPropagate = wx.CheckBox(self, -1, label="Allow state to propagate?")
        self.elements_epoch_box.Add(self.lblelements_reference_epoch)
        self.elements_epoch_box.AddSpacer(5)
        self.elements_epoch_box.Add(self.txtelements_reference_epoch)
        self.elements_epoch_box.AddSpacer(5)
        self.elements_epoch_box.Add(self.chkAllowJourneyFreePointToPropagate)

        empty_cell = wx.StaticText(self, -1, "")
        self.lblvaryelements = wx.StaticText(self, -1, "Vary?")
        self.lblelementsvalue = wx.StaticText(self, -1, "Value")
        self.lblelementslower = wx.StaticText(self, -1, "Lower bound")
        self.lblelementsupper = wx.StaticText(self, -1, "Upper bound")
        self.lblSMA = wx.StaticText(self, -1, "SMA (km)")
        self.lblECC = wx.StaticText(self, -1, "ECC")
        self.lblINC = wx.StaticText(self, -1, "INC (degrees)")
        self.lblRAAN = wx.StaticText(self, -1, "RAAN (degrees)")
        self.lblAOP = wx.StaticText(self, -1, "AOP (degrees)")
        self.lblMA = wx.StaticText(self, -1, "MA (degrees)")
        self.chkSMA = wx.CheckBox(self, -1)
        self.chkECC = wx.CheckBox(self, -1)
        self.chkINC = wx.CheckBox(self, -1)
        self.chkRAAN = wx.CheckBox(self, -1)
        self.chkAOP = wx.CheckBox(self, -1)
        self.chkMA = wx.CheckBox(self, -1)
        self.txtSMA = wx.TextCtrl(self, -1, "SMA_val")
        self.txtSMA0 = wx.TextCtrl(self, -1, "SMA_val0")
        self.txtSMA1 = wx.TextCtrl(self, -1, "SMA_val1")
        self.txtECC = wx.TextCtrl(self, -1, "ECC_val")
        self.txtECC0 = wx.TextCtrl(self, -1, "ECC_val0")
        self.txtECC1 = wx.TextCtrl(self, -1, "ECC_val1")
        self.txtINC = wx.TextCtrl(self, -1, "INC_val")
        self.txtINC0 = wx.TextCtrl(self, -1, "INC_val0")
        self.txtINC1 = wx.TextCtrl(self, -1, "INC_val1")
        self.txtRAAN = wx.TextCtrl(self, -1, "RAAN_val")
        self.txtRAAN0 = wx.TextCtrl(self, -1, "RAAN_val0")
        self.txtRAAN1 = wx.TextCtrl(self, -1, "RAAN_val1")
        self.txtAOP = wx.TextCtrl(self, -1, "AOP_val")
        self.txtAOP0 = wx.TextCtrl(self, -1, "AOP_val0")
        self.txtAOP1 = wx.TextCtrl(self, -1, "AOP_val1")
        self.txtMA = wx.TextCtrl(self, -1, "MA_val")
        self.txtMA0 = wx.TextCtrl(self, -1, "MA_val0")
        self.txtMA1 = wx.TextCtrl(self, -1, "MA_val1")
        self.ElementsSizer = wx.FlexGridSizer(14,5,5,5)
        self.ElementsSizer.AddMany([empty_cell, self.lblvaryelements, self.lblelementsvalue, self.lblelementslower, self.lblelementsupper, 
                                            self.lblSMA, self.chkSMA, self.txtSMA, self.txtSMA0, self.txtSMA1,
                                            self.lblECC, self.chkECC, self.txtECC, self.txtECC0, self.txtECC1,
                                            self.lblINC, self.chkINC, self.txtINC, self.txtINC0, self.txtINC1,
                                            self.lblRAAN, self.chkRAAN, self.txtRAAN, self.txtRAAN0, self.txtRAAN1,
                                            self.lblAOP, self.chkAOP, self.txtAOP, self.txtAOP0, self.txtAOP1,
                                            self.lblMA, self.chkMA, self.txtMA, self.txtMA0, self.txtMA1])

    def update():
        raise NotImplementedError("Subclass must implement abstract method")