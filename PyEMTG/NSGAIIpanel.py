import wx
import NSGAIIpopulation
import NSGAIIFilter
import numpy
import copy

import matplotlib
#matplotlib.use('WXAgg') 
from matplotlib.backends.backend_agg import FigureCanvas as FigureCanvas
from matplotlib.figure import Figure

class NSGAIIPlotOptions:
    def __init__(self):
        self.UpperBounds = []
        self.LowerBounds = []
        self.TimeUnit = 0
        self.EpochUnit = 0
        self.BaseMarkerSize = 20.0
        self.FontSize = 10.0

class NSGAIIpanel(wx.Panel):
    def __init__(self, parent, Population):
        wx.Panel.__init__(self, parent)

        self.NSGAIIpopulation = Population
        self.plotoptions = NSGAIIPlotOptions()
        self.Xobjective = 0
        self.Yobjective = 1
        self.Zobjective = 2
        self.Cobjective = 3
        self.Sobjective = 4
        self.filter = []

        self.plotoptions.LowerBounds = numpy.zeros(5)
        self.plotoptions.UpperBounds = 1.0e+50 * numpy.ones(5)

        #first we want an array of [Objective # label, combobox to select objective, lowerbound for objective display, upperbound for objective display]
        #we want one of these rows for every objective in the population file
        #each combobox after the second should have an option for "do not display"
        #put the array of objective selectors in a frame
        self.AxisOptionsBox = wx.StaticBox(self, -1, "Axis options")
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.AxisOptionsBox.SetFont(font)
        AxisOptionsBoxSizer = wx.StaticBoxSizer(self.AxisOptionsBox, wx.VERTICAL)

        self.objective_selectors = []
        self.objective_upperbound_fields = []
        self.objective_lowerbound_fields = []
        self.objective_row_sizers = []
        
        #x axis
        xaxislabel = wx.StaticText(self, -1, "X axis", size=(100, -1))
        self.objective_selectors.append(wx.ComboBox(self, -1, choices = self.NSGAIIpopulation.objective_column_headers, style = wx.CB_READONLY, size=(200, -1)))
        self.objective_selectors[-1].SetSelection(self.Xobjective)
        self.objective_lowerbound_fields.append(wx.TextCtrl(self, -1, str(0.0)))
        self.objective_upperbound_fields.append(wx.TextCtrl(self, -1, str(1.0e+50)))
        self.objective_selectors[-1].Bind(wx.EVT_COMBOBOX, self.ChangeXObjective)
        self.objective_lowerbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeXLowerBound)
        self.objective_upperbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeXUpperBound)
        XaxisSizer = wx.BoxSizer(wx.HORIZONTAL)
        XaxisSizer.AddMany([xaxislabel, self.objective_selectors[-1], self.objective_lowerbound_fields[-1], self.objective_upperbound_fields[-1]])
        AxisOptionsBoxSizer.Add(XaxisSizer)

        #y axis
        yaxislabel = wx.StaticText(self, -1, "Y axis", size=(100, -1))
        self.objective_selectors.append(wx.ComboBox(self, -1, choices = self.NSGAIIpopulation.objective_column_headers, style = wx.CB_READONLY, size=(200, -1)))
        self.objective_selectors[-1].SetSelection(self.Yobjective)
        self.objective_lowerbound_fields.append(wx.TextCtrl(self, -1, str(0.0)))
        self.objective_upperbound_fields.append(wx.TextCtrl(self, -1, str(1.0e+50)))
        self.objective_selectors[-1].Bind(wx.EVT_COMBOBOX, self.ChangeYObjective)
        self.objective_lowerbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeYLowerBound)
        self.objective_upperbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeYUpperBound)
        YaxisSizer = wx.BoxSizer(wx.HORIZONTAL)
        YaxisSizer.AddMany([yaxislabel, self.objective_selectors[-1], self.objective_lowerbound_fields[-1], self.objective_upperbound_fields[-1]])
        AxisOptionsBoxSizer.Add(YaxisSizer)

        #z axis
        if len(self.NSGAIIpopulation.objective_column_headers) > 2:
            zaxislabel = wx.StaticText(self, -1, "Z axis", size=(100, -1))
            zaxischoices = copy.deepcopy(self.NSGAIIpopulation.objective_column_headers)
            zaxischoices.append('do not display')
            self.objective_selectors.append(wx.ComboBox(self, -1, choices = zaxischoices, style = wx.CB_READONLY, size=(200, -1)))
            self.objective_selectors[-1].SetSelection(self.Zobjective)
            self.objective_lowerbound_fields.append(wx.TextCtrl(self, -1, str(0.0)))
            self.objective_upperbound_fields.append(wx.TextCtrl(self, -1, str(1.0e+50)))
            self.objective_selectors[-1].Bind(wx.EVT_COMBOBOX, self.ChangeZObjective)
            self.objective_lowerbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeZLowerBound)
            self.objective_upperbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeZUpperBound)
            ZaxisSizer = wx.BoxSizer(wx.HORIZONTAL)
            ZaxisSizer.AddMany([zaxislabel, self.objective_selectors[-1], self.objective_lowerbound_fields[-1], self.objective_upperbound_fields[-1]])
            AxisOptionsBoxSizer.Add(ZaxisSizer)

        #color axis
        if len(self.NSGAIIpopulation.objective_column_headers) > 3:
            caxislabel = wx.StaticText(self, -1, "Color axis", size=(100, -1))
            caxischoices = copy.deepcopy(self.NSGAIIpopulation.objective_column_headers)
            caxischoices.append('do not display')
            self.objective_selectors.append(wx.ComboBox(self, -1, choices = caxischoices, style = wx.CB_READONLY, size=(200, -1)))
            self.objective_selectors[-1].SetSelection(self.Cobjective)
            self.objective_lowerbound_fields.append(wx.TextCtrl(self, -1, str(0.0)))
            self.objective_upperbound_fields.append(wx.TextCtrl(self, -1, str(1.0e+50)))
            self.objective_selectors[-1].Bind(wx.EVT_COMBOBOX, self.ChangeCObjective)
            self.objective_lowerbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeCLowerBound)
            self.objective_upperbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeCUpperBound)
            CaxisSizer = wx.BoxSizer(wx.HORIZONTAL)
            CaxisSizer.AddMany([caxislabel, self.objective_selectors[-1], self.objective_lowerbound_fields[-1], self.objective_upperbound_fields[-1]])
            AxisOptionsBoxSizer.Add(CaxisSizer)

        #size axis
        if len(self.NSGAIIpopulation.objective_column_headers) > 4:
            saxislabel = wx.StaticText(self, -1, "Size axis", size=(100, -1))
            saxischoices = copy.deepcopy(self.NSGAIIpopulation.objective_column_headers)
            saxischoices.append('do not display')
            self.objective_selectors.append(wx.ComboBox(self, -1, choices = saxischoices, style = wx.CB_READONLY, size=(200, -1)))
            self.objective_selectors[-1].SetSelection(self.Sobjective)
            self.objective_lowerbound_fields.append(wx.TextCtrl(self, -1, str(0.0)))
            self.objective_upperbound_fields.append(wx.TextCtrl(self, -1, str(1.0e+50)))
            self.objective_selectors[-1].Bind(wx.EVT_COMBOBOX, self.ChangeSObjective)
            self.objective_lowerbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeSLowerBound)
            self.objective_upperbound_fields[-1].Bind(wx.EVT_KILL_FOCUS, self.ChangeSUpperBound)
            SaxisSizer = wx.BoxSizer(wx.HORIZONTAL)
            SaxisSizer.AddMany([saxislabel, self.objective_selectors[-1], self.objective_lowerbound_fields[-1], self.objective_upperbound_fields[-1]])
            AxisOptionsBoxSizer.Add(SaxisSizer)

        #next we want checkboxes for any other plot options
        self.lblTimeUnit = wx.StaticText(self, -1, "Display time unit", size=(200,-1))
        TimeUnitChoices = ['years','days']
        self.cmbTimeUnit = wx.ComboBox(self, -1, choices=TimeUnitChoices, style = wx.CB_READONLY)
        self.cmbTimeUnit.SetSelection(self.plotoptions.TimeUnit)
        self.cmbTimeUnit.Bind(wx.EVT_COMBOBOX, self.ChangeTimeUnit)
        TimeUnitSizer = wx.BoxSizer(wx.HORIZONTAL)
        TimeUnitSizer.AddMany([self.lblTimeUnit, self.cmbTimeUnit])

        self.lblEpochUnit = wx.StaticText(self, -1, "Display epoch unit", size=(200,-1))
        EpochUnitChoices = ['TDB Gregorian','TDB MJD']
        self.cmbEpochUnit = wx.ComboBox(self, -1, choices=EpochUnitChoices, style = wx.CB_READONLY)
        self.cmbEpochUnit.SetSelection(self.plotoptions.EpochUnit)
        self.cmbEpochUnit.Bind(wx.EVT_COMBOBOX, self.ChangeEpochUnit)
        EpochUnitSizer = wx.BoxSizer(wx.HORIZONTAL)
        EpochUnitSizer.AddMany([self.lblEpochUnit, self.cmbEpochUnit])

        #marker size control
        self.lblMarkerSize = wx.StaticText(self, -1, "Base marker size")
        self.spnctrlMarkerSizeControl = wx.SpinCtrl(self, -1, min=1, max=10000000, initial=20, name="Marker size")
        self.spnctrlMarkerSizeControl.Bind(wx.EVT_SPINCTRL, self.ChangeMarkerSize)

        #font size control
        self.lblFontSize = wx.StaticText(self, -1, "Font size")
        self.spnctrlFontSizeControl = wx.SpinCtrl(self, -1, min=1, max=100, initial=10, name="Font size")
        self.spnctrlFontSizeControl.Bind(wx.EVT_SPINCTRL, self.ChangeFontSize)
        
        
        FormatBoxSizer = wx.FlexGridSizer(2,2,5,5)
        FormatBoxSizer.AddMany([self.lblMarkerSize, self.spnctrlMarkerSizeControl,
                                self.lblFontSize, self.spnctrlFontSizeControl])


        self.PlotOptionsBox = wx.StaticBox(self, -1, "Plot Options", size = (300, 300))
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.PlotOptionsBox.SetFont(font)
        PlotOptionsSizer = wx.StaticBoxSizer(self.PlotOptionsBox, wx.VERTICAL)
        PlotOptionsSizer.AddMany([TimeUnitSizer, EpochUnitSizer, FormatBoxSizer])

        
        
        #button to make the plot
        self.btnPlotPopulation = wx.Button(self, -1, "Plot Population")
        self.btnPlotPopulation.Bind(wx.EVT_BUTTON, self.ClickPlotPopulation)

        #output window
        self.lblOutputTextWindow = wx.StaticText(self, -1, "Output window")
        self.lblOutputTextWindow.SetFont(font)
        self.txtOutputTextWindow = wx.TextCtrl(self, -1, size=(800, 200), style=wx.TE_MULTILINE)
        self.btnClearOutputTextWindow = wx.Button(self, -1, "Clear output")
        self.btnClearOutputTextWindow.Bind(wx.EVT_BUTTON, self.ClickClearOutputTextWindow)

        #filter panel
        ObjectFilterBox = wx.StaticBox(self, -1, "NSGAII object filter")
        ObjectFilterBox.SetFont(font)
        ObjectFilterBoxSizer = wx.StaticBoxSizer(ObjectFilterBox, wx.VERTICAL)

        self.lblObjectFilterFile = wx.StaticText(self, -1, "NSGAII object filter file")
        self.txtObjectFilterFile = wx.TextCtrl(self, -1, "ObjectFilterFile", size=(600, -1))
        self.btnObjectFilterFile = wx.Button(self, -1, "...")
        ObjectFilterFileSizer = wx.BoxSizer(wx.HORIZONTAL)
        ObjectFilterFileSizer.AddMany([self.lblObjectFilterFile, self.txtObjectFilterFile, self.btnObjectFilterFile])

    

        self.txtObjectFilterFile.Bind(wx.EVT_KILL_FOCUS, self.ChangeObjectFilterFile)
        self.btnObjectFilterFile.Bind(wx.EVT_BUTTON, self.ClickObjectFilterFileButton)

        self.lbxObjectFilterList = wx.ListBox(self, -1, size=(600,200), style=wx.LB_EXTENDED)
        #self.lbxObjectFilterList.Bind(wx.EVT_LISTBOX, self.ChangeFilterListSelections)

        ObjectFilterBoxSizer.AddMany([ObjectFilterFileSizer, self.lbxObjectFilterList])

        #plot panel
        self.PopulationFigure = Figure(facecolor = 'white',figsize=(10,10))

        self.PlotCanvas = FigureCanvas(self, -1, self.PopulationFigure)

        if len(self.NSGAIIpopulation.ordered_list_of_objectives) == 2:
            self.PopulationAxes = self.PopulationFigure.add_axes([0,0,1,1])
        else:
            self.PopulationAxes = self.PopulationFigure.add_axes([0,0,1,1], projection='3d')

        #put everything in a big sizer

        leftsizer = wx.BoxSizer(wx.VERTICAL)
        leftsizer.AddMany([AxisOptionsBoxSizer, PlotOptionsSizer, self.btnPlotPopulation, self.lblOutputTextWindow, self.txtOutputTextWindow, self.btnClearOutputTextWindow, ObjectFilterBoxSizer])

        mainsizer = wx.BoxSizer(wx.HORIZONTAL)
        mainsizer.AddMany([leftsizer, (self.PlotCanvas, wx.SHAPED|wx.ALL)])
        self.SetSizer(mainsizer)

        
    #methods    
    def ChangeXObjective(self, event):
        self.Xobjective = self.objective_selectors[0].GetSelection()

    def ChangeYObjective(self, event):
        self.Yobjective = self.objective_selectors[1].GetSelection()

    def ChangeZObjective(self, event):
        self.Zobjective = self.objective_selectors[2].GetSelection()

    def ChangeCObjective(self, event):
        self.Cobjective = self.objective_selectors[3].GetSelection()
    
    def ChangeSObjective(self, event):
        self.Sobjective = self.objective_selectors[4].GetSelection()

    def ChangeXLowerBound(self, event):
        event.Skip()
        self.plotoptions.LowerBounds[0] = float(eval(self.objective_lowerbound_fields[0].GetValue()))
        self.objective_lowerbound_fields[0].SetValue(str(self.plotoptions.LowerBounds[0]))

    def ChangeXUpperBound(self, event):
        event.Skip()
        self.plotoptions.UpperBounds[0] = float(eval(self.objective_upperbound_fields[0].GetValue()))
        self.objective_upperbound_fields[0].SetValue(str(self.plotoptions.UpperBounds[0]))

    def ChangeYLowerBound(self, event):
        event.Skip()
        self.plotoptions.LowerBounds[1] = float(eval(self.objective_lowerbound_fields[1].GetValue()))
        self.objective_lowerbound_fields[1].SetValue(str(self.plotoptions.LowerBounds[1]))

    def ChangeYUpperBound(self, event):
        event.Skip()
        self.plotoptions.UpperBounds[1] = float(eval(self.objective_upperbound_fields[1].GetValue()))
        self.objective_upperbound_fields[1].SetValue(str(self.plotoptions.UpperBounds[1]))

    def ChangeZLowerBound(self, event):
        event.Skip()
        self.plotoptions.LowerBounds[2] = float(eval(self.objective_lowerbound_fields[2].GetValue()))
        self.objective_lowerbound_fields[2].SetValue(str(self.plotoptions.LowerBounds[2]))

    def ChangeZUpperBound(self, event):
        event.Skip()
        self.plotoptions.UpperBounds[2] = float(eval(self.objective_upperbound_fields[2].GetValue()))
        self.objective_upperbound_fields[2].SetValue(str(self.plotoptions.UpperBounds[2]))

    def ChangeCLowerBound(self, event):
        event.Skip()
        self.plotoptions.LowerBounds[3] = float(eval(self.objective_lowerbound_fields[3].GetValue()))
        self.objective_lowerbound_fields[3].SetValue(str(self.plotoptions.LowerBounds[3]))

    def ChangeCUpperBound(self, event):
        event.Skip()
        self.plotoptions.UpperBounds[3] = float(eval(self.objective_upperbound_fields[3].GetValue()))
        self.objective_upperbound_fields[3].SetValue(str(self.plotoptions.UpperBounds[3]))
        
    def ChangeSLowerBound(self, event):
        event.Skip()
        self.plotoptions.LowerBounds[4] = float(eval(self.objective_lowerbound_fields[4].GetValue()))
        self.objective_lowerbound_fields[4].SetValue(str(self.plotoptions.LowerBounds[4]))

    def ChangeSUpperBound(self, event):
        event.Skip()
        self.plotoptions.UpperBounds[4] = float(eval(self.objective_upperbound_fields[4].GetValue()))
        self.objective_upperbound_fields[4].SetValue(str(self.plotoptions.UpperBounds[4]))

    def ChangeTimeUnit(self, event):
        self.plotoptions.TimeUnit = self.cmbTimeUnit.GetSelection()

    def ChangeEpochUnit(self, event):
        self.plotoptions.EpochUnit = self.cmbEpochUnit.GetSelection()

    def ChangeMarkerSize(self, event):
        self.plotoptions.BaseMarkerSize = self.spnctrlMarkerSizeControl.GetValue()

    def ChangeFontSize(self, e):
        self.plotoptions.FontSize = self.spnctrlFontSizeControl.GetValue()

    def ClickPlotPopulation(self, event):
        #first assemble the ordered list of objectives
        #note that if C is set but not Z, throw an error

        if self.Cobjective < len(self.NSGAIIpopulation.objective_column_headers) - 1 and self.Zobjective == len(self.NSGAIIpopulation.objective_column_headers):
            errordlg = wx.MessageDialog(self, "You cannot set the color axis without setting the Z axis first", "EMTG Error", wx.OK)
            errordlg.ShowModal()
            errordlg.Destroy()

        if self.Sobjective < len(self.NSGAIIpopulation.objective_column_headers) - 1 and self.Cobjective == len(self.NSGAIIpopulation.objective_column_headers):
            errordlg = wx.MessageDialog(self, "You cannot set the size axis without setting the color axis first", "EMTG Error", wx.OK)
            errordlg.ShowModal()
            errordlg.Destroy()

        else:
            ordered_list_of_objectives = [self.Xobjective, self.Yobjective]

            if self.Zobjective < len(self.NSGAIIpopulation.objective_column_headers):
                ordered_list_of_objectives.append(self.Zobjective)

            if self.Cobjective < len(self.NSGAIIpopulation.objective_column_headers):
                ordered_list_of_objectives.append(self.Cobjective)

            if self.Sobjective < len(self.NSGAIIpopulation.objective_column_headers):
                ordered_list_of_objectives.append(self.Sobjective)

            #check for duplicate objectives. If present, throw an error
            s = set()
            if any(obj in s or s.add(obj) for obj in ordered_list_of_objectives):
                errordlg = wx.MessageDialog(self, "Objective axes must be unique", "EMTG Error", wx.OK)
                errordlg.ShowModal()
                errordlg.Destroy()
            
            else:
                if not self.filter == []:
                    filter_indices = self.lbxObjectFilterList.GetSelections()
                    self.filter.filter_list = []
                    for index in filter_indices:
                        self.filter.filter_list.append(self.filter.object_list[index])
                #now make the plot
                self.NSGAIIpopulation.plot_population(ordered_list_of_objectives, LowerBounds = self.plotoptions.LowerBounds, UpperBounds = self.plotoptions.UpperBounds, TimeUnit = self.plotoptions.TimeUnit, EpochUnit = self.plotoptions.EpochUnit, FontSize = self.plotoptions.FontSize, BaseMarkerSize = self.plotoptions.BaseMarkerSize, OutputWindow = self.txtOutputTextWindow, PopulationFigure = self.PopulationFigure, PopulationAxes = self.PopulationAxes, Filter = self.filter)
                self.PlotCanvas.draw()

    def ClickClearOutputTextWindow(self, e):
        e.Skip()
        self.txtOutputTextWindow.Clear()


    def ChangeObjectFilterFile(self, e):
        e.Skip()
        self.ObjectFilterFile = self.txtObjectFilterFile.GetValue()
        self.filter = NSGAIIFilter.NSGAIIFilter(ObjectFilterFile)
        self.lbxObjectFilterList.Clear()
        self.lbxObjectFilterList.InsertItems(self.filter.object_list, 0)
        
    def ClickObjectFilterFileButton(self, e):
        dlg = wx.FileDialog(self, "Select an NSGAII object filter file", self.GetParent().dirname, "", '*.NSGAIIfilter', wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.ObjectFilterFile = dlg.GetPath()
            self.txtObjectFilterFile.SetValue(self.ObjectFilterFile)
            self.filter = NSGAIIFilter.NSGAIIFilter(self.ObjectFilterFile)
            self.lbxObjectFilterList.Clear()
            self.lbxObjectFilterList.InsertItems(self.filter.object_list, 0)
        dlg.Destroy()
        
