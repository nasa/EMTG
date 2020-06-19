import MissionOptions as MO
import Mission
import Universe
import BodyPicker
import Archive
import ArchivePanel
import NSGAIIpopulation
import NSGAIIpanel
import numpy as np
import OptionsNotebook
import UniverseNotebook
import MissionPanel
import webbrowser
import os
import subprocess
import platform
import wx
import wx.adv
import copy

class PyEMTG_interface(wx.Frame):
    
    def __init__(self, *args, **kwargs):
        super(PyEMTG_interface, self).__init__(*args, **kwargs)

        
        self.homedir = os.path.dirname(__file__)
        self.emtgpath = ''
        self.filename = ''
        self.dirname = ''
        self.mode = ''
        self.de_file = ''
        self.leapseconds_file = ''
        
        self.default_universe_path = ""

        self.read_PyEMTG_options_file()
        
        icon = wx.Icon("clemonaut.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        self.initialize_GUI()
        self.Maximize()

    def read_PyEMTG_options_file(self):
        input_file_name = os.path.join(self.homedir, "PyEMTG.options")
        self.default_thruster_file = " "
        self.default_small_bodies_file = None
        if os.path.isfile(input_file_name):
            inputfile = open(input_file_name, "r")
        else:
            print("Unable to open", input_file_name, "EMTG Error")
            return

        for line in inputfile:
            linecell = line.split(' ')

            if linecell[0] == "EMTG_path":
                self.emtgpath = linecell[1]
                if len(linecell) > 2:
                    for term in linecell[2:]:
                        self.emtgpath += ' ' + term
                self.emtgpath = self.emtgpath.strip('\n')

            elif linecell[0] == "default_universe_path":
                self.default_universe_path = linecell[1]
                if len(linecell) > 2:
                    for term in linecell[2:]:
                        self.default_universe_path += ' ' + term
                self.default_universe_path = self.default_universe_path.strip('\n')

            elif linecell[0] == "de_file":
                self.de_file = linecell[1]
                if len(linecell) > 2:
                    for term in linecell[2:]:
                        self.de_file += ' ' + term
                self.de_file = self.de_file.strip('\n')                

            elif linecell[0] == "leapseconds_file":
                self.leapseconds_file = linecell[1]
                if len(linecell) > 2:
                    for term in linecell[2:]:
                        self.leapseconds_file += ' ' + term
                self.leapseconds_file = self.leapseconds_file.strip('\n')

            elif linecell[0] == "default_small_bodies_file":
                self.default_small_bodies_file = linecell[1]
                if len(linecell) > 2:
                    for term in linecell[2:]:
                        self.default_small_bodies_file += ' ' + term
                self.default_small_bodies_file = self.default_small_bodies_file.strip('\n')

            elif linecell[0] == "default_HardwarePath":
                self.default_HardwarePath = linecell[1]
                if len(linecell) > 2:
                    for term in linecell[2:]:
                        self.default_HardwarePath += ' ' + term
                self.default_HardwarePath = self.default_HardwarePath.strip('\n')

            elif linecell[0] == "default_ThrottleTableFile":
                self.default_ThrottleTableFile = linecell[1]
                if len(linecell) > 2:
                    for term in linecell[2:]:
                        self.default_ThrottleTableFile += ' ' + term
                self.default_ThrottleTableFile = self.default_ThrottleTableFile.strip('\n')
                            
            elif linecell[0] == "default_LaunchVehicleLibraryFile":
                self.default_LaunchVehicleLibraryFile = linecell[1]
                if len(linecell) > 2:
                    for term in linecell[2:]:
                        self.default_LaunchVehicleLibraryFile += ' ' + term
                self.default_LaunchVehicleLibraryFile = self.default_LaunchVehicleLibraryFile.strip('\n')
                            
            elif linecell[0] == "default_PowerSystemsLibraryFile":
                self.default_PowerSystemsLibraryFile = linecell[1]
                if len(linecell) > 2:
                    for term in linecell[2:]:
                        self.default_PowerSystemsLibraryFile += ' ' + term
                self.default_PowerSystemsLibraryFile = self.default_PowerSystemsLibraryFile.strip('\n')
                            
            elif linecell[0] == "default_PropulsionSystemsLibraryFile":
                self.default_PropulsionSystemsLibraryFile = linecell[1]
                if len(linecell) > 2:
                    for term in linecell[2:]:
                        self.default_PropulsionSystemsLibraryFile += ' ' + term
                self.default_PropulsionSystemsLibraryFile = self.default_PropulsionSystemsLibraryFile.strip('\n')
                            
            elif linecell[0] == "default_SpacecraftOptionsFile":
                self.default_SpacecraftOptionsFile = linecell[1]
                if len(linecell) > 2:
                    for term in linecell[2:]:
                        self.default_SpacecraftOptionsFile += ' ' + term
                self.default_SpacecraftOptionsFile = self.default_SpacecraftOptionsFile.strip('\n')

        inputfile.close()
    
    def reset_file_menu(self):
        self.menubar = wx.MenuBar()
        self.fileMenu = wx.Menu()
        
        newMenu = wx.Menu()
        fnewmission = newMenu.Append(wx.ID_NEW, '&Mission\tCtrl+m')
        fnewuniverse = newMenu.Append(wx.ID_ANY, '&Universe\tCtrl+u')
        
        self.fileMenu.Append(wx.ID_ANY, '&New', newMenu)
        fopen = self.fileMenu.Append(wx.ID_OPEN, '&Open\tCtrl+o')
        fsave = self.fileMenu.Append(wx.ID_SAVE, '&Save\tCtrl+s')
        self.fileMenu.AppendSeparator()
        frun = self.fileMenu.Append(wx.ID_ANY, '&Run\tCtrl+r')
        self.fileMenu.AppendSeparator()
        fedit = self.fileMenu.Append(wx.ID_EDIT, 'Open file in &Editor\tCtrl+e')        
        self.fileMenu.AppendSeparator()
        fexit = self.fileMenu.Append(wx.ID_EXIT, 'E&xit\tCtrl+q')
        self.menubar.Append(self.fileMenu, '&File')
        self.SetMenuBar(self.menubar)
        
        self.Bind(wx.EVT_MENU, self.OnNewMission, fnewmission, id=wx.ID_NEW)
        self.Bind(wx.EVT_MENU, self.OnNewUniverse, fnewuniverse, id=wx.ID_NEW)

        self.Bind(wx.EVT_MENU, self.OnOpen, fopen, id=wx.ID_OPEN)
        self.Bind(wx.EVT_MENU, self.OnSave, fsave, id=wx.ID_SAVE)
        
        self.Bind(wx.EVT_MENU, self.OnRun, frun, id=wx.ID_ANY)

        self.Bind(wx.EVT_MENU, self.OnEdit, fedit, id=wx.ID_EDIT)
        
        self.Bind(wx.EVT_MENU, self.OnExit, fexit, id=wx.ID_EXIT)
        self.Bind(wx.EVT_CLOSE, self.OnExit)

        self.Bind(wx.EVT_SIZE, self.OnResize)
        self.Bind(wx.EVT_MAXIMIZE, self.OnResize)
                        
    def initialize_GUI(self):
        import random

        self.reset_file_menu()
        
        #Disable various menu options at program start
        self.fileMenu.Enable(wx.ID_SAVE, False)
        self.fileMenu.Enable(wx.ID_EDIT, False)
        
        #Otherwise load welcome message
        messages = ["Are thy wings plumed indeed for such far flights? \n -Walt Whitman",
                    "Up from Earth's Centre, through the Seventh Gate I rose, and on the Throne of Saturn sate. \n -Rubáiyát of Omar Khayyám",
                    "I was thinking this globe enough till there sprang out \n so noiseless around me myriads of other globes. \n -Walt Whitman",
                    "Every generation has the obligation to free men's minds for a look at new worlds...\n to look out from a higher plateau than the last generation. \n -Ellison Onizuka",
                    "Ah, who shall soothe these feverish children? \n Who will justify these restless explorations? \n -Walt Whitman",
                    "I'm not insane sir. I have a finely calibrated sense of acceptable risk. \n -John Scalzi",
                    "To go faster, slow down. Everybody who knows about orbital mechanics understands that. \n -Scott Cherf",
                    "When I’m working on a problem, I never think about beauty. I think only how to solve the problem. \n But when I have finished, if the solution is not beautiful, I know it is wrong. \n -Freeman Dyson",
                    "I am putting myself to the fullest possible use, \n which is all I think that any conscious entity can ever hope to do. \n -EMTG (HAL 9000)",
                    "I was taught that the way of progress was neither swift nor easy. \n - Marie Curie",
                    "Num num num. \n - Curtis Englander",
                    "The world is my toilet. \n - Cinnamon Englander",
                    "I am Willis, destroyer of worlds. Or at least your furniture. \n - Willis Englander",
                    "I have come to reap your socks. \n - The Socknal Reaper"]
        self.lblWelcome = wx.StaticText(self, -1, random.choice(messages), style = wx.ALIGN_CENTER)
        font = self.GetFont()
        font.SetPointSize(20)
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.lblWelcome.SetFont(font)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.lblWelcome, 1, flag = wx.CENTER)
        self.SetSizer(sizer)

        #self.SetSize((800,600))
        self.SetTitle("EMTG Python Interface")

        self.Show()

    def OnResize(self, e):
        self.Layout()
        MySize = self.GetClientSize()
        if self.mode == "options":
            self.optionsnotebook.SetSize((MySize.x, MySize.y))
            self.optionsnotebook.Layout()

        elif self.mode == "mission":
            self.missionpanel.SetSize((MySize.x, MySize.y))
            self.missionpanel.Layout()

        elif self.mode == "universe":
            self.universenotebook.SetSize((MySize.x, MySize.y))
            self.universenotebook.Layout()

        elif self.mode == "archive":
            self.archivepanel.SetSize((MySize.x, MySize.y))
            self.archivepanel.Layout()

        elif self.mode == "NSGAII":
            self.NSGAIIpanel.SetSize((MySize.x, MySize.y))
            self.NSGAIIpanel.Layout()
            
    def OnNewMission(self, e):
        #If the GUI has not yet loaded anything then create a new mission. Otherwise as for permission first.
        if self.mode == "":
            self.dirname = self.homedir
            self.filename = "default.emtgopt"
            self.missionoptions = MO.MissionOptions(os.path.join(self.dirname, self.filename))
            
            self.mode = "options"
            self.fileMenu.Enable(wx.ID_SAVE, True)
            self.fileMenu.Enable(wx.ID_EDIT, True)
            self.lblWelcome.Show(False)
            self.InitializeMissionOptionsEditor()
        else:
            dlg = wx.MessageDialog(self,
                    "Do you really want to create a new mission? This will clear your current GUI settings.",
                    "Confirm New", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
            result = dlg.ShowModal()
            dlg.Destroy()
            if result == wx.ID_OK:
                #destroy previous options memory block if you have one
                if self.mode == "options":
                    self.mode = ""
                    self.missionoptions = []
                    self.optionsnotebook.Destroy()

                elif self.mode == "mission":
                    self.mode = ""
                    self.mission = []
                    self.missionpanel.Destroy()

                elif self.mode == "universe":
                    self.mode = ""
                    self.universe = []
                    self.universenotebook.Destroy()

                elif self.mode == "archive":
                    self.mode = ""
                    self.archive = []
                    self.archivepanel.Destroy()

                elif self.mode == "NSGAII":
                    self.mode = ""
                    self.NSGAII = []
                    self.NSGAIIpanel.Destroy()

                #attempt to open a new options file
                self.dirname = self.homedir
                self.filename = "default.emtgopt"
                self.missionoptions = MO.MissionOptions()
                
                if self.missionoptions.success == 1:
                    self.mode = "options"
                    self.fileMenu.Enable(wx.ID_SAVE, True)
                    self.fileMenu.Enable(wx.ID_EDIT, True)
                    self.lblWelcome.Show(False)
                    self.InitializeMissionOptionsEditor()
                    
    def OnNewUniverse(self, e):
        if self.mode == "":
            self.dirname = self.homedir
            self.filename = "default.emtg_universe"
            self.universe = Universe.Universe(os.path.join(self.dirname, self.filename))
            
            if self.universe.success == 1:
                self.mode = "universe"
                self.fileMenu.Enable(wx.ID_SAVE, True)
                self.fileMenu.Enable(wx.ID_EDIT, True)
                self.lblWelcome.Show(False)
                self.InitializeUniverseOptionsEditor()
        else:
            dlg = wx.MessageDialog(self,
                    "Do you really want to create a new universe? This will clear your current GUI settings.",
                    "Confirm New", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
            result = dlg.ShowModal()
            dlg.Destroy()
            if result == wx.ID_OK:
                #destroy previous options memory block if you have one
                if self.mode == "options":
                    self.mode = ""
                    self.missionoptions = []
                    self.optionsnotebook.Destroy()

                elif self.mode == "mission":
                    self.mode = ""
                    self.mission = []
                    self.missionpanel.Destroy()

                elif self.mode == "universe":
                    self.mode = ""
                    self.universe = []
                    self.universenotebook.Destroy()

                elif self.mode == "archive":
                    self.mode = ""
                    self.archive = []
                    self.archivepanel.Destroy()

                elif self.mode == "NSGAII":
                    self.mode = ""
                    self.NSGAII = []
                    self.NSGAIIpanel.Destroy()

                #attempt to open a new universe file
                self.dirname = self.homedir
                self.filename = "default.emtg_universe"
                self.universe = Universe.Universe(os.path.join(self.dirname, self.filename))
                
                if self.universe.success == 1:
                    self.mode = "universe"
                    self.fileMenu.Enable(wx.ID_SAVE, True)
                    self.fileMenu.Enable(wx.ID_EDIT, True)
                    self.lblWelcome.Show(False)
                    self.InitializeUniverseOptionsEditor()
        
        
    def OnOpen(self, e):
        self.OpenFile(e)
        
    def OnSave(self, e):
        self.SetFocus()
        self.SaveFile(e)

    def OnRun(self, e):
        self.SetFocus()
        saved = self.SaveFile(e)

        if saved:
            if platform.system() == 'Windows':
                os.system('start ' + self.emtgpath.strip() + ' ' + os.path.join(self.dirname, self.filename))
            else:
                pyemtg_dir = os.path.split( os.path.realpath(__file__) )[0]
                emtg_dir = pyemtg_dir.replace('/PyEMTG'.replace('/',os.sep),'')
                emtg_exec = os.path.join(emtg_dir,'emtg')
                full_file_name = os.path.join(self.dirname,self.filename)
                str_cmd = emtg_exec + ' ' + full_file_name
                print(str_cmd)
                proc = subprocess.call(str_cmd,shell=True)
                if proc != 0:
                    print("call failed")
        
    def OnEdit(self, e):
        webbrowser.open(os.path.join(self.dirname, self.filename))
        
    def OnExit(self, e):
        self.SetFocus()
        if self.mode == "options" or self.mode == "universe":
            dlg = wx.MessageDialog(self,
            "Save changes?",
            "Confirm Exit", wx.YES|wx.NO|wx.CANCEL|wx.ICON_QUESTION)
            result = dlg.ShowModal()
            dlg.Destroy()
            if result == wx.ID_NO:
                self.Destroy()
            elif result == wx.ID_YES:
                saved = self.SaveFile(e)
                if saved:
                    self.Destroy()
        
        else:
            dlg = wx.MessageDialog(self,
                "Do you really want to exit?",
                "Confirm Exit", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
            result = dlg.ShowModal()
            dlg.Destroy()
            if result == wx.ID_OK:
                self.Destroy()
        
    
    def OpenFile(self, e):
        dlg = wx.FileDialog(self, message="Open an EMTG file", defaultDir=self.dirname, defaultFile="",
                            wildcard="*.emtgopt;*.emtg_universe;*.emtg;*.NSGAII", style=wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()

            fileparts = self.filename.split(".")
        

            #before we actually open the new file, we need to clear memory associated with whatever file we currently have open
            if self.mode == "options":
                self.missionoptions = []
                self.optionsnotebook.Destroy()
                
            elif self.mode == "mission":
                self.mission = []
                self.missionpanel.Destroy()
                

            elif self.mode == "universe":
                self.universe = []
                self.universenotebook.Destroy()

            elif self.mode == "archive":
                self.archive = []
                self.archivepanel.Destroy()

            elif self.mode == "NSGAII":
                self.NSGAII = []
                self.NSGAIIpanel.Destroy()

            self.mode = ""

            #next open the new file
            if fileparts[1] == "emtgopt":
                

                import sys
                import inspect
                currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
                sys.path.append(currentdir + "/" + 'Converters')
                from Convert_emtgopt_v1_to_v2 import Convert_emtgopt_v1_to_v2

                self.missionoptions = Convert_emtgopt_v1_to_v2(os.path.join(self.dirname, self.filename))

                if self.missionoptions.success == 1:
                    self.mode = "options"
                    self.lblWelcome.Show(False)
                    self.InitializeMissionOptionsEditor()
                    self.fileMenu.Enable(wx.ID_SAVE, True)
                    self.fileMenu.Enable(wx.ID_EDIT, True)

            elif fileparts[1] == "emtg":
                self.mission = Mission.Mission(os.path.join(self.dirname, self.filename))
                if self.mission.success == 1:
                    self.mode = "mission"
                    self.lblWelcome.Show(False)
                    self.missionpanel = MissionPanel.MissionPanel(self, self.mission)
                    self.missionpanel.SetSize(self.GetSize())
                    self.fileMenu.Enable(wx.ID_EDIT, True)
        
            elif fileparts[1] == "emtg_universe":
                self.universe = Universe.Universe(os.path.join(self.dirname, self.filename))
                if self.universe.success == 1:
                    self.mode = "universe"
                    self.lblWelcome.Show(False)
                    self.InitializeUniverseOptionsEditor()
                    self.fileMenu.Enable(wx.ID_SAVE, True)
                    self.fileMenu.Enable(wx.ID_EDIT, True)

            elif fileparts[1] == "emtg_archive":
                self.archive = Archive.Archive(os.path.join(self.dirname, self.filename))
                if self.archive.success == 1:
                    self.mode = "archive"
                    self.lblWelcome.Show(False)
                    self.InitializeArchiveProcessor()
                    self.fileMenu.Enable(wx.ID_EDIT, True)
                    self.fileMenu.Enable(wx.ID_SAVE, False)

            elif fileparts[1] == "NSGAII":
                self.NSGAII = NSGAIIpopulation.NSGAII_outerloop_population(os.path.join(self.dirname, self.filename))
                if self.NSGAII.success == 1:
                    self.mode = "NSGAII"
                    self.lblWelcome.Show(False)
                    self.InitializeNSGAIIPanel()
                    self.fileMenu.Enable(wx.ID_EDIT, True)
                    self.fileMenu.Enable(wx.ID_SAVE, False)

            else:
                errordlg = wx.MessageDialog(self, "Unrecognized file type.", "EMTG Error", wx.OK)
                errordlg.ShowModal()
                errordlg.Destroy()

        dlg.Destroy()


    
    def SaveFile(self, e):
        extension = []
        if self.mode == "options":
            dialogtitle = self.missionoptions.mission_name
            ext = ".emtgopt"
            extension = "*" + ext
        elif self.mode == "universe":
            dialogtitle = self.universe.central_body_name
            ext = ".emtg_universe"
            extension = "*" + ext

        dlg = wx.FileDialog(self, "Save", self.dirname, dialogtitle + ext, extension, wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            saved = True
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            
            if self.mode == "options":
                self.missionoptions.write_options_file(os.path.join(self.dirname, self.filename), not self.missionoptions.print_only_non_default_options)
            if self.mode == "universe":
                self.universe.write_universe_file(os.path.join(self.dirname, self.filename))
        else:
            saved = False
        
        dlg.Destroy()

        return saved

    def InitializeUniverseOptionsEditor(self):
        #create and size a universe notebook object
        self.universenotebook = UniverseNotebook.UniverseNotebook(self, self.universe)
        self.universenotebook.SetSize(self.GetClientSize())

    def InitializeArchiveProcessor(self):
        #create and size an Archive panel object
        self.archivepanel = ArchivePanel.ArchivePanel(self, self.archive)
        self.archivepanel.SetSize(self.GetSize())

    def InitializeNSGAIIPanel(self):
        self.NSGAIIpanel = NSGAIIpanel.NSGAIIpanel(self, self.NSGAII)
        self.NSGAIIpanel.SetSize(self.GetSize())
        
    def InitializeMissionOptionsEditor(self):
        #create an options notebook object
        self.optionsnotebook = OptionsNotebook.OptionsBook(self, self.missionoptions)
        
        #update the latest information from the missionoptions object into the GUI
        self.optionsnotebook.update()
        self.optionsnotebook.SetSize(self.GetClientSize())