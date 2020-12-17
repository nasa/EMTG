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

class SolverOptionsPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, missionoptions):
        self.missionoptions = missionoptions
        self.parent = parent

        self.only_call_once = False

        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)
        
        innerloopgrid = wx.FlexGridSizer(50,2,5,5)
        
        self.lblInnerLoopSolver = wx.StaticText(self, -1, "Inner-loop Solver Mode")
        innerloopsolvertypes = ['Evaluate trialX', 'Monotonic Basin Hopping',
                                'Adaptive Constrained Differential Evolution (NOT IMPLEMNTED)','NLP with initial guess',
                                'Filament Walker (experimental)']
        self.cmbInnerLoopSolver = wx.ComboBox(self, -1, choices = innerloopsolvertypes, style=wx.CB_READONLY)

        self.lblNLP_solver_type = wx.StaticText(self, -1, "NLP solver")
        NLP_solver_types = ['SNOPT','WORHP']
        self.cmbNLP_solver_type = wx.ComboBox(self, -1, choices = NLP_solver_types, style=wx.CB_READONLY)

        self.lblNLP_solver_mode = wx.StaticText(self, -1, "NLP solver mode")
        NLP_solver_modes = ['Feasible point','Optimize','Satisfy equality constraints']
        self.cmbNLP_solver_mode = wx.ComboBox(self, -1, choices = NLP_solver_modes, style = wx.CB_READONLY)

        self.lblquiet_NLP = wx.StaticText(self, -1, "Quiet NLP solver?")
        self.chkquiet_NLP = wx.CheckBox(self, -1)

        self.lblenable_NLP_chaperone = wx.StaticText(self, -1, "Enable NLP chaperone?")
        self.chkenable_NLP_chaperone = wx.CheckBox(self, -1)
        
        self.lblNLP_stop_on_goal_attain = wx.StaticText(self, -1, "Stop NLP upon attaining goal?")
        self.chkNLP_stop_on_goal_attain = wx.CheckBox(self, -1)

        self.lblNLP_objective_goal = wx.StaticText(self, -1, "NLP objective goal")
        self.txtNLP_objective_goal = wx.TextCtrl(self, -1, "NLP_objective_goal")

        self.lblquiet_MBH = wx.StaticText(self, -1, "Quiet MBH solver?")
        self.chkquiet_MBH = wx.CheckBox(self, -1)
        
        self.lblACE_feasible_point_finder = wx.StaticText(self, -1, "Enable ACE feasible point finder?")
        self.chkACE_feasible_point_finder = wx.CheckBox(self, -1)
        
        self.lblMBH_max_not_improve = wx.StaticText(self, -1, "MBH Impatience")
        self.txtMBH_max_not_improve = wx.TextCtrl(self, -1, "MBH_max_not_improve")
        
        self.lblMBH_max_trials = wx.StaticText(self, -1, "Maximum number of innerloop trials")
        self.txtMBH_max_trials = wx.TextCtrl(self, -1, "MBH_max_trials")
        
        self.lblMBH_max_run_time = wx.StaticText(self, -1, "Maximum run-time (s)")
        self.txtMBH_max_run_time = wx.TextCtrl(self, -1, "MBH_max_run_time")

        self.lblMBH_hop_distribution = wx.StaticText(self, -1, "MBH hop probability distribution")
        hop_distribution_choices = ["Uniform","Cauchy","Pareto","Gaussian"]
        self.cmbMBH_hop_distribution = wx.ComboBox(self, -1, choices = hop_distribution_choices, style = wx.CB_READONLY)
        
        self.lblMBH_max_step_size = wx.StaticText(self, -1, "MBH maximum perturbation size")
        self.txtMBH_max_step_size = wx.TextCtrl(self, -1, "MBH_max_step_size")

        self.lblMBH_Pareto_alpha = wx.StaticText(self, -1, "MBH Pareto distribution alpha")
        self.txtMBH_Pareto_alpha = wx.TextCtrl(self, -1, "MBH_Pareto_alpha")

        self.lblMBH_time_hop_probability = wx.StaticText(self, -1, "Probability of MBH time hop")
        self.txtMBH_time_hop_probability = wx.TextCtrl(self, -1, "MBH_time_hop_probability")

        self.lblMBH_always_write_archive = wx.StaticText(self, -1, "Always write MBH archive file?")
        self.chkMBH_always_write_archive = wx.CheckBox(self, -1)

        self.lblMBH_write_every_improvement = wx.StaticText(self, -1, "Write output file for all MBH improvements? (for later animation)")
        self.chkMBH_write_every_improvement = wx.CheckBox(self, -1)

        self.lblprint_NLP_movie_frames = wx.StaticText(self, -1, "Print NLP movie frames at every major iteration?")
        self.chkprint_NLP_movie_frames = wx.CheckBox(self, -1)
        
        self.lblsnopt_feasibility_tolerance = wx.StaticText(self, -1, "Feasibility tolerance")
        self.txtsnopt_feasibility_tolerance = wx.TextCtrl(self, -1, "snopt_feasibility_tolerance")
        
        self.lblsnopt_optimality_tolerance = wx.StaticText(self, -1, "Optimality tolerance")
        self.txtsnopt_optimality_tolerance = wx.TextCtrl(self, -1, "snopt_optimality_tolerance")

        self.lblNLP_max_step = wx.StaticText(self, -1, "NLP max step")
        self.txtNLP_max_step = wx.TextCtrl(self, -1, "NLP_max_step")

        self.lblsnopt_major_iterations = wx.StaticText(self, -1, "SNOPT major iterations limit")
        self.txtsnopt_major_iterations = wx.TextCtrl(self, -1, "snopt_major_iterations")
        self.lblsnopt_minor_iterations = wx.StaticText(self, -1, "SNOPT minor iterations limit")
        self.txtsnopt_minor_iterations = wx.TextCtrl(self, -1, "snopt_minor_iterations")
        
        self.lblsnopt_max_run_time = wx.StaticText(self, -1, "SNOPT maximum run time (s)")
        self.txtsnopt_max_run_time = wx.TextCtrl(self, -1, "snopt_max_run_time")
        
        self.lblcheck_derivatives = wx.StaticText(self, -1, "Check derivatives via finite differencing?")
        self.chkcheck_derivatives = wx.CheckBox(self, -1)
        
        self.lblseed_MBH = wx.StaticText(self, -1, "Seed MBH?")
        self.chkseed_MBH = wx.CheckBox(self, -1)
        self.lblskip_first_nlp_run = wx.StaticText(self, -1, "Skip first NLP run?")
        self.chkskip_first_nlp_run = wx.CheckBox(self, -1)
                
        innerloopgrid.AddMany(  [self.lblInnerLoopSolver, self.cmbInnerLoopSolver,
                                self.lblNLP_solver_type, self.cmbNLP_solver_type,
                                self.lblNLP_solver_mode, self.cmbNLP_solver_mode,
                                self.lblquiet_NLP, self.chkquiet_NLP,
                                self.lblenable_NLP_chaperone, self.chkenable_NLP_chaperone,
                                self.lblNLP_stop_on_goal_attain, self.chkNLP_stop_on_goal_attain,
                                self.lblNLP_objective_goal, self.txtNLP_objective_goal,
                                self.lblquiet_MBH, self.chkquiet_MBH,
                                self.lblACE_feasible_point_finder, self.chkACE_feasible_point_finder,
                                self.lblMBH_max_not_improve, self.txtMBH_max_not_improve,
                                self.lblMBH_max_trials, self.txtMBH_max_trials,
                                self.lblMBH_max_run_time, self.txtMBH_max_run_time,
                                self.lblMBH_hop_distribution, self.cmbMBH_hop_distribution,
                                self.lblMBH_max_step_size, self.txtMBH_max_step_size,
                                self.lblMBH_Pareto_alpha, self.txtMBH_Pareto_alpha,
                                self.lblMBH_time_hop_probability, self.txtMBH_time_hop_probability,
                                self.lblMBH_always_write_archive, self.chkMBH_always_write_archive,
                                self.lblMBH_write_every_improvement, self.chkMBH_write_every_improvement,
                                self.lblprint_NLP_movie_frames, self.chkprint_NLP_movie_frames,
                                self.lblsnopt_feasibility_tolerance, self.txtsnopt_feasibility_tolerance,
                                self.lblsnopt_optimality_tolerance, self.txtsnopt_optimality_tolerance,
                                self.lblNLP_max_step, self.txtNLP_max_step,
                                self.lblsnopt_major_iterations, self.txtsnopt_major_iterations,
                                self.lblsnopt_minor_iterations, self.txtsnopt_minor_iterations,
                                self.lblsnopt_max_run_time, self.txtsnopt_max_run_time,
                                self.lblcheck_derivatives, self.chkcheck_derivatives,
                                self.lblseed_MBH, self.chkseed_MBH,
                                self.lblskip_first_nlp_run, self.chkskip_first_nlp_run])
                                

                                
        vboxleft = wx.BoxSizer(wx.VERTICAL)
        lblLeftTitle = wx.StaticText(self, -1, "Inner-Loop Solver Parameters")
        vboxleft.Add(lblLeftTitle)
        vboxleft.Add(innerloopgrid)
        
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        lblLeftTitle.SetFont(font)

        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        hbox.Add(vboxleft)
        
        self.lbltrialX = wx.StaticText(self, -1, "Trial decision vector or initial guess")
        self.btntrialX = wx.Button(self, -1, "...")
    
        trialbox = wx.GridSizer(3,2,0,0)
        trialbox.AddMany([self.lbltrialX,self.btntrialX])
        
        self.mainbox = wx.BoxSizer(wx.VERTICAL)
        self.mainbox.AddMany([hbox, trialbox])

        self.SetSizer(self.mainbox)
        self.SetupScrolling()

        #bindings
        self.cmbInnerLoopSolver.Bind(wx.EVT_COMBOBOX, self.ChangeInnerLoopSolver)
        self.cmbNLP_solver_type.Bind(wx.EVT_COMBOBOX, self.ChangeNLP_solver_type)
        self.cmbNLP_solver_mode.Bind(wx.EVT_COMBOBOX, self.ChangeNLP_solver_mode)
        self.chkquiet_NLP.Bind(wx.EVT_CHECKBOX, self.Changequiet_NLP)
        self.chkenable_NLP_chaperone.Bind(wx.EVT_CHECKBOX, self.Changeenable_NLP_chaperone)
        self.chkNLP_stop_on_goal_attain.Bind(wx.EVT_CHECKBOX, self.ChangeNLP_stop_on_goal_attain)
        self.txtNLP_objective_goal.Bind(wx.EVT_KILL_FOCUS, self.ChangeNLP_objective_goal)
        self.chkquiet_MBH.Bind(wx.EVT_CHECKBOX, self.Changequiet_MBH)
        self.chkACE_feasible_point_finder.Bind(wx.EVT_CHECKBOX, self.ChangeACE_feasible_point_finder)
        self.txtMBH_max_not_improve.Bind(wx.EVT_KILL_FOCUS, self.ChangeMBH_max_not_improve)
        self.txtMBH_max_trials.Bind(wx.EVT_KILL_FOCUS, self.ChangeMBH_max_trials)
        self.txtMBH_max_run_time.Bind(wx.EVT_KILL_FOCUS, self.ChangeMBH_max_run_time)
        self.txtMBH_max_step_size.Bind(wx.EVT_KILL_FOCUS, self.ChangeMBH_max_step_size)
        self.cmbMBH_hop_distribution.Bind(wx.EVT_COMBOBOX, self.ChangeMBH_hop_distribution)
        self.txtMBH_Pareto_alpha.Bind(wx.EVT_KILL_FOCUS, self.ChangeMBH_Pareto_alpha)
        self.txtMBH_time_hop_probability.Bind(wx.EVT_KILL_FOCUS, self.ChangeMBH_time_hop_probability)
        self.chkMBH_always_write_archive.Bind(wx.EVT_CHECKBOX, self.ChangeMBH_always_write_archive)
        self.chkMBH_write_every_improvement.Bind(wx.EVT_CHECKBOX, self.ChangeMBH_write_every_improvement)
        self.chkprint_NLP_movie_frames.Bind(wx.EVT_CHECKBOX, self.Changeprint_NLP_movie_frames)
        self.txtsnopt_feasibility_tolerance.Bind(wx.EVT_KILL_FOCUS, self.Changesnopt_feasibility_tolerance)
        self.txtsnopt_optimality_tolerance.Bind(wx.EVT_KILL_FOCUS, self.Changesnopt_optimality_tolerance)
        self.txtNLP_max_step.Bind(wx.EVT_KILL_FOCUS, self.ChangeNLP_max_step)
        self.txtsnopt_major_iterations.Bind(wx.EVT_KILL_FOCUS, self.Changesnopt_major_iterations)
        self.txtsnopt_minor_iterations.Bind(wx.EVT_KILL_FOCUS, self.Changesnopt_minor_iterations)
        self.txtsnopt_max_run_time.Bind(wx.EVT_KILL_FOCUS, self.Changesnopt_max_run_time)
        self.chkcheck_derivatives.Bind(wx.EVT_CHECKBOX, self.ChangeCheckDerivatives)
        self.chkseed_MBH.Bind(wx.EVT_CHECKBOX, self.ChangeSeedMBH)
        self.chkskip_first_nlp_run.Bind(wx.EVT_CHECKBOX, self.Changeskip_first_nlp_run)
        
        self.btntrialX.Bind(wx.EVT_BUTTON, self.ClickTrialXButton)

    def update(self):
        import wx

        if self.missionoptions.snopt_max_run_time > self.missionoptions.MBH_max_run_time:
            self.missionoptions.snopt_max_run_time = self.missionoptions.MBH_max_run_time - 1

        #inner-loop solver options
        self.cmbInnerLoopSolver.SetSelection(self.missionoptions.run_inner_loop)
        self.cmbNLP_solver_type.SetSelection(self.missionoptions.NLP_solver_type)
        self.cmbNLP_solver_mode.SetSelection(self.missionoptions.NLP_solver_mode)
        self.chkquiet_NLP.SetValue(self.missionoptions.quiet_NLP)
        self.chkenable_NLP_chaperone.SetValue(self.missionoptions.enable_NLP_chaperone)
        self.chkNLP_stop_on_goal_attain.SetValue(self.missionoptions.NLP_stop_on_goal_attain)
        self.txtNLP_objective_goal.SetValue(str(self.missionoptions.NLP_objective_goal))
        self.chkquiet_MBH.SetValue(self.missionoptions.quiet_basinhopping)
        self.chkACE_feasible_point_finder.SetValue(self.missionoptions.ACE_feasible_point_finder)
        self.txtMBH_max_not_improve.SetValue(str(self.missionoptions.MBH_max_not_improve))
        self.txtMBH_max_trials.SetValue(str(self.missionoptions.MBH_max_trials))
        self.txtMBH_max_run_time.SetValue(str(self.missionoptions.MBH_max_run_time))
        self.txtMBH_max_step_size.SetValue(str(self.missionoptions.MBH_max_step_size))
        self.cmbMBH_hop_distribution.SetSelection(self.missionoptions.MBH_hop_distribution)
        self.txtMBH_Pareto_alpha.SetValue(str(self.missionoptions.MBH_Pareto_alpha))
        self.txtMBH_time_hop_probability.SetValue(str(self.missionoptions.MBH_time_hop_probability))
        self.chkMBH_always_write_archive.SetValue(self.missionoptions.MBH_always_write_archive)
        self.chkMBH_write_every_improvement.SetValue(self.missionoptions.MBH_write_every_improvement)
        self.chkprint_NLP_movie_frames.SetValue(self.missionoptions.print_NLP_movie_frames)
        self.txtsnopt_feasibility_tolerance.SetValue(str(self.missionoptions.snopt_feasibility_tolerance))
        self.txtsnopt_optimality_tolerance.SetValue(str(self.missionoptions.snopt_optimality_tolerance))
        self.txtNLP_max_step.SetValue(str(self.missionoptions.NLP_max_step))
        self.txtsnopt_major_iterations.SetValue(str(self.missionoptions.snopt_major_iterations))
        self.txtsnopt_minor_iterations.SetValue(str(self.missionoptions.snopt_minor_iterations))
        self.txtsnopt_max_run_time.SetValue(str(self.missionoptions.snopt_max_run_time))
        self.chkcheck_derivatives.SetValue(self.missionoptions.check_derivatives)
        self.chkseed_MBH.SetValue(self.missionoptions.seed_MBH)
        self.chkskip_first_nlp_run.SetValue(self.missionoptions.skip_first_nlp_run)

        if self.missionoptions.run_inner_loop == 0: #trialX
            self.lblMBH_max_not_improve.Show(False)
            self.lblMBH_max_trials.Show(False)
            self.lblMBH_max_run_time.Show(False)
            self.lblMBH_max_step_size.Show(False)
            self.lblMBH_time_hop_probability.Show(False)
            self.lblsnopt_feasibility_tolerance.Show(False)
            self.lblsnopt_optimality_tolerance.Show(False)
            self.lblsnopt_major_iterations.Show(False)
            self.lblsnopt_minor_iterations.Show(False)
            self.lblsnopt_max_run_time.Show(False)
            self.lblcheck_derivatives.Show(True)
            self.lblseed_MBH.Show(False)
            self.lblMBH_hop_distribution.Show(False)
            self.lblNLP_solver_type.Show(False)
            self.lblNLP_solver_mode.Show(False)
            self.lblquiet_NLP.Show(False)
            self.lblenable_NLP_chaperone.Show(False)
            self.lblNLP_stop_on_goal_attain.Show(False)
            self.lblNLP_objective_goal.Show(False)
            self.lblACE_feasible_point_finder.Show(False)

            self.txtMBH_max_not_improve.Show(False)
            self.txtMBH_max_trials.Show(False)
            self.txtMBH_max_run_time.Show(False)
            self.txtMBH_max_step_size.Show(False)
            self.txtMBH_time_hop_probability.Show(False)
            self.txtsnopt_feasibility_tolerance.Show(False)
            self.txtsnopt_optimality_tolerance.Show(False)
            self.txtsnopt_major_iterations.Show(False)
            self.txtsnopt_minor_iterations.Show(False)
            self.txtsnopt_max_run_time.Show(False)
            self.chkcheck_derivatives.Show(True)
            self.chkseed_MBH.Show(False)
            self.cmbMBH_hop_distribution.Show(False)
            self.cmbNLP_solver_type.Show(False)
            self.cmbNLP_solver_mode.Show(False)
            self.chkquiet_NLP.Show(False)
            self.chkenable_NLP_chaperone.Show(False)
            self.chkNLP_stop_on_goal_attain.Show(False)
            self.txtNLP_objective_goal.Show(False)
            self.chkACE_feasible_point_finder.Show(False)

            self.lblMBH_Pareto_alpha.Show(False)
            self.txtMBH_Pareto_alpha.Show(False)
            self.lblMBH_always_write_archive.Show(False)
            self.chkMBH_always_write_archive.Show(False)
            self.lblMBH_write_every_improvement.Show(False)
            self.chkMBH_write_every_improvement.Show(False)
            self.lblprint_NLP_movie_frames.Show(False)
            self.chkprint_NLP_movie_frames.Show(False)
            self.lblNLP_max_step.Show(False)
            self.txtNLP_max_step.Show(False)

            self.lblquiet_MBH.Show(False)
            self.chkquiet_MBH.Show(False)

            self.lblskip_first_nlp_run.Show(False)
            self.chkskip_first_nlp_run.Show(False)

        elif self.missionoptions.run_inner_loop == 1: #MBH
            self.lblMBH_max_not_improve.Show(True)
            self.lblMBH_max_trials.Show(True)
            self.lblMBH_max_run_time.Show(True)
            self.lblMBH_max_step_size.Show(True)
            self.lblMBH_time_hop_probability.Show(True)
            self.lblsnopt_feasibility_tolerance.Show(True)
            self.lblsnopt_optimality_tolerance.Show(True)
            self.lblsnopt_major_iterations.Show(True)
            self.lblsnopt_minor_iterations.Show(True)
            self.lblsnopt_max_run_time.Show(True)
            self.lblcheck_derivatives.Show(True)
            self.lblseed_MBH.Show(True)
            self.lblMBH_hop_distribution.Show(True)
            self.lblNLP_solver_type.Show(True)
            self.lblNLP_solver_mode.Show(True)
            self.lblquiet_NLP.Show(True)
            self.lblenable_NLP_chaperone.Show(True)
            self.lblACE_feasible_point_finder.Show(True)

            self.txtMBH_max_not_improve.Show(True)
            self.txtMBH_max_trials.Show(True)
            self.txtMBH_max_run_time.Show(True)
            self.txtMBH_max_step_size.Show(True)
            self.txtMBH_time_hop_probability.Show(True)
            self.txtsnopt_feasibility_tolerance.Show(True)
            self.txtsnopt_optimality_tolerance.Show(True)
            self.txtsnopt_major_iterations.Show(True)
            self.txtsnopt_minor_iterations.Show(True)
            self.txtsnopt_max_run_time.Show(True)
            self.chkcheck_derivatives.Show(True)
            self.chkseed_MBH.Show(True)
            self.cmbMBH_hop_distribution.Show(True)
            self.cmbNLP_solver_type.Show(True)
            self.cmbNLP_solver_mode.Show(True)
            self.chkquiet_NLP.Show(True)
            self.chkenable_NLP_chaperone.Show(True)
            self.chkACE_feasible_point_finder.Show(True)

            self.lblMBH_always_write_archive.Show(True)
            self.chkMBH_always_write_archive.Show(True)
            self.lblMBH_write_every_improvement.Show(True)
            self.chkMBH_write_every_improvement.Show(True)
            self.lblprint_NLP_movie_frames.Show(True)
            self.chkprint_NLP_movie_frames.Show(True)
            self.lblNLP_max_step.Show(True)
            self.txtNLP_max_step.Show(True)

            self.lblquiet_MBH.Show(True)
            self.chkquiet_MBH.Show(True)

            self.lblskip_first_nlp_run.Show(True)
            self.chkskip_first_nlp_run.Show(True)

            if self.missionoptions.enable_NLP_chaperone:
                self.lblNLP_stop_on_goal_attain.Show(True)
                self.chkNLP_stop_on_goal_attain.Show(True)
                if self.missionoptions.NLP_stop_on_goal_attain:
                    self.lblNLP_objective_goal.Show(True)
                    self.txtNLP_objective_goal.Show(True)
                else:
                    self.lblNLP_objective_goal.Show(False)
                    self.txtNLP_objective_goal.Show(False)
            else:                
                self.lblNLP_stop_on_goal_attain.Show(False)
                self.lblNLP_objective_goal.Show(False)
                self.chkNLP_stop_on_goal_attain.Show(False)
                self.txtNLP_objective_goal.Show(False)

            #change the available parameters and labels based on which distribution is selected
            if self.missionoptions.MBH_hop_distribution == 0: #uniform
                self.lblMBH_max_step_size.SetLabel("MBH uniform hop ball size")
                self.lblMBH_Pareto_alpha.Show(False)
                self.txtMBH_Pareto_alpha.Show(False)
            elif self.missionoptions.MBH_hop_distribution == 1: #Cauchy
                self.lblMBH_max_step_size.SetLabel("MBH hop scale factor")
                self.lblMBH_Pareto_alpha.Show(False)
                self.txtMBH_Pareto_alpha.Show(False)
            elif self.missionoptions.MBH_hop_distribution == 2: #Pareto
                self.lblMBH_max_step_size.SetLabel("MBH hop scale factor")
                self.lblMBH_Pareto_alpha.Show(True)
                self.txtMBH_Pareto_alpha.Show(True)
            elif self.missionoptions.MBH_hop_distribution == 3: #Gaussian
                self.lblMBH_max_step_size.SetLabel("MBH hop standard deviation")
                self.lblMBH_Pareto_alpha.Show(False)
                self.txtMBH_Pareto_alpha.Show(False)

        elif self.missionoptions.run_inner_loop == 2: #ACDE
            self.lblMBH_max_not_improve.Show(True)
            self.lblMBH_max_trials.Show(True)
            self.lblMBH_max_run_time.Show(True)
            self.lblMBH_max_step_size.Show(False)
            self.lblMBH_time_hop_probability.Show(False)
            self.lblsnopt_feasibility_tolerance.Show(True)
            self.lblsnopt_optimality_tolerance.Show(True)
            self.lblsnopt_major_iterations.Show(False)
            self.lblsnopt_minor_iterations.Show(False)
            self.lblsnopt_max_run_time.Show(False)
            self.lblcheck_derivatives.Show(False)
            self.lblseed_MBH.Show(False)
            self.lblMBH_hop_distribution.Show(False)
            self.lblNLP_solver_type.Show(False)
            self.lblNLP_solver_mode.Show(False)
            self.lblquiet_NLP.Show(False)
            self.lblenable_NLP_chaperone.Show(False)
            self.lblNLP_stop_on_goal_attain.Show(False)
            self.lblNLP_objective_goal.Show(False)
            self.lblACE_feasible_point_finder.Show(False)

            self.txtMBH_max_not_improve.Show(True)
            self.txtMBH_max_trials.Show(True)
            self.txtMBH_max_run_time.Show(True)
            self.txtMBH_max_step_size.Show(False)
            self.txtMBH_time_hop_probability.Show(False)
            self.txtsnopt_feasibility_tolerance.Show(True)
            self.txtsnopt_optimality_tolerance.Show(True)
            self.txtsnopt_major_iterations.Show(False)
            self.txtsnopt_minor_iterations.Show(False)
            self.txtsnopt_max_run_time.Show(False)
            self.chkcheck_derivatives.Show(False)
            self.chkseed_MBH.Show(False)
            self.cmbMBH_hop_distribution.Show(False)
            self.cmbNLP_solver_type.Show(False)
            self.cmbNLP_solver_mode.Show(False)
            self.chkquiet_NLP.Show(False)
            self.chkenable_NLP_chaperone.Show(False)
            self.chkNLP_stop_on_goal_attain.Show(False)
            self.txtNLP_objective_goal.Show(False)
            self.chkACE_feasible_point_finder.Show(False)

            self.lblMBH_Pareto_alpha.Show(False)
            self.txtMBH_Pareto_alpha.Show(False)
            self.lblMBH_write_every_improvement.Show(False)
            self.chkMBH_write_every_improvement.Show(False)
            self.lblprint_NLP_movie_frames.Show(False)
            self.chkprint_NLP_movie_frames.Show(False)
            self.lblNLP_max_step.Show(False)
            self.txtNLP_max_step.Show(False)

            self.lblquiet_MBH.Show(False)
            self.chkquiet_MBH.Show(False)

            self.lblskip_first_nlp_run.Show(False)
            self.chkskip_first_nlp_run.Show(False)

        elif self.missionoptions.run_inner_loop == 3: #NLP
            self.lblMBH_max_not_improve.Show(False)
            self.lblMBH_max_trials.Show(False)
            self.lblMBH_max_run_time.Show(False)
            self.lblMBH_max_step_size.Show(False)
            self.lblMBH_time_hop_probability.Show(False)
            self.lblsnopt_feasibility_tolerance.Show(True)
            self.lblsnopt_optimality_tolerance.Show(True)
            self.lblsnopt_major_iterations.Show(True)
            self.lblsnopt_minor_iterations.Show(True)
            self.lblsnopt_max_run_time.Show(True)
            self.lblcheck_derivatives.Show(True)
            self.lblseed_MBH.Show(False)
            self.lblMBH_hop_distribution.Show(False)
            self.lblNLP_solver_type.Show(True)
            self.lblNLP_solver_mode.Show(True)
            self.lblquiet_NLP.Show(True)
            self.lblenable_NLP_chaperone.Show(True)
            self.lblACE_feasible_point_finder.Show(False)

            self.txtMBH_max_not_improve.Show(False)
            self.txtMBH_max_trials.Show(False)
            self.txtMBH_max_run_time.Show(False)
            self.txtMBH_max_step_size.Show(False)
            self.txtMBH_time_hop_probability.Show(False)
            self.txtsnopt_feasibility_tolerance.Show(True)
            self.txtsnopt_optimality_tolerance.Show(True)
            self.txtsnopt_major_iterations.Show(True)
            self.txtsnopt_minor_iterations.Show(True)
            self.txtsnopt_max_run_time.Show(True)
            self.chkcheck_derivatives.Show(True)
            self.chkseed_MBH.Show(False)
            self.cmbMBH_hop_distribution.Show(False)
            self.cmbNLP_solver_type.Show(True)
            self.cmbNLP_solver_mode.Show(True)
            self.chkquiet_NLP.Show(True)
            self.chkenable_NLP_chaperone.Show(True)
            self.chkACE_feasible_point_finder.Show(False)

            self.lblMBH_Pareto_alpha.Show(False)
            self.txtMBH_Pareto_alpha.Show(False)
            self.lblMBH_write_every_improvement.Show(False)
            self.chkMBH_write_every_improvement.Show(False)
            self.lblprint_NLP_movie_frames.Show(True)
            self.chkprint_NLP_movie_frames.Show(True)
            self.lblNLP_max_step.Show(True)
            self.txtNLP_max_step.Show(True)

            self.lblquiet_MBH.Show(False)
            self.chkquiet_MBH.Show(False)

            self.lblskip_first_nlp_run.Show(False)
            self.chkskip_first_nlp_run.Show(False)
            
            if self.missionoptions.enable_NLP_chaperone:
                self.lblNLP_stop_on_goal_attain.Show(True)
                self.chkNLP_stop_on_goal_attain.Show(True)
                if self.missionoptions.NLP_stop_on_goal_attain:
                    self.lblNLP_objective_goal.Show(True)
                    self.txtNLP_objective_goal.Show(True)
                else:
                    self.lblNLP_objective_goal.Show(False)
                    self.txtNLP_objective_goal.Show(False)
            else:                
                self.lblNLP_stop_on_goal_attain.Show(False)
                self.lblNLP_objective_goal.Show(False)
                self.chkNLP_stop_on_goal_attain.Show(False)
                self.txtNLP_objective_goal.Show(False)


        if (self.missionoptions.run_inner_loop == 1 and self.missionoptions.seed_MBH == 1) or self.missionoptions.run_inner_loop in [0, 3]:
            self.lbltrialX.Show(True)
            self.btntrialX.Show(True)
        else:
            self.lbltrialX.Show(False)
            self.btntrialX.Show(False)

        #re-size the panel
        self.Layout()
        if platform.system() == 'Windows':
            self.SetupScrolling(scrollToTop=False)
            
    #event handlers for solver options
    def ChangeInnerLoopSolver(self, e):
        e.Skip()
        self.missionoptions.run_inner_loop = self.cmbInnerLoopSolver.GetSelection()
        self.update()

    def ChangeNLP_solver_type(self, e):
        e.Skip()
        self.missionoptions.NLP_solver_type = self.cmbNLP_solver_type.GetSelection()
        self.update()

    def ChangeNLP_solver_mode(self, e):
        e.Skip()
        self.missionoptions.NLP_solver_mode = self.cmbNLP_solver_mode.GetSelection()
        self.update()

    def Changequiet_NLP(self, e):
        e.Skip()
        self.missionoptions.quiet_NLP = int(self.chkquiet_NLP.GetValue())
        self.update()

    def Changeenable_NLP_chaperone(self, e):
        e.Skip()
        self.missionoptions.enable_NLP_chaperone = int(self.chkenable_NLP_chaperone.GetValue())
        self.update()

        
    def ChangeNLP_stop_on_goal_attain(self, e):
        e.Skip()
        self.missionoptions.NLP_stop_on_goal_attain = int(self.chkNLP_stop_on_goal_attain.GetValue())
        self.update()
        
    def ChangeNLP_objective_goal(self, e):
        e.Skip()
        self.missionoptions.NLP_objective_goal = eval(self.txtNLP_objective_goal.GetValue())
        self.update()

    def Changequiet_MBH(self, e):
        e.Skip()
        self.missionoptions.quiet_basinhopping = int(self.chkquiet_MBH.GetValue())
        self.update()

    def ChangeACE_feasible_point_finder(self, e):
        e.Skip()
        self.missionoptions.ACE_feasible_point_finder = int(self.chkACE_feasible_point_finder.GetValue())
        self.update()
        
    def ChangeMBH_max_not_improve(self, e):
        e.Skip()
        self.missionoptions.MBH_max_not_improve = eval(self.txtMBH_max_not_improve.GetValue())
        self.update()

    def ChangeMBH_max_trials(self, e):
        e.Skip()
        self.missionoptions.MBH_max_trials = eval(self.txtMBH_max_trials.GetValue())
        self.update()
        
    def ChangeMBH_max_run_time(self, e):
        e.Skip()
        self.missionoptions.MBH_max_run_time = eval(self.txtMBH_max_run_time.GetValue())
        self.update()
        
    def ChangeMBH_max_step_size(self, e):
        e.Skip()
        self.missionoptions.MBH_max_step_size = eval(self.txtMBH_max_step_size.GetValue())
        self.update()

    def ChangeMBH_hop_distribution(self, e):
        e.Skip()
        self.missionoptions.MBH_hop_distribution = self.cmbMBH_hop_distribution.GetSelection()
        self.update()

    def ChangeMBH_Pareto_alpha(self, e):
        e.Skip()
        self.missionoptions.MBH_Pareto_alpha = self.txtMBH_Pareto_alpha.GetValue()
        self.update()

    def ChangeMBH_time_hop_probability(self, e):
        e.Skip()
        self.missionoptions.MBH_time_hop_probability = eval(self.txtMBH_time_hop_probability.GetValue())
        self.update()

    def ChangeMBH_always_write_archive(self, e):
        e.Skip()
        self.missionoptions.MBH_always_write_archive = int(self.chkMBH_always_write_archive.GetValue())
        self.update()
     
    def ChangeMBH_write_every_improvement(self, e):
        e.Skip()
        self.missionoptions.MBH_write_every_improvement = int(self.chkMBH_write_every_improvement.GetValue())
        self.update()
     
    def Changeprint_NLP_movie_frames(self, e):
        e.Skip()
        self.missionoptions.print_NLP_movie_frames = int(self.chkprint_NLP_movie_frames.GetValue())
        self.update()

    def Changesnopt_feasibility_tolerance(self, e):
        e.Skip()
        self.missionoptions.snopt_feasibility_tolerance = eval(self.txtsnopt_feasibility_tolerance.GetValue())
        self.update()

    def Changesnopt_optimality_tolerance(self, e):
        e.Skip()
        self.missionoptions.snopt_optimality_tolerance = eval(self.txtsnopt_optimality_tolerance.GetValue())
        self.update()

    def ChangeNLP_max_step(self, e):
        e.Skip()
        self.missionoptions.NLP_max_step = eval(self.txtNLP_max_step.GetValue())
        self.update()
                        
    def Changesnopt_major_iterations(self, e):
        e.Skip()
        self.missionoptions.snopt_major_iterations = eval(self.txtsnopt_major_iterations.GetValue())
        self.update()
                        
    def Changesnopt_minor_iterations(self, e):
        e.Skip()
        self.missionoptions.snopt_minor_iterations = eval(self.txtsnopt_minor_iterations.GetValue())
        self.update()
        
    def Changesnopt_max_run_time(self, e):
        e.Skip()
        self.missionoptions.snopt_max_run_time = eval(self.txtsnopt_max_run_time.GetValue())
        self.update()

    def ChangeCheckDerivatives(self, e):
        e.Skip()
        self.missionoptions.check_derivatives = int(self.chkcheck_derivatives.GetValue())
        
    def ChangeSeedMBH(self, e):
        self.missionoptions.seed_MBH = int(self.chkseed_MBH.GetValue())
        self.update()

    def Changeskip_first_nlp_run(self, e):
        self.missionoptions.skip_first_nlp_run = int(self.chkskip_first_nlp_run.GetValue())
        self.update()

    def ClickTrialXButton(self, e):
        e.Skip()
        import Mission
        import os
        self.dirname = self.parent.Parent.dirname
        dlg = wx.FileDialog(self, "Choose a mission file", self.dirname, "", "*.emtg; *.emtg_initialguess", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()

            fileparts = self.filename.split(".")

            if fileparts[1] == "emtg":
                #read the file into a mission object and extract the trialX entry
                tempmission = Mission.Mission(os.path.join(self.dirname, self.filename))

                self.missionoptions.trialX = []
                for Xindex in range(0, len(tempmission.DecisionVector)):
                    self.missionoptions.trialX.append([tempmission.Xdescriptions[Xindex].replace('\n','').replace('\r',''), tempmission.DecisionVector[Xindex]])

                self.missionoptions.DisassembleMasterDecisionVector()
                self.missionoptions.ConvertDecisionVector()
                self.missionoptions.AssembleMasterDecisionVector()


                del tempmission


        self.update()