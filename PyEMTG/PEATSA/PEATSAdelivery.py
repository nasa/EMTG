#EMTG: Evolutionary Mission Trajectory Generator
#An open-source global optimization tool for preliminary mission design
#Provided by NASA Goddard Space Flight Center
#
#Copyright (c) 2014 - 2024 United States Government as represented by the
#Administrator of the National Aeronautics and Space Administration.
#All Other Rights Reserved.
#
#Licensed under the NASA Open Source License (the "License"); 
#You may not use this file except in compliance with the License. 
#You may obtain a copy of the License at:
#https://opensource.org/license/nasa1-3-php
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
#express or implied.   See the License for the specific language
#governing permissions and limitations under the License.

"""
PEATSAdelivery.py
==================

File holds the PEATSAdelivery class, which contains methods that make plots of PEATSA results.

"""

class PEATSAdelivery(object):
    """
    Holds the PEATSAdelivery class, which contains methods that make plots of PEATSA results.
    
    Parameters
    ----------
    None.

    Returns
    -------
    None.
    
    """
    
    def jd2datelist(self,dateval,ifMJD = 0):
        import math
        
        JD = dateval
        if ifMJD:
            JD += 2400000.5
        Q = JD+0.5
        Z = math.floor(Q)
        W = math.floor((Z - 1867216.25)/36524.25)
        X = math.floor(W/4.0)
        A = Z+1+W-X
        B = A+1524
        C = math.floor((B-122.1)/365.25)
        D = math.floor(365.25*C)
        E = math.floor((B-D)/30.6001)
        F = math.floor(30.6001*E)
        day  = B-D-F+(Q-Z)
        month = E-1
        if month >= 13:
            month -=12
        if month == 1 or month == 2:
            year = C-4715
        else:
            year = C-4716
            
        return int(year),int(month),int(day)
        
    def post_process(self,PEATSAorder,PEATSAboxes):
        import logging

        # Call the default plots if requested
        if PEATSAorder.generate_default_plots:
            if PEATSAorder.PEATSA_type == 1:
                logging.info("Generating default Missed Thrust plot")
                self.MissedThrustPlotter(PEATSAorder,PEATSAboxes)
            elif PEATSAorder.PEATSA_type == 2:
                logging.info("Generating default Trade Study plot")
                self.TradeStudyPlotter(PEATSAorder,PEATSAboxes)
            else:
                logging.info("Default plot for this peatsa type not created yet. Sorry.")
        
        # Call the built-in plotter if requested
        if len(PEATSAorder.built_in_plotter_files):
            for filename in PEATSAorder.built_in_plotter_files:
                self.BuiltInPlotterFromFile(PEATSAorder,PEATSAboxes,filename)
    
        # Call the custom post processing script if so desired    
        if PEATSAorder.custom_post_processing_script != "":
            self.Custom(PEATSAorder,PEATSAboxes)
        
        if PEATSAorder.move_results != "":
            self.moveFiles(PEATSAorder)
                        
    def TradeStudyPlotter(self,PEATSAorder,PEATSAboxes):
        import PEATSAchef
        
        # Create a peatsa chef
        chef = PEATSAchef.PEATSAchef()
        
        # Loop through all of the options files
        for options_file_tuple in PEATSAorder.trade_study_options_files:
            
            # Check how this options file is formatted
            if options_file_tuple[1] == 0:
                # Call the type zero trade study peatsa chef
                options_list = chef.TradeStudyTypeZero(PEATSAorder,options_file_tuple[0],1)
                
                # Call the type zero plotter
                self.TradeStudyTypeZeroPlotter(PEATSAorder,PEATSAboxes,options_list)
            elif options_file_tuple[1] == 1:
                logging.info("Sorry, the built in plotter not setup yet for this trade study type")
            elif options_file_tuple[1] == 2:
                logging.info("Sorry, the built in plotter not setup yet for this trade study type")
    
    def convert_variable_to_label(self,string):
        
        # Remove the Mission and or mission options 
        string = string.lstrip("MO").lstrip("M")
        
        # Check if the variable is on a journey or journey options object
        if string.startswith("Journeys["):
            
            # Get the journey number
            journey_no = string.split("[")[1].split("]")[0]
            
            # Check if the journey_number is negative
            if journey_no == "-1":
                journey_str = "Final_Journey"
            else:
                journey_str = "Journey_" + journey_no
            
            # replace the journey bracketing
            string = journey_str + "_" + string.split("eys[" + journey_no + "]")[1]
        
        # Check if the variable has a mission event object
        if ".missionevents[" in string:
            
            # Get the mission event
            mission_event_no = string.split("missionevents[")[1].split("]")[0]
            
            # Check if the mission event is the first or last
            if mission_event_no == "0":
                mission_event_string = ",_Initial"
            elif mission_event_no == "-1":
                mission_event_string = ",_Final"
            else:
                mission_event_string = ",_Mission_Event_" + mission_event_no
            
            # Replace the mission event bracketing
            string = string.split("_.missionevents[")[0] + mission_event_string + "_" + string.split("events[" + mission_event_no + "]")[1]
        
        # Convert underscores to spaces
        string = string.replace("_"," ")
        
        # Capitalize each word
        string = self.capitalize_each_word(string)
    
        # Return the formatted label
        return string
        
    def capitalize_each_word(self,string):
        # Split the string by spaces
        stringsplit = string.split(" ")
        
        # Start a new string and capitalize the first word
        new_string = string[0].capitalize()
        
        # Loop through the rest of the words
        for string_part in stringsplit:
            # Capitalize each and add to the new string
            new_string += " " + string_part.capitalize()
            
        return new_string
    
    def TradeStudyTypeZeroPlotter(self,PEATSAorder,PEATSAboxes,options_list):
            
        # Initialize a counter
        plot_ctr = 0
        
        # Loop through all of the options in the options list
        for option in options_list.keys():
            
            # Initialize a new list of boxes
            new_boxes = []
            
            # Loop through all of the boxes
            for box in PEATSAboxes:
                
                # Assume this box belongs on the plot
                box_belongs_on_plot = True
                
                # Loop through the options again
                for option2 in options_list.keys():
                    
                    # Check if the second option is the same as the first
                    if option == option2:
                        # Don't care
                        continue
                    
                    # Check if this box's value matches the default value
                    if eval("box.PEATSAdough." + option[0].lstrip("MO")) != options_list[option][0]:
                        # It does not
                        box_belongs_on_plot = False
                        # No need to keep looping
                        break
                        
                # Check if the box belongs on the plot
                if box_belongs_on_plot:
                    # It does, add it to the list
                    new_boxes.append(box)

            # Check if this is the first iteration or not
            if PEATSAorder.iteration == 0:
                # This is the first iteration
                
                # Initialize a set of plot options
                PO = PEATSAmenu.PEATSAinstagram()

                # Set the y_variable
                PO.y_variable = "box.objective_value"
                # Set the x_variable
                PO.x_variable = "box.PEATSAdough." + option[0].lstrip("MO").lstrip(".")
                # set the y label
                PO.ylabel = self.convert_variable_to_label(PEATSAorder.objective_formula)
                # set the x label
                PO.xlabel = self.convert_variable_to_label(option[0])
                # Set the title
                PO.title = "Trade Study"
                # Set the file name string
                PO.file_name_str = "default_" + option[0].lstrip("MO").lstrip(".") + "_sensitivity_plot"
                # Tell the plotter to use lines not dots
                PO.lines_or_dots = 'lines'
            
                # Check if the x variable is a date and set the flag if so
                if "date" in PO.xlabel:
                    PO.x_is_date = 1
                
                # Check if the x variable is a date and set the flag if so
                if "date" in PO.ylabel:
                    PO.y_is_date = 1
                
                # Write the options to file   
                PO.write_options(PEATSAorder.working_directory + "Default_TradeStudy_Plot" + str(plot_ctr) + "Options.py")
            else:
                # It is not the first iteration, so this plot should have midbake options saved to file. Load from there in case 
                # the user modified them
                
                # Create a peatsa plot object from the mid bake file
                PO = PEATSAmenu.PEATSAinstagram(PEATSAorder.working_directory + "Default_TradeStudy_Plot" + str(plot_ctr) + "Options.py")
            
            
            # Clear the plot
            plt.clf()
        
            # Call the plotter
            self.SingleDataSetPlotter(PEATSAorder,PO,PEATSAboxes)
            
            # Update the counter
            plot_ctr += 1
                
    def MissedThrustPlotter(self,PEATSAorder,PEATSAboxes):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import datetime
        import matplotlib.dates as mdates       
        import PEATSAmenu
        
        # Check if this is the first iteration or not
        if PEATSAorder.iteration == 0:
            # It is. 
        
            # Create a plot options object
            PO = PEATSAmenu.PEATSAinstagram()
            
            # Set the xlabel for the plot
            PO.xlabel = "Missed Thrust Event Date"
            # Set the ylabel for the plot
            PO.ylabel = "Longest Survivable Coast, days"
            # Set the title
            PO.title = "Missed Thrust Analysis"
            # Set the ylimits
            for box in PEATSAboxes:
                if "FAILURE" not in box.PEATSAdough_path and box.objective_value != -1e100:
                    PO.ylim = [-5.0,box.PEATSAdough.Journeys[0].wait_time_bounds[1]+5.0]
                    break
            # Set the file name string
            PO.file_name_str = "default_missed_thrust_plot"
            # Tell the plotter that the x axis is dates
            PO.x_is_date = 1
            # Tell the plotter that the y axis is not dates
            PO.y_is_date = 0
            # Tell the plotter to use lines not dots
            PO.lines_or_dots = 'dots'
            # Set the x_variable
            PO.x_variable = "box.seed_criteria[PEATSAorder.seed_criteria[0][0]]"
            # Set the y_variable
            PO.y_variable = "box.objective_value"
            
            PO.write_options(PEATSAorder.working_directory + "Default_MissedThrust_PlotOptions.py")
        else:
            # It is not the first iteration, so this plot should have midbake options saved to file. Load from there in case 
            # the user modified them
            
            # Create a plot options object from the midbake plot options file
            PO = PEATSAmenu.PEATSAinstagram(PEATSAorder.working_directory + "Default_MissedThrust_PlotOptions.py")
        
        # Clear the plot
        plt.clf()
        
        # Call the plotter
        self.SingleDataSetPlotter(PEATSAorder,PO,PEATSAboxes)
        
        if PEATSAorder.iteration > 0:
            self.MakeIterationHistory(PEATSAorder,"default_missed_thrust_plot")

    def BuiltInPlotterFromFile(self,PEATSAorder,PEATSAboxes,filename):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import PEATSAmenu
        
        # Load the options from file
        PO = PEATSAmenu.PEATSAinstagram(filename)
        
        # Clear the current figure
        plt.clf()
        
        # Check if the user wants multiple groupings on one plot
        if len(PO.groupings):
            # The user does.
            
            # Loop through all the boxes
            for box in PEATSAboxes:
                
                # Re-initialize its plot fingerprint
                box.plot_fingerprint = []
                
                # Loop through the different entries in the groupings list
                for fingerprint in PO.groupings:
                    
                    # Eval this fingerprint value and add it to the net fingerprint
                    box.plot_fingerprint.append(eval(fingerprint[0]))
                    
                box.plot_fingerprint.append(eval(PO.x_variable))
                    
            # Sort the plot boxes by the newly created fingerprint        
            plot_boxes = sorted(PEATSAboxes, key = lambda box: box.plot_fingerprint)
            
            # Create a list of everytime we reach a new fingerprint
            switches = [0]
            
            # Create a flag for the last good switch
            last_switch = 0
            
            # Create a list of whether to include this fingerprint in the plot
            includes = []
            
            # Get the initial fingerprint
            previous_fingerprint = plot_boxes[0].plot_fingerprint[0:-1]
            
            # Initialize the legend entries
            legend_entries = []
        
            # Loop through all of the boxes
            for box_idx in range(1,len(plot_boxes)):
            
                # Check if we are on a new plot fingerprint
                if plot_boxes[box_idx].plot_fingerprint[0:-1] != previous_fingerprint:
                    # We are.
                    
                    # check if we need to skip this fingerprint
                    if len(PO.grouping_requirements):
                        # Assume that this grouping will work
                        include = False
                        
                        # Loop through all requirements
                        for req in PO.grouping_requirements:
                            
                            # Assume this requirement will be met
                            local_include = True
                            
                            # Loop through each entry
                            for req_idx in range(len(previous_fingerprint)):
                                
                                # Check if we care about this entry                                
                                if req[req_idx] != "NA":
                                    # We care
                                    
                                    # Check if this requirement matches
                                    if req[req_idx] != previous_fingerprint[req_idx]:
                                        # It does not
                                        local_include = False
                                        # No need to check the rest of the requirements
                                        break
                                        
                            # Check if we made it through all of the fingerprint matching
                            if local_include == True:
                                # We did! Huzzah. Include this in the plot
                                include = True
                                # No need to check the rest of the requirements
                                break
                        # Add the flag to include this fingerprint or not
                        includes.append(include)
                        
                        # Update the last switch flag
                        if include:
                            last_switch = len(includes) - 1
                    else:
                        # Update the last switch flag
                        last_switch += 1
                        # Add the flag to include this fingerprint or not
                        includes.append(True)
                
                    # Add this to the switch
                    switches.append(box_idx)
                
                    # Update the fingerprint
                    previous_fingerprint = plot_boxes[box_idx].plot_fingerprint[0:-1]
                
                    # Get the legend entry for this fingerprint
                    box = plot_boxes[box_idx - 1]
                    legend_str = eval(PO.groupings[0][1])
                    # Loop through the entries
                    for fingerprint in PO.groupings[1:]:
                        legend_str += ", " + eval(fingerprint[1])
                        
                    if includes[-1]:
                        # Add this legend string to the list
                        legend_entries.append(legend_str)
             
            # Add the final index       
            switches.append(len(plot_boxes))
            
            # Get the legend entry for the final fingerprint
            box = plot_boxes[-1]
            
            legend_str = eval(PO.groupings[0][1])
            # Loop through the entries
            for fingerprint in PO.groupings[1:]:
                legend_str += ", " + eval(fingerprint[1])
            
            # check if we need to skip this fingerprint
            if len(PO.grouping_requirements):
                # Assume that this grouping will work
                include = False
                
                # Loop through all requirements
                for req in PO.grouping_requirements:
                    
                    # Assume this requirement will be met
                    local_include = True
                    
                    # Loop through each entry
                    for req_idx in range(len(plot_boxes[-1].plot_fingerprint[0:-1])):
                        
                        # Check if we care about this entry                                
                        if req[req_idx] != "NA":
                            # We care
                            
                            # Check if this requirement matches
                            if req[req_idx] != plot_boxes[-1].plot_fingerprint[req_idx]:
                                # It does not
                                local_include = False
                                # No need to check the rest of the requirements
                                break
                                
                    # Check if we made it through all of the fingerprint matching
                    if local_include == True:
                        # We did! Huzzah. Include this in the plot
                        include = True
                        # No need to check the rest of the requirements
                        break
                # Add the flag to include this fingerprint or not
                includes.append(include)
                        
                # Update the last switch flag
                if include:
                    last_switch = len(includes) - 1
            else:
                # Update the last switch flag
                last_switch += 1
                # Add the flag to include this fingerprint or not
                includes.append(True)
                
            if includes[-1]:    
                # Add this legend string to the list
                legend_entries.append(legend_str)
                
            # Loop through all of the switches
            for switch_idx in range(0,len(switches)-1):
    
                # Check if we are including this data
                if includes[switch_idx]:
                
                    # Call the plotting routine either with the legend entries or not (if it is the final set of data)
                    if switch_idx == last_switch and PO.legend == True:
                        self.SingleDataSetPlotter(PEATSAorder,PO,plot_boxes[switches[switch_idx]:switches[switch_idx+1]],legend_entries)
                    else:
                        self.SingleDataSetPlotter(PEATSAorder,PO,plot_boxes[switches[switch_idx]:switches[switch_idx+1]])
                                        
        else:
        
            # Call the plotting routine
            self.SingleDataSetPlotter(PEATSAorder,PO,PEATSAboxes)
        
    def SingleDataSetPlotter(self,PEATSAorder,PO,PEATSAboxes,label_list = []):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import datetime
        import matplotlib.dates as mdates       
        import itertools
        import PEATSAbusboy
        
        # Create a busboy for sorting and filtering
        busboy = PEATSAbusboy.PEATSAbusboy()

        # Get the y data
        if PO.y_is_date:
            # Check if the date is MJD or JD by finding the first converged case and seeing
            # if it is greater than 2400000
            for box in PEATSAboxes:
                # Check if this case converged
                if abs(box.objective_value) != 1e100 and "FAILURE" not in box.mission_name:
                    # It did
                    
                    # Evaluate the y variable and see if it is MJD or JD range
                    if eval(PO.y_variable) > 2400000.5:
                        # Variable is already JD
                        MJD_to_JD = 0
                    else:
                        # Variable needs to be converted
                        MJD_to_JD = 2400000.5
                        
                    # No need to keep looping
                    break

            # Loop through all cases and get the data for the cases that converged
            ydata = [datetime.date(*self.jd2datelist(eval(PO.y_variable) + MJD_to_JD)) for box in PEATSAboxes if abs(box.objective_value) != 1e100 and "FAILURE" not in box.mission_name]    
        else:
            # Loop through all cases and get the data for the cases that converged
            ydata = [eval(PO.y_variable) for box in PEATSAboxes if abs(box.objective_value) != 1e100 and "FAILURE" not in box.mission_name]   
       
        # Get the x data
        if PO.x_is_date:

            # Check if the date is MJD or JD by finding the first converged case and seeing
            # if it is greater than 2400000
            for box in PEATSAboxes:
                # Check if this case converged
                if abs(box.objective_value) != 1e100 and "FAILURE" not in box.mission_name:
                    # It did
                    
                    # Evaluate the y variable and see if it is MJD or JD range
                    if eval(PO.x_variable) > 2400000.5:
                        # Variable is already JD
                        MJD_to_JD = 0
                    else:
                        # Variable needs to be converted
                        MJD_to_JD = 2400000.5
                        
                    # No need to keep looping
                    break
            # Loop through all cases and get the data for the cases that converged        
            xdata = [datetime.date(*self.jd2datelist(eval(PO.x_variable) + 2400000.5)) for box in PEATSAboxes if abs(box.objective_value) != 1e100 and "FAILURE" not in box.mission_name]    
        else:
            # Loop through all cases and get the data for the cases that converged
            xdata = [eval(PO.x_variable) for box in PEATSAboxes if abs(box.objective_value) != 1e100 and "FAILURE" not in box.mission_name]   
        
        if PO.lines_or_dots == "dots":
            # Create a scatter plot
            plt.scatter(xdata,ydata)
        elif PO.lines_or_dots == 'lines':
            # Create a line plot
            plt.plot(xdata,ydata)
            
        # Set the x label
        if PO.xlabel != "":
            plt.xlabel(PO.xlabel)
    
        # Set the y label
        if PO.ylabel != "":
            plt.ylabel(PO.ylabel)
    
        # Set the title of the plot
        if PO.title != "":
            plt.title(PO.title)
    
        # Set the ylimits of the plot
        if PO.ylim[1] != PO.ylim[0]:
            plt.ylim(PO.ylim)
    
        # Set the xlimits of the plot
        if PO.xlim[1] != PO.xlim[0]:
            plt.xlim(PO.xlim)
    
        # Check if the x-axis is dates
        if PO.x_is_date:
            # It is.
            
            import math
            
            if len(xdata) > 1:
                
                # Get the range of dates
                x_range = max(xdata) - min(xdata)            
                        
                # In order to set the labels correctly, dtermine the date range
                if x_range.days < 180:
                    # Fewer than 6 months, write out the day and month
                    fmt = mdates.DateFormatter('%b %d')
                
                    # Set the interval to be 1/6th of the range in days
                    days = mdates.DayLocator(range(1,366),interval=int(math.ceil(x_range.days/6.0)))
            
                    # Get handles to the axes and figure
                    ax = plt.gca()
                    fig = plt.gcf()
                
                    # Set the tick marks to be at the 6 locations
                    ax.xaxis.set_major_locator(days)
                    # Apply the format
                    ax.xaxis.set_major_formatter(fmt)
                elif x_range.days < 365*6:
                    # More than 6 months, but less than 6 years
                    # Write out the month and year
                    fmt = mdates.DateFormatter('%b %Y')
                
                
                    # Set the interval to be 1/6th of the range in months
                    months = mdates.MonthLocator(range(1,13),bymonthday=1,interval=int(math.ceil(x_range.days/30.0/6.0)))
        
                    # Get handles to the axes and figure
                    ax = plt.gca()
                    fig = plt.gcf()
            
                    # Set the tick marks to be at the 6 locations
                    ax.xaxis.set_major_locator(months)
                    # Apply the format
                    ax.xaxis.set_major_formatter(fmt)
                else:
                    # More than 6 years. just put out the year
                    fmt = mdates.DateFormatter('%Y')
                
                    # Set the interval to be 1/6th of the range in years
                    years = mdates.YearLocator(int(math.ceil(x_range.days/365.25/6.0)))
        
                    # Get handles to the axes and figure
                    ax = plt.gca()
                    fig = plt.gcf()
            
                    # Set the tick marks to be at the 6 locations
                    ax.xaxis.set_major_locator(years)
                    # Apply the format
                    ax.xaxis.set_major_formatter(fmt)
                
            # Check if the y-axis is dates
            if PO.y_is_date:
                # It is.
                
                import math
            
                if len(ydata) > 1:
                    
                    # Get the range of dates
                    y_range = max(ydata) - min(ydata)
            
                    # In order to set the labels correctly, dtermine the date range
                    if y_range.days < 180:
                        # Fewer than 6 months, write out the day and month
                        fmt = mdates.DateFormatter('%b %d')
                
                        # Set the interval to be 1/6th of the range in days
                        days = mdates.DayLocator(range(1,366),interval=int(math.ceil(y_range.days/6.0)))
            
                        # Get handles to the axes and figure
                        ax = plt.gca()
                        fig = plt.gcf()
                
                        # Set the tick marks to be at the 6 locations
                        ax.yaxis.set_major_locator(days)
                        # Apply the format
                        ax.yaxis.set_major_formatter(fmt)
                    elif y_range.days < 365*6:
                        # More than 6 months, but less than 6 years
                        # Write out the month and year
                        fmt = mdates.DateFormatter('%b %Y')
                
                        # Set the interval to be 1/6th of the range in months
                        months = mdates.MonthLocator(range(1,13),bymonthday=1,interval=int(math.ceil(y_range.days/30.0/6.0)))
        
                        # Get handles to the axes and figure
                        ax = plt.gca()
                        fig = plt.gcf()
            
                        # Set the tick marks to be at the 6 locations
                        ax.yaxis.set_major_locator(months)
                        # Apply the format
                        ax.yaxis.set_major_formatter(fmt)
                    else:
                        # More than 6 years. just put out the year
                        fmt = mdates.DateFormatter('%Y')
                
                        # Set the interval to be 1/6th of the range in years
                        years = mdates.YearLocator(int(math.ceil(y_range.days/365.25/6.0)))
        
                        # Get handles to the axes and figure
                        ax = plt.gca()
                        fig = plt.gcf()
            
                        # Set the tick marks to be at the 6 locations
                        ax.yaxis.set_major_locator(years)
                        # Apply the format
                        ax.yaxis.set_major_formatter(fmt)
        
        # Apply any extra plot commands
        for command in PO.plot_commands:
            exec(command)
        
        # Check if we need to add a legend
        if label_list != []:
            # we do
            lgd = plt.gca().legend(label_list, loc = 'center left', bbox_to_anchor = (1.0, 0.5))
        
            # Save the image
            plt.savefig(PEATSAorder.images_dir + "Iteration" + str(PEATSAorder.iteration) + "_" + PO.file_name_str + ".png",bbox_etra_artists=(lgd,),bbox_inches='tight')
        else:
        
            # Save the image
            plt.savefig(PEATSAorder.images_dir + "Iteration" + str(PEATSAorder.iteration) + "_" + PO.file_name_str + ".png")
    
        # Check if an iteration history gif should be made
        if PEATSAorder.iteration != "Impatient" and PEATSAorder.iteration > 0:
            # It should. Call the method
            self.MakeIterationHistory(PEATSAorder,PO.file_name_str)

    # Method to make a iteration history gif
    def MakeIterationHistory(self,PEATSAorder,title_string):
        import logging
        import os
        import matplotlib
        from matplotlib.font_manager import findfont, FontProperties
        
        # Make sure imageio and PIL are actually installed
        try:
            import imageio
            from PIL import Image, ImageFont, ImageDraw
            import numpy
            if int(numpy.__version__.split(".")[1]) < 9:
                raise Exception("Version of numpy is old and incompatible with imageio")
        except:
            # They arent, cant do this
            logging.info("imageio and/or PIL not installed, so iteration history will not be generated")
            return
            
        
        # If we are here, then we have imageio and PIL and can make a gif!
        
        # First load the raw images
        raw_images = [file for file in os.listdir(PEATSAorder.images_dir) if file.endswith(title_string + ".png")]

        # Sort the images based on iteration
        raw_images = sorted(raw_images,key=lambda image: int(image.lstrip("Iteration").split("_")[0]))

        # Create a list to store the imageio objects
        images = []
        
        # Loop through all of the images to draw their iteration number on the plot
        for raw_image in raw_images:
            
            # First get the iteration number
            Iteration_no = raw_image.lstrip("Iteration").split("_")[0]
            
            # Second, open the raw image
            img = Image.open(PEATSAorder.images_dir + raw_image).convert('RGBA')
            
            # Create the draw object on the image
            draw = ImageDraw.Draw(img,'RGBA')
            
            try:
                # Get the default matplotlib font name
                font_name = findfont(FontProperties(family=[matplotlib.rcParams['font.family'][0]]))
            
                # Setup the matplotlib default font
                fnt = ImageFont.truetype(font_name,20)
            except:
                logging.info("Cant find matplotlib default font, using default and probably small font for gif iteration tag")
                
                # Set up the default font
                fnt = ImageFont.load_default()
            
            # Draw the iteration number
            draw.text((30,10),"Iteration " + Iteration_no,"Black",font=fnt)
            
            # Save the revised image
            img.save(PEATSAorder.images_dir + raw_image.rstrip("png").rstrip(".") + "_with_text.png")
            
            # Add the revised image to the list for imageio to work with
            images.append(imageio.imread(PEATSAorder.images_dir + raw_image.rstrip("png").rstrip(".") + "_with_text.png"))

        # Set the duration of each frame
        kargs = {"duration":1}
        
        # Create the gif
        imageio.mimsave(PEATSAorder.images_dir + "Iteration_history_" + title_string + ".gif",images,"GIF",**kargs)

        # Clean up the working directory
        raw_images = [file for file in os.listdir(PEATSAorder.images_dir) if file.endswith("_with_text.png")]
        
        # Loop through and delete all of the revised images
        for img in raw_images:
            os.remove(PEATSAorder.images_dir + img)
                             
    def Custom(self,PEATSAorder,PEATSAboxes):
        import sys
        sys.path.append(PEATSAorder.emtg_root_directory + "/PyEMTG")
        import MissionOptions
        import logging
        import importlib
                
        logging.info("Calling custom PEATSA post processing script")
        
        # Import the custom script
        custom_module = importlib.import_module(PEATSAorder.custom_post_processing_script)
        
        # Calling the custom script
        custom_module.PEATSAPostProcess(PEATSAorder,PEATSAboxes)
        
        logging.info("Done running custom post processing script.")

        # # Create an iter for the colors, based on how many lines are to be drawn.
        # # These colors were selected from here: http://colorbrewer2.org/
        # if PEATSAorder.plot_nx[idx] == 2:
        #     colorIter = itertools.cycle(["#d8b365","#5ab4ac"])
        # elif PEATSAorder.plot_nx[idx] == 3:
        #     colorIter = itertools.cycle(["#d8b365","#000000","#5ab4ac"])
        # elif PEATSAorder.plot_nx[idx] == 4:
        #     colorIter = itertools.cycle(["#a6611a","#dfc27d","#80cdc1","#018571"])
        # elif PEATSAorder.plot_nx[idx] == 5:
        #     colorIter = itertools.cycle(["#a6611a","#dfc27d","#000000","#80cdc1","#018571"])
        # elif PEATSAorder.plot_nx[idx] == 6:
        #     colorIter = itertools.cycle(["#8c510a","#d8b365","#f6e8c3","#c7eae5","#5ab4ac","#01665e"])
        # elif PEATSAorder.plot_nx[idx] == 7:
        #     colorIter = itertools.cycle(["#8c510a","#d8b365","#f6e8c3","#000000","#c7eae5","#5ab4ac","#01665e"])
        # elif PEATSAorder.plot_nx[idx] == 8:
        #     colorIter = itertools.cycle(["#8c510a","#bf812d","#dfc27d","#f6e8c3","#c7eae5","#80cdc1","#35978f","#01665e"])
        # elif PEATSAorder.plot_nx[idx] == 9:
        #     colorIter = itertools.cycle(["#8c510a","#bf812d","#dfc27d","#f6e8c3","#000000","#c7eae5","#80cdc1","#35978f","#01665e"])
        # elif PEATSAorder.plot_nx[idx] == 10:
        #     colorIter = itertools.cycle(["#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"])
        # else:
        #     colorIter = itertools.cycle(["#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#000000","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"])
            
            # # Create a locator for years, months and every three months
            # years = mdates.YearLocator()   # every year
            # months = mdates.MonthLocator()  # every month
            # threemonths = mdates.MonthLocator(range(1,13),bymonthday=1,interval=3)
            # monthsFmt = mdates.DateFormatter("%b '%y")
       
            # # Sort and filter the results based on this plot's fingerprint

    def CopyBestPeatsas(self,working_directory,destination_dir,iteration=0):
        import shutil
        import os
        
        if iteration == -1:
            csvfiles = [file for file in os.listdir(working_directory + "/docs") if file.endswith(".csv") and file.startswith("Iteration")]
            for csvfile in csvfiles:
                itno = int(csvfile.split("Iteration")[1].rstrip("csv").rstrip("."))
                if itno > iteration:
                    iteration = itno
        
        # Open the latest iteration file
        csv_handle = open(working_directory + "/docs/Iteration" + str(iteration) + ".csv",'r')

        # Loop through all lines in the file
        for line in csv_handle.readlines()[1:]:
    
            # Split the line by commas
            linesplit = line.split(",")
    
            # Make sure this isnt a fake failure file
            if "FAILURE" not in linesplit[1]:
    
                # Copy the file to the destination
                shutil.copyfile(linesplit[0] + "/" + linesplit[1],destination_dir + "/" + linesplit[1])
                # Copy the file to the destination
                shutil.copyfile(linesplit[0] + "/" + linesplit[1] + "opt",destination_dir + "/" + linesplit[1] + "opt")
            else:
                shutil.copyfile(linesplit[0] + "/" + linesplit[2] + ".emtgopt",destination_dir + "/" + linesplit[2] + ".emtgopt")

    def moveFiles(self,PEATSAorder):
        import os
        
        results_base = PEATSAorder.move_results + "/latest_" + PEATSAorder.PEATSAorder.run_name + "_results"
        
        # Check if the destination already has a csv file
        if os.path.isfile(results_base + ".csv"):
            # It does, remove it
            os.remove(results_base + ".csv")
        # Copy the latest result file to the destination
        shutil.copyfile(PEATSAorder.working_directory + "/docs/Iteration" + str(PEATSAorder.iteration) + ".csv",results_base + ".csv")
        
        # Check if the defination has a history csv file
        if os.path.isfile(results_base + "_history.csv"):
            # It does, remove it
            os.remove(results_base + "_history.csv")
        # Copy the latest result history file to the destination
        shutil.copyfile(PEATSAorder.working_directory + "/docs/PEATSAhistory.csv",results_base + "_history.csv")
        
        # Handle gif files
        giffiles = [file for file in os.listdir(PEATSAorder.working_directory + "/images") if file.endswith(".gif")]
        if len(giffiles):
            for idx,giffile in enumerate(giffiles):
                # Check if the destination already has this gif file
                if os.path.isfile(results_base + "_img" + str(idx) + ".gif"):
                    # It does, remove it
                    os.remove(results_base + "_img" + str(idx) + ".gif")
                shutil.copyfile(PEATSAorder.working_directory + "/images/" + giffile,results_base + "_img" + str(idx) + ".gif")
      
        # Handle images
        img_files = [file for file in os.listdir(PEATSAorder.working_directory + "/images") if file.endswith(".png") and file.startswith("Iteration" + str(PEATSAorder.iteration))]
        for idx,img in enumerate(img_files):
            if os.path.isfile(results_base + "_img" + str(idx) + ".png"):
                os.remove(results_base + "_img" + str(idx) + ".png")
            shutil.copyfile(PEATSAorder.working_directory + "/images/" + img,results_base + "_img" + str(idx) + ".png")

        # Get the most recent results
        if os.path.isfile(results_base + "_bestResults.tar.gz"):
            os.remove(results_base + "_bestResults.tar.gz")
        os.makedirs(results_base + "_bestResults")
        self.CopyBestPeatsas(PEATSAorder.working_directory,results_base + "_bestResults",PEATSAorder.iteration)
        os.system("tar czf " + results_base + "_bestResults.tar.gz " + results_base + "_bestResults")
        shutil.rmtree(results_base + "_bestResults")
        
        # Write out a file to update when the iteration results were written
        handle = open(results_base + ".txt",'w')
        handle.write("Iteration " + str(latest_iteration) + " written at " + datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
        handle.close()