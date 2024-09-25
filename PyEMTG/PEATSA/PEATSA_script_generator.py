if __name__ == "__main__":

    # Define what options must be set, and acceptable responses
    needed_options = ["start_type","PEATSA_type","objective_type","generate_default_plots"]
    acceptable_options = {"PEATSA_type"            : [0,1,2,3,4,5],
                          "start_type"             : ["Warm","Fresh","Hot"],
                          "objective_type"         : [0,1],
                          "generate_default_plots" : [0,1],
                          "initial_guess_type"     : ["x_days","n_steps","copy"]}
    
    # -----------------------------------------------------------------
    # Shouldn't need to change anything below here
    import sys
    import PEATSAmenu
    
    # Create a peatsa menu with default options
    opts = PEATSAmenu.PEATSAmenu()
    
    # Check what version of python is running
    if (sys.version_info > (3, 0)):
        # Python 3
        input_function = input
    else:
        # Python 2
        input_function = raw_input
                  
    # Loop through the options
    for var in needed_options:
        
        # Make sure this option is still valid
        if not (eval(opts.requirements[var].replace("self","opts"))):
            continue
            
        # Print the description of the start_types     
        print(opts.descriptions[var])   
        # Loop until we get a valid answer 
        while True:
            # Ask the user to enter the start type
            user_input = input_function("Enter " + var + " = ")
            
            # Check if the start type is valid
            # First, Check the type of variable and then apply it
            if isinstance(eval("opts." + var),str):
                if user_input in acceptable_options[var]:
                    # It was, break the while loop
                    break
                else:
                    # It was not valid, try again.
                    print("Invalid input. Try again.")
            else:
                try: 
                    user_input_val = int(user_input)
                except:
                    print("Invalid input. Try again.")
                    continue
                    
                if user_input_val in acceptable_options[var]:
                    # It was, break the while loop
                    break
                else:
                    # It was not valid, try again.
                    print("Invalid input. Try again.")
                
                    
               
        # Check the type of variable and then apply it
        if isinstance(eval("opts." + var),str):
            # Put the user input into the options object
            exec("opts." + var + " = '" + user_input + "'")
        else:
            # Put the user input into the options object
            exec("opts." + var + " = " + user_input)
    
    user_input = input_function("Enter filename: ")
    
    if user_input != "":
        filename = user_input
    else:
        filename = "PEATSAsample.py"
        
    if not filename.endswith(".py"):
        filename += ".py"    
        
    if opts.start_type == "Fresh" and opts.PEATSA_type == 1:
        opts.objective_formula = 'M.Journeys[0].missionevents[0].JulianDate - 2400000.5 - MO.launch_window_open_date'
        opts.peatsa_goal_objective = 59.0
        opts.seed_criteria = [("MO.launch_window_open_date",1,-1,2)]
    
    if opts.PEATSA_type == 4:
        opts.seed_criteria = [(1,0,2)]
        opts.seed_from_cases_that_havent_met_target = 1
        
    if opts.PEATSA_type ==5 and opts.start_type == "Fresh":
        opts.seed_criteria = [("MO.user_data['MonteCarloSampleIndex']",1,-1,0)]
        opts.fingerprint = ["MO.user_data['MonteCarloDepthIndex']"]
    
    # Output the options file
    opts.conditionally_write_to_file(filename)
    
    
    


