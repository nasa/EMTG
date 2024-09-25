import sys

if len(sys.argv) != 3:
    print("Syntax:")
    print("python impatiently_run_post_process_now.py csv_and_plot_destination_folder /full/path/to/PEATSA_run")
else:
    import PEATSAmenu
    import PEATSAbusboy
    import PEATSAdelivery
    import os
    import logging
    
    # Set up a log file script
    try:
        logging.basicConfig(filename="impatient.log",filemode='w',level=logging.INFO,format='%(asctime)s %(message)s')
        logging.info("Successfully opened logfile")
        print("Successfully opened logfile, all status updates will now be logged instead of printed to console.")
    except:
        raise Exception("Error opening logfile")

    # Define a custom error handler so that the errors go to the logfile, not the console
    def myError(excType, excValue, traceback):
        logging.error("*************** Oh Noes! A PEATSA Error! ****************",
                     exc_info=(excType, excValue, traceback))
    sys.excepthook = myError
    
    dir_path = os.path.dirname(os.path.realpath(__file__))
    sys.path.append(dir_path + "/..")

    opts = PEATSAmenu.PEATSAmenu()

    opts.load_options_from_file(sys.argv[2] + "/PEATSA_ExecutedOptions.py")
    
    original_plot_files = opts.built_in_plotter_files
    
    opts.load_options_from_file(sys.argv[2] + "/PEATSA_MidBake_Options.py")
    
    opts.built_in_plotter_files = original_plot_files
    
    busboy = PEATSAbusboy.PEATSAbusboy()
    
    delivery = PEATSAdelivery.PEATSAdelivery()
    
    PEATSAorder.restart_run_root_directory = [sys.argv[2] + "/results"]
    PEATSAorder.csv_file = sys.argv[1] + "/Impatient_results.csv"
    PEATSAorder.iteration = "Impatient"
    PEATSAorder.images_dir = sys.argv[1]
    
    # Call the waiter to parse the results
    boxes = busboy.fresh_parse_results(opts)

    # Sort the results and eliminate lesser cases
    boxes = busboy.sort_and_filter(opts,boxes)
    
    # Write the 0th iteration results to file
    busboy.write_to_csv(opts,boxes)
    
    # Call the delivery
    delivery.post_process(opts,boxes)
    
        
    