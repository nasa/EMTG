import psutil
import sys
import os
import time

if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise Exception('Syntax: python kill_iteration.py /path/to/peatsa/folder')
    elif len(sys.argv) > 2:
        raise Exception('Too many options specified.\nSyntax: python kill_iteration.py /path/to/peatsa/folder')
     
    # expand username if a relative path was used
    source_dir = os.path.expanduser(sys.argv[1])
                     
    # Grab all summary files
    csv_summary_files = [csvfile for csvfile in os.listdir(source_dir + "/docs/") if csvfile.endswith(".csv")]
    
    # Assume the latest iteration is negative 1 to start
    latest_iteration = -1
    
    # Loop through all csv files to get their iteration number
    for csvfile in csv_summary_files:
        
        # Get the iteration number
        it_number = int(csvfile.lstrip("Iteration").rstrip(".csv\r\n "))
        
        # Check if this is later than the last iteration number found
        if it_number > latest_iteration:
            # It is. store it
            latest_iteration = it_number
    
    # Check if the latest iteration number was updated    
    if latest_iteration == -1:
        # It was not.
        raise Exception("No csv files were found in that peatsa folder.\n\
                         Syntax: python kill_iteration.py /path/to/peatsa/folder")
        
    # Loop until the next iteration folder exists                     
    while "/Iteration" + str(latest_iteration+1) + ".csv" not in csv_summary_files:
        
        # Loop through the process list
        for proc in psutil.process_iter():
            
            # Check if this is emtg
            if proc.name() == 'EMTGv9' or proc.name() == 'emtg':
                
                # It does. Lets kill it
                proc.terminate()
        
        # Wait 5 seconds
        time.sleep(5)
        
        # Grab all summary files
        csv_summary_files = [csvfile for csvfile in os.listdir(source_dir + "/docs/") if csvfile.endswith(".csv")]
    