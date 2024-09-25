import shutil
import sys
import os
import glob

if __name__ == "__main__":
    # Ensure the correct number of command line options were provided
    if len(sys.argv) != 3:
        raise Exception("Unknown number of command line options!\n\
                         Syntax: python get_most_recent_doc.py peatsa_doc_folder destination_folder")
                         
    source_dir = os.path.expanduser(sys.argv[1])
    destination_dir = os.path.expanduser(sys.argv[2])
    
    list_of_files = glob.glob(source_dir + '/Iteration*') # * means all if need specific format then *.csv
    latest_iteration_file = max(list_of_files, key=os.path.getctime)
                         
    shutil.copyfile(latest_iteration_file,destination_dir + '/' + latest_iteration_file.split('/')[-1])
    shutil.copyfile(source_dir + '/PEATSAhistory.csv',destination_dir + '/PEATSAhistory.csv')