import shutil
import sys
import os
import PEATSAdelivery

if __name__ == "__main__":
    # Ensure the correct number of command line options were provided
    if len(sys.argv) != 3:
        raise Exception("Unknown number of command line options!\n\
                         Syntax: python grab_best_peatsa_results.py peatsa_folder destination_folder")
    
    # expand username if a relative path was used
    source_dir = os.path.expanduser(sys.argv[1])
    destination_dir = os.path.expanduser(sys.argv[2])
                     
    delivery = PEATSAdelivery.PEATSAdelivery()
    
    delivery.CopyBestPeatsas(source_dir,destination_dir,-1)