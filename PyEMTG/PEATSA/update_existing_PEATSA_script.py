"""
update_existing_PEATSA_script.py
=================================

This script is run to take an existing PEATSA script and update it to include add to the script any options
that were not included in the input script. Calling sequence is

python update_existing_PEATSA_script.py <path/to/existing/peatsa/script.py> [optional: if_mid_bake; default = False]

<path/to/existing/peatsa/script.py> is the script you wish to update. If if_mid_bake is included and True (or 1),
then the script is updated as a midbake script rather than a base script.

"""

import PEATSAmenu
import sys

def main():
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print("Syntax:")
        print("python update_existing_PEATSA_script.py existing_script.py [if_mid_bake, default=False]")
    else:
    
        opts = PEATSAmenu.PEATSAmenu(sys.argv[1])
    
        if len(sys.argv) == 2:
            if_mid_bake = False
        else:
            if_mid_bake = int(sys.argv[2])
    
        opts.write_to_file(sys.argv[1],True,if_mid_bake)
        
if __name__ == "__main__":
    main()