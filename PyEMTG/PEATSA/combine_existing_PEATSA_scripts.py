import PEATSAmenu
import sys

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Syntax:")
        print("python combine_existing_PEATSA_scripts.py new_script.py existing_script1.py ... existing_script#.py")
        print("or")
        print("python combine_existing_PEATSA_scripts.py new_script.py PEATSA_path (this will combine the executed options and the mid bake options for the peatsa run)")
    else:
        opts = PEATSAmenu.PEATSAmenu()
    
        if len(sys.argv) > 3:
            for script in sys.argv[2:]:
                opts.load_options_from_file(script)  
        
            opts.conditionally_write_to_file(sys.argv[1])
        else:
            opts.load_options_from_file(sys.argv[2] + "/PEATSA_ExecutedOptions.py")
            opts.load_options_from_file(sys.argv[2] + "/PEATSA_MidBake_Options.py")
            opts.conditionally_write_to_file(sys.argv[1])
