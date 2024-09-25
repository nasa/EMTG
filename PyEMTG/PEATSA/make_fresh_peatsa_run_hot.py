import PEATSAmenu
import sys

if len(sys.argv) != 3:
    print("Syntax:")
    print("python make_fresh_peatsa_run_hot.py new_script.py /full/path/to/PEATSA_run")
else:

    opts = PEATSAmenu.PEATSAmenu()
    opts.load_options_from_file(sys.argv[2] + "/PEATSA_ExecutedOptions.py")
    opts.load_options_from_file(sys.argv[2] + "/PEATSA_MidBake_Options.py")
    opts.start_type = "Hot"
    opts.restart_run_root_directory = [(sys.argv[2] + "/results/","Full")]
    opts.working_directory = opts.working_directory.rstrip("/").rstrip(opts.run_name + "_1234567890")
    opts.conditionally_write_to_file(sys.argv[1])