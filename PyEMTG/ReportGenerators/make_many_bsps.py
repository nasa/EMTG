#script to make .bsp files out of everything in a ``best'' folder

working_directory = ''

def create_ephemeris_file(options_file, destination_dir, EMTGpath):
    import MissionOptions
    import Mission
    
    #does a seed exist?
    seedfilename = os.path.join(working_directory, options_file.replace('.emtgopt','.emtg'))
    if not os.path.exists(seedfilename):
        print('Seed file for case ' + options_file + ' does not exist.')
        return
        
    #load
    myOptions = MissionOptions.MissionOptions(os.path.join(working_directory, options_file))
        
    #seed
    seedMission = Mission.Mission(seedfilename)

    myOptions.trialX = []
    for Xindex in range(0, len(seedMission.DecisionVector)):
        myOptions.trialX.append([seedMission.Xdescriptions[Xindex].replace('\n','').replace('\r',''), seedMission.DecisionVector[Xindex]])

    myOptions.DisassembleMasterDecisionVector()
    myOptions.ConvertDecisionVector()    
        
    #redirect output
    myOptions.forced_working_directory = destination_dir + '/'
    myOptions.override_mission_subfolder = 0
        
    #enable ephemeris printing
    myOptions.generate_forward_integrated_ephemeris = 1

    #set up SPICE writer
    myOptions.spice_utilities_path = '/Utilities/cspice/exe/'
    myOptions.spice_utility_extension ='' #Linux
    myOptions.pyemtg_path = '~/emtg/PyEMTG/'

    #trialX mode
    myOptions.run_inner_loop = 0
        
    #save
    myOptions.write_options_file(destination_dir + '/' + options_file)
        
    #run!
    os.system(EMTGpath + ' ' + destination_dir + '/' + options_file)

def create_bsp(bsp_to_make):
    os.system('python ' + bsp_to_make)

def make_bsps(args):
    import os
    import time
    nproc = 192

    #Windows needs this for some reason
    #if os.name == "nt":
    #    multiprocessing.freeze_support()

    #myPool = multiprocessing.Pool(multiprocessing.cpu_count())
    
    working_directory = args[0]

    if working_directory == '.':
        from os import getcwd
        working_directory = getcwd()

    destination_dir = working_directory + '/ephemeris'
    EMTGpath = args[1]
    
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)
    
    options_files = [f for f in os.listdir(working_directory) if '.emtgopt' in f]
    
    #myPool.map(create_ephemeris_file, options_files)
    counter = 1
    for options_file in options_files:
        os.system('python ~/emtg/PyEMTG/ReportGenerators/print_ephemeris.py ' + options_file + ' ' + working_directory + ' ' + destination_dir + ' ' + EMTGpath + ' &')
        counter += 1
        if counter > nproc:
            counter = 1
            print('resetting counter')
            time.sleep(30) #sleep for 30 seconds to catch up
            print('resuming counter')
        #create_ephemeris_file(options_file, destination_dir, EMTGpath)
        
    #now walk through and run every single bspwriter.py script
    bsps_to_make = []
    for root,dirs,files in os.walk(destination_dir):
        for file in files:
            if 'bspwriter.py' in file:
                bsps_to_make.append(os.path.join(root, file))
                
    counter = 1
    for bsp_to_make in bsps_to_make:
        sys.path.append(os.path.dirname(bsp_to_make))
        os.system('python ' + bsp_to_make + ' &')
        counter += 1
        if counter > nproc:
            counter = 1
            time.sleep(30) #sleep for 30 seconds to catch up
    #for bsp_to_make in bsps_to_make:
    #    create_bsp(bsp_to_make)
    #myPool.map(create_bsp, bsps_to_make)

if __name__ == '__main__':
    import sys
    import os

    if len(sys.argv) != 3:
        raise Exception("Unknown number of command line options!\n\
                         Syntax: python make_many_bsps.py working_directory path_to_EMTG_executable")
    
    thisdir = os.path.dirname(os.path.realpath(sys.argv[0]))
    sys.path.append(thisdir)
    sys.path.append(thisdir + '/../')
    sys.path.append(thisdir + '/../SpiceyPy_Utilities')
    
    make_bsps(sys.argv[1:])
