'''
This program allows the user to run any or all of the 
EMTG regression tests in the test suite. 

--------------------------------------- USAGE -------------------------------------------
1) Open anaconda prompt & change directory to the testatron folder 
   e.g. > cd C:\emtg\testatron
   
2) Use 'python testatron.py -h' to see help on command-line args
'''

'''
Import utilities
'''
from os import makedirs, getcwd, listdir, system, walk, path
from time import strftime
import ast
import pdb # Python debugger; use q in command line to quit
import sys
import argparse

# Can't use a relative path because we use importlib later
test_directory = getcwd().replace('\\','/') + '/tests/' 


'''
PARSE USER INPUT ------------------------------------------------------------------------
'''
test_cases = [];

# use argparse package
parser = argparse.ArgumentParser(description="EMTG test system driver. Select only ONE of -c, -f, --failure, -u, -m, -a, --update_truths. If none is provided, -a is the default")

# arguments

# path to emtg executable
parser.add_argument('-e', '--emtg', dest = 'emtgPath', 
    default = 'c:/emtg/bin/EMTGv9.exe',
    help = 'Set the path to the EMTG executable (default: c:/emtg/bin/EMTGv9.exe)')
    
# path to pyemtg
parser.add_argument('-p', '--pyemtg', dest = 'pyemtgPath',
    default = 'c:/emtg/PyEMTG/',
    help = 'Set the path to PyEMTG (default: c:/emtg/PyEMTG/)')
    
# run specific cases
parser.add_argument('-c', '--cases', dest = 'runCases', nargs = '*',
    help = 'Specify full path to one or more individual cases to run. Separate multiple cases with spaces only. DO NOT include file extension. DO NOT put in quotes. DO NOT put in brackets. DO NOT use commas.')

# run specific folders of cases
parser.add_argument('-f', '--folders', dest = 'runFolders', nargs = '*',
    help = 'Run all cases in one or more folders. Folders must be in testatron/tests. End folder names with a /. Separate multiple folders with spaces only. DO NOT prepend folder name with /path/to/testatron/tests. DO NOT put in quotes. DO NOT put in brackets. DO NOT use commas.')

# run failure cases from csv
parser.add_argument('--failure', dest = 'runFailure', nargs = 1,
    help = 'Run all tests that failed, as per a failed_tests.csv file. Give path to folder that contains failed_tests.csv, ending with a /. DO NOT include the "failed_tests.csv" file in the path to the file.')

# run unit tests (ie, all except failures)
parser.add_argument('-u', '--unit', dest = 'runUnit', nargs = '?',
    const = 1, default = 0, type = int,
    help = 'If active, run unit tests (currently, this means all non-failing tests)')
    
# run mission tests only
parser.add_argument('-m', '--mission', dest = 'runMission', nargs = '?',
    const = 1, default = 0, type = int,
    help = 'If active, run mission tests')
    
# run all cases
parser.add_argument('-a', '--all', dest = 'runAll', nargs = '?',
    const = 1, default = 0, type = int,
    help = 'If active, run all cases. This is also the default behavior if no other specific behavior is requested.')
    
# update truth .emtg files
parser.add_argument('--update_truths', dest = 'updateTruths', nargs = '?',
    const = 1, default = 0, type = int,
    help = 'If selected, the test system IS NOT RUN. Instead, all test .emtgopt files are executed, and all truth .emtg files are replaced with the results of executing the .emtgopt files. This option is useful if something has changed that breaks every test. For example, if a new attribute has been added to the PyEMTG Mission class.')

# attributes to ignore when running comparatron
parser.add_argument('--ignore', dest = 'attributes_to_ignore', nargs = '*',
	default = [],
    help = 'List attributes to ignore when running Comparatron. Start Mission attributes with M., Journey attributes with J., and MissionEvent attributes with E.. DO NOT put in quotes. DO NOT put in brackets. DO NOT use commas.')


args = parser.parse_args()

EMTG_path = args.emtgPath
PyEMTG_path = args.pyemtgPath

# we want run all to be the default behavior
runCases = 0
runFolders = 0
runFailure = 0
runUnit = 0
runMission = 0
updateTruths = 0
runAll = args.runAll
testCases = None
testFolders = None
failure = None
if (args.runCases == None and args.runFolders == None and args.runFailure == None and args.runUnit == 0 and args.runMission == 0 and args.updateTruths == 0):
    runAll = 1
    run_type = 'all'
    
if (args.runCases != None):
    runCases = 1
    testCases = args.runCases
    run_type = 'cases'
elif (args.runFolders != None):
    runFolders = 1
    testFolders = args.runFolders
    run_type = 'folders'
elif (args.runFailure != None):
    runFailure = 1
    failure = args.runFailure
    run_type = 'failed'
elif (args.runUnit == 1):
    runUnit = 1
    run_type = 'folders'
elif (args.runMission == 1):
    runMission = 1
    run_type = 'folders'
elif (args.updateTruths == 1):
    updateTruths = 1
    run_type = 'all'

if runCases == 1:
    print('Running user cases')
    for casse in testCases:
        test_cases.append(casse) 
        
elif runFolders == 1:
    print('Running all cases in user-specified folders')
    for folderr in testFolders:
        test_cases.append(folderr) 
        
elif runFailure == 1:
    print('Running failed tests')
    test_cases.append(failure[0]) # Needs to be a *.csv file from the testatron\output dir
    
elif runUnit == 1:
    print('Running unit (feature) tests')
    run_type = 'folders'
    test_cases = ['global_mission_options','journey_options','mission_tests', 'output_options', 'physics_options', 'script_constraint_tests', 'solver_options','spacecraft_options','state_representation_tests','transcription_tests']
    
elif runMission == 1:
    print('Running mission tests')
    run_type = 'folders'
    test_cases = ['mission_tests']
    
elif updateTruths == 1:
    print('Updating truth .emtg files')
    run_type = 'all'
    test_cases = []
    
else: # Else run all b/c can't use 'cases' without actual cases
    print('Running all tests')
    run_type   = 'all'
    test_cases = []




'''
METHOD: MAKE THE LIST OF TESTS ----------------------------------------------------------
'''
def MakeTestsList(test_cases):
    tests_list = [] # Initialize
    test_folders = []
    if run_type == 'all':
        test_folders = next(walk(test_directory))[1] 
        test_folders.append('')
    
    elif run_type == 'folders':
        # Uses only the user-provided subfolders
        test_folders = test_cases 
    
    elif run_type == 'cases':
        for casse in test_cases:
            tests_list.append(casse.replace('\\','/'))
        return tests_list
    
    elif run_type == 'failed':
        # Can give testatron multiple fail files
        for failpath in test_cases: 
            with open(failpath+'failed_tests.csv','r') as f:
                # Reads the text in the fail file
                lines = f.readlines()  
                # Converts the list inside the string in the last line of 
                # the file into a usable list
                listOfFailedFlag = False
                
                for line in lines:
                    if 'List of failed runs:' in line:
                        listOfFailedFlag = True
                    elif listOfFailedFlag == True:
                        tests_list.append(line.strip('\n'))
        # Removes tests that appear more than once (in case the same
        # test appears in multiple fail files)
        tests_list = list(set(tests_list)) 
        return tests_list      
    
    # Convert list of folders to list of tests in those folders
    for test in test_folders:
        # Creates a list of all files in each subfolder
        #pdb.set_trace()
        filenames = listdir(test_directory+test) 
        for file in filenames:
            if file.endswith('.emtgopt'):
                # Appends all .emtgopt files to the list of tests
                tests_list.append(test_directory+test+'/'+file.replace('.emtgopt','')) 

    return tests_list
'''
------------------------------------------------------------------------------------------
'''

# Create tests list
tests_to_run = MakeTestsList(test_cases)

print(tests_to_run)

# create directories for output if we are actually running the test system
if updateTruths == 0:
    # Get the current epoch
    now = strftime('%c')
    now_formatted = now.replace(' ','_').replace(':','')

    # Create an output directory and summary file
    outputdir = getcwd() + '/output/' + now_formatted
    print('OUTPUT_DIRECTORY: '+outputdir)
    makedirs(outputdir)

    # Overall summary file with all run tests & their ?success status
    summaryFile = open(outputdir + '/test_results.csv','w') 
    summaryFile.write('Beginning test run ' + now + '\n\n')

    # File that only records when a test fails, records the test 
    # and all mission objects that failed
    failFile = open(outputdir + '/failed_tests.csv','w') 
    failFile.write('Beginning test run ' + now + '\n\n')

'''
Run the tests ---------------------------------------------------------------------------

'''
import os
import sys
sys.path.append(test_directory)
sys.path.append(PyEMTG_path)
import importlib

import Mission
import MissionOptions

failedTests      = 0  # Continuous counter of the number of failed tests
failedTests_list = [] # List of file paths for failed tests

for test in tests_to_run:
    test_name=test.split('/')[-1]
    if updateTruths == 1:
        print('Updating test "' + test + '"...')
    else:
        summaryFile.write('Beginning test "' + test + '"...')

    if path.isfile(test+'.emtgopt'): #True: #try:
        testOptions = MissionOptions.MissionOptions(test + '.emtgopt') # Emtg ops object

        # Override the output, universe, and HardwareModels paths
        testOptions.override_working_directory = 1
        testOptions.short_output_file_names    = 1
        testOptions.background_mode            = 1
        testOptions.override_mission_subfolder = 1
        testOptions.forced_mission_subfolder   = '.'
        testOptions.universe_folder            = test_directory.replace('tests/',\
                                                                        'universe/')
        testOptions.HardwarePath               = test_directory.replace('tests/',\
                                                                'HardwareModels/')
                                                                
        if updateTruths == 1:
            # if updating truths, want the output to go to the same directory as the input
            testOptions.forced_working_directory   = os.path.dirname(test)
        else:
            # if running tests, want all output to go to the output directory
            testOptions.forced_working_directory   = outputdir
        
                                                                
        # update all gravity file paths to use tests/universe/gravity_files
        for journey in testOptions.Journeys:
            grav_file_name = journey.central_body_gravity_file.replace('\\\\','/').replace('\\','/').split('/')[-1]
            journey.central_body_gravity_file = test_directory.replace('tests/','universe/gravity_files/') + grav_file_name

        # Save updated emtg options file and always write out all of the options to the file
        if updateTruths == 1:
            # if updating truths, overwrite the .emtgopt files
            testOptions.write_options_file(test + '.emtgopt', not testOptions.print_only_non_default_options)
        else:
            # if running tests, save to output directory
            testOptions.write_options_file(outputdir + '/' + test_name + '.emtgopt', not testOptions.print_only_non_default_options)

        try:
            # Run EMTG
            if updateTruths == 1:
                system(EMTG_path + ' ' + test + '.emtgopt')
            else:
                system(EMTG_path + ' ' + outputdir + '/' + test_name + '.emtgopt')
        except:
            if updateTruths == 1:
                print('\nFAILURE to run "' + test + '" options file.\n\n')
            else:
                summaryFile.write('\nFAILURE to run "' + test + '" options file.\n\n')
                failFile.write('\nFAILURE to run "' + test + '" options file.\n\n')
                failedTests += 1
                failedTests_list.append(test)
            continue
        
        if updateTruths == 0:
            try:
                # Post-process
                testMission = Mission.Mission(outputdir + '/' + test_name + '.emtg')
            except:
                summaryFile.write('\nFAILURE to parse output for "' + test + '" options file.\n\n')
                failFile.write('\nFAILURE to parse output for "' + test + '" options file.\n\n')
                failedTests += 1
                failedTests_list.append(test)
                continue
                    
            print("\nRunning comparator...\n")
            # Run standard comparator that checks every mission event
            try:
                success, output = testMission.Comparatron(baseline_path = test + '.emtg',\
                          csv_file_name = outputdir +'/' + test_name + '_comparison.csv',\
                                                                       full_output=False,\
                                                                       tolerance_dict={},\
                                                              default_tolerance = 1.0e-10,\
                                                              attributes_to_ignore = args.attributes_to_ignore)                                                        

                if success:
                    summaryFile.write('successful\n\n')
                else:
                    summaryFile.write('failed\n\n')
                    failFile.write('\nTest "' + test + '" failed for parameters:\n')
                    failFile.close() 

                    # Writes output to the end of the failFile
                    output.loc[output['Match']==False].to_csv(outputdir + '/failed_tests.csv',\
                                                                         mode='a', index=False)

                    # Reopens the failFile so that the driver can write the next test
                    failFile = open(outputdir + '/failed_tests.csv','a') 
                    failFile.write('\n')
                    failedTests += 1
                    failedTests_list.append(test)
            except:
                summaryFile.write('\nFAILURE to compare "' + test + '" options file.\n\n')
                failFile.write('\nFAILURE to compare "' + test + '" options file.\n\n')
                failedTests += 1
                failedTests_list.append(test)
    else:
        if updateTruths == 0:
            summaryFile.write('\nFAILURE to load "' + test + '" options file.\n\n')
            failFile.write('\nFAILURE to load "' + test + '" options file.\n\n')
            failedTests += 1
            failedTests_list.append(test)

if updateTruths == 1:
    # clean up: don't keep bspwriter.py, XFfile.csv, mission_maneuver_spec, mission_target_spec, cmd, or emtg_spacecraftopt files
    test_folders = []
    test_folders = next(walk(test_directory))[1] 
    test_folders.append('')
    for folder in test_folders:
        filesToDelete = [] # list of files to delete
        filesInDirectory = listdir(test_directory + folder) # all files in this directory
        
        # gather up all the files to delete
        
        files = [file for file in filesInDirectory if file.endswith(".py")] # get rid of bspwriter.py
        filesToDelete.extend(files)
        
        files = [file for file in filesInDirectory if file.endswith("XFfile.csv")] # get rid of XFfile.csv
        filesToDelete.extend(files)
        
        files = [file for file in filesInDirectory if file.endswith(".mission_maneuver_spec")] # get rid of .mission_maneuver_spec
        filesToDelete.extend(files)
        
        files = [file for file in filesInDirectory if file.endswith(".mission_target_spec")] # get rid of .mission_target_spec
        filesToDelete.extend(files)
        
        files = [file for file in filesInDirectory if file.endswith(".cmd")] # get rid of .cmd
        filesToDelete.extend(files)
        
        files = [file for file in filesInDirectory if file.endswith(".emtg_spacecraftopt")] # get rid of .emtg_spacecraftopt
        filesToDelete.extend(files)
        
        # loop through the files we found
        for file in filesToDelete:
            pathToFile = os.path.join(test_directory + folder, file) # add the path to the file
            #print("Want to delete " + pathToFile)
            os.remove(pathToFile) # delete the file

        
    
    print("Finished updating truth files")
else:
    # Find and write the time that the tests finished running
    end = strftime('%c')
    summaryFile.write('Finished test run ' + end)
    failFile.write('Finished test run ' + end + '\n')

    # Add a list of failed cases to the end of the failFile. 
    # This can be used for future testatron runs
    failFile.write('\nList of failed runs:\n')
    for failTest in failedTests_list:
        failFile.write(failTest + '\n')

    print("\nAll tests completed.\n")
    print("Failed " + str(failedTests) + " test(s).")
    for failTest in failedTests_list:
        print(failTest)

    summaryFile.close()
    failFile.close()




