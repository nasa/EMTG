'''
This program allows the user to run any or all of the 
EMTG regression tests in the test suite. 

--------------------------------------- USAGE -------------------------------------------
1) Open anaconda prompt & change directory to the testatron folder 
   e.g. > cd C:\emtg\testatron
   
2) Optionally specify a run type & list of cases. Syntax:
   
   *run_type = all or <empty> (default; run the whole suite))
       > python testatron.py all 
       OR 
       > python testatron.py
       
   *run_type = cases (run all specific test files)
       > python testatron.py cases "relative\path\to\case\nofileextension"
       *EXAMPLE: > python testatron.py cases tests\spacecraft_options\spacecraftoptions_varinitmass
       OR
       > python testatron.py cases
       ^ Runs all the tests
       
   *run_type = folders (run all tests in given folder(s)
       > python testatron.py folders folder_name_1 ... folder_name_N
       *EXAMPLE: > python testatron.py folders journey_options output_options
       ^ DO NOT PREFIX YOUR FOLDER NAME WITH testatron\tests\ . This gets appended later.
       ^EDIT: This syntax assumes your test folder is in testatron\tests\
       
   *run_type = failed (run all tests that failed, as per a failed_tests.csv file)
       > python testatron.py failed "\path\to\failed_tests_CSV
       ^ Do not include the "failed_tests.csv" file in the path.       
       
3) Run...
'''

'''
Import utilities
'''
from os import makedirs, getcwd, listdir, system, walk, path
from time import strftime
import ast
import pdb # Python debugger; use q in command line to quit
import sys

# Can't use a relative path because we use importlib later
test_directory = getcwd().replace('\\','/') + '/tests/' 
EMTG_path      = 'c:/emtg/bin/EMTGv9.exe'
PyEMTG_path    = 'c:/emtg/PyEMTG/'

'''
PARSE USER INPUT ------------------------------------------------------------------------
'''
test_cases = [];

if len(sys.argv) > 1:
    run_type = sys.argv[1] # In cmd, type all OR cases OR folders OR failed-no ''   
    # If run_type is 'cases' & user specifies a case, set them
    if sys.argv[1] == 'cases' and len(sys.argv) >=3: # Type which case you want to run
        print('Running user cases')
        for casse in sys.argv[2:len(sys.argv)]:
            test_cases.append(casse) 
            
    elif sys.argv[1] == 'folders' or sys.argv[1] == 'folder':
        print('Running all cases in user-specified folders')
        for folderr in sys.argv[2:len(sys.argv)]:
            test_cases.append(folderr) 
            
    elif 'fail' in sys.argv[1]: #this way the user could accidentally type "fail" or "failure" and it will still work
        print('Running failed tests')
        test_cases.append(sys.argv[2]) # Needs to be a *.csv file from the testatron\output dir
        
    elif sys.argv[1] == 'unit':
        print('Running unit (feature) tests')
        run_type = 'folders'
        test_cases = ['global_mission_options','journey_options','physics_options','output_options','physics_options','solver_options','spacecraft_options','transcription_tests','script_constraint_tests']
        
    elif sys.argv[1] == 'mission':
        print('Running mission tests')
        run_type = 'folders'
        test_cases = ['mission_tests']
        
    else: # Else run all b/c can't use 'cases' without actual cases
        print('Running all tests')
        run_type   = 'all'
        test_cases = []
        
else: # Else run all by default
    print('Running all tests');
    run_type = 'all'




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
    summaryFile.write('Beginning test "' + test + '"...')

    if path.isfile(test+'.emtgopt'): #True: #try:
        testOptions = MissionOptions.MissionOptions(test + '.emtgopt') # Emtg ops object

        # Override the output, universe, and HardwareModels paths
        testOptions.override_working_directory = 1
        testOptions.forced_working_directory   = outputdir
        testOptions.short_output_file_names    = 1
        testOptions.background_mode            = 1
        testOptions.override_mission_subfolder = 1
        testOptions.forced_mission_subfolder   = '.'
        testOptions.universe_folder            = test_directory.replace('tests/',\
                                                                        'universe/')
        testOptions.HardwarePath               = test_directory.replace('tests/',\
                                                                'HardwareModels/')
        testOptions.LaunchVehicleLibraryFile = "default.emtg_launchvehicleopt" #LaunchVehicleLibraryFile
        testOptions.PowerSystemsLibraryFile = "default.emtg_powersystemsopt" #PowerSystemsLibraryFile
        testOptions.PropulsionSystemsLibraryFile = "default.emtg_propulsionsystemopt" #PropulsionSystemsLibraryFile

        # Save updated emtg options file
        testOptions.write_options_file(outputdir + '/' + test_name + '.emtgopt') 

        try:
            # Run EMTG
            system(EMTG_path + ' ' + outputdir + '/' + test_name + '.emtgopt')
        except:
            summaryFile.write('\nFAILURE to run "' + test + '" options file.\n\n')
            failFile.write('\nFAILURE to run "' + test + '" options file.\n\n')
            failedTests += 1
            failedTests_list.append(test)
            continue

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
                                                          default_tolerance = 1.0e-10)                                                        
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
        summaryFile.write('\nFAILURE to load "' + test + '" options file.\n\n')
        failFile.write('\nFAILURE to load "' + test + '" options file.\n\n')
        failedTests += 1
        failedTests_list.append(test)

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




