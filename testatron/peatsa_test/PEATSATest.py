

def PEATSATest(path_to_peatsa_input, path_to_pyemtg_folder):

    path_to_peatsa_input = path_to_peatsa_input.replace('\\', '/')
    path_in_pieces = path_to_peatsa_input.split('/')

    input_name = path_in_pieces[-1]
    input_dir = '/'.join(path_in_pieces[:-1]) + '/'


    import os
    import sys
    path_to_pyemtg_folder = path_to_pyemtg_folder.replace('\\', '/')
    sys.path.append(path_to_pyemtg_folder)

    from MissionOptions import MissionOptions
    from Mission import Mission

    import pandas
    import itertools
    #csv_file = pandas.read_csv('peatsa_test.csv', header = 2)
    exec(open(input_dir + input_name, 'r').read())

    test_dir = input_dir + run_name + '/'

    #region  check case generation (count & fingerprints)
    n_cases = 0

    #need to loop through each options file (may be more than 1 csv informning the PEATSA run)
    #need to obtain the number of expected cases as well as the list of expected fingerprints
    all_cases = []
    all_fingerprints = []
    case_fingerprint_map = {}
    case_input_map = {}
    for tup in trade_study_options_files:
        file = tup[0] #options file is the first item in the tuple
        type = tup[1] #options file type is the second item (all combinations vs. explicitly provided sets)

        file_name = file.split('/')[-1]
        temp = pandas.read_csv(input_dir + file_name, header = 2) #header is the 3rd row of the csv data

        fingerprint_attrs = list(temp.columns) #assuming this is consistent across csvs (should be unless problem isn't set up correctly)

        if type == 2: #user explicitly provided the case list
            start_case = n_cases
            n_cases += len(temp) #add the number of rows in the csv to the number of cases

            case_list = ['Case' + str(n) for n in range(start_case, n_cases)]
            fingerprints = list(zip(*[temp[col] for col in temp]))

            all_cases = all_cases + case_list
            all_fingerprints = all_fingerprints + fingerprints
            case_fingerprint_map = {**case_fingerprint_map, **dict(zip(case_list, fingerprints))}

        elif type == 1: #need to perform the combinatorics
            cols_wout_nans = []
            for col in temp.columns:
                #NOTE: The columns won't all be the same length but they will all have the same number of rows in the pandas
                #      df. Need to filter out the nan values in each column.
                non_nans = temp.loc[~temp[col].isna()]
                cols_wout_nans.append(non_nans[col].tolist())

            fingerprints = list(itertools.product(*cols_wout_nans))

            start_case = n_cases
            n_cases += len(fingerprints)
            case_list = ['Case' + str(n) for n in range(start_case, n_cases)]

            all_cases = all_cases + case_list
            all_fingerprints = all_fingerprints + fingerprints
            case_fingerprint_map = {**case_fingerprint_map, **dict(zip(case_list, fingerprints))}

    #now have the total number of cases stored in the n_cases variable
    #have to create a list of unique Case numbers from the cases folder (can just use Iteration0 for this)
    #need to keep in mind that a Case may appear more than once in the folder if it has multiple seeds
    unique_cases = []

    #create a list of all the input files in Iteration0
    all_files = [file for file in os.listdir(test_dir + 'cases/Iteration0/') if file.endswith('.emtgopt')]

    #loop through the input files and parse out unique cases
    for file in all_files:
        foo = file.split('.emtgopt')[0]
        bar = foo.split('_')
        case = bar[1]

        if case not in unique_cases:
            unique_cases.append(case)
            case_input_map[case] = file #just store the first file for each case since the information we care about will be the same across all case instances

    unique_cases = set(unique_cases)
    expected_cases = set(['Case' + str(n) for n in range(n_cases)]) #also need a set of the cases we expect to be there

    #now figure out if any cases are missing
    missing_cases = list(expected_cases - unique_cases)
    excess_cases = list(unique_cases - expected_cases)


    def GetFingerprint(emtgopt, attributes):
        MO = MissionOptions(emtgopt)

        fingerprint = tuple()
        for attr in attributes:
            new_fingerprint_val = eval(attr)
            fingerprint = fingerprint + (new_fingerprint_val,)

        return fingerprint


    #now validate the fingerprints of each case
    mismatched_fingerprints = []
    for case in case_fingerprint_map:
        full_path = test_dir + 'cases/Iteration0/' + case_input_map[case]

        case_fingerprint = GetFingerprint(full_path, fingerprint_attrs)

        if case_fingerprint != case_fingerprint_map[case]:
            mismatched_fingerprints.append((case, case_fingerprint_map[case], case_fingerprint))

    #endregion


    #region  emtg running sucessfully
    iterations_list = os.listdir(test_dir + 'results/')

    #make sure there are the same number of outputs as inputs
    #make sure no output file starts with the 'FAILURE_' flag
    failed_outputs = {}
    missing_output = {}
    excess_output   = {}
    for iteration in iterations_list:
        n_inputs = [file.split('.emtgopt')[0] for file in os.listdir(test_dir + 'cases/' + iteration) if file.endswith('.emtgopt')]
        n_outputs = [file.split('.emtg')[0] for file in os.listdir(test_dir + 'results/' + iteration) if file.endswith('.emtg')]

        #first check if there were duplicated input files (would have overwrite issues in the outputs)
        duplicated_inputs = [file for file in n_inputs if n_inputs.count(file) > 1]

        #also check if any outputs failed
        failed_outputs[iteration] = []
        fixed_output = [] #need a list without any 'FAILURE_' flags to compare against the inputs
        for file in n_outputs:
            if file.startswith('FAILURE_'):
                failed_outputs[iteration].append(file)
                fixed_output.append(file.split('FAILURE_')[-1])
            else:
                fixed_output.append(file)

        #now see if any inputs failed to generate outputs
        missing_output[iteration] = list(set(n_inputs) - set(fixed_output))
        #and see if any outputs show up that don't belong
        excess_output[iteration] = list(set(fixed_output) - set(n_inputs))

    #endregion


    #region  seeding correctly

    #NOTE: Neighbors won't change per iteration, only need to determine neighbors for original cases
    #first unpack the seed_criteria information into lists
    seed_attrs  = [] #attributes to be used in determining neighbors (e.g. 'MO.launch_window_open_date')
    seed_bounds = {} #bounds that define what a neighboring case is, stored in a tuple (lb, ub)
    seed_types  = {} #seed selection information (e.g. seed from all, seed from best, etc.)
    for criteria in seed_criteria:
        seed_attrs.append(criteria[0])
        seed_bounds[criteria[0]] = (criteria[1], criteria[2])
        seed_types[criteria[0]] = criteria[3]

    #now collect the seed attribute information into a dataframe for each case
    df_cols = []
    for attr in seed_attrs:
        df_cols = df_cols + [attr, attr + ' neighbors']

    seed_dat = pandas.DataFrame(columns = df_cols)
    for case in case_fingerprint_map:
        MO = MissionOptions(test_dir + 'cases/Iteration0/TradeStudy_' + case + '.emtgopt')
        for attr in seed_attrs:
            seed_dat.at[case, attr] = eval(attr)

    #now that the attribute values are all collected, we can locate the neighbors
    #NOTE: Need to find neighbors per-criteria (depends on user inputs but will collect it this way to accommodate all runs)
    for case in case_fingerprint_map:
        for attr in seed_attrs:
            val = seed_dat.at[case, attr]
            lb = val + seed_bounds[attr][0]
            ub = val + seed_bounds[attr][1]

            neighbor_list = list(seed_dat.loc[(lb <= seed_dat[attr]) & (seed_dat[attr] <= ub)].index)
            if not allow_cases_to_seed_themselves:
                neighbor_list.remove(case)
            seed_dat.at[case, attr + ' neighbors'] = neighbor_list

    #if only_one_seed_from_best_in_all_seed_directions then finding seeds is simpler, just need to identify the best neighbor from all possible lists
    neighbor_cols = [col for col in seed_dat.columns if 'neighbors' in col]
    seed_dat['all neighbors'] = seed_dat[neighbor_cols].apply(lambda row: sum(row, []), axis = 1)

    #now we know who each case's neighbors are; the seeds themselves can only be determined by the data in each iteration
    #look at the results for each iteration, determine which case(s) should have seeded each, then check the input file names
    # in the next iteration to determine if the correct seed(s) was used
    def FindSeeds(case, iteration, neighbors_list, selection_type):
        #for either selection_type == 1 or selection_type == 2, need to know the best result for each neighboring case
        best_result_for_each_neighbor = []
        best_objective_val_for_each_neighbor = []
        for neighbor in neighbors_list:
            neighbor_results = [file for file in os.listdir(test_dir + 'results/' + iteration)
                                if (file.startswith('TradeStudy_' + neighbor + '_')
                                or file.startswith('TradeStudy_' + neighbor + '.'))
                                and file.endswith('.emtg')]
            best_case = None
            best_objective = 1e20
            for result in neighbor_results:
                M = Mission(test_dir + 'results/' + iteration + '/' + result)
                if M.objective_value < best_objective:
                    best_case = result
                    best_objective = M.objective_value
            best_result_for_each_neighbor.append(best_case)
            best_objective_val_for_each_neighbor.append(best_objective)

        if selection_type == 1: #use all seeds, just return the best result from each possible neighbor
            return best_result_for_each_neighbor

        elif selection_type == 2: #use the best seed, need to check results for each potential neighbor
            best_result_ix = best_objective_val_for_each_neighbor.index(min(best_objective_val_for_each_neighbor))
            best_neighbor = best_result_for_each_neighbor[best_result_ix]
            return [best_neighbor]


    missing_seeds = {} #first key is the iteration
    excess_seeds  = {} #first key is the iteration
    for ix, iteration in enumerate(iterations_list[:-1]):
        emtg_output = [file for file in os.listdir(test_dir + 'results/' + iteration) if file.endswith('.emtg')]
        missing_seeds[iteration] = {} #second key is the case
        excess_seeds[iteration]  = {} #second key is the case

        for case in case_fingerprint_map:
            missing_seeds[iteration][case] = []
            excess_seeds[iteration][case]  = []

            if only_one_seed_from_best_in_all_seed_directions:
                neighbors = seed_dat.at[case, 'all neighbors']
                seeds = FindSeeds(case, iteration, neighbors, 2)

            else:
                seeds = []
                for attr in seed_attrs:
                    neighbors = seed_dat.at[case, attr + ' neighbors']
                    new_seeds = FindSeeds(case, iteration, neighbors, seed_types[attr])

                    seeds = seeds + new_seeds

            seeds = list(set(seeds)) #remove duplicates
            altered_seeds_list = [seed.split('.emtg')[0] for seed in seeds]
            expected_cases = ['TradeStudy_' + case + '_seeded_by_TradeStudy_' + seed.split('_')[1] + '.emtgopt'
                              for seed in altered_seeds_list]
            #the list of cases in the seeds list is all of the neighbors taht our current case should be seeded by in the next iteration
            #now check the next iteration's cases/ folder and see if all seeds are accounted for
            list_of_seeded_cases = [file for file in os.listdir(test_dir + 'cases/' + iterations_list[ix + 1])
                                    if file.endswith('.emtgopt') and file.startswith('TradeStudy_' + case + '_')]

            missing_seeds[iteration][case] = missing_seeds[iteration][case] + list(set(expected_cases) - set(list_of_seeded_cases))
            excess_seeds[iteration][case] = excess_seeds[iteration][case] + list(set(list_of_seeded_cases) - set(expected_cases))

    #endregion


    #region  csv generated correctly

    #final check is to make sure the csv column names are correct and that data appears for each case in the csv
    csv_header_cols = ['Folder', 'Filename', 'Mission Name', 'Instances Created', 'Instances Ran', 'Instances Finished',
                       'Instances Converged', 'First Feasible', 'Improved', 'Last Improved Iteration', 'Seed History',
                       'SeedCriteria 0', 'Final Mass [kg]', 'Flight Time [years]', 'Feasibility', 'PEATSA objective',
                       'First NLP Solve Feasible', 'First NLP Solve Objective Function Value', 'Number of Solution Attempts',
                       'Best Solution Attempt']
    csv_user_cols = [tup[0] for tup in extra_csv_column_definitions]
    csv_seed_cols = ['S0F' + str(num) for num in range(len(seed_attrs))]

    expected_csv_cols = csv_header_cols + csv_user_cols + csv_seed_cols

    missing_csv_cols = {}; excess_csv_cols  = {}
    missing_csv_data = {}; excess_csv_data  = {}
    #now load in each csv file and check it
    for iteration in iterations_list:
        csv_dat = pandas.read_csv(test_dir + 'docs/' + iteration + '.csv')
        csv_cols = [col for col in csv_dat.columns if 'Unnamed: ' not in col]
        #NOTE: Sometimes pandas reads in blank columns; depends on how the csv is written. Have to removed "Unnamed" cols.

        #check that the csv has the correct data columns
        missing_csv_cols[iteration] = list(set(expected_csv_cols) - set(csv_cols))
        excess_csv_cols[iteration] = list(set(csv_cols) - set(expected_csv_cols))

        #check that each case is represented in the data
        output_files = [file for file in os.listdir(test_dir + 'results/' + iteration) if file.endswith('.emtg')]
        missing_csv_data[iteration] = list(set(output_files) - set(csv_dat['Filename']))
        excess_csv_data[iteration] = list(set(csv_dat['Filename']) - set(output_files))

    #endregion


    #region  collect and report

    report = pandas.DataFrame(columns = [''])
    passed = True

    #are any cases missing/did any unexpected cases show up?
    if len(missing_cases) > 0:
        report.at['Cases that should have been created but were not: ', ''] = missing_cases
        passed = False
    if len(excess_cases) > 0:
        report.at['Cases that should not have been created but were: ', ''] = excess_cases
        passed = False

    #are the fingerprints correct for each case?
    if len(mismatched_fingerprints) > 0:
        #(case, expected fingerprint, actual fingerprint)
        report.at['Cases with incorrect fingerprints (case, expected fp, actual fp): ', ''] = mismatched_fingerprints
        passed = False

    #did any cases fail?
    for iteration in failed_outputs:
        if len(failed_outputs[iteration]) > 0:
            report.at['Cases failed in ' + iteration + ': ', ''] = failed_outputs[iteration]
            passed = False

    #are any results missing/did any unexpected results show up?
    for iteration in missing_output:
        if len(missing_output[iteration]) > 0:
            report.at['Cases missing output files in ' + iteration + ': ', ''] = missing_output[iteration]
            passed = False
        if len(excess_output[iteration]) > 0:
            report.at['Cases with unexpected output files in ' + iteration + ': ', ''] = excess_output[iteration]
            passed = False

    #are the cases being seeded correctly?
    for ix, iteration in enumerate(missing_seeds):
        # want to report the seeds for the iteration they're impacting, not the iteration they're from
        report_itr = 'Iteration' + str(ix + 1)
        for case in missing_seeds[iteration]:
            if len(missing_seeds[iteration][case]) > 0:
                report.at['Epected seeding(s) for ' + case + ' missing in ' + report_itr + ': ', ''] = missing_seeds[iteration][case]
                passed = False
            if len(excess_seeds[iteration][case]) > 0:
                report.at['Unexpected seeding(s) for ' + case + ' in ' + report_itr + ': ', ''] = excess_seeds[iteration][case]
                passed = False

    #are the correct columns appearing in the csv output?
    for iteration in missing_csv_cols:
        if len(missing_csv_cols[iteration]) > 0:
            report.at['Missing csv column(s) in ' + iteration + ': ', ''] = missing_csv_cols[iteration]
            passed = False
        if len(excess_csv_cols[iteration]) > 0:
            report.at['Unexpected csv column(s) in ' + iteration + ': ', ''] = excess_csv_cols[iteration]
            passed = False
        if len(missing_csv_data[iteration]) > 0:
            report.at['Missing data row in ' + iteration + ' for case(s): ', ''] = missing_csv_data[iteration]
            passed = False
        if len(excess_csv_data[iteration]) > 0:
            report.at['Unexpected data row in ' + iteration + ' for case(s): ', ''] = excess_csv_data[iteration]
            passed = False

    return passed, report

    #endregion