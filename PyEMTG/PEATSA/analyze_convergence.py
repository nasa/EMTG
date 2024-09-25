"""
analyze_convergence.py
========================

Standlone executable for producing the equivalent of a PEATSA csv summary file.

Calling sequence is

python analyze_convergence.py <path/to/peatsa/output/folder> <'max' or 'min'>

Max or min sets whether the PEATSA run is attempting to maximize or minimize the objective function.

"""

import os
import sys

def main():
    if len(sys.argv) !=3:
        print("Syntax:")
        print("python analyze_convergence.py /path/to/peatsa/folder max_or_min")
    else:
        
        sys.path.append("/Utilities/emtg/PyEMTG")
        import MissionOptions
        import Mission
        
        if sys.argv[2] == "max":
            ifMax = 1
        elif sys.argv[2] == "min":
            ifMax = 0
        else:
            raise Exception("Enter 'max' or 'min' as second command line option")
        
        root = sys.argv[1]
        case_root = root + "/cases/"
        result_root = root + "/results/"
                
        docs = [root + "/docs/" + file for file in os.listdir(root + "/docs") if file.endswith(".csv")]
            
        docs = sorted(docs,key = lambda doc: int(doc.split("/")[-1].lstrip("Iteration").rstrip(".csv")))
        
        final_data = {}
        
        data = {}
        
        has_iter0_cases = False
        
        for doc in docs:
            print(doc) # assume python3
            #if (sys.version_info > (3,0)):
            #    print(doc)
            #else:
            #    print doc
            
            iter_no = int(doc.split("/")[-1].lstrip("Iteration").rstrip(".csv"))
            
            data.update({iter_no:{}})
            
            cases = [file for file in os.listdir(case_root + "Iteration" + str(iter_no)) if file.endswith(".emtgopt")]
            result_cases = [file for file in os.listdir(result_root + "Iteration" + str(iter_no)) if file.endswith(".emtgopt")]
            
            if iter_no == 0 and len(cases):
                has_iter0_cases = True
                
            docH = open(doc,'r')
                    
            for line in docH.readlines():
            
                linesplit = line.split(",")
                
                if line.startswith("Folder"):
                    
                    for idx,entry in enumerate(linesplit):
                        
                        if entry == "PEATSA objective":
                            objcolumn = idx
                        if entry == "Mission Name":
                            namecolumn = idx
                        if entry == "Filename":
                            filecolumn = idx
                                        
                    continue
                    
                name = linesplit[namecolumn].split("_seeded_by")[0]
                objval = float(linesplit[objcolumn])
                filenames = [file for file in cases if file.startswith(name + "_") or file.startswith(name + ".")]
                result_filenames = [file for file in result_cases if (file.startswith(name + "_") or file.startswith(name + ".")) and file.endswith(".emtgopt")]
                emtg_files = [file for file in result_cases if (file.startswith(name) or file.startswith("FAILURE_" + name)) and file.endswith(".emtg")]
                converged_files = [file for file in emtg_files if "FAILURE" not in file]
                seeds = [file.split("_seeded_by_")[1].rstrip("emtgopt").rstrip(".") for file in filenames if len(file.split("_seeded_by"))>1]
                numCases = len(filenames)
                numRans = len(result_filenames)
                numFinished = len(emtg_files)
                numConverged = len(converged_files)
                if "FAIL" in linesplit[filecolumn]:
                    ifConverge = 0
                else:
                    ifConverge = 1
                
                data[iter_no].update({name:[objval,ifConverge,seeds,numCases,numRans,numFinished,numConverged]})
                        
                
            for name in data[iter_no].keys():
    
                if iter_no == 0:
                    final_data.update({name:[]})
                        
                numCases = data[iter_no][name][3]
                numRans = data[iter_no][name][4]
                numFin = data[iter_no][name][5]
                numConv = data[iter_no][name][6]
                
                if data[iter_no][name][1] == 0:
                    ifConverged = 0
                    ifImprove = 0
                    ifFirstFeasible = 0
                else:
                    ifConverged = 1
                    
                
                if iter_no > 0:                
                                    
                    if name not in prev_iteration_data.keys():
                        raise Exception("Case " + name + " is in iteration " + str(iter_no) + " but not in previous.")
                        
                    if data[iter_no][name][1] == 0 and prev_iteration_data[name][1] == 1:
                        raise Exception("Iteration " + str(iter_no) + " has no solution for case " + name + ", but the previous iteration does")
                        
                    
                    seed_code = ""
                    if len(data[iter_no][name][2]) != 0:
                        
                        for seed in data[iter_no][name][2]:
                            
                            if "External" not in seed:
                                if prev_iteration_data[seed][1] == 0: 
                                    seed_code += "F"
                                
                                if seed not in prev_iteration_data[name][2]:
                                    seed_code += "1"
                                else:
                                    iter_no_2 = iter_no
                                                            
                                    while True:
                                        if seed in data[iter_no_2-1][name][2] and data[iter_no_2][seed][0] == data[iter_no_2-1][seed][0]:
                                            iter_no_2 -= 1
                                            if iter_no_2 == 0:
                                                break
                                        else:
                                            break
                                    seed_code += str(iter_no-iter_no_2+1)
                            else:
                                seed_code += "E"
                    
                    if data[iter_no][name][1] != 0:
                        
                        if prev_iteration_data[name][1] == 0:
                            ifFirstFeasible = 1
                            ifImprove = 0
                        else: 
                            ifFirstFeasible = 0
                               
                            prev_obj = prev_iteration_data[name][0]
                            if ifMax:
                                if data[iter_no][name][0] > prev_obj:
                                    ifImprove = 1
                                else:
                                    ifImprove = 0
                            else:
                                if data[iter_no][name][0] < prev_obj:
                                    ifImprove = 1
                                else:
                                    ifImprove = 0
                                            
                    # final_data[name].append(str(ifConverged) + "/" + str(ifFirstFeasible) + "/" + str(ifImprove) + "/" + str(ifImprovedWithoutSeed) + "/" + str(ifRepeatedSeed))
                    final_data[name].append(str(numCases) + "/" + str(numRans) + "/" + str(numFin) + "/" + str(numConv) + "/" + str(ifConverged) + "/" + str(ifFirstFeasible) + "/" + str(ifImprove) + "/" + str(seed_code))
               
                    
               
                else:
                    final_data[name].append(str(numCases) + "/" + str(numRans) + "/" + str(numFin) + "/" + str(numConv) + "/" + str(ifConverged) + "/" + str(ifFirstFeasible) + "//")                
            
            if iter_no > 0:
                if len(data[iter_no].keys()) != len(prev_iteration_data.keys()):
                    for name in prev_iteration_data,keys():
                       if name not in data[iter_no].keys():
                           raise Exception("Case " + name + " is in iteration " + str(iter_no-1) + " but not in the next.")    
            prev_iteration_data = data[iter_no]
            
        
        outfile = open("ConvergenceSummary.csv",'w')
        
        outfile_strings = [""]
        
        new_feasible = []
        improvements = []
        converged = []
        improved_unseeded = []
        ran_unseeded = []
        used_infeasible_ig = []
        cases_created = []
        cases_run = []
        waiting_for_seeds = []
        iterations_before_improvement = []
        converged_this_iteration = []
        cases_finished = []
        
        ifFirst = 1
        for name in final_data.keys():
            outfile_strings.append(name)
            
            ctr = 0
            
            for iter_result in final_data[name]:
                
                if ifFirst:
                    outfile_strings[0] += ",Iteration" + str(ctr)
                    new_feasible.append(0)
                    improvements.append(0)
                    converged.append(0)
                    improved_unseeded.append(0)
                    ran_unseeded.append(0)
                    used_infeasible_ig.append(0)
                    cases_run.append(0)
                    cases_created.append(0)
                    waiting_for_seeds.append(0)
                    iterations_before_improvement.append(0)
                    converged_this_iteration.append(0)
                    cases_finished.append(0)
                
                outfile_strings[-1] += "," + iter_result
                
                result_split = iter_result.split("/")
                
                if ctr!= 0 or has_iter0_cases:                
                    if result_split[0] == "0":
                        waiting_for_seeds[ctr] += 1 
              
                cases_created[ctr] += int(result_split[0])                
                            
                cases_run[ctr] += int(result_split[1])
                
                cases_finished[ctr] += int(result_split[2])
                
                converged_this_iteration[ctr] += int(result_split[3])
                   
                if result_split[4] == "1":
                    converged[ctr] += 1
                    
                if result_split[5] == "1":
                    new_feasible[ctr] += 1  
                    # TODO: This wont work if the case has multiple seeds. in that case we need to figure out which seed fixed it by parsing the seedcode string
                    iterations_before_improvement[int(result_split[7])] += 1          
                
                if result_split[6] == "1":
                    improvements[ctr] += 1
                    # TODO: This wont work if the case has multiple seeds. in that case we need to figure out which seed fixed it by parsing the seedcode string
                    iterations_before_improvement[int(result_split[7])] += 1
                
                if result_split[7] == "" and result_split[0] != "0":
                    if ctr != 0:
                        ran_unseeded[ctr] += 1
                    if result_split[5] == "1" or result_split[6] == "1":
                        improved_unseeded[ctr] += 1
                elif result_split[7] != "":
                    # TODO: This wont work if the case has multiple seeds. in that case we need to figure out which seed fixed it by parsing the seedcode string
                    if result_split[7][0] == "F":
                        used_infeasible_ig[ctr] += 1 
                
                ctr += 1
        
            ifFirst = 0
        
            outfile_strings[-1] += "\n"
        
        outfile_strings[0] += "\n"
        
        outfile.write(outfile_strings[0])
        
        outfile.write("CasesCreated")    
        for im in cases_created:
            outfile.write("," + str(im))
    
        outfile.write("\n")
        
        outfile.write("CasesRun")    
        for im in cases_run:
            outfile.write("," + str(im))
    
        outfile.write("\n")
        
        outfile.write("CasesFinished")    
        for im in cases_finished:
            outfile.write("," + str(im))
    
        outfile.write("\n")
        
        outfile.write("ConvergedThisIteration")    
        for im in converged_this_iteration:
            outfile.write("," + str(im))
    
        outfile.write("\n")
        
        outfile.write("ConvergedTotal")    
        for im in converged:
            outfile.write("," + str(im))
    
        outfile.write("\n")
        
        outfile.write("FirstConverged")
        for im in new_feasible:
            outfile.write("," + str(im))
    
        outfile.write("\n")
            
        outfile.write("Improvements")    
        for im in improvements:
            outfile.write("," + str(im))
            
        outfile.write("\n")
            
        outfile.write("RanUnseeded")    
        for im in ran_unseeded:
            outfile.write("," + str(im))
                
        outfile.write("\n")
            
        outfile.write("UnseededImprovements")    
        for im in improved_unseeded:
            outfile.write("," + str(im))
                
        outfile.write("\n")
        
        outfile.write("InfeasibleIG")    
        for im in used_infeasible_ig:
            outfile.write("," + str(im))
                
        outfile.write("\n")
        
        outfile.write("TimesUsingSeedBeforeImprovement")    
        for im in iterations_before_improvement:
            outfile.write("," + str(im))
                
        outfile.write("\n")
        
        for string in outfile_strings[1:]:
            outfile.write(string)    
            
        outfile.close()
        return
    

if __name__ == "__main__":
    main()