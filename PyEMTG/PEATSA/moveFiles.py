import sys
import os
import shutil
import time
import datetime

if len(sys.argv) != 4:
    print("Syntax:")
    print("python moveFiles.py /full/path/to/destination results_name /full/path/to/PEATSA_run")
else:
    
    peatsa_dir = sys.argv[3]

    results_name = sys.argv[2]

    results_base = sys.argv[1] + "/" + results_name

    latest_iteration = -1
    while True:
        csvfiles = [file for file in os.listdir(peatsa_dir + "/docs") if file.endswith(".csv")]
        allimgfiles = [file for file in os.listdir(peatsa_dir + "/images")]


        iteration = -1


        for csvfile in csvfiles:
            if csvfile.startswith("PEATSA"):
                continue
            iterno = int(csvfile.lstrip("Iteration").rstrip(".csv"))

            if iterno > iteration:
                iteration = iterno

        if iteration == -1:
            continue

        if iteration != latest_iteration:
            need_to_update = True
        else:
            need_to_update = False


        if need_to_update:
            latest_iteration = iteration
            if os.path.isfile(results_base + ".csv"):
                os.remove(results_base + ".csv")
            shutil.copyfile(peatsa_dir + "/docs/Iteration" + str(latest_iteration) + ".csv",results_base + ".csv")
            if os.path.isfile(results_base + "_history.csv"):
                os.remove(results_base + "_history.csv")
            shutil.copyfile(peatsa_dir + "/docs/PEATSAhistory.csv",results_base + "_history.csv")
            if os.path.isfile(results_base + ".gif"):
                os.remove(results_base + ".gif")
            giffiles = [file for file in os.listdir(peatsa_dir + "/images") if file.endswith(".gif")]
            if len(giffiles):
                gifcounter = 0
                for giffile in giffiles:
                    shutil.copyfile(peatsa_dir + "/images/" + giffile,results_base + "_img" + str(gifcounter) + ".gif")
                    gifcounter+=1
            img_files = [file for file in os.listdir(peatsa_dir + "/images") if file.endswith(".png")]
            img_ctr = 1
            for img in img_files:
                if "Iteration" + str(latest_iteration) in img:
                    if os.path.isfile(results_base + "_img" + str(img_ctr) + ".png"):
                        os.remove(results_base + "_img" + str(img_ctr) + ".png")
                    shutil.copyfile(peatsa_dir + "/images/" + img,results_base + "_img" + str(img_ctr) + ".png")
                    img_ctr += 1
            os.makedirs(results_name)
            os.system("python /workdir/jeremymknittel/EMTG/PyEMTG/PEATSA/grab_best_peatsa_results.py " + peatsa_dir + " " + results_name)
            os.system("tar czf " + results_name + ".tar.gz " + results_name)
            shutil.rmtree(results_name)
            if os.path.isfile(results_base + ".tar.gz"):
                os.remove(results_base + ".tar.gz")
            shutil.move(results_name + ".tar.gz",results_base + ".tar.gz")
            handle = open(results_base + ".txt",'w')
            handle.write("Iteration " + str(latest_iteration) + " written at " + datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
            handle.close()


        time.sleep(60)




