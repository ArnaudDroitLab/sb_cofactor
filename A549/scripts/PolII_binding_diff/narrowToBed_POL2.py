import subprocess
import os
import time

shs = ["shCTRL-1", "shCTRL-2", "shNIPBL-3", "shNIPBL-5"]
conditions = ["Dex", "EtOH"]
targets = ["POL2"]

pwd = subprocess.check_output("pwd", shell = True)
i = pwd.index("A549")
wd = pwd[0:i+4]

peak_path = os.path.join(wd, "output/chip-pipeline-PolII-GRCh38/peak_call")

i = 1
for sh in shs:
    for condition in conditions:
        for target in targets:
            start = time.time()
            basename = target + "_" + sh + "_" + condition + "N"
            
            narrowfile = basename + "_peaks.narrowPeak"
            bedfile = narrowfile + ".bed"
            input_narrow = os.path.join(peak_path, basename, narrowfile)
            output_bed = os.path.join(peak_path, basename, bedfile)

            sub = ["cut", "-f", "1-6", input_narrow, ">", output_bed]
                   
            sub2 = " ".join(sub)
            print("##### " + str(i) + " / 8")
            print(sub2)
            
            subprocess.call(sub2, shell=True)
            
            end = time.time()
            timer = end - start
            print("Conversion from narrowPeak to BED of " + basename + " :\nDONE in " + str(timer) + " seconds")
            print("\n#################################################\n")
            
            i += 1