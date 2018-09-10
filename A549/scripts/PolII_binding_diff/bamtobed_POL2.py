import subprocess
import os
import time

shs = ["shCTRL-1", "shCTRL-2", "shNIPBL-3", "shNIPBL-5"]
conditions = ["Dex", "EtOH"]
targets = ["POL2"]

pwd = subprocess.check_output("pwd", shell = True)
i = pwd.index("A549")
wd = pwd[0:i+4]

alignment_path = os.path.join(wd, "output/chip-pipeline-PolII-GRCh38/alignment")


i = 1
for sh in shs:
    for condition in conditions:
        for target in targets:
            start = time.time()
            basename = target + "_" + sh + "_" + condition + "_" + "Rep1"
            
            bamfile = basename + ".sorted.dup.bam"
            bedfile = basename + ".sorted.dup.bed"
            input_bam = os.path.join(alignment_path, basename, bamfile)
            output_bed = os.path.join(alignment_path, basename, bedfile)
            
            sub = ["bedtools", "bamtobed",
                   "-i", input_bam,
                   ">", output_bed]
                   
            sub2 = " ".join(sub)
            print("##### " + str(i) + " / 8")
            print(sub2)
            
            # subprocess.call(sub2, shell=True)
            
            end = time.time()
            timer = end - start
            print("Conversion from BAM to BED of " + basename + " :\nDONE in " + str(timer) + " seconds")
            print("\n#################################################\n")
            
            i += 1
