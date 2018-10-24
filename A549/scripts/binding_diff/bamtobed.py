import subprocess
import os
import time

cofactors = ["BRD4", "CDK9", "MED1", "NIPBL", "SMC1A"]
conditions = ["CTRL", "DEX"]
alignment_path = "/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment"

i = 1
for cofactor in cofactors:
    for condition in conditions:
        start = time.time()
        basename = "A549" + "_" + condition + "_" + cofactor + "_" + "rep1"
        bamfile = basename + ".sorted.dup.bam"
        bedfile = basename + ".sorted.dup.bed"
        input_bam = os.path.join(alignment_path, basename, bamfile)
        output_bed = os.path.join(alignment_path, basename, bedfile)
        
        sub = ["bedtools", "bamtobed",
               "-i", input_bam,
               ">", output_bed]
               
        sub2 = " ".join(sub)
        print("##### " + str(i) + " / 10")
        print(sub2)
        
        subprocess.call(sub2, shell=True)
        
        end = time.time()
        timer = end - start
        print("Conversion from BAM to BED of " + basename + " :\nDONE in " + str(timer) + " seconds")
        print("\n#################################################\n")
        
        i += 1