import subprocess
import os
import time

# cofactors = ["BRD4", "CDK9", "MED1", "NIPBL", "SMC1A"]
cofactors = ["SMC1A"]

manorm = "/usr/local/bin/manorm"
read_path = "/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment"
peak_path = "/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/peak_call"
output_path = "/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/binding_diff"

i = 1
for cofactor in cofactors:
    start = time.time()
    basename_dex = "A549_DEX_" + cofactor + "_rep1"
    basename_ctrl = "A549_CTRL_" + cofactor + "_rep1"
    extension_peak = "_peaks.narrowPeak.bed"
    extension_read = ".sorted.dup.bed"
    
    p1 = os.path.join(peak_path, basename_dex, basename_dex + extension_peak)
    p2 = os.path.join(peak_path, basename_ctrl, basename_ctrl + extension_peak)
    r1 = os.path.join(read_path, basename_dex, basename_dex + extension_read)
    r2 = os.path.join(read_path, basename_ctrl, basename_ctrl + extension_read)
    output = os.path.join(output_path, "A549_" + cofactor )
        
    sub = [manorm,
           "--p1", p1,
           "--p2", p2,
           "--r1", r1,
           "--r2", r2,
           "-o", output]
               
    sub2 = " ".join(sub)
    print("##### " + str(i) + " / 5")
    print(sub2)
        
    subprocess.call(sub2, shell=True)
        
    end = time.time()
    timer = end - start
    print("\nDifferential binding on " + output + " samples has been done")
    print("Processed time : " + str(timer) + " seconds")
    print("\n#################################################\n")
    
    i += 1


"""
/usr/local/bin/manorm
--p1 A549_DEX_BRD4_rep1_peaks.narrowPeak.bed
--p2 A549_CTRL_BRD4_rep1_peaks.narrowPeak.bed
--r1 A549_DEX_BRD4_rep1.sorted.dup.bed
--r2 A549_CTRL_BRD4_rep1.sorted.dup.bed
-o A549_BRD4
"""
