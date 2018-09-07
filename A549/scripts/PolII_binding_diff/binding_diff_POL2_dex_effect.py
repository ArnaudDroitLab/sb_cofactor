import subprocess
import os
import time

shs = ["shCTRL-1", "shCTRL-2", "shNIPBL-3", "shNIPBL-5"]

pwd = subprocess.check_output("pwd", shell = True)
i = pwd.index("A549")
wd = pwd[0:i+4]

read_path = os.path.join(wd, "output/chip-pipeline-PolII-GRCh38/alignment")
peak_path = os.path.join(wd, "output/chip-pipeline-PolII-GRCh38/peak_call")
output_path = os.path.join(wd, "output/chip-pipeline-PolII-GRCh38/binding_diff")

i = 1
for sh in shs:
    start = time.time()
    basename_dex = "POL2_" + sh + "_Dex"
    basename_ctrl = "POL2_" + sh + "_EtOH"
    extension_peak = "_peaks.narrowPeak.bed"
    extension_read = ".sorted.dup.bed"
    
    p1 = os.path.join(peak_path, basename_dex + "N", basename_dex + "N" + extension_peak)
    p2 = os.path.join(peak_path, basename_ctrl + "N", basename_ctrl + "N" + extension_peak)
    r1 = os.path.join(read_path, basename_dex + "_Rep1", basename_dex + "_Rep1" + extension_read)
    r2 = os.path.join(read_path, basename_ctrl + "_Rep1", basename_ctrl + "_Rep1" + extension_read)
    
    output = os.path.join(output_path, "POL2_Dex_vs_EtOH_" + sh)
    mkdir_sub = ["mkdir", "-v", output]
    mkdir_sub2 = " ".join(mkdir_sub)
    subprocess.call(mkdir_sub2, shell=True)
        
    sub = ["manorm",
           "--p1", p1,
           "--p2", p2,
           "--r1", r1,
           "--r2", r2,
           "-o", output]
               
    sub2 = " ".join(sub)
    print("##### " + str(i) + " / 4")
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
