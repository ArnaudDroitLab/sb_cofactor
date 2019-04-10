import subprocess
import os
import time

cofactors = ["BRD4", "CDK9", "MED1", "NIPBL", "SMC1A"]
bd_path = "/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/binding_diff"

i = 1
for cofactor in cofactors:
    start = time.time()
    
    A_cofactor = "A549_" + cofactor
    basename_dex = "A549_DEX_" + cofactor + "_rep1"
    basename_ctrl = "A549_CTRL_" + cofactor + "_rep1"
    peaks_path = os.path.join(bd_path, A_cofactor, "output_filters")
    
    up_peaks = os.path.join(peaks_path, basename_dex + "_peaks.narrowPeak_M_above_0.8_biased_peaks.bed")
    down_peaks = os.path.join(peaks_path, basename_ctrl + "_peaks.narrowPeak_M_below_-0.8_biased_peaks.bed")
    unbiased_peaks = os.path.join(peaks_path, A_cofactor + "_unbiased_peaks.bed")
    
    basename_up_peaks = A_cofactor + "_up_peaks"
    basename_down_peaks = A_cofactor + "_down_peaks"
    basename_unbiased_peaks = A_cofactor + "_unbiased_peaks"
    
    ###
    sub = ["sort -k1,1 -k2,2n", up_peaks, "|", "cut -f1-4", ">", os.path.join(peaks_path, basename_up_peaks + ".sorted.bed")]
    sub2 = " ".join(sub)
    print(sub2)
    subprocess.call(sub2, shell=True)
    
    sub = ["sort -k1,1 -k2,2n", down_peaks, "|", "cut -f1-4", ">", os.path.join(peaks_path, basename_down_peaks + ".sorted.bed")]
    sub2 = " ".join(sub)
    print(sub2)
    subprocess.call(sub2, shell=True)
    
    sub = ["sort -k1,1 -k2,2n", unbiased_peaks, "|", "cut -f1-4", ">", os.path.join(peaks_path, basename_unbiased_peaks + ".sorted.bed")]
    sub2 = " ".join(sub)
    print(sub2)
    subprocess.call(sub2, shell=True)
    
    ###
    bedToBidBed = "/home/chris/Documents/bedToBigBed"
    hg38_chromsizes = "/home/chris/Documents/hg38.chrom.sizes"
    
    sub = [bedToBidBed,
           os.path.join(peaks_path, basename_up_peaks + ".sorted.bed"),
           hg38_chromsizes,
           os.path.join(peaks_path, basename_up_peaks + ".sorted.bb")]
    sub2 = " ".join(sub)
    print(sub2)
    subprocess.call(sub2, shell=True)
    
    sub = [bedToBidBed,
           os.path.join(peaks_path, basename_down_peaks + ".sorted.bed"),
           hg38_chromsizes,
           os.path.join(peaks_path, basename_down_peaks + ".sorted.bb")]
    sub2 = " ".join(sub)
    print(sub2)
    subprocess.call(sub2, shell=True)

    sub = [bedToBidBed,
           os.path.join(peaks_path, basename_unbiased_peaks + ".sorted.bed"),
           hg38_chromsizes,
           os.path.join(peaks_path, basename_unbiased_peaks + ".sorted.bb")]
    sub2 = " ".join(sub)
    print(sub2)
    subprocess.call(sub2, shell=True)
        
    end = time.time()
    timer = end - start
    print("\nDifferential binding on " + A_cofactor + " samples has been done")
    print("Processed time : " + str(timer) + " seconds")
    print("\n#################################################\n")