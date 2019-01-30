setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")
source("scripts/chris/countRead_overtime/load_Reddy_peaks_consensus.R")
source("scripts/chris/countRead_overtime/countReads.utils.R")
library(Rsamtools)
library(knitr)
library(dplyr)

all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"

##### List of bam SMC3
report_smc3_bam <- make_report_bam(target_name = "SMC3", all_chip_bam)
bamPath_smc3 <- generate_bamPath_from_report(report_smc3_bam, "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/")
indexBam(bamPath_smc3) # Index SMC3 bam files (to run only one time is sufficient)

##################################################
#     SMC3 Peaks
##################################################
##### Gather all SMC3 peaks all over the time frame
make_ENCODE_Reddy_ChIP_experiments_file("SMC3")
smc3_regions <- load_reddy_binding_consensus(target = "SMC3")
all_smc3_regions_reduced <- get_reduced_peaks_from_list(smc3_regions)
summary(width(all_smc3_regions_reduced))
peaks_smc3_coordVector <- generate_coordVector(all_smc3_regions_reduced)

##### Count reads in all SMC3 peaks all over the time frame
count_total_SMC3 <- countRead(all_smc3_regions_reduced, peaks_smc3_coordVector, bamPath_smc3, report_smc3_bam)

save(count_total_SMC3, file = file.path(output_path, "count_total_SMC3.RData"))
write.table(count_total_SMC3, file = file.path(output_path, "count_total_SMC3.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

