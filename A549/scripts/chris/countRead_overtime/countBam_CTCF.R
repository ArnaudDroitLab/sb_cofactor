setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")
source("scripts/chris/countRead_overtime/load_Reddy_peaks_consensus.R")
source("scripts/chris/countRead_overtime/countReads.utils.R")
library(Rsamtools)
library(knitr)
library(dplyr)

all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"

##### List of bam CTCF
report_ctcf_bam <- make_report_bam(target_name = "CTCF", all_chip_bam)
bamPath_ctcf <- generate_bamPath_from_report(report_ctcf_bam, "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/")
indexBam(bamPath_ctcf) # Index CTCF bam files (to run only one time is sufficient)

##################################################
#     CTCF Peaks
##################################################
##### Gather all CTCF peaks all over the time frame
make_ENCODE_Reddy_ChIP_experiments_file("CTCF")
ctcf_regions <- load_reddy_binding_consensus(target = "CTCF")
all_ctcf_regions_reduced <- get_reduced_peaks_from_list(ctcf_regions)
summary(width(all_ctcf_regions_reduced))
peaks_ctcf_coordVector <- generate_coordVector(all_ctcf_regions_reduced)

##### Count reads in all CTCF peaks all over the time frame
count_total_CTCF <- countRead(all_ctcf_regions_reduced, peaks_ctcf_coordVector, bamPath_ctcf, report_ctcf_bam)

save(count_total_CTCF, file = file.path(output_path, "count_total_CTCF.RData"))
write.table(count_total_CTCF, file = file.path(output_path, "count_total_CTCF.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)