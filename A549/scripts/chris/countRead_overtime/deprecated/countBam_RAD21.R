setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")
source("scripts/chris/countRead_overtime/load_Reddy_peaks_consensus.R")
source("scripts/chris/countRead_overtime/countReads.utils.R")
library(Rsamtools)
library(knitr)
library(dplyr)

all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"

##### List of bam RAD21
report_rad21_bam <- make_report_bam(target_name = "RAD21", all_chip_bam)
bamPath_rad21 <- generate_bamPath_from_report(report_rad21_bam, "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/")
indexBam(bamPath_rad21) # Index RAD21 bam files (to run only one time is sufficient)

##################################################
#     RAD21 Peaks
##################################################
##### Gather all RAD21 peaks all over the time frame
make_ENCODE_Reddy_ChIP_experiments_file("RAD21")
rad21_regions <- load_reddy_binding_consensus(target = "RAD21")
all_rad21_regions_reduced <- get_reduced_peaks_from_list(rad21_regions)
summary(width(all_rad21_regions_reduced))
peaks_rad21_coordVector <- generate_coordVector(all_rad21_regions_reduced)

##### Count reads in all RAD21 peaks all over the time frame
count_total_RAD21 <- countRead(all_rad21_regions_reduced, peaks_rad21_coordVector, bamPath_rad21, report_rad21_bam)

save(count_total_RAD21, file = file.path(output_path, "count_total_RAD21.RData"))
write.table(count_total_RAD21, file = file.path(output_path, "count_total_RAD21.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

