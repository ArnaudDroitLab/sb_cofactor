setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")
source("scripts/chris/countRead_overtime/load_Reddy_peaks_consensus.R")
source("scripts/chris/countRead_overtime/countReads.utils.R")
library(Rsamtools)
library(knitr)
library(dplyr)

all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"

##### List of bam FOSL2
report_fosl2_bam <- make_report_bam(target_name = "FOSL2", all_chip_bam)
bamPath_fosl2 <- generate_bamPath_from_report(report_fosl2_bam, "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/")
indexBam(bamPath_fosl2) # Index FOSL2 bam files (to run only one time is sufficient)

##################################################
#     FOSL2 Peaks
##################################################
##### Gather all FOSL2 peaks all over the time frame
make_ENCODE_Reddy_ChIP_experiments_file("FOSL2")
fosl2_regions <- load_reddy_binding_consensus(target = "FOSL2")
all_fosl2_regions_reduced <- get_reduced_peaks_from_list(fosl2_regions)
summary(width(all_fosl2_regions_reduced))
peaks_fosl2_coordVector <- generate_coordVector(all_fosl2_regions_reduced)

##### Count reads in all FOSL2 peaks all over the time frame
count_total_FOSL2 <- countRead(all_fosl2_regions_reduced, peaks_fosl2_coordVector, bamPath_fosl2, report_fosl2_bam)

save(count_total_FOSL2, file = file.path(output_path, "count_total_FOSL2.RData"))
write.table(count_total_FOSL2, file = file.path(output_path, "count_total_FOSL2.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)