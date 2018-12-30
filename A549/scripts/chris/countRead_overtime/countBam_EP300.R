setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")
source("scripts/chris/countRead_overtime/load_Reddy_EP300_consensus.R")
source("scripts/chris/countRead_overtime/countReads.utils.R")
library(Rsamtools)
library(knitr)
library(dplyr)

##### List of bam EP300
all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
report_ep300_bam <- make_report_bam(target_name = "EP300", all_chip_bam)
bamPath_ep300 <- generate_bamPath_from_report(report_ep300_bam, "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/")

##### Index EP300 bam files (to run only one time is sufficient)
# index_bam(bamPath_ep300)

##################################################
#     EP Peaks
##################################################
##### Gather all EP300 peaks all over the time frame
ep300_regions <- load_reddy_ep300_binding_consensus()
all_ep300_regions_reduced <- get_reduced_peaks_from_list(ep300_regions)
summary(width(all_ep300_regions_reduced))
peaks_ep300_coordVector <- generate_coordVector(all_ep300_regions_reduced)

##### Count reads in all EP300 peaks all over the time frame
count_total_EP300 <- countRead(all_ep300_regions_reduced, peaks_ep300_coordVector, bamPath_ep300, report_ep300_bam)

output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"
save(count_total_EP300, file = file.path(output_path, "count_total_EP300.RData"))
write.table(count_total_EP300, file = file.path(output_path, "count_total_EP300.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

##################################################
#     EP300 Background
##################################################
##### Shift all EP300 peaks to 10000 to get the background
background_ep300_regions_reduced <- shift(all_ep300_regions_reduced, 10000)
peaks_EP300_background_coordVector <- generate_coordVector(background_ep300_regions_reduced)
##### Count reads in EP300 background
count_total_EP300_background <- countRead(background_ep300_regions_reduced, peaks_EP300_background_coordVector, bamPath_ep300, report_ep300_bam)

output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"
save(count_total_EP300, file = file.path(output_path, "count_total_EP300_background.RData"))
write.table(count_total_EP300, file = file.path(output_path, "count_total_EP300_background.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

##################################################
#     EP300 Gap Background
##################################################
##### Get the background: all regions where EP300 are not binding on the genome
gapbackground_ep300_regions_reduced <- gaps(all_ep300_regions_reduced)
summary(width(gapbackground_ep300_regions_reduced))
length(gapbackground_ep300_regions_reduced)
peaks_EP300_gapbackground_coordVector <- generate_coordVector(gapbackground_ep300_regions_reduced)
##### Count reads in EP300 background
count_total_EP300_gapbackground <- countRead(gapbackground_ep300_regions_reduced, peaks_EP300_gapbackground_coordVector, bamPath_ep300, report_ep300_bam)

output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"
save(count_total_EP300, file = file.path(output_path, "count_total_EP300_gaps_background.RData"))
write.table(count_total_EP300, file = file.path(output_path, "count_total_EP300_gaps_background.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)