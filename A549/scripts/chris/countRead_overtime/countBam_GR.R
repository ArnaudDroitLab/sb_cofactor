setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")
source("scripts/chris/countRead_overtime/countReads.utils.R")
library(Rsamtools)
library(knitr)
library(dplyr)

all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"

##### List of bam GR
report_gr_bam <- make_report_bam(target_name = "NR3C1", all_chip_bam)
bamPath_gr <- generate_bamPath_from_report(report_gr_bam, "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/")
# index_bam(bamPath_gr) # Index GR bam files (to run only one time is sufficient)

##### List of bam GR WCE
report_gr_wce_bam <- make_report_WCE_bam(report_gr_bam, all_chip_bam)
bamPath_gr_wce <- generate_bamPath_from_report(report_gr_wce_bam, "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/")
# index_bam(bamPath_gr_wce) # Index GR WCE bam files (to run only one time is sufficient)

##################################################
#     GR Peaks
##################################################
##### Gather all GR peaks all over the time frame
gr_regions <- load_reddy_gr_binding_consensus()
all_gr_regions_reduced <- get_reduced_peaks_from_list(gr_regions)
summary(width(all_gr_regions_reduced))
peaks_GR_coordVector <- generate_coordVector(all_gr_regions_reduced)

##### Count reads in all GR peaks all over the time frame in GR bam
count_total_GR <- countRead(all_gr_regions_reduced, peaks_GR_coordVector, bamPath_gr, report_gr_bam)

save(count_total_GR, file = file.path(output_path, "count_total_GR2.RData"))
write.table(count_total_GR, file = file.path(output_path, "count_total_GR2.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

##### Count reads in all GR peaks all over the time frame in GR WCE
count_total_GR_wce <- countRead(all_gr_regions_reduced, peaks_GR_coordVector, bamPath_gr_wce, report_gr_wce_bam)

save(count_total_GR_wce, file = file.path(output_path, "count_total_GR_wce.RData"))
write.table(count_total_GR_wce, file = file.path(output_path, "count_total_GR_wce.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

##################################################
#     GR Background
##################################################
##### Shift all GR peaks to 10000 to get the background
background_gr_regions_reduced <- shift(all_gr_regions_reduced, 10000)
peaks_GR_background_coordVector <- generate_coordVector(background_gr_regions_reduced)
##### Count reads in GR background
count_total_GR_background <- countRead(background_gr_regions_reduced, peaks_GR_background_coordVector, bamPath_gr, report_gr_bam)

save(count_total_GR_background, file = file.path(output_path, "count_total_GR_background.RData"))
write.table(count_total_GR_background, file = file.path(output_path, "count_total_GR_background.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

##################################################
#     GR Gap Background
##################################################
##### Get the background: all regions where GR are not binding on the genome
gapbackground_gr_regions_reduced <- gaps(all_gr_regions_reduced)
summary(width(gapbackground_gr_regions_reduced))
length(gapbackground_gr_regions_reduced)
peaks_GR_gapbackground_coordVector <- generate_coordVector(gapbackground_gr_regions_reduced)
##### Count reads in GR background
count_total_GR_gapbackground <- countRead(gapbackground_gr_regions_reduced, peaks_GR_gapbackground_coordVector, bamPath_gr, report_gr_bam)

save(count_total_GR, file = file.path(output_path, "count_total_GR_gaps_background.RData"))
write.table(count_total_GR, file = file.path(output_path, "count_total_GR_gaps_background.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)