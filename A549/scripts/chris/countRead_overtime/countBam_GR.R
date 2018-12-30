setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")
source("scripts/chris/countRead_overtime/countReads.utils.R")
library(Rsamtools)
library(knitr)
library(dplyr)

##### List of bam GR
all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
report_gr_bam <- make_report_bam(target_name = "NR3C1", all_chip_bam)
bamPath_gr <- generate_bamPath_from_report(report_gr_bam, "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/")

##################################################
#     GR Peaks
##################################################
##### Gather all GR peaks all over the time frame
gr_regions <- load_reddy_gr_binding_consensus()
all_gr_regions_reduced <- get_reduced_peaks_from_list(gr_regions)
summary(width(all_gr_regions_reduced))
peaks_GR_coordVector <- generate_coordVector(all_gr_regions_reduced)

##### Count reads in all GR peaks all over the time frame
count_total_GR <- countRead(all_gr_regions_reduced, peaks_GR_coordVector, bamPath_gr, report_gr_bam)

output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"
save(count_total_GR, file = file.path(output_path, "count_total_GR2.RData"))
write.table(count_total_GR, file = file.path(output_path, "count_total_GR2.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

##################################################
#     GR Background
##################################################
##### Shift all GR peaks to 10000 to get the background
background_gr_regions_reduced <- shift(all_gr_regions_reduced, 10000)
peaks_GR_background_coordVector <- generate_coordVector(background_gr_regions_reduced)
##### Count reads in GR background
count_total_GR_background <- countRead(background_gr_regions_reduced, peaks_GR_background_coordVector, bamPath_gr, report_gr_bam)

output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"
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

output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"
save(count_total_GR, file = file.path(output_path, "count_total_GR_gaps_background.RData"))
write.table(count_total_GR, file = file.path(output_path, "count_total_GR_gaps_background.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ##### Count reads
# which <- all_gr_regions_reduced
# what <- c("rname", "strand", "pos", "qwidth", "seq")
# param <- ScanBamParam(which=which, what=what)
# 
# count_total_GR <- data.frame(peaks_GR_coordVector)
# for (bam in bam_path) {
#   message(bam)
#   index_bam <- gsub("\\.bam", "", bam)
#   count_bySample <- countBam(bam, index = index_bam, param = param)
#   records <- count_bySample$records
#   count_total_GR <- data.frame(count_total_GR, records)
# }
# 
# names(count_total_GR) <- c("Coordinates", paste0(report_gr_bam$target, "_", report_gr_bam$treatment_duration, report_gr_bam$treatment_duration_unit,
#                                                  "_rep", report_gr_bam$biological_replicates,
#                                                  "_", report_gr_bam$file_accession))