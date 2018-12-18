setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")
library(Rsamtools)
library(knitr)
library(dplyr)

####
generate_coordVector <- function(granges) {
  df <- as.data.frame(granges)
  coord_Vector <- paste0(df$seqnames, ":", df$start, "-", df$end)
  return (coord_Vector)
}

#### Gather all GR peaks all over the time frame
gr_regions <- load_reddy_gr_binding_consensus()

all_gr_regions <- GRanges()
for (name in names(gr_regions)) {
  all_gr_regions <- c(all_gr_regions, gr_regions[[name]])
  message(length(gr_regions[[name]]))
}

all_gr_regions_reduced <- reduce(all_gr_regions)
summary(width(all_gr_regions_reduced))

peaks_GR_coordVector <- generate_coordVector(all_gr_regions_reduced)

##### List of bam GR
all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
gr_bam <- all_chip_bam %>% filter(target == "NR3C1", assembly == "GRCh38", lab == "Tim Reddy, Duke")
report_gr_bam <- gr_bam %>% select(accession, file_accession, submitted_by, file_format, target, treatment, treatment_duration, treatment_duration_unit, biological_replicates, controls)
report_gr_bam$treatment_duration[is.na(report_gr_bam$treatment_duration)] <- 0
report_gr_bam$treatment_duration_unit[is.na(report_gr_bam$treatment_duration_unit)] <- "minute"
report_gr_bam <- report_gr_bam %>% arrange(desc(treatment_duration_unit), treatment_duration, submitted_by, biological_replicates)
kable(report_gr_bam)

bam_path <- paste0("/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/",
                   paste0(report_gr_bam$target, "_", report_gr_bam$treatment_duration, report_gr_bam$treatment_duration_unit,
                          "_rep", report_gr_bam$biological_replicates,
                          "_", report_gr_bam$file_accession, ".bam"))

##### Count reads
which <- all_gr_regions_reduced
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)

count_total_GR <- data.frame(peaks_GR_coordVector)
for (bam in bam_path) {
  message(bam)
  index_bam <- gsub("\\.bam", "", bam)
  count_bySample <- countBam(bam, index = index_bam, param = param)
  records <- count_bySample$records
  count_total_GR <- data.frame(count_total_GR, records)
}

names(count_total_GR) <- c("Coordinates", paste0(report_gr_bam$target, "_", report_gr_bam$treatment_duration, report_gr_bam$treatment_duration_unit,
                                                 "_rep", report_gr_bam$biological_replicates,
                                                 "_", report_gr_bam$file_accession))

output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_GR_overtime"
save(count_total_GR, file = file.path(output_path, "count_total_GR2.RData"))
write.table(count_total_GR, file = file.path(output_path, "count_total_GR2.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
