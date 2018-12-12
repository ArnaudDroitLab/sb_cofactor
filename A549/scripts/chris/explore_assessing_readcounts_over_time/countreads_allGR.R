setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")
library(ENCODExplorer)
library(dplyr)
library(knitr)

### regions to assess
gr_regions <- load_reddy_gr_binding_consensus()

all_gr_regions <- GRanges()
for (name in names(gr_regions)) {
  all_gr_regions <- c(all_gr_regions, gr_regions[[name]])
  message(length(gr_regions[[name]]))
}

all_gr_regions_reduced <- reduce(all_gr_regions)
summary(width(all_gr_regions_reduced))

###
download_bam_from_ENCODE <- function(encode_accession) {
  chip_bw_dir = "input/ENCODE/A549/GRCh38/chip-seq/bam"
  dir.create(chip_bw_dir, recursive=TRUE, showWarnings=FALSE)
  downloaded_file = ENCODExplorer::downloadEncode(file_acc = encode_accession, dir=chip_bam_dir, force=FALSE)
}

### explore bam to download
all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
gr_bam <- all_chip_bam %>% filter(target == "NR3C1", assembly == "GRCh38", lab == "Tim Reddy, Duke", submitted_by == "Ian McDowell")
report_bam <- gr_bam %>% select(accession, file_accession, file_format, target, treatment, treatment_duration, treatment_duration_unit, biological_replicates, controls)
kable(report_bam)

