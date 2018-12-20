setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
library(ENCODExplorer)
library(dplyr)
library(knitr)

###
chip_bam_dir = "input/ENCODE/A549/GRCh38/chip-seq/bam"

###
download_bam_from_ENCODE <- function(encode_accession) {
  chip_bam_dir = "input/ENCODE/A549/GRCh38/chip-seq/bam"
  dir.create(chip_bam_dir, recursive=TRUE, showWarnings=FALSE)
  downloaded_file = ENCODExplorer::downloadEncode(file_acc = encode_accession, dir=chip_bam_dir, force=FALSE)
}

### Download bam GR
all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")

gr_bam <- all_chip_bam %>% filter(target == "NR3C1", assembly == "GRCh38", lab == "Tim Reddy, Duke")
report_gr_bam <- gr_bam %>% select(accession, file_accession, submitted_by, file_format, target, treatment, treatment_duration, treatment_duration_unit, biological_replicates, controls)
report_gr_bam$treatment_duration[is.na(report_gr_bam$treatment_duration)] <- 0
report_gr_bam$treatment_duration_unit[is.na(report_gr_bam$treatment_duration_unit)] <- "minute"
kable(report_gr_bam)

accession_list_gr <- unique(report_gr_bam$file_accession)
cptmax <- length(accession_list_gr)

cpt <- 0
for (acc in accession_list_gr) {
  cpt <- cpt + 1
  message("#########################################################################")
  target <- report_gr_bam %>% filter(file_accession == acc) %>% select(target)
  rep <- report_gr_bam %>% filter(file_accession == acc) %>% select(biological_replicates)
  td <- report_gr_bam %>% filter(file_accession == acc) %>% select(treatment_duration) # td = treatment_duration
  tdu <- report_gr_bam %>% filter(file_accession == acc) %>% select(treatment_duration_unit) # tdu = treatment_duration_unit
  
  message("####### ", acc, "\t", target, "\trep", rep, "\t", cpt, "\\", cptmax)
  message("####### ", td, " ", tdu)
  # download_bam_from_ENCODE(acc)
  
  oldfilename <- file.path(chip_bam_dir, paste0(acc, ".bam"))
  newfilename <- file.path(chip_bam_dir, paste0(target, "_", td, tdu, "_rep", rep, "_", acc, ".bam"))
  
  cmd <- paste0("mv", " ", oldfilename, " ", newfilename)
  message(cmd)
  # system(cmd)
}

### Download bam EP300
ep300_bam <- all_chip_bam %>% filter(target == "EP300", assembly == "GRCh38", lab == "Tim Reddy, Duke")
report_ep300_bam <- ep300_bam %>% select(accession, file_accession, submitted_by, file_format, target, treatment, treatment_duration, treatment_duration_unit, biological_replicates, controls)
report_ep300_bam$treatment_duration[is.na(report_ep300_bam$treatment_duration)] <- 0
report_ep300_bam$treatment_duration_unit[is.na(report_ep300_bam$treatment_duration_unit)] <- "minute"
kable(report_ep300_bam)

accession_list_ep300 <- unique(report_ep300_bam$file_accession)
cptmax <- length(accession_list_ep300)

cpt <- 0
for (acc in accession_list_ep300) {
  cpt <- cpt + 1
  message("#########################################################################")
  target <- report_ep300_bam %>% filter(file_accession == acc) %>% select(target)
  rep <- report_ep300_bam %>% filter(file_accession == acc) %>% select(biological_replicates)
  td <- report_ep300_bam %>% filter(file_accession == acc) %>% select(treatment_duration) # td = treatment_duration
  tdu <- report_ep300_bam %>% filter(file_accession == acc) %>% select(treatment_duration_unit) # tdu = treatment_duration_unit
  
  message("####### ", acc, "\t", target, "\trep", rep, "\t", cpt, "\\", cptmax)
  message("####### ", td, " ", tdu)
  # download_bam_from_ENCODE(acc)
  
  oldfilename <- file.path(chip_bam_dir, paste0(acc, ".bam"))
  newfilename <- file.path(chip_bam_dir, paste0(target, "_", td, tdu, "_rep", rep, "_", acc, ".bam"))
  
  cmd <- paste0("mv", " ", oldfilename, " ", newfilename)
  message(cmd)
  # system(cmd)
}
