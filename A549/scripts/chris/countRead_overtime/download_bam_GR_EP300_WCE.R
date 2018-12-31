setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/chris/countRead_overtime/countReads.utils.R")
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

download_bam_in_AccessionList <- function(accession_list, report_bam) {
  cpt <- 0
  cptmax <- length(accession_list)
  for (acc in accession_list) {
    cpt <- cpt + 1
    message("#########################################################################")
      target <- report_bam %>% filter(file_accession == acc) %>% select(target)
      rep <- report_bam %>% filter(file_accession == acc) %>% select(biological_replicates)
      td <- report_bam %>% filter(file_accession == acc) %>% select(treatment_duration) # td = treatment_duration
      tdu <- report_bam %>% filter(file_accession == acc) %>% select(treatment_duration_unit) # tdu = treatment_duration_unit
      
        message("####### ", acc, "\t", target, "\trep", rep, "\t", cpt, "\\", cptmax)
        message("####### ", td, " ", tdu)
        # download_bam_from_ENCODE(acc)

        oldfilename <- file.path(chip_bam_dir, paste0(acc, ".bam"))
        newfilename <- file.path(chip_bam_dir, paste0(target, "_", td, tdu, "_rep", rep, "_", acc, ".bam"))

        cmd <- paste0("mv", " ", oldfilename, " ", newfilename)
        message(cmd)
        system(cmd)
  }
}

### All bam files associated with ChIP-seq experiments in A549
all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")

### Download bam GR
report_gr_bam <- make_report_bam(target_name = "NR3C1", all_chip_bam)
accession_list_gr <- unique(report_gr_bam$file_accession)
download_bam_in_AccessionList(accession_list_gr, report_gr_bam)

### Download bam EP300
report_ep300_bam <- make_report_bam(target_name = "EP300", all_chip_bam)
accession_list_ep300 <- unique(report_ep300_bam$file_accession)
download_bam_in_AccessionList(accession_list_ep300, report_ep300_bam)

### Download bam WCE GR
report_gr_wce_bam <- make_report_WCE_bam(report_gr_bam, all_chip_bam)
accession_list_gr_wce <- unique(report_gr_wce_bam$file_accession)
download_bam_in_AccessionList(accession_list_gr_wce, report_gr_wce_bam)

### Download bam WCE GR
report_ep300_wce_bam <- make_report_WCE_bam(report_ep300_bam, all_chip_bam)
accession_list_ep300_wce <- unique(report_ep300_wce_bam$file_accession)
download_bam_in_AccessionList(accession_list_ep300_wce, report_ep300_wce_bam)
