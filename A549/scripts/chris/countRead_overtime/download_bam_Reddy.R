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
      size <- report_bam %>% filter(file_accession == acc) %>% select(file_size)
      td <- report_bam %>% filter(file_accession == acc) %>% select(treatment_duration) # td = treatment_duration
      tdu <- report_bam %>% filter(file_accession == acc) %>% select(treatment_duration_unit) # tdu = treatment_duration_unit
      
        message("####### ", acc, "\t", target, "\trep", rep, "\t",
                td, " ", tdu, "\t",
                "| ", size, "\t",
                cpt, "\\", cptmax)
        download_bam_from_ENCODE(acc)

        oldfilename <- file.path(chip_bam_dir, paste0(acc, ".bam"))
        newfilename <- file.path(chip_bam_dir, paste0(target, "_", td, tdu, "_rep", rep, "_", acc, ".bam"))

        cmd <- paste0("mv", " ", oldfilename, " ", newfilename)
        message(cmd)
        system(cmd)
  }
}

download_Reddy_ChIP <- function(protein, WCE = FALSE) {
  report_bam <- make_report_bam(target_name = protein, all_chip_bam)
  if (WCE == TRUE) {
    report_wce_bam <- make_report_WCE_bam(report_bam, all_chip_bam)
    accession_list_wce <- unique(report_wce_bam$file_accession)
    download_bam_in_AccessionList(accession_list_wce, report_wce_bam)
  } else {
    accession_list <- unique(report_bam$file_accession)
    download_bam_in_AccessionList(accession_list, report_bam)
  }
}

### All bam files associated with ChIP-seq experiments in A549
all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")

### Download Bam
download_Reddy_ChIP("NR3C1")
download_Reddy_ChIP("EP300")
download_Reddy_ChIP("NR3C1", WCE = TRUE)
download_Reddy_ChIP("EP300", WCE = TRUE)
download_Reddy_ChIP("CTCF")
download_Reddy_ChIP("SMC3")
download_Reddy_ChIP("RAD21")
download_Reddy_ChIP("FOSL2")
download_Reddy_ChIP("BCL3")
download_Reddy_ChIP("JUN")
download_Reddy_ChIP("JUNB")
download_Reddy_ChIP("HES2")
download_Reddy_ChIP("CEPBP")