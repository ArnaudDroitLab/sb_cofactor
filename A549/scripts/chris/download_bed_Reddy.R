# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(ENCODExplorer)
library(tidyverse)
library(knitr)

###
make_report_bed <- function(target_name, all_chip_bed) {
  bed <- all_chip_bed %>% filter(target == target_name, assembly == "GRCh38", lab == "Tim Reddy, Duke")
  report_bed <- bed %>% dplyr::select(accession, file_accession, file_size, submitted_by, target, treatment_duration, treatment_duration_unit, biological_replicates, controls)
  report_bed$treatment_duration[is.na(report_bed$treatment_duration)] <- 0
  report_bed$treatment_duration_unit[is.na(report_bed$treatment_duration_unit) &
                                     report_bed$submitted_by == "Ian McDowell"] <- "hour"
  report_bed$treatment_duration_unit[is.na(report_bed$treatment_duration_unit) &
                                       report_bed$submitted_by == "Alejandro Barrera"] <- "minute"
  report_bed <- report_bed %>% arrange(desc(treatment_duration_unit), treatment_duration, submitted_by, biological_replicates)
  print(kable(report_bed))
  return(report_bed)
}

###
download_bed_from_ENCODE <- function(encode_accession) {
  chip_bed_dir = "input/ENCODE/A549/GRCh38/chip-seq/narrow"
  dir.create(chip_bed_dir, recursive=TRUE, showWarnings=FALSE)
  downloaded_file = ENCODExplorer::downloadEncode(file_acc = encode_accession, dir=chip_bed_dir, force=FALSE)
}

###
download_bed_in_AccessionList <- function(accession_list, report_bed) {
  cpt <- 0
  cptmax <- length(accession_list)
  for (acc in accession_list) {
    cpt <- cpt + 1
    message("#########################################################################")
    target <- report_bed %>% filter(file_accession == acc) %>% select(target)
    rep <- report_bed %>% filter(file_accession == acc) %>% select(biological_replicates)
    size <- report_bed %>% filter(file_accession == acc) %>% select(file_size)
    td <- report_bed %>% filter(file_accession == acc) %>% select(treatment_duration) # td = treatment_duration
    tdu <- report_bed %>% filter(file_accession == acc) %>% select(treatment_duration_unit) # tdu = treatment_duration_unit
    
    message("####### ", acc, "\t", target, "\trep", rep, "\t",
            td, " ", tdu, "\t",
            "| ", size, "\t",
            cpt, "\\", cptmax)
    download_bed_from_ENCODE(acc)
    
    oldfilename <- file.path(chip_bed_dir, paste0(acc, ".bed.gz"))
    newfilename <- file.path(chip_bed_dir, paste0(target, "_", td, tdu, "_rep", rep, "_", acc, ".bed.gz"))
    
    cmd <- paste0("mv", " ", oldfilename, " ", newfilename)
    message(cmd)
    system(cmd)
  }
}

###
download_Reddy_Peaks <- function(protein, WCE = FALSE) {
  report_bed <- make_report_bed(target_name = protein, all_chip_bed)
  if (WCE == TRUE) {
    report_wce_bed <- make_report_WCE_bed(report_bed, all_chip_bed)
    accession_list_wce <- unique(report_wce_bed$file_accession)
    download_bed_in_AccessionList(accession_list_wce, report_wce_bed)
  } else {
    accession_list <- unique(report_bed$file_accession)
    message("Number of files : ", length(accession_list))
    download_bed_in_AccessionList(accession_list, report_bed)
  }
}

### Output directory
chip_bed_dir = "input/ENCODE/A549/GRCh38/chip-seq/narrow"

### All bed files associated with ChIP-seq experiments in A549
all_chip_bed <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bed", assay="ChIP-seq")

# download_Reddy_Peaks("NR3C1")
download_Reddy_Peaks("EP300")
