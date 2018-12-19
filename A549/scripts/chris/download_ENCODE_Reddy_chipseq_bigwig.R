# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(dplyr)
library(knitr)
library(ENCODExplorer)

###
download_bigwig_from_ENCODE <- function(encode_accession) {
  chip_bw_dir = "input/ENCODE/A549/GRCh38/chip-seq/bigWig"
  dir.create(chip_bw_dir, recursive=TRUE, showWarnings=FALSE)
  downloaded_file = ENCODExplorer::downloadEncode(file_acc = encode_accession, dir=chip_bw_dir, force=FALSE)
}

#####
all_chip_bw <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bigWig", assay="ChIP-seq")
chip_bw_dir = "input/ENCODE/A549/GRCh38/chip-seq/bigWig"

##### DEX Condition (1 hour)
df_bw <- all_chip_bw %>% filter(assembly == "GRCh38", lab == "Tim Reddy, Duke",
                                treatment == "dexamethasone", treatment_duration == "1")

report_bw <- df_bw %>% select(accession, file_accession, file_format, target, treatment, treatment_duration, treatment_duration_unit, biological_replicates, controls)
kable(report_bw)

accession_list <- unique(df_bw$file_accession)
cptmax <- length(accession_list)

cpt <- 0
for (acc_dex in accession_list) {
  cpt <- cpt + 1
  message("#########################################################################")
  target <- report_bw %>% filter(file_accession == acc_dex) %>% select(target)
  rep <- report_bw %>% filter(file_accession == acc_dex) %>% select(biological_replicates)
  
  message("####### ", acc_dex, "\t", target, "\trep", rep, "\t", cpt, "\\", cptmax)
  # download_bigwig_from_ENCODE(acc_dex)
  
  oldfilename <- file.path(chip_bw_dir, paste0(acc_dex, ".bigWig"))
  newfilename <- file.path(chip_bw_dir, paste0(target, "_dex1h", "_rep", rep, ".bigWig"))
  
  cmd <- paste0("mv", " ", oldfilename, " ", newfilename)
  message(cmd)
  # system(cmd)
}

##### EtOH condition
df_bw_ctrl <- all_chip_bw %>% filter(assembly == "GRCh38", lab == "Tim Reddy, Duke")
df_bw_ctrl <- df_bw_ctrl[is.na(df_bw_ctrl$treatment), ]

report_bw_ctrl <- df_bw_ctrl %>% select(accession, file_accession, file_format, target, treatment, treatment_duration, treatment_duration_unit, biological_replicates, controls)
kable(report_bw_ctrl)

accession_list_ctrl <- unique(report_bw_ctrl$file_accession)
cptmax <- length(accession_list_ctrl)

chip_bw_dir = "input/ENCODE/A549/GRCh38/chip-seq/bigWig"

cpt <- 0
for (acc_ctrl in accession_list_ctrl) {
  cpt <- cpt + 1
  message("#########################################################################")
  target <- report_bw_ctrl %>% filter(file_accession == acc_ctrl) %>% select(target)
  rep <- report_bw_ctrl %>% filter(file_accession == acc_ctrl) %>% select(biological_replicates)
  
  message("####### ", acc_ctrl, "\t", target, "\trep", rep, "\t", cpt, "\\", cptmax)
  #download_bigwig_from_ENCODE(acc_ctrl)
  
  oldfilename <- file.path(chip_bw_dir, paste0(acc_ctrl, ".bigWig"))
  newfilename <- file.path(chip_bw_dir, paste0(target, "_etoh", "_rep", rep, ".bigWig"))
  
  cmd <- paste0("mv", " ", oldfilename, " ", newfilename)
  message(cmd)
  # system(cmd)
}
