# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(tidyverse)
library(knitr)
library(DiffBind)

#####
build_sSheet <- function(target, bam_folder, bed_folder) {
  # BAM
  bam_pattern <- paste0("^", target, "_([0-9]+minute)_rep(.)_(.*\\.bam$)")
  bam_files <- list.files(path = bam_folder, pattern = bam_pattern, full.names = TRUE)
  nb_bam <- length(bam_files)
  
  SampleID <- basename(bam_files) %>% gsub(pattern = bam_pattern, replacement = paste0(target, "_", "\\1", "_rep\\2")) %>%
    gsub(pattern = "minute", replacement = "m")
  Timepoint <- SampleID %>% gsub(pattern = paste0("^", target, "_([0-9]+m)_rep(.)"), replacement = "\\1")
  Replicate <- SampleID %>% gsub(pattern = paste0("^", target, "_([0-9]+m)_rep(.)"), replacement = "\\2")
  Treatment <- Timepoint
  Treatment_bis <- ifelse(Timepoint == "0m", "EtOH", "DEX")
  Tissue <- rep("A549", nb_bam)
  Antibody <- rep(target, nb_bam)
  bamReads <- bam_files
  
  # BED
  bed_pattern <- paste0("^", target, "_([0-9]+minute)_rep(.)_(.*\\.bed.gz$)")
  bed_files <- list.files(path = bed_folder, pattern = bed_pattern, full.names = TRUE)
  Peaks <- bed_files
  PeakCaller <- rep("bed", nb_bam)
  
  # sSheet
  sSheet <- data.frame(SampleID, Tissue, Antibody,
                       Treatment, Treatment_bis,
                       Timepoint, Replicate, bamReads,
                       Peaks, PeakCaller)
  return(sSheet)
}

#####
perform_diffbind <- function(sSheet, tp1, tp2) {
  message("### ", tp2, " vs ", tp1)
  contrast_name <- paste0(tp2, "VS", tp1)
  sSheet_filtered <- sSheet %>% dplyr::filter(Timepoint %in% c(tp1, tp2))
  
  dba <- dba(sampleSheet = sSheet_filtered)
  
  message("Counting...")
  count <- dba.count(dba)
  print(count)
  
  category = DBA_TREATMENT
  print(category)
  
  contrast <- dba.contrast(count, categories = category, minMembers = 2)
  print(contrast)
  
  message("Analyzing...")
  analyze <- dba.analyze(contrast)
  pritn(analyse)
  
  report_pval <- dba.report(analyze, bCounts = T, bUsePval = TRUE)
  df_filename_pval <- paste0("diffbind_", contrast_name, "_pval.txt")
  output_path <- file.path("output/analyses/GR_diffbind", df_filename_pval)
  write.table(report_pval, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE)
  message(" > Differential binding saved in", output_path)
  
  report <- dba.report(analyze, bCounts = T)
  df_filename <- paste0("diffbind_", contrast_name, "_fdr.txt")
  output_path <- file.path("output/analyses/GR_diffbind", df_filename)
  write.table(report, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE)
  message(" > Differential binding saved in", output_path)
}

#####
open_diffBind <- function(tp1, tp2, pval = FALSE) {
  message("### ", tp2, "VS", tp1)
  contrast_name <- paste0(tp2, "VS", tp1)
  filename <- paste0("diffbind_", contrast_name)
  
  if (pval) {
    df_filename <- paste(filename, "pval.txt", sep = "_")
  } else {
    df_filename <- paste(filename, "fdr.txt", sep = "_")
  }
  
  output_path <- file.path("output/analyses/GR_diffbind", df_filename)
  message("   # >>> ", output_path)
  
  report <- read.delim(output_path, sep = "\t" , header = TRUE)
  report$Coord <- paste0(report$seqnames, ":", report$start, "-", report$end)
  report_up <- report %>% filter(Fold > 0)
  report_down <- report %>% filter(Fold < 0)
  
  message("    Number of differential regions : ", nrow(report))
  message("       Increased signal : ", nrow(report_up))
  message("       Decreased signal : ", nrow(report_down))  
  return(report)
}

#####
bam_folder <- "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam"
bed_folder <- "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/narrow"
sSheet_GR <- build_sSheet("NR3C1", bam_folder, bed_folder)

timepoint <- c("0m", "5m", "10m", "15m", "20m", "25m")
ltp <- length(timepoint)
for (i in 1:(ltp-1)) {
  for (j in (i+1):ltp) {
    tp1 <- timepoint[i]
    tp2 <- timepoint[j]
    perform_diffbind(sSheet_GR, tp1, tp2)
  }
}
  
for (i in 1:(ltp-1)) {
  for (j in (i+1):ltp) {
    for (TP in c(FALSE, TRUE)) {
      tp1 <- timepoint[i]
      tp2 <- timepoint[j]
      open_diffBind(tp1, tp2, pval = TF)
    }
  }
}

#####
open_diffBind("5m", "25m", pval = FALSE) %>% filter(Fold < 0)
open_diffBind("5m", "25m", pval = TRUE)