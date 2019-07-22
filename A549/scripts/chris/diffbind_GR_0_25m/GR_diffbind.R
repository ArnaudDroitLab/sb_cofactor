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
  Timepoint <- gsub(pattern = paste0("^", target, "_([0-9]+m)_rep(.)"), replacement = "\\1", SampleID)
  Replicate <- gsub(pattern = paste0("^", target, "_([0-9]+m)_rep(.)"), replacement = "\\2", SampleID)
  Treatment <- Timepoint
  Treatment_bis <- ifelse(Timepoint == "0m", "EtOH", "DEX")
  Tissue <- "A549"
  Antibody <- target
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
perform_diffbind <- function(target, sSheet, tp1, tp2, output_dir) {
  message("#######################")
  message("### ", tp2, " vs ", tp1)
  message("#######################")
  contrast_name <- paste0(tp2, "VS", tp1)
  sSheet_filtered <- sSheet %>% dplyr::filter(Timepoint %in% c(tp1, tp2))
  
  message("  > Reading...")
  dba <- dba(sampleSheet = sSheet_filtered)
  print(dba)
  
  message("  > Counting...")
  count <- dba.count(dba)
  print(count)
  
  message("  > Define category...")
  category = DBA_TREATMENT
  print(category)
  
  message("  > Define contrast...")
  contrast <- dba.contrast(count, categories = category, minMembers = 2,
                           group1 = count$mask[[tp2]], group2 = count$mask[[tp1]],
                           name1 = tp2, name2 = tp1)
  print(contrast)
  
  message("  > Analyzing...")
  analyze <- dba.analyze(contrast)
  print(analyze)
  
  message("  > Reporting...")
  report_pval <- dba.report(analyze, bCounts = T, bUsePval = TRUE)
  df_filename_pval <- paste0("diffbind_", target, "_", contrast_name, "_pval.txt")
  output_path <- file.path(output_dir, df_filename_pval)
  write.table(report_pval, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE)
  message("     > Differential binding saved in ", output_path)

  report <- dba.report(analyze, bCounts = T)
  df_filename <- paste0("diffbind_", target, "_", contrast_name, "_fdr.txt")
  output_path <- file.path(output_dir, df_filename)
  write.table(report, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE)
  message("     > Differential binding saved in ", output_path)
}

#####
open_diffBind <- function(target, tp1, tp2, pval = FALSE, output_dir) {
  message("#######################")
  message("### ", tp2, " VS ", tp1)
  message("#######################")
  contrast_name <- paste0(tp2, "VS", tp1)
  filename <- paste0("diffbind_", target, "_", contrast_name)
  
  if (pval) {
    df_filename <- paste(filename, "pval.txt", sep = "_")
  } else {
    df_filename <- paste(filename, "fdr.txt", sep = "_")
  }
  
  output_path <- file.path(output_dir, df_filename)
  message("  # >>> ", output_path)
  
  report <- try(read.delim(output_path, sep = "\t" , header = TRUE))
  if (!inherits(report, 'try-error')) {
    report$Coord <- paste0(report$seqnames, ":", report$start, "-", report$end)
    report_up <- report %>% filter(Fold > 0)
    report_down <- report %>% filter(Fold < 0)
    
    if (pval) {
      message("  ### PVAL")
    } else {
      message("  ### FDR")
    }
    
    message("      Number of differential regions : ", nrow(report))
    message("         Increased signals : ", nrow(report_up))
    message("         Decreased signals : ", nrow(report_down))  
    return(report)
  } else {
    message("  ### Empty file")
  }
}

#####
bam_folder <- "input/ENCODE/A549/GRCh38/chip-seq/bam"
bed_folder <- "input/ENCODE/A549/GRCh38/chip-seq/narrow"
sSheet_GR <- build_sSheet("NR3C1", bam_folder, bed_folder)
sSheet_EP300 <- build_sSheet("EP300", bam_folder, bed_folder)

timepoint <- c("0m", "5m", "10m", "15m", "20m", "25m")
ltp <- length(timepoint)
for (i in 1:(ltp-1)) {
  for (j in (i+1):ltp) {
    tp1 <- timepoint[i]
    tp2 <- timepoint[j]
    # perform_diffbind("GR", sSheet_GR, tp1, tp2, output_dir = "output/analyses/GR_diffbind")
  }
}
  
for (i in 1:(ltp-1)) {
  for (j in (i+1):ltp) {
    for (TF in c(TRUE)) {
      tp1 <- timepoint[i]
      tp2 <- timepoint[j]
      report <- open_diffBind("GR", tp1, tp2, pval = TF, output_dir = "output/analyses/GR_diffbind")
    }
  }
}
