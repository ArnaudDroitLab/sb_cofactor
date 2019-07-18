# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(tidyverse)
library(knitr)

build_sSheet <- function(target, bam_folder) {
  # BAM
  bam_pattern <- paste0("^", target, "_([0-9]+minute)_rep(.)_(.*\\.bam$)")
  bam_files <- list.files(path = bam_folder, pattern = bam_pattern, full.names = TRUE)
  nb_bam <- length(bam_files)
  
  SampleID <- basename(bam_files) %>% gsub(pattern = bam_pattern, replacement = paste0(target, "_", "\\1", "_rep\\2")) %>%
    gsub(pattern = "minute", replacement = "m")
  Timepoint <- SampleID %>% gsub(pattern = paste0("^", target, "_([0-9]+m)_rep(.)"), replacement = "\\1")
  Replicate <- SampleID %>% gsub(pattern = paste0("^", target, "_([0-9]+m)_rep(.)"), replacement = "\\2")
  Treatment <- ifelse(Timepoint == "0m", "EtOH", "DEX")
  Tissue <- rep("A549", nb_bam)
  Antibody <- rep(target, nb_bam)
  bamReads <- bam_files
  
  # BED
  bed_pattern <- paste0("^", target, "_([0-9]+minute)_rep(.)_(.*\\.bed$)")
  bed_files <- list.files(path = bed_folder, pattern = bed_pattern, full.names = TRUE)
  
}

bam_folder <- "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam"
bed_folder <- "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/narrow"
sSheet_GR <- build_sSheet("NR3C1")

ControlID <- # check with ENCODEExplorer
bamControl <- # check with ENCODEExplorer
Peaks <- # use bed_pattern <- paste0("(", target, ").*\\.bed$")
PeakCaller <- rep("bed", nb_bam)

sSheet <- data.frame(SampleID, Tissue, Antibody, Treatment,
                     Timepoint, Replicate, bamReads, 
                     ControlID, bamControl,
                     Peaks, PeakCaller)
# colonnes:
# SampleID
# Tissue
# Antibody
# Treatment
# Replicate
# bamReads
# ControlID
# bamControl
# Peaks
# PeakCaller

perform_diffbind <- function(sSheet, tp1, tp2) {
  message("### ", tp2, "VS", tp1)
  contrast_name <- paste0(tp2, "VS", tp1)
  sSheet_filtered <- sSheet %>% dplyr::filter(Timepoint %in% c(tp1, tp2))
  
  dba <- dba(sampleSheet = sSheet_filtered)
  
  message("Counting...")
  count <- dba.count(dba)
  
  category = Timepoint
  
  contrast <- dba.contrast(count, categories = category, minMembers = 2)
  
  analyze <- dba.analyze(contrast)
  
  report_pval <- dba.report(analyze, bCounts = T, bUsePval = pval)
  df_filename_pval <- paste0("diffbind_", contrast_name, "_pval.txt")
  output_path <- file.path("output/analyses/GR_diffbind", df_filename_pval)
  write.table(report, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE)
  message(" > Differential binding saved in", output_path)
  
  report <- dba.report(analyze, bCounts = T)
  df_filename <- paste0("diffbind_", contrast_name, "_fdr.txt")
  output_path <- file.path("output/analyses/GR_diffbind", df_filename)
  write.table(report, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE)
  message(" > Differential binding saved in", output_path)
}

timepoint <- c("10m", "15m", "20m", "25m")
for (tp in timepoint) {
    perform_diffbind(sSheet_GR, tp, "5m")
}

stats_diffBind <- function(tp1, tp2, pval = FALSE) {
  message("### ", tp2, "VS", tp1)
  contrast_name <- paste0(tp2, "VS", tp1)
  filename <- paste0("diffbind_", contrast_name)
  
  if (pval) {
    df_filename <- paste(filename, "pval.txt", sep = "_")
  } else {
    df_filename <- paste(filename, "fdr.txt", sep = "_")
  }
  
  output_path <- file.path("output/analyses/GR_diffbind", df_filename)
  #   message("   # >>> ", output_path)
  #   
  #   annodf <- read.delim(output_path, sep = "\t" , header = TRUE)
  #   annodf$geneId <- as.character(annodf$geneId)
  #   annodf$SYMBOL <- as.character(annodf$SYMBOL)
  #   
  #   annodf_up <- annodf %>% filter(Fold > 0)
  #   annodf_down <- annodf %>% filter(Fold < 0)
  #   
  #   message("    Number of differential regions : ", nrow(annodf))
  #   message("       Increased signal : ", nrow(annodf_up), " (including ", nrow(annodf_up %>% filter(Annot == "Promoter")), " regions at promoters)")
  #   message("       Decreased signal : ", nrow(annodf_down), " (including ", nrow(annodf_down %>% filter(Annot == "Promoter")), " regions at promoters)")  
}
# perform_diffbind <- function(pol, cst, effect, peak, pval) {
#   message("######\t", pol, " | ", cst, " | ", effect, " effect | ", peak, "Peak")
#   
#   basename <- paste0(pol, "_", cst, "_", effect, "_effect_", peak)
#   
#   filename = paste0("sSheet_", basename, "Peak.csv")
#   message("   # >>> ", filename)
#   
#   sSheet <- read.table(file.path("scripts/chris/diffbind_pol2/sampleSheet", filename), sep= "," , header = TRUE)
#   
#   dba <- dba(sampleSheet = sSheet)
#   
#   message("Counting...")
#   count <- dba.count(dba)
#   
#   if (effect == "DEX") {
#     category = DBA_TREATMENT
#   } else if (effect == "shNIPBL") {
#     category = DBA_CONDITION
#   }
#   message("effect : ", effect , " | ", category)
#   
#   contrast <- dba.contrast(count, categories = category, minMembers = 2)
#   
#   analyze <- dba.analyze(contrast)
#   
#   if (pval) {
#     report <- dba.report(analyze, bCounts = T, bUsePval = pval)
#   } else {
#     report <- dba.report(analyze, bCounts = T)
#   }
#   
#   annodf <- annotatePeaks(report, output = "df", tss = 5000)
#   # print(nrow(annodf)); print(sum(annodf$Annot == "Promoter")); print(sort(unique(annodf$SYMBOL)))
#   
#   if (pval) {
#     annodf_filename <- paste0("diffbind_", basename, "_pval.txt")
#   } else {
#     annodf_filename <- paste0("diffbind_", basename, "_fdr.txt")
#   }
#   
#   output_path <- file.path("output/analyses/diffbind_pol2", annodf_filename)
#   write.table(annodf, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE)
# }
# 
# ###
# stats_diffbind <- function(pol, cst, effect, peak, pval, geneList = FALSE) {
#   message("######\t", pol, " | ", cst, " | ", effect, " effect | ", peak, "Peak")
#   
#   basename <- paste0(pol, "_", cst, "_", effect, "_effect_", peak)
#   
#   if (pval) {
#     annodf_filename <- paste0("diffbind_", basename, "_pval.txt")
#   } else {
#     annodf_filename <- paste0("diffbind_", basename, "_fdr.txt")
#   }
#   
#   output_path <- file.path("output/analyses/diffbind_pol2", annodf_filename)
#   message("   # >>> ", output_path)
#   
#   annodf <- read.delim(output_path, sep = "\t" , header = TRUE)
#   annodf$geneId <- as.character(annodf$geneId)
#   annodf$SYMBOL <- as.character(annodf$SYMBOL)
#   
#   annodf_up <- annodf %>% filter(Fold > 0)
#   annodf_down <- annodf %>% filter(Fold < 0)
#   
#   message("    Number of differential regions : ", nrow(annodf))
#   message("       Increased signal : ", nrow(annodf_up), " (including ", nrow(annodf_up %>% filter(Annot == "Promoter")), " regions at promoters)")
#   message("       Decreased signal : ", nrow(annodf_down), " (including ", nrow(annodf_down %>% filter(Annot == "Promoter")), " regions at promoters)")
#   
#   if (geneList == "increase") {
#     geneList <- annodf_up %>% filter(Annot == "Promoter") %>% pull(geneId) %>% unique
#     message("   Length of geneList (increase) : ", length(geneList))
#     return(geneList)
#   } else if (geneList == "decrease") {
#     geneList <- annodf_down %>% filter(Annot == "Promoter") %>% pull(geneId) %>% unique
#     message("   Length of geneList (decrease) : ", length(geneList))
#     return(geneList)
#   }
# }