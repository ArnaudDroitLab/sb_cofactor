# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(tidyverse)
library(knitr)
library(DiffBind)
library(ComplexHeatmap)
source("scripts/chris/metagene2_Reddy.utils.R")
source("scripts/ckn_utils.R")

#####
build_sSheet <- function(target, bam_folder, bed_folder, reps = "123") {
  # BAM
  bam_pattern <- paste0("^", target, "_([0-9]+.*)_rep([", reps, "])_(.*\\.bam$)")
  bam_files <- list.files(path = bam_folder, pattern = bam_pattern, full.names = TRUE)
  nb_bam <- length(bam_files)
  
  SampleID <- basename(bam_files) %>% gsub(pattern = bam_pattern, replacement = paste0(target, "_", "\\1", "_rep\\2")) %>%
    gsub(pattern = "minute", replacement = "m") %>% gsub(pattern = "hour", replacement = "h")
  Timepoint <- gsub(pattern = paste0("^", target, "_([0-9]+(m|h))_rep(.)"), replacement = "\\1", SampleID)
  Replicate <- gsub(pattern = paste0("^", target, "_([0-9]+(m|h))_rep(.)"), replacement = "\\3", SampleID)
  Treatment <- Timepoint
  Treatment_bis <- ifelse(Timepoint == "0m" | Timepoint == "0h", "EtOH", "DEX")
  Tissue <- "A549"
  Antibody <- target
  bamReads <- bam_files
  
  # BED
  bed_pattern <- paste0("^", target, "_([0-9]+.*)_rep([", reps, "])_(.*\\.bed.gz$)")
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
perform_diffbind <- function(target, sSheet, tp1, tp2, reps = "123", output_dir) {
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
  df_filename_pval <- paste0("diffbind_", target, "_", contrast_name, "_reps", reps, "_pval.txt")
  output_path <- file.path(output_dir, df_filename_pval)
  write.table(report_pval, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE)
  message("     > Differential binding saved in ", output_path)

  report <- dba.report(analyze, bCounts = T)
  df_filename <- paste0("diffbind_", target, "_", contrast_name, "_reps", reps, "_fdr.txt")
  output_path <- file.path(output_dir, df_filename)
  write.table(report, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE)
  message("     > Differential binding saved in ", output_path)
}

#####
open_diffBind <- function(target, tp1, tp2, reps = "123", pval = FALSE, output_dir) {
  message("#######################")
  message("### ", tp2, " VS ", tp1)
  message("#######################")
  contrast_name <- paste0(tp2, "VS", tp1)
  filename <- paste0("diffbind_", target, "_", contrast_name, "_reps", reps)
  
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
# sSheet_GR <- build_sSheet("NR3C1", bam_folder, bed_folder, reps = "12")
sSheet_EP300 <- build_sSheet("EP300", bam_folder, bed_folder, reps = "123")

timepoint <- c("0m", "5m", "10m", "15m", "20m", "25m")
ltp <- length(timepoint)
for (i in 1:(ltp-1)) {
  for (j in (i+1):ltp) {
    tp1 <- timepoint[i]
    tp2 <- timepoint[j]
    # perform_diffbind("GR", sSheet_GR, tp1, tp2, reps = "12", output_dir = "output/analyses/GR_diffbind")
    perform_diffbind("EP300", sSheet_EP300, tp1, tp2, reps = "123", output_dir = "output/analyses/EP300_diffbind")
    # open_diffBind("EP300", tp1, tp2, pval = TRUE, reps = "123", output_dir = "output/analyses/EP300_diffbind")
  }
}

## perform diffbind for EP300 0h 1h
# perform_diffbind("EP300", sSheet_EP300, "0h", "1h", reps = "123", output_dir = "output/analyses/EP300_diffbind")
EP300_DOWN <- open_diffBind("EP300", tp1 = "0h", tp2 = "1h", pval = TRUE, reps = "123", output_dir = "output/analyses/EP300_diffbind") %>%
  dplyr::filter(Fold < 0) %>% makeGRangesFromDataFrame

## perform diffbind for EP300 0h 0m
# perform_diffbind("EP300", sSheet_EP300, "0h", "0m", reps = "123", output_dir = "output/analyses/EP300_diffbind")
open_diffBind("EP300", "0h", "0m", pval = TRUE, reps = "123", output_dir = "output/analyses/EP300_diffbind")

# explore EP300 diffbind 
EP300_diffbind_downreg <- list()
for (i in 1:(ltp-1)) {
  for (j in (i+1):ltp) {
    for (TF in c(TRUE)) {
      tp1 <- timepoint[i]
      tp2 <- timepoint[j]
      contrast_name <- paste0(tp2, "VS", tp1)
      report <- open_diffBind("EP300", tp1, tp2, pval = TF, reps = "123", output_dir = "output/analyses/EP300_diffbind")
      
      if (!is.null(report)) {
        upreg <- report %>% dplyr::filter(Fold > 0)
        downreg <- report %>% dplyr::filter(Fold < 0)
        if (nrow(downreg) != 0) {
          EP300_diffbind_downreg[[paste0(contrast_name, "_downreg")]] <- makeGRangesFromDataFrame(downreg)
        }
      }
    }
  }
}

names(EP300_diffbind_downreg)
sapply(EP300_diffbind_downreg, length)

# Overlaps Cofactors_down and GR_DOWN
cofactors <- load_diffbind_cofactors_peaks()
cofactors <- cofactors[grep("_DOWN", names(cofactors))]
cofactors <- append(cofactors, GRangesList("EP300_DOWN" = EP300_DOWN))
names(cofactors)

set_EP300_list_with_cofactors <- c(cofactors, EP300_diffbind_downreg)
names(set_EP300_list_with_cofactors)
sapply(set_EP300_list_with_cofactors, length)

inter_cofactors <- build_intersect(set_EP300_list_with_cofactors)
matrix_cofactors <- inter_cofactors$Matrix
sum(matrix_cofactors > 1)
matrix_cofactors[matrix_cofactors > 1] <- 1
sum(matrix_cofactors > 1)
colnames(matrix_cofactors)

m4 <- make_comb_mat(matrix_cofactors, remove_empty_comb_set = TRUE)
m4 <- m4[comb_size(m4) >= 50]
UpSet(m4)
comb_size(m4)

annot_top <- HeatmapAnnotation("Intersection\nsize" = anno_barplot(comb_size(m4), 
                                                                   border = FALSE, gp = gpar(fill = "black"), height = unit(3, "cm")), 
                               annotation_name_side = "left", annotation_name_rot = 0,
                               "Size" = anno_text(comb_size(m4), rot = 0, just = "center", location = 0.25))
annot_right <- rowAnnotation("Set size" = anno_barplot(set_size(m4), 
                                                       border = FALSE, 
                                                       gp = gpar(fill = "black"), 
                                                       width = unit(2, "cm")),
                             "Size" = anno_text(set_size(m4))
)

UpSet(m4, top_annotation = annot_top, right_annotation = annot_right)

# Let's extract the intersection set
idToName <- function(id, set_names) {
  name <- c()
  for (i in 1:nchar(id)) {
    if (substring(id, i, i) == "1") {
      good_set <- strsplit(set_names[i], "_")[[1]][1]
      name <- paste(name, good_set, sep = "+")
    }
  }
  return(substring(name, 2))
}

cs <- comb_size(m4)
de <- comb_degree(m4)
set_GR_list <- list()
for (id in names(cs)){
  i <- extract_comb(m4, id)
  setname <- idToName(id, colnames(matrix_cofactors))
  message("##### ", id, " | degree ", de[id], " | ", setname)
  regions <- inter_cofactors$Regions[i]
  
  set_GR_list[[setname]] <- regions
}

names(set_GR_list)
sapply(set_GR_list, length)

# for (i in c(35, 30, 36, 28, 38, 33, 31, 27, 21)) {
#   contrast <- names(set_GR_list)[i]
#   message("##### ", contrast)
#   df_set_GR_list <- make_df_metagene_Reddy(chip_target = c("GR", "EP300"), peaks = set_GR_list[[i]], merge_replicates = TRUE, reps = "12")
#   title_group <- paste(contrast, paste(length(set_GR_list[[i]]), "regions"), sep = " | ")
#   p <- plot_metagene_Reddy(df_set_GR_list, title = title_group)
#   saveMetagene(metagene_plot = p,
#                output_dir = "output/analyses/GR_diffbind_downreg_setlist_2",
#                output_file = paste("GR_diffbind_downreg_GR_EP300", contrast, "reps12", sep = "_"),
#                width = 20, height = 9)
# }

t <- set_GR_list[["NIPBL+BRD4+CDK9+MED1+SMC1A+25mVS10m"]] %>% as.data.frame
t$seqnames <- as.character(t$seqnames)
t <- t %>% arrange(seqnames) %>%
  mutate(coord = paste0(seqnames, ":", start, "-", end),
         index = 1:length(coord))
kable(t)
