# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(GenomicRanges)
library(org.Hs.eg.db)
library(plotly)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
source("scripts/load_reddy.R")

most_expressed = load_most_expressed_transcripts()
most_expressed_TxDb = load_most_expressed_TxDb()
seqlevelsStyle(most_expressed_TxDb) <- "UCSC"
txdb.hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene

#####
load_cofactor_peaks <- function(cofactors = c("NIPBL", "BRD4", "CDK9", "MED1", "SMC1A")) {
  peaks_dir <- "output/chip-pipeline-GRCh38/peak_call"
  cofactors_peaks <- GRangesList()
  name_cofactors_peaks <- c()
  for (cofactor in cofactors) {
    for (condition in c("CTRL", "DEX")) {
      message("####\t", cofactor, " | ", condition)
      basename <- paste0("A549_", condition, "_", cofactor, "_rep1")
      peaks_path <- file.path(peaks_dir, basename, paste0(basename, "_peaks.narrowPeak.bed"))
      message(peaks_path)
      peaks <- rtracklayer::import(peaks_path)
      message("Number of regions : ", length(peaks))
      cofactors_peaks <- append(cofactors_peaks, GRangesList(peaks))
      name_cofactors_peaks <- c(name_cofactors_peaks, paste0(cofactor, "_", condition))
    }
  }
  names(cofactors_peaks) <- name_cofactors_peaks
  message("#####################################")
  message("Available set of regions: ")
  print(names(cofactors_peaks))
  return(cofactors_peaks)
}

#####
keepStdChr <- function(gr) {
  message("With all chromosomes, including contigs : ", length(gr), " regions")
  stdChr <- paste0("chr", c(seq(1:22), "X", "Y"))
  gr_StdChr <- keepSeqlevels(gr, stdChr[stdChr %in% seqlevels(gr)], pruning.mode = "coarse")
  message("Keeping standard chromosomes : ", length(gr_StdChr), " regions")
  message("\t--> ", length(gr) - length(gr_StdChr), " regions removed")
  return(gr_StdChr)
}

#####
load_cofactor_stdchr_peaks <- function(cofactors = c("NIPBL", "BRD4", "CDK9", "MED1","SMC1A")) {
  peaks_dir <- "output/chip-pipeline-GRCh38/peak_call"
  cofactors_peaks <- GRangesList()
  name_cofactors_peaks <- c()
  for (cofactor in cofactors) {
    for (condition in c("CTRL", "DEX")) {
      message("####\t", cofactor, " | ", condition)
      basename <- paste0("A549_", condition, "_", cofactor, "_rep1")
      peaks_path <- file.path(peaks_dir, basename, paste0(basename, "_peaks.narrowPeak.stdchr.bed"))
      message(" > ", peaks_path)
      peaks <- rtracklayer::import(peaks_path)
      message("   > Number of regions : ", length(peaks))
      cofactors_peaks <- append(cofactors_peaks, GRangesList(peaks))
      name_cofactors_peaks <- c(name_cofactors_peaks, paste0(cofactor, "_", condition))
    }
  }
  names(cofactors_peaks) <- name_cofactors_peaks
  message("#####################################")
  message("Available set of regions: ")
  print(names(cofactors_peaks))
  return(cofactors_peaks)
}

#####
load_diffbind_cofactors_peaks <- function(cofactors = c("NIPBL", "BRD4", "CDK9", "MED1", "SMC1A")) {
  peaks_dir <- "output/chip-pipeline-GRCh38/binding_diff"
  diffbind_cofactors_peaks <- GRangesList()
  name_diffbind_cofactors_peaks <- c()
  for (cofactor in cofactors) {
    message("####\t", cofactor)
    cofactor_folder <- file.path(paste("A549", cofactor, sep = "_"), "output_filters")
    
    # up
    peaks_up_path <- file.path(peaks_dir, cofactor_folder, paste("A549_DEX", cofactor, "rep1_peaks.narrowPeak_M_above_0.5_biased_peaks.bed", sep = "_"))
    message(" > ",peaks_up_path)
    peaks_up <- rtracklayer::import(peaks_up_path)
    message("   > Number of regions : ", length(peaks_up))
      
    # down
    peaks_down_path <- file.path(peaks_dir, cofactor_folder, paste("A549_CTRL", cofactor, "rep1_peaks.narrowPeak_M_below_-0.5_biased_peaks.bed", sep = "_"))
    message(" > ",peaks_down_path)
    peaks_down <- rtracklayer::import(peaks_down_path)
    message("   > Number of regions : ", length(peaks_down))
    
    # unbiased
    peaks_unbiased_path <- file.path(peaks_dir, cofactor_folder, paste("A549", cofactor, "unbiased_peaks.bed", sep = "_"))
    message(" > ",peaks_unbiased_path)
    peaks_unbiased <- rtracklayer::import(peaks_unbiased_path)
    message("   > Number of regions : ", length(peaks_unbiased))
    
    # gather everything
    diffbind_cofactors_peaks <- append(diffbind_cofactors_peaks, GRangesList(peaks_up, peaks_down, peaks_unbiased))
    name_diffbind_cofactors_peaks <- c(name_diffbind_cofactors_peaks, c(paste(cofactor, c("UP", "DOWN", "UNBIASED"), sep = "_")))
  }
  names(diffbind_cofactors_peaks) <- name_diffbind_cofactors_peaks
  message("#####################################")
  message("Available set of regions: ")
  print(names(diffbind_cofactors_peaks))
  return(diffbind_cofactors_peaks)
}

#####
annotatePeaks <- function(gr, output = "df", tss = 3000, TxDb = most_expressed_TxDb) {
  # difference between txdb and most expressed txdb???
  # output can be: "df" or "anno"
  # TxDb can be: most_expressed_TxDb or txdb.hg38
  gr_anno <- ChIPseeker::annotatePeak(gr, tssRegion = c(-tss, tss), TxDb=most_expressed_TxDb, annoDb = "org.Hs.eg.db")
  if (output == "anno") {
    message("Return a csAnno object")
    return(gr_anno)
  } else if (output == "df") {
    gr_anno_df <- as.data.frame(gr_anno)
    gr_anno_df$Annot <- gr_anno_df$annotation
    gr_anno_df$Annot <- gsub(" \\(.*\\)", "", gr_anno_df$Annot)
    gr_anno_df$Annot <- as.factor(gr_anno_df$Annot)
    gr_anno_df$Symbol <- mapIds(org.Hs.eg.db, gr_anno_df$geneId, "SYMBOL", "ENSEMBL")
    gr_anno_df$Coordinates <- paste0(gr_anno_df$seqnames, ":", gr_anno_df$start, "-", gr_anno_df$end)
    message("Return a data.frame object")
    return(gr_anno_df)
  }
}

#####
plotAnnotation <- function(anno_df) {
  data <- as.data.frame(table(anno_df$Annot))
  colnames(data) <- c("Annotation", "Freq")
  message("Number of regions: ", nrow(anno_df))
  print(data)
  p <- plot_ly(data, labels = ~Annotation, values = ~Freq,
               textinfo = "label+percent",
               insidetextfont = list(color = "#FFFFFF")) %>%
    add_pie(hole = 0.4)
  return(p)
}

#####
plotVenn <- function(gr_list, labels = TRUE) {
  overlaps <- build_intersect(gr_list)
  mat <- overlaps$Matrix > 0
  fit <- euler(mat, shape = "ellipse")
  p <- plot(fit, quantities = TRUE, labels = labels,
            fills = list(fill = c("#FFB027", "#2B70AB", "#CD3301", "#449B2B")),
            edges = list(alpha = 0))
  return(p)
}

##### load POLR2A peaks from Myers Lab (EtOH and DEX)
# nb: import_files_into_grl function from ef.utils
load_POLR2A_peaks_Myers <- function() {
  peaks_dir <- "output/chip-POLR2A-Myers/peaks"
  ctrl <- "A549_CTRL_POLR2A_ENCFF664KTN.bed"
  dex <- "A549_DEX_POLR2A_ENCFF915LKZ.bed"
  
  bed_path <- paste(peaks_dir, c(ctrl, dex), sep="/")
  names(bed_path) <- c("POLR2A_CTRL", "POLR2A_DEX")
  
  peaks <- import_files_into_grl(bed_path,
                                 file.format = "narrow",
                                 file.ext = "")
  
  message("#####################################")
  message("Available set of regions: ")
  print(names(peaks))
  return(peaks)
}

##### load POLR2A diffbind peaks from Myers Lab (EtOH and DEX)
# nb: import_files_into_grl function from ef.utils
load_diffbind_POLR2A_peaks_Myers <- function() {
  peaks_dir <- "output/chip-POLR2A-Myers/peaks"
  up_fdr <- "POLR2A_up_peaks_fdr.bed"
  down_fdr <- "POLR2A_down_peaks_fdr.bed"
  up_pval <- "POLR2A_up_peaks_pval.bed"
  down_pval <- "POLR2A_down_peaks_pval.bed"
  
  bed_path <- paste(peaks_dir, c(up_fdr, down_fdr, up_pval, down_pval), sep="/")
  names(bed_path) <- c("POLR2A_UP_FDR", "POLR2A_DOWN_FDR", "POLR2A_UP_PVAL", "POLR2A_DOWN_PVAL")
  
  peaks <- lapply(bed_path, rtracklayer::import)
  
  message("#####################################")
  message("Available set of regions: ")
  print(names(peaks))
  return(peaks)
}

##### make UpSet plot
# 
displayUpSet <- function(combMat, threshold = 1, customSetOrder = FALSE) {
  combMat <- combMat[comb_size(combMat) >= threshold]
  annot_top <- HeatmapAnnotation("Intersection\nsize" = anno_barplot(comb_size(combMat), 
                                                                     border = FALSE,
                                                                     gp = gpar(fill = "black"),
                                                                     height = unit(3, "cm")), 
                                 "Size" = anno_text(comb_size(combMat),
                                                    rot = 0,
                                                    just = "center",
                                                    location = 0.25),
                                 annotation_name_side = "left", annotation_name_rot = 0)
  annot_right <- rowAnnotation("Set size" = anno_barplot(set_size(combMat), 
                                                         border = FALSE, 
                                                         gp = gpar(fill = "black"), 
                                                         width = unit(2, "cm")),
                               "Size" = anno_text(set_size(combMat)))
  
  if (customSetOrder != FALSE) {
    UpSet(combMat, top_annotation = annot_top, right_annotation = annot_right,
          set_order = customSetOrder)
  } else {
    UpSet(combMat, top_annotation = annot_top, right_annotation = annot_right)
  }
}