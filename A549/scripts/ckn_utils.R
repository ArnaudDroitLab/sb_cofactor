# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(GenomicRanges)
library(org.Hs.eg.db)
source("scripts/load_reddy.R")

most_expressed = load_most_expressed_transcripts()
most_expressed_TxDb = load_most_expressed_TxDb()
seqlevelsStyle(most_expressed_TxDb) <- "UCSC"

#####
load_cofactor_peaks <- function(cofactors = c("NIPBL", "BRD4", "CDK9", "MED1","SMC1A")) {
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

annotatePeaks <- function(gr, output = "df") {
  # difference between txdb and most expressed txdb???
  gr_anno <- ChIPseeker::annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb=most_expressed_TxDb, annoDb = "org.Hs.eg.db")
  if (output == "anno") {
    message("Return a csAnno object")
    return(gr_anno)
  } else if (output == "df") {
    gr_anno_df <- as.data.frame(gr_anno)
    gr_anno_df$Annot <- gr_anno_df$annotation
    gr_anno_df$Annot <- gsub(" \\(.*\\)", "", gr_anno_df$Annot)
    gr_anno_df$Annot <- as.factor(gr_anno_df$Annot)
    # gr_anno_df$Symbol <- mapIds(org.Hs.eg.db, gr_anno_df$geneId, "SYMBOL", "ENSEMBL")
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
plotVenn <- function(gr_list) {
  overlaps <- build_intersect(gr_list)
  mat <- overlaps$Matrix > 0
  fit <- euler(mat, shape = "ellipse")
  p <- plot(fit, quantities = TRUE,
            fills = list(fill = c("#FFB027", "#2B70AB", "#CD3301", "#449B2B")),
            edges = list(alpha = 0))
  return(p)
}