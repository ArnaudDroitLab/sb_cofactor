# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(ef.utils)
library(plotly)
library(ChIPseeker)
source("scripts/ckn_utils.R")

nbc_peaks_dir <- file.path("output/chip-pipeline-GRCh38/peak_call/A549_NBC")

###
nbc_spectrl_path <- file.path(nbc_peaks_dir, "A549_NBC_CTRL_specific.bed")
nbc_common_path <- file.path(nbc_peaks_dir, "A549_NBC_common.bed")
nbc_spedex_path <- file.path(nbc_peaks_dir, "A549_NBC_DEX_specific.bed")

nbc_spectrl <- rtracklayer::import(nbc_spectrl_path) # 4176
nbc_common <- rtracklayer::import(nbc_common_path) # 1376
nbc_spedex <- rtracklayer::import(nbc_spedex_path) # 808

###
nbc_spectrl_anno <- annotatePeaks(nbc_spectrl, output = "anno")
nbc_common_anno <- annotatePeaks(nbc_common, output = "anno")
nbc_spedex_anno <- annotatePeaks(nbc_spedex, output = "anno")

NBC_anno_list <- list(nbc_spectrl_anno, nbc_common_anno, nbc_spedex_anno)
names(NBC_anno_list) <- c("Specific NBC (CTRL)", "Common NBC", "Specific NBC (DEX)")
plotAnnoBar(NBC_anno_list)
plotDistToTSS(NBC_anno_list)

###
nbc_spectrl_df <- annotatePeaks(nbc_spectrl, output = "df")
nbc_common_df <- annotatePeaks(nbc_common, output = "df")
nbc_spedex_df <- annotatePeaks(nbc_spedex, output = "df")

plotAnnotation(nbc_spectrl_df)
plotAnnotation(nbc_common_df)
plotAnnotation(nbc_spedex_df)

###
get_unique_genes <- function(anno_df) {
  cat("Number of regions : ", nrow(anno_df), "\n")
  ugenes <- unique(anno_df$SYMBOL)
  cat("Number of unique genes : ", length(ugenes), "\n")
  return(ugenes)
}

ugenes_spectrl <- get_unique_genes(nbc_spectrl_df)
ugenes_common <- get_unique_genes(nbc_common_df)
ugenes_spedex <- get_unique_genes(nbc_spedex_df)

##### Identify NBC regions binding to promoter
nbc_spedex <- tbl_df(nbc_spedex_df) %>% filter(Annot == "Promoter")
nbc_spedex_prom <- nbc_spedex %>% select(annotation, distanceToTSS, geneId, ENTREZID, SYMBOL, GENENAME, Annot) %>% arrange(distanceToTSS)
kable(nbc_spedex_prom)
