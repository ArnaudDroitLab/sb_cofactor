# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(tidyverse)
library(purrr)
library(knitr) 
library(DESeq2)
library(EnsDb.Hsapiens.v86)

# Explore GSE104714
raw <- read_tsv("input/GSE104714/GSE104714_DEX_EtOH_expression.tsv") %>% dplyr::select(-X2)
dim(raw)
names(raw)

# Make matrix for DESeq2
mat <- raw %>% dplyr::select(-transcript) %>% as.matrix
rownames(mat) <- raw$transcript

# Make design for DESeq2
samples <- names(raw)[2:length(raw)]
condition <- strsplit(samples, split = "_") %>% purrr::map(1) %>% unlist
time_point <- strsplit(samples, split = "_") %>% purrr::map(2) %>% unlist
replicate <- strsplit(samples, split = "_") %>% purrr::map(3) %>% unlist
design <- data.frame(samples, condition, time_point, replicate)

# Map transcript_id EST to SYMBOL
edb <- EnsDb.Hsapiens.v86
gene2symbol <- genes(edb, columns = "gene_name", return.type = "DataFrame") %>% as.data.frame
transcript2gene <- transcripts(edb, columns = "gene_id", return.type = "DataFrame") %>% as.data.frame
transcript2symbol <- left_join(transcript2gene, gene2symbol, by = "gene_id")
transcript2symbol <- transcript2symbol[, c(2,1,3)]
colnames(transcript2symbol) <- c("tx_id", "gene_id", "symbol")

# DESeq2 functions
get_res <- function(x) {
  message("###\t", x)
  current_design <- dplyr::filter(design, time_point == x)
  current_matrix <- mat[, current_design$samples]
  dds <- DESeqDataSetFromMatrix(current_matrix, current_design, ~ condition)
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  res <- results(dds, c("condition", "DEX", "EtOH"))
}

get_df <- function(x) {
  as.data.frame(x) %>%
    mutate(transcript_id = rownames(x),
           tx_id = strsplit(transcript_id, split = "\\.") %>% purrr::map(1) %>% unlist) %>%
    left_join(transcript2symbol, by = "tx_id") %>%
    dplyr::select(transcript_id, tx_id, gene_id, symbol, everything()) %>%
    arrange(padj)
}

# Set comparisons to do
comparisons <- design$time_point %>% unique %>% as.vector
names(comparisons) <- comparisons
comparisons <- comparisons[paste0(c(1, 3, 5, 7, 9, 11), "hr")]

# DESeq2 analysis
de_res <- map(comparisons, get_res)
de_df <- map(de_res, get_df)

# Control
map_int(de_res, ~ dplyr::filter(as.data.frame(.x), padj <= 0.1) %>% nrow)
map_int(de_df, ~ dplyr::filter(.x, padj <= 0.1) %>% nrow)

# DEG output files
walk(names(de_df), ~ write_csv(de_df[[.x]], paste0("results/a549_dex_time_points_1_11hr/",.x)))