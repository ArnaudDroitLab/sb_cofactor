library(tidyverse)
library(DESeq2)
library(EnsDb.Hsapiens.v86)

# Prepare inputs
design <- read_csv("input/sample_sheet_ENCSR897XFT.csv") %>%
    dplyr::select(file_accession, time_point)

raw_counts <- read_csv("results/a549_dex_time_points/raw_counts.csv")
m <- dplyr::select(raw_counts, -gene_id) %>% as.matrix
rownames(m) <- raw_counts$gene_id

edb <- EnsDb.Hsapiens.v86
gene2symbol <- genes(edb, columns = "gene_name", return.type = "DataFrame") %>%
    as.data.frame
gene2symbol <- gene2symbol[,2:1]
colnames(gene2symbol) <- c("gene_id", "symbol")

# DESeq2 analysis
get_res <- function(x) {
    message(x)
    current_design <- filter(design, time_point %in% c("0h", x))
    current_matrix <- m[,current_design$file_accession]
    dds <- DESeqDataSetFromMatrix(current_matrix, current_design, ~ time_point)
    dds <- dds[ rowSums(counts(dds)) > 1, ]
    dds <- DESeq(dds)

    res <- results(dds, c("time_point", "0h", x))
}

get_df <- function(x) {
    as.data.frame(x) %>%
        mutate(gene_id = rownames(x)) %>%
        left_join(gene2symbol, by = "gene_id") %>%
        dplyr::select(gene_id, symbol, everything()) %>%
        arrange(padj)
}

comparisons <- filter(design, time_point != "0h")$time_point %>% unique
names(comparisons) <- comparisons
comparisons <- comparisons[paste0(c(0.5,1:8,10,12), "h")]

de_res <- map(comparisons, get_res)
de_df <- map(de_res, get_df)

map_int(de_res, ~ filter(as.data.frame(.x), padj <= 0.05) %>% nrow)
map_int(de_df, ~ filter(.x, padj <= 0.05) %>% nrow)

walk(names(de_df), ~ write_csv(de_df[[.x]], paste0("results/a549_dex_time_points/",.x)))
