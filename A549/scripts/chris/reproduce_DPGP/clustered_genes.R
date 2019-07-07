# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(tidyverse)
library(knitr)
library(purrr)
library(EnsDb.Hsapiens.v86)

# load clusters
clusters <- read_tsv("output/analyses/DPGP_on_a549_dex_1_11hr/alabama_optimal_clustering.txt")
clusters$gene <- strsplit(clusters$gene, split = "_") %>% purrr::map(1) %>% unlist

# Get map transcript_id EST to SYMBOL
edb <- EnsDb.Hsapiens.v86
gene2symbol <- genes(edb, columns = "gene_name", return.type = "DataFrame") %>% as.data.frame
transcript2gene <- transcripts(edb, columns = "gene_id", return.type = "DataFrame") %>% as.data.frame
transcript2symbol <- left_join(transcript2gene, gene2symbol, by = "gene_id")
transcript2symbol <- transcript2symbol[, c(2,1,3)]
colnames(transcript2symbol) <- c("tx_id", "gene_id", "symbol")

# Mapping
clusters_symbol <- left_join(clusters, transcript2symbol, by = c("gene" = "tx_id"))

# How many genes per clusters?
clusters_genes <- clusters_symbol %>% group_by(cluster)%>%
  summarise(nb_elements = n(), nb_genes = length(unique(symbol)))

# clusters_only_symbol
clusters_only_symbol <- clusters_symbol %>% distinct(cluster, gene_id, .keep_all = TRUE)
length(clusters_only_symbol$symbol)
length(unique(clusters_only_symbol$symbol))

# where is my gene?
clusters_symbol %>% dplyr::filter(symbol == "ANGPTL4")
clusters_symbol %>% dplyr::filter(symbol == "IL11")
clusters_symbol %>% dplyr::filter(symbol == "PER1")
clusters_symbol %>% dplyr::filter(symbol == "SOCS1")
clusters_symbol %>% dplyr::filter(symbol == "GPR1")
clusters_symbol %>% dplyr::filter(symbol == "CSF3")
clusters_symbol %>% dplyr::filter(symbol == "KLF6")

clusters_symbol %>% dplyr::filter(symbol == "SNORD41")
clusters_symbol %>% dplyr::filter(symbol == "SNORD99")
clusters_symbol %>% dplyr::filter(symbol == "DRD1")
clusters_symbol %>% dplyr::filter(symbol == "RASSF10")
clusters_symbol %>% dplyr::filter(symbol == "IER2")
clusters_symbol %>% dplyr::filter(symbol == "CXCL8")
clusters_symbol %>% dplyr::filter(symbol == "GDF15")
clusters_symbol %>% dplyr::filter(symbol == "MAFK")
clusters_symbol %>% dplyr::filter(symbol == "AREG")

# how many NAs?
# Check that website to associate trasncript with gene
# https://biobankengine.stanford.edu/transcript/ENST00000260630
for (i in 1:10) {
  message("### cluster ", i)
  nas <- clusters_symbol %>% dplyr::filter(is.na(symbol), cluster == i)
  message(nrow(nas))
}
clusters_symbol %>% dplyr::filter(is.na(symbol), cluster == 10)

