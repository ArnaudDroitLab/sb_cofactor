setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(tidyverse)
library(knitr)
library(purrr)
library(EnsDb.Hsapiens.v86)

# load clusters
clusters <- read_tsv("output/analyses/DPGP_on_a549_dex_0_12hr/bruxelles/bruxelles_FC1_optimal_clustering.txt")

# Get map transcript_id EST to SYMBOL
edb <- EnsDb.Hsapiens.v86
gene2symbol <- genes(edb, columns = "gene_name", return.type = "DataFrame") %>% as.data.frame
gene2symbol <- gene2symbol[, 2:1]
colnames(gene2symbol) <- c("gene_id", "symbol")

# Mapping
clusters_symbol <- left_join(clusters, gene2symbol, by = c("gene" = "gene_id"))
colnames(clusters_symbol) <- c("cluster", "gene_id", "symbol")
write.table(clusters_symbol, file = "output/analyses/DPGP_on_a549_dex_0_12hr/bruxelles/bruxelles_FC1_clusters_gene_symbols.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

# How many genes per clusters?
clusters_genes <- clusters_symbol %>% group_by(cluster) %>% summarise(nb_genes = length(unique(symbol)))
kable(clusters_genes)

# where is my gene?
myGenes <- c("ANGPTL4", "IL11", "PER1", "SOCS1", "GPR1", "CSF3", "KLF6",
             "SNORD41", "SNORD99", "DRD1", "RASSF10", "IER2" ,"CXCL8", "GDF15", "MAFK", "AREG")
where_are_myGenes <- clusters_symbol %>% dplyr::filter(symbol %in% myGenes) %>% arrange(cluster, symbol)
kable(where_are_myGenes)