setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(tidyverse)
library(knitr)
library(purrr)
library(EnsDb.Hsapiens.v86)

print_cluster_n <- function(clusters_symbol, n) {
  cluster_members <- clusters_symbol %>% dplyr::filter(cluster == n) %>% pull(symbol)
  message("Cluster ", n, " has ", length(cluster_members), " members")
  print(cluster_members)
}

# load clusters
clusters <- read_tsv("output/analyses/DPGP_on_a549_dex_0_6hr/evry_FC2/evry_FC2_optimal_clustering.txt")

# Get map transcript_id EST to SYMBOL
edb <- EnsDb.Hsapiens.v86
gene2symbol <- genes(edb, columns = "gene_name", return.type = "DataFrame") %>% as.data.frame
gene2symbol <- gene2symbol[, 2:1]
colnames(gene2symbol) <- c("gene_id", "symbol")

# Mapping
clusters_symbol <- left_join(clusters, gene2symbol, by = c("gene" = "gene_id"))
colnames(clusters_symbol) <- c("cluster", "gene_id", "symbol")
write.table(clusters_symbol, file = "output/analyses/DPGP_on_a549_dex_0_6hr/evry_FC2/evry_FC2_clusters_gene_symbols.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

# How many genes per clusters?
clusters_genes <- clusters_symbol %>% group_by(cluster) %>% summarise(nb_genes = length(unique(symbol)))
kable(clusters_genes)
mean(clusters_genes$nb_genes)
sd(clusters_genes$nb_genes)

# where is my gene?
# TODO
myGenes_activated <- c("PER1", "ZFP36",
                       "TFCP2L1", "IGFBP1", "ANGPTL4", "BIRC3", "ENTPD2", "CIDEC", "TIPARP", "BCL6", "SLC19A2",
                       "NFKBIA", "CEBPD", "SDPR", "PAMCI", "C9orf150", "KLF9", "KLF6")

myGenes_repressed <- c("MAFK", "IL11", "EDN1", "ED3", "FZD2", "BDKRB2", "MIDN", "GDF15")

# "SOCS1", "GPR1", "CSF3", "SNORD41", "SNORD99", "DRD1", "RASSF10", "IER2" ,"CXCL8", "GDF15", "MAFK", "AREG")
where_are_myGenes_activated <- clusters_symbol %>% dplyr::filter(symbol %in% myGenes_activated) %>% arrange(cluster, symbol)
kable(where_are_myGenes_activated)

where_are_myGenes_repressed <- clusters_symbol %>% dplyr::filter(symbol %in% myGenes_repressed) %>% arrange(cluster, symbol)
kable(where_are_myGenes_repressed)

#
print_cluster_n(clusters_symbol, 3)

