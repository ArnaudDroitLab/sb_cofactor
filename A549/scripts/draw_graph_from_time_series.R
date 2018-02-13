# setwd("/home/chris/Bureau/sb_cofactor_hr/A549/scripts")

raw_counts <- read.csv("../results/a549_dex_time_points/raw_counts.csv", sep=",", header=TRUE, row.names=1)
nb_genes <- nrow(raw_counts)
nb_sample <- ncol(raw_counts)
# 60483 genes
# 46 samples (corresponding to several time points)

total_expression <- rowSums(raw_counts)
no_expression <- sum(total_expression == 0)
# 22569 genes do not expressed during the 12 hours (because sum=0)

who_no_expression <- which(total_expression == 0)

# Control
expressed_raw_counts <- raw_counts[-who_no_expression,]
expressed_total_expression <- rowSums(expressed_raw_counts)
expressed_no_expression <- sum(expressed_total_expression == 0)
nb_expressed_genes <- dim(expressed_raw_counts)[1]
# Genes that do not expressed have been removed

total_reads <- (colSums(raw_counts)/1000000)
expressed_total_reads <- (colSums(expressed_raw_counts)/1000000)
identical(total_reads, expressed_total_reads)

expressed_counts_norm <- data.frame(t(t(expressed_raw_counts) / expressed_total_reads))
# write.csv(expressed_counts_norm, file="expressed_counts_norm.csv")

########################

map_time <- read.csv("../input/sample_sheet_ENCSR897XFT_without_h.csv", sep=",", header=TRUE)[,-3]

hours <- c(0)
for (i in 1:nb_sample) {
  exp <- colnames(raw_counts)[i]
  ind <- which(exp == map_time$file_accession )
  hours[i] <- map_time$time_point[ind]
}

colnames(expressed_counts_norm) <- hours

########################
map_ens_hgnc <- read.csv("../input/map_ensembl_hgnc.csv", header=TRUE)
map_ens_hgnc$ensembl_gene_id <- as.character(map_ens_hgnc$ensembl_gene_id)
map_ens_hgnc$hgnc_symbol <- as.character(map_ens_hgnc$hgnc_symbol)

# "ENSG00000231709" not in the database

##############################
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)

### draw_a549_dex_time_points function
# genes_list: vector of genes
# counts_mat: counts matrix (row: genes, col: time)
# mapping: correspondance between ensembl_id and symbol_hgnc_id

draw_a549_dex_time_points <- function(genes_list, counts_mat, mapping) {
  inds <- c(NULL)
  for (gene in genes_list) {
    ind_hgnc <- which(gene == mapping$hgnc_symbol)
    id_ensembl <- mapping$ensembl_gene_id[ind_hgnc]
    inds <- c(inds, which(id_ensembl == rownames(counts_mat)))
  }
  
  mat_for_graph <- counts_mat[inds,]
  rownames(mat_for_graph) <- genes_list
  
  # Data
  hours <- as.numeric(colnames(mat_for_graph))
  
  # select(hours, as.factor(rownames(mat_for_graph))) %>%
  df <- data.frame(hours, t(mat_for_graph)) %>% gather(key = "genes", value = "counts", -hours)
  
  # Calculation of mean
  df_mean <- aggregate(log2(df$counts + 1), by=list(df$hours, df$genes), mean)
  colnames(df_mean) <- c("hours", "genes", "mean")

  # Make the plot
  ggplot(df, aes(x = hours, y = log2(counts + 1))) +
    geom_point(aes(color = genes), size = 1) +
    geom_line(data=df_mean, aes(x=hours, y=mean, group=genes, color=genes)) +
    scale_x_continuous(name="Hours", labels=c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12), breaks=c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12)) +
    scale_y_continuous(name="log2(RPM)", labels=comma)
}

### TO RUN
test_genes <- c("BRD4", "HEXIM1", "NELFA", "NCOR1", "NCOR2", "TBL1XR1", "NR3C1", "NIPBL", "WNT7B", "FZD8", "HACD3", "FZD5", "DVL2")
draw_a549_dex_time_points(test_genes, expressed_counts_norm, map_ens_hgnc)

