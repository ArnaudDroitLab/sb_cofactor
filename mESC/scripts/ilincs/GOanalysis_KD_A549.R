# setwd("/home/chris/Bureau/sb_cofactor_hr/mESC")
setwd("/Users/chris/Desktop/sb_cofactor_hr/mESC")

library(dplyr)
library(httr)
library(jsonlite)
source("scripts/ilincs/lincs.utils.R")

A549_ILINCs_KD_matrix_path <- "output/analysis/ILINCS/A549_ILINCs_KD_matrix.rds"
A549_ILINCs_KD_matrix <- readRDS(A549_ILINCs_KD_matrix_path)

genes_A549 <- A549_ILINCs_KD_matrix %>% pull(Name_GeneSymbol)

# Obtain gene-to-GO mappings
library(biomaRt)
listMarts() 
listMarts(host="www.ensembl.org")
m <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(m)
m <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
listAttributes(m)
go <- getBM(attributes = c("hgnc_symbol", "go_id", "name_1006", "namespace_1003"), mart = m)

go[, 4] != 