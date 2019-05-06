# setwd("/home/chris/Bureau/sb_cofactor_hr/mESC")
setwd("/Users/chris/Desktop/sb_cofactor_hr/mESC")

library(dplyr)
library(httr)
library(jsonlite)
library(systemPipeR)
source("scripts/ilincs/lincs.utils.R")

A549_ILINCs_KD_matrix_path <- "output/analysis/ILINCS/A549_ILINCs_KD_matrix.rds"
A549_ILINCs_KD_matrix <- readRDS(A549_ILINCs_KD_matrix_path)

genes_A549 <- A549_ILINCs_KD_matrix %>% pull(Name_GeneSymbol)
a549 <- c()
a549$genes <- genes_A549

# Obtain gene-to-GO mappings
library(biomaRt)
listMarts() 
listMarts(host="www.ensembl.org")
m <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(m)
m <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
listAttributes(m)
go <- getBM(attributes = c("hgnc_symbol", "go_id", "name_1006", "namespace_1003"), mart = m)
go <- go[go[, 4] != "",]
go[, 4] <- as.character(go[, 4])
go[go[, 4] == "molecular_function", 4] <- "F"
go[go[, 4] == "biological_process", 4] <- "P"
go[go[, 4]=="cellular_component", 4] <- "C"
go[1:4, ]
dir.create("./input/GO")
write.table(go, "input/GO/GOannotationsBiomart_mod.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
catdb <- makeCATdb(myfile = "input/GO/GOannotationsBiomart_mod.txt", lib = NULL, org = "", colno = c(2,1,4), idconv = NULL)
save(catdb, file="input/GO/catdb.RData")

# Batch GO termn enrichment analysis
BatchResult <- GOCluster_Report(catdb=catdb, setlist=a549, method="all", id_type="gene", CLSZ=2, cutoff=0.9, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
