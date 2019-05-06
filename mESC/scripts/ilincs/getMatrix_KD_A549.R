# setwd("/home/chris/Bureau/sb_cofactor_hr/mESC")
setwd("/Users/chris/Desktop/sb_cofactor_hr/mESC")

library(dplyr)
library(httr)
library(jsonlite)
source("scripts/ilincs/lincs.utils.R")

all_signatures <- read.csv("input/LINCS/LINCS_consensus_gene_KD_Signatures_all_37275.xls", sep = "\t")

all_cell_lines <- table(all_signatures$CellLine)

A549_cell_line <- all_signatures %>% filter(CellLine == "A549")
unique(A549_cell_line$Time)

signIds_A549 <- get_signIds(A549_cell_line, cell_line = "A549")
targets_A549 <- get_target(A549_cell_line, cell_line = "A549")
sum(mutated_cofactors %in% targets_A549)

###
A549_ILINCs_KD_matrix <- downloadSignatureInBatch(signIds_A549, targets_A549)
saveRDS(A549_ILINCs_KD_matrix, file = "output/analysis/ILINCS/A549_ILINCs_KD_matrix.rds")
