# setwd("/home/chris/Bureau/sb_cofactor_hr/mESC")
setwd("/Users/chris/Desktop/sb_cofactor_hr/mESC")

library(dplyr)
library(knitr)
source("scripts/ilincs/lincs.utils.R")

all_signatures <- read.csv("input/LINCS/LINCS_consensus_gene_KD_Signatures_all_37275.xls", sep = "\t")

all_cell_lines <- table(all_signatures$CellLine)

nb_sign <- c()
nb_mutated <- c()
for (cell_line in names(all_cell_lines)) {
  df <- all_signatures %>% filter(CellLine == cell_line)
  
  nb_sign <- c(nb_sign, nrow(df))
  
  tmp_mutated <- sum(mutated_cofactors %in% df$TargetGene)
  nb_mutated<- c(nb_mutated, tmp_mutated)
}

overview_KD <- data.frame("CellLine" = names(all_cell_lines),
                          "NumberOfKD" = nb_sign,
                          "MutatedCofactor" = nb_mutated)

kable(overview_KD)