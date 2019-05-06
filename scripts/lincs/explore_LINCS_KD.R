# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(knitr)
source("scripts/lincs/lincs.utils.R")

# Load mutated_cofactors gene symbol
MUTATED_COFACTORS <- get_mutated_cofactors()

# Load all available gene KD signatures from LINCS database
all_signatures <- read.csv("input/lincs/LINCS_consensus_gene_KD_Signatures_all_37275.xls", sep = "\t")

all_cell_lines <- table(all_signatures$CellLine)

overview_KD <- all_signatures %>% group_by(CellLine) %>%
  summarize(NumberOfKD = n(), MutatedCofactor = sum(TargetGene %in% MUTATED_COFACTORS))
  
kable(overview_KD)