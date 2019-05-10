# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

source("scripts/lincs/lincs.utils.R")

# Load mutated_cofactors gene symbol
MUTATED_COFACTORS <- get_mutated_cofactors()

# Load all available gene KD signatures from LINCS database
all_signatures <- read.csv("input/lincs/LINCS_consensus_gene_KD_Signatures_all_37275.xls", sep = "\t")

# Explore LINCS KD database
all_cell_lines <- table(all_signatures$CellLine)
overview_KD <- all_signatures %>% group_by(CellLine) %>%
  summarize(NumberOfKD = n(), MutatedCofactor = sum(TargetGene %in% MUTATED_COFACTORS))
kable(overview_KD)


A549_cell_line <- all_signatures %>% filter(CellLine == "A549")
unique(A549_cell_line$Time)

signIds_A549 <- get_signIds(A549_cell_line, cell_line = "A549")
targets_A549 <- get_target(A549_cell_line, cell_line = "A549")
sum(MUTATED_COFACTORS %in% targets_A549)

###
A549_ILINCs_KD_matrix <- downloadSignatureInBatch(signIds_A549, targets_A549)
saveRDS(A549_ILINCs_KD_matrix, file = "output/analysis/lincs/A549_LINCS_KD_matrix.rds")