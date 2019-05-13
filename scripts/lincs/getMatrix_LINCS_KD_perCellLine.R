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

# cell_lines <- c("A549", "HA1E")
cell_lines <- c("HA1E")

for (CL in cell_lines) {
  message("#####\t", CL)
  CL_cell_line <- all_signatures %>% filter(CellLine == CL)

  # Are samples associated with an unique timepoint?
  print(length(unique(CL_cell_line$Time)) == 1)
  
  signIds_CL <- get_signIds(CL_cell_line, cell_line = CL)
  targets_CL <- get_target(CL_cell_line, cell_line = CL)
  sum(MUTATED_COFACTORS %in% targets_CL)
  
  CL_ILINCs_KD_matrix <- downloadSignatureInBatch(signIds_CL, targets_CL)
  saveRDS(CL_ILINCs_KD_matrix, file = paste0("output/analysis/lincs/", CL, "_LINCS_KD_matrix.rds"))
}