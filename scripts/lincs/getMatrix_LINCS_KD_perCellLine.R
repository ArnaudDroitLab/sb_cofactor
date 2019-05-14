# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(knitr)
source("scripts/lincs/lincs.utils.R")

# Load mutated_cofactors gene symbol
MUTATED_COFACTORS <- get_mutated_cofactors()

# Load all available gene KD signatures from LINCS database
all_signatures <- read.csv("input/lincs/LINCS_consensus_gene_KD_Signatures_all_37275.xls", sep = "\t")

#### Explore LINCS KD database
all_cell_lines <- table(all_signatures$CellLine)
overview_KD <- all_signatures %>% group_by(CellLine, Time)

# All time points
overview_KD_allTime <- overview_KD %>% summarize(NumberOfKD = n(), MutatedCofactor = sum(TargetGene %in% MUTATED_COFACTORS), Time = "TO COMPLETE")
kable(overview_KD_allTime)

# Only 96h
overview_KD_96h <- overview_KD %>% filter(Time == "96 h") %>% summarize(NumberOfKD = n(), MutatedCofactor = sum(TargetGene %in% MUTATED_COFACTORS))
kable(overview_KD_96h)
kable(overview_KD_96h %>% filter(MutatedCofactor != 0))
our_cell_lines <- overview_KD_96h %>% filter(MutatedCofactor != 0) %>% pull(CellLine)
our_cell_lines

#### Downloading matrix
cell_lines_toDL <- c("A375")
for (CL in cell_lines_toDL) {
  downloadSignature_KD_CellLine(CellLine = CL, time = "96 h", output_path = "output/analysis/lincs")
}