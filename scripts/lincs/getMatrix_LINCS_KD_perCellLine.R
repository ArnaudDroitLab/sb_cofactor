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
overview_KD_allTime <- overview_KD %>% summarize(NumberOfKD = n(), MutatedCofactor = sum(TargetGene %in% MUTATED_COFACTORS)) %>%
  arrange(Time, CellLine) %>% select(Time, CellLine, everything())
kable(overview_KD_allTime)

# Remove CellLine where MutatedCofactor == 0
overview_KD_allTime_with_MutCof <- overview_KD_allTime %>% filter(MutatedCofactor != 0) %>% arrange(Time)
kable(overview_KD_allTime_with_MutCof)

# Only 96h
overview_KD_96h <- overview_KD_allTime_with_MutCof %>% filter(Time == "96 h")
kable(overview_KD_96h)
our_cell_lines <- overview_KD_96h %>% filter(MutatedCofactor != 0) %>% pull(CellLine) %>% as.character
our_cell_lines

#### Downloading matrix
cell_lines_toDL <- c("A375", "ASC")
for (CL in cell_lines_toDL) {
  downloadSignature_KD_CellLine(CellLine = CL, time = "96 h", output_path = "output/analysis/lincs")
}