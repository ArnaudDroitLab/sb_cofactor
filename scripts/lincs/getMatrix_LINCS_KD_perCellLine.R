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
overview_KD <- all_signatures %>% group_by(CellLine)

# All time points
overview_KD_allTime <- overview_KD %>% summarize(NumberOfKD = n(), MutatedCofactor = sum(TargetGene %in% MUTATED_COFACTORS))
kable(overview_KD_allTime)

# Only 96h: time point available in all our favorites cell lines
overview_KD_96h <- overview_KD %>% filter(Time == "96 h") %>% summarize(NumberOfKD = n(), MutatedCofactor = sum(TargetGene %in% MUTATED_COFACTORS))
kable(overview_KD_96h)
kable(overview_KD_96h %>% filter(MutatedCofactor != 0))

# Our cell lines at 96 h
cell_lines <- c("A375", "A549", "ASC", "HA1E", "HCC515", "HEKTE", "HEPG2", "HT29", "MCF7", "NPC", "PC3")
overview_KD_96h_OurCellLines <- overview_KD_96h %>% filter(CellLine %in% cell_lines)
kable(overview_KD_96h_OurCellLines)

#### Downloading matrix
cell_lines_toDL <- c("HEKTE")
for (CL in cell_lines_toDL) {
  downloadSignature_KD_CellLine(CellLine = CL, time = "96 h", output_path = "output/analysis/lincs")
}