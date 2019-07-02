# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(knitr)
source("scripts/lincs/lincs.utils.R")

# Load mutated_cofactors gene symbol
MUTATED_COFACTORS <- get_mutated_cofactors()

# Load all available gene OE signatures from LINCS database
all_signatures <- read.csv("input/lincs/LINCS_gene_overexpression_Signatures_all_9291.xls", sep = "\t")

#### Explore LINCS OE database
all_cell_lines <- table(all_signatures$CellLine)
overview_OE <- all_signatures %>% group_by(CellLine, Time)

# All time points
overview_OE_allTime <- overview_OE %>% summarize(NumberOfOE = n(), MutatedCofactor = sum(TargetGene %in% MUTATED_COFACTORS)) %>%
  arrange(Time, CellLine) %>% select(Time, CellLine, everything())
kable(overview_OE_allTime)

# Remove CellLine where MutatedCofactor == 0
overview_OE_allTime_with_MutCof <- overview_OE_allTime %>% filter(MutatedCofactor != 0) %>% arrange(Time)
kable(overview_OE_allTime_with_MutCof)

# Only 96h
overview_OE_96h <- overview_OE_allTime_with_MutCof %>% filter(Time == "96 h")
kable(overview_OE_96h)
our_cell_lines <- overview_OE_96h %>% filter(MutatedCofactor != 0) %>% pull(CellLine) %>% as.character
our_cell_lines