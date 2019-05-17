# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(knitr)
source("scripts/lincs/lincs.utils.R")

all_signatures <- read.csv("input/lincs/LINCS_consensus_gene_KD_Signatures_all_37275.xls", sep = "\t")

KD_HEPG2 <- get_target(all_signatures, cell_line = "HEPG2", time = "96 h")

OurFavoritesGenes <-c("ZNF707",
                      "UBP1", "SF3A1", "MARS", "EDC4", "USP7", "RBM6", "SNW1", "WTAP",
                      "ZNP707", "KAP1", "RPN2",
                      "PFN1", "LBR",
                      "SULT2A1", "RPL36", "COPG", "ALDH18A1", "MAP4", "MYO1B", "ACADSB", "TMEM214")

available_KD_HEPG2 <- OurFavoritesGenes %in% KD_HEPG2
names(available_KD_HEPG2) <- OurFavoritesGenes
kable(as.data.frame(available_KD_HEPG2))