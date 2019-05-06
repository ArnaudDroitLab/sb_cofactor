# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(knitr)
source("scripts/lincs/lincs.utils.R")

all_signatures <- read.csv("input/LINCS/LINCS_consensus_gene_KD_Signatures_all_37275.xls", sep = "\t")

hsf1 <- all_signatures %>% filter(TargetGene == "HSF1")

kable(hsf1)