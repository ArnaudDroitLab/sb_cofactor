# setwd("/home/chris/Bureau/sb_cofactor_hr/mESC")
setwd("/Users/chris/Desktop/sb_cofactor_hr/mESC")

library(dplyr)
library(httr)
library(jsonlite)
library(knitr)
source("scripts/ilincs/lincs.utils.R")

all_signatures <- read.csv("input/LINCS/LINCS_consensus_gene_KD_Signatures_all_37275.xls", sep = "\t")

hsf1 <- all_signatures %>% filter(TargetGene == "HSF1")

kable(hsf1)
