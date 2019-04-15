setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(dplyr)

df <- read.csv("scripts/chris/iLINCS/iLINCs_KD_A549.csv", sep = "\t")

A549 <- df %>% filter(CellLine == "A549")

TargetGene <- A549 %>% pull(TargetGene) %>% as.character %>% unique

cofactors <- c("NIPBL", "SMC1A", "SMC3", "RAD21", "HDAC8", "BRD4",
               "CREBBP", "EP300",
               "KAT6A", "KAT6B", "SRCAP", "BRPF1",
               "RAI1",
               "MBD5", "EHMT1", "HDAC4",
               "ARID1A", "ARID1B", "SMARCA4", "SMARCB1", "SMARCE1",
               "ANKRD11",
               "KMT2A",
               "SMARCA2",
               "AFF4",
               "KMT2D", "KDM6A",
               "TAF6",
               "ESCO2",
               "CHD7",
               "DDX11",
               "STAG1",
               "STAG2",
               "SETD5",
               "MED12",
               "MED13L",
               "MED17",
               "MED23",
               "MED25") # 39

there <- cofactors %in% TargetGene
names(there) <- cofactors
there
sum(there) # 14

A549_bis <- A549 %>% filter(TargetGene %in% cofactors)
A549_bis
