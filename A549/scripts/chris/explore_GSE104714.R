setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(tidyverse)

raw <- read_tsv("input/GSE104714/GSE104714_DEX_EtOH_expression.tsv") %>% select(-X2)
names(raw)
