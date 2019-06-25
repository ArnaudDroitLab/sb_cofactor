setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(tidyverse)

deg_dir <- "results/a549_dex_time_points"
time_point <- c("1h", "2h", "3h", "4h")

deg_file <- read_csv("results/a549_dex_time_points/1h")
upreg <- deg_file %>% filter(padj <= 0.05, log2FoldChange <= -2) %>% arrange(log2FoldChange)
upreg
downreg <- deg_file %>% filter(padj <= 0.05, log2FoldChange >= 1.5) %>% arrange(desc(log2FoldChange))
downreg

"MAFK" %in% downreg$symbol

    