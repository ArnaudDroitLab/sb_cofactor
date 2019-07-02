setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(tidyverse)

deg_dir <- "results/a549_dex_time_points"

### 0-12h
time_point_0_12h <- c("0.5h", "1h", "2h", "3h", "4h", "5h", "6h", "7h", "8h", "10h", "12h")
upreg_0_12h <- list()
for (timepoint in time_point_0_12h) {
  deg_file <- read_csv(file.path("results/a549_dex_time_points", timepoint))
  upreg <- deg_file %>% filter(padj <= 0.05, log2FoldChange <= -2) %>% arrange(log2FoldChange) 
  upreg_0_12h[[timepoint]] <-  upreg
}

names(upreg_0_12h)
sapply(upreg_0_12h, nrow)

downreg_0_12h <- list()
for (timepoint in time_point_0_12h) {
  deg_file <- read_csv(file.path("results/a549_dex_time_points", timepoint))
  downreg <- deg_file %>% filter(padj <= 0.05, log2FoldChange >= 1.5) %>% arrange(desc(log2FoldChange)) 
  downreg_0_12h[[timepoint]] <- downreg
}

names(downreg_0_12h)
sapply(downreg_0_12h, nrow)

# "MAFK" %in% downreg$symbol

### 0-25m
time_point_0_25m <- c("5m", "10m", "15m", "20m", "25m")
upreg_0_25m <- list()
for (timepoint in time_point_0_25m) {
  deg_file <- read_csv(file.path("results/a549_dex_time_points", timepoint))
  upreg <- deg_file %>% filter(padj <= 0.05, log2FoldChange >= 1) %>% arrange(log2FoldChange) 
  upreg_0_25m[[timepoint]] <-  upreg
}

names(upreg_0_25m)
sapply(upreg_0_25m, nrow)

downreg_0_25m <- list()
for (timepoint in time_point_0_25m) {
  deg_file <- read_csv(file.path("results/a549_dex_time_points", timepoint))
  downreg <- deg_file %>% filter(padj <= 0.05, log2FoldChange <= -1) %>% arrange(desc(log2FoldChange)) 
  downreg_0_25m[[timepoint]] <- downreg
}

names(downreg_0_25m)
sapply(downreg_0_25m, nrow)

    