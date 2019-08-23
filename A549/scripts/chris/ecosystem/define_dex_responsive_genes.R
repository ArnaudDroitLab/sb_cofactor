# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(knitr)
library(tidyverse)
options(readr.num_columns = 0)

deg_dir <- "results/a549_dex_time_points"
time_point <- paste0(c(0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12), "h")

deg <- list()
deg_numbers <- data.frame(0)
for (time in time_point) {
  message("##### ", time)
  deg_file <- read_csv(file.path(deg_dir, time), progress = FALSE)
  
  significant <- deg_file %>% dplyr::filter(padj <= 0.1)
  deg[[time]]$fdr0p1 <- significant
  
  numbers_per_time <- c(nrow(deg_file), nrow(significant))
  for (FC in c(0, 0.5, 1, 1.5, 2, 2.5, 3)) {
    message("    # FC : ", FC)
    upreg <- significant %>% dplyr::filter(log2FoldChange >= FC) %>% arrange(desc(log2FoldChange))
    downreg <- significant %>% dplyr::filter(log2FoldChange <= -FC) %>% arrange(log2FoldChange)
    message("      # Number of upreg genes : ", nrow(upreg))
    message("      # Number of downreg genes : ", nrow(downreg))
    numbers_per_time <- c(numbers_per_time, nrow(upreg), nrow(downreg))
    
    FCname <- paste0("FC", FC)
    FCname <- gsub(pattern = "\\.", "p", FCname)
    
    deg[[time]][[FCname]]$upreg <- upreg
    deg[[time]][[FCname]]$downreg <- downreg
  }
  deg_numbers <- cbind(deg_numbers, numbers_per_time)
}

deg_numbers <- deg_numbers <- deg_numbers[, 2:ncol(deg_numbers)]
colnames(deg_numbers) <- time_point
rownames(deg_numbers) <- c("# of genes", "# of significant genes",
                           "# of upregulated (FC > 0)", "# of downregulated (FC < 0)",
                           "# of upregulated (FC > 0.5)", "# of downregulated (FC < -0.5)",
                           "# of upregulated (FC > 1)", "# of downregulated (FC < -1)",
                           "# of upregulated (FC > 1.5)", "# of downregulated (FC < -1.5)",
                           "# of upregulated (FC > 2)", "# of downregulated (FC < -2)",
                           "# of upregulated (FC > 2.5)", "# of downregulated (FC < -2.5)",
                           "# of upregulated (FC > 3)", "# of downregulated (FC < -3)")
kable(deg_numbers)

# How many unique differentially expressed genes over time depending on the FC?
deg_overtime <- c(NA, NA)
deg_overtime_unique <- c(NA, NA)

for (FC in paste0("FC", c("0", "0p5", "1", "1p5", "2", "2p5", "3"))) {
  message("##### ", FC)
  upreg_list <- c()
  downreg_list <- c()
  is_deg <- read_tsv(file.path(deg_dir, "raw_counts_with_colnames.txt")) %>% dplyr::select(gene_id)
  for (time in time_point) {
    # message("     # ", time)
    upreg <- deg[[time]][[FC]]$upreg
    downreg <- deg[[time]][[FC]]$downreg
    # message("      # Number of upreg genes : ", nrow(upreg))
    # message("      # Number of downreg genes : ", nrow(downreg))
    
    upreg_list <- c(upreg_list, upreg$gene_id)
    downreg_list <- c(downreg_list, downreg$gene_id)
    inter_updown <- length(intersect(upreg_list, downreg_list))
    message(" Length of intersection (", time, ") : ", inter_updown)
    if (inter_updown > 0) {
      print(intersect(upreg_list, downreg_list))
    }
    
    vna <- rep(0, length(is_deg$gene_id))
    vna[is_deg$gene_id %in% upreg$gene_id] = 2
    vna[is_deg$gene_id %in% downreg$gene_id] = 1
    is_deg <- cbind(is_deg, vna)
  }
  message("     # Number of differentially upregulated genes : ", length(upreg_list))
  message("     # Number of differentially downregulated genes : ", length(downreg_list))
  message("     # Number of differentially upregulated genes (unique): ", length(unique(upreg_list)))
  message("     # Number of differentially downregulated genes (unique): ", length(unique(downreg_list)))
  
  deg_overtime <- c(deg_overtime, length(upreg_list), length(downreg_list))
  deg_overtime_unique <- c(deg_overtime_unique, length(unique(upreg_list)), length(unique(downreg_list)))
  
  deg$gene_list[[FC]]$upreg <- unique(upreg_list)
  deg$gene_list[[FC]]$downreg <- unique(downreg_list)
  
  colnames(is_deg) <- c("gene_id", time_point)
  deg$is_deg[[FC]] <- is_deg
}

deg_numbers <- cbind(deg_numbers, deg_overtime, deg_overtime_unique)
colnames(deg_numbers) <- c(time_point, "total", "total_unique")
kable(deg_numbers)

deg$deg_numbers <- deg_numbers

# saveRDS(deg, file = "output/analyses/deg.rds")


