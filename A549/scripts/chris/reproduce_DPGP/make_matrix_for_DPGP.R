  setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
  # setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
  
  library(tidyverse)
  library(knitr)
  library(ComplexHeatmap)
  library(circlize)
  library(EnsDb.Hsapiens.v86)
  
  #### Save heatmap
  saveHeatmap <- function(heatmap_obj, output_dir, output_file, width_val = 25, height_val = 22, format = "pdf") {
    output_filepath <- file.path(output_dir, paste0(output_file, ".", format))
    pdf(file = output_filepath, width = width_val, height = height_val)
    print(heatmap_obj)
    dev.off()
    message(" > Heatmap saved in ", output_filepath)
  }
  
  #### Parameters
  # output_dir <- "output/analyses/DPGP_on_a549_dex_0_12hr"
  output_dir <- "output/analyses/DPGP_on_a549_dex_0_6hr"
  FC_val <- 2
  
  ####
  deg_dir <- "results/a549_dex_time_points"
  #time_point <- paste0(c(0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12), "h")
  time_point <- paste0(c(0.5, 1, 2, 3, 4, 5, 6), "h")
  
  deg <- list()
  deg_numbers <- data.frame(0)
  for (time in time_point) {
    message("##### ", time)
    deg_file <- read_csv(file.path(deg_dir, time), progress = FALSE)
    deg[[time]]$DEG <- deg_file
    
    significant <- deg_file %>% dplyr::filter(padj <= 0.1)
    deg[[time]]$fdr0p1 <- significant
    
    upreg_0 <- significant %>% dplyr::filter(log2FoldChange >= 0) %>% arrange(desc(log2FoldChange))
    downreg_0 <- significant %>% dplyr::filter(log2FoldChange <= -0) %>% arrange(log2FoldChange)
    deg[[time]]$upreg_0 <- upreg_0
    deg[[time]]$downreg_0 <- downreg_0
    
    upreg_0p5 <- significant %>% dplyr::filter(log2FoldChange >= 0.5) %>% arrange(desc(log2FoldChange))
    downreg_0p5 <- significant %>% dplyr::filter(log2FoldChange <= -0.5) %>% arrange(log2FoldChange)
    deg[[time]]$upreg_0p5 <- upreg_0p5
    deg[[time]]$downreg_0p5 <- downreg_0p5
    
    upreg_1 <- significant %>% dplyr::filter(log2FoldChange >= 1) %>% arrange(desc(log2FoldChange))
    downreg_1 <- significant %>% dplyr::filter(log2FoldChange <= -1) %>% arrange(log2FoldChange)
    deg[[time]]$upreg_1 <- upreg_1
    deg[[time]]$downreg_1 <- downreg_1
    
    upreg_1p5 <- significant %>% dplyr::filter(log2FoldChange >= 1.5) %>% arrange(desc(log2FoldChange))
    downreg_1p5 <- significant %>% dplyr::filter(log2FoldChange <= -1.5) %>% arrange(log2FoldChange)
    deg[[time]]$upreg_1p5 <- upreg_1p5
    deg[[time]]$downreg_1p5 <- downreg_1p5
    
    upreg_2 <- significant %>% dplyr::filter(log2FoldChange >= 2) %>% arrange(desc(log2FoldChange))
    downreg_2 <- significant %>% dplyr::filter(log2FoldChange <= -2) %>% arrange(log2FoldChange)
    deg[[time]]$upreg_2 <- upreg_2
    deg[[time]]$downreg_2 <- downreg_2
    
    upreg_2p5 <- significant %>% dplyr::filter(log2FoldChange >= 2.5) %>% arrange(desc(log2FoldChange))
    downreg_2p5 <- significant %>% dplyr::filter(log2FoldChange <= -2.5) %>% arrange(log2FoldChange)
    deg[[time]]$upreg_2p5 <- upreg_2p5
    deg[[time]]$downreg_2p5 <- downreg_2p5
    
    deg_numbers <- cbind(deg_numbers, c(nrow(deg_file), nrow(significant),
                                        nrow(upreg_0), nrow(downreg_0),
                                        nrow(upreg_0p5), nrow(downreg_0p5),
                                        nrow(upreg_1), nrow(downreg_1),
                                        nrow(upreg_1p5), nrow(downreg_1p5),
                                        nrow(upreg_2), nrow(downreg_2),
                                        nrow(upreg_2p5), nrow(downreg_2p5)))
  }
  
  deg_numbers <- deg_numbers[, 2:ncol(deg_numbers)]
  colnames(deg_numbers) <- time_point
  rownames(deg_numbers) <- c("# of genes", "# of significant genes",
                             "# of upregulated (FC > 0)", "# of downregulated (FC < 0)",
                             "# of upregulated (FC > 0.5)", "# of downregulated (FC < -0.5)",
                             "# of upregulated (FC > 1)", "# of downregulated (FC < -1)",
                             "# of upregulated (FC > 1.5)", "# of downregulated (FC < -1.5)",
                             "# of upregulated (FC > 2)", "# of downregulated (FC < -2)",
                             "# of upregulated (FC > 2.5)", "# of downregulated (FC < -2.5)")
  kable(deg_numbers)
  
  # Keep only genes that are differentially expressed at two consecutive timepoints
  de_at_two_consecutive_timepoint <- c()
  for (i in 1:(length(time_point)-1)) {
    t1 <- time_point[i]
    t2 <- time_point[i+1]
    message("### ", t1, " vs ", t2)
    de_gene_t1 <- deg[[t1]]$fdr0p1 %>% dplyr::filter(abs(log2FoldChange) >= FC_val) %>% pull(gene_id)
    de_gene_t2 <- deg[[t2]]$fdr0p1 %>% dplyr::filter(abs(log2FoldChange) >= FC_val) %>% pull(gene_id)
    inter_t1_t2 <- intersect(de_gene_t1, de_gene_t2)
    
    message(" # at ", t1, " : ", length(de_gene_t1), " DEGs")
    message(" # at ", t2, " : ", length(de_gene_t2), " DEGs")
    message(" # at ", t1, " and ", t2, " : ", length(inter_t1_t2), " DEGs")
    de_at_two_consecutive_timepoint <- c(de_at_two_consecutive_timepoint, inter_t1_t2)
  }
  
  length(de_at_two_consecutive_timepoint)
  unique_DEGs <- unique(de_at_two_consecutive_timepoint)
  length(unique_DEGs)

# Get mapping ENSG to SYMBOL
edb <- EnsDb.Hsapiens.v86
gene2symbol <- genes(edb, columns = "gene_name", return.type = "DataFrame") %>% as.data.frame
gene2symbol <- gene2symbol[, 2:1]
colnames(gene2symbol) <- c("gene_id", "symbol")

# Left_join
unique_DEGS_wSymbol <- dplyr::left_join(data.frame(gene_id = unique_DEGs), gene2symbol, by = "gene_id")

unique_DEGS_wSymbol %>% dplyr::filter(is.na(symbol))
sum(is.na(unique_DEGS_wSymbol$symbol))

# without SYMBOL
# ENSG00000278309 _ RP11-102K13.5 _ lincRNA
# ENSG00000255414 _ not in the current EnsEMBL _ LINC01059 _ long intergenic non-protein coding RNA
# ENSG00000207770 _ hsa-mir-568 _ miRNA
# ENSG00000263470 _ RP11-160O5.1 _ lincRNA
# ENSG00000254667 _ not in the current EnsEMBL _ LOC100128242 _ uncharacterized
# ENSG00000188206 _ HNRNPU-AS1 _ HNRNPU antisense RNA 1
# ENSG00000259264 _ RP11-60L3.1 _ antisense

# Make matrix for DPGP
raw <- read_tsv("results/a549_dex_time_points/raw_counts_with_colnames.txt")

#
# matfiltered <- raw %>% dplyr::filter(gene_id %in% unique_DEGs)
matfiltered <- raw %>% dplyr::select(matches("gene_id|^0h|^0.5h|^1h|^2h|^3h|^4h|^5h|^6h")) %>%
  dplyr::filter(gene_id %in% unique_DEGs)
names(matfiltered)

mat <- matfiltered %>% dplyr::select(-gene_id) %>% as.matrix
rownames(mat) <- matfiltered$gene_id

##### Perform correlation matrix
# mat <- mat[rowSums(mat) > 100,]

# Annotation for heatmap > time_point
annot_timepoint <- strsplit(colnames(mat), split = "_") %>% purrr::map(1) %>% unlist
customColors <- c("0h" = "#f44336", "0.5h" = "#e91e63", "1h" = "#9c27b0", "2h" = "#673ab7",
                  "3h" = "#3f51b5", "4h" = "#2196f3", "5h" = "#009688", "6h" = "#4caf50",
                  "7h" = "#8bc34a", "8h" = "#ffc107", "10h" = "#ff5722", "12h" = "#607d8b")
ha <- HeatmapAnnotation(Timepoint = annot_timepoint,
                        col = list(Timepoint = customColors))
rowha = rowAnnotation(Timepoint = annot_timepoint,
                      col = list(Timepoint = c(customColors)),
                      show_legend = FALSE)

# Correlation analysis : Pearson method
cor_pearson <- cor(mat, method = "pearson")
min(cor_pearson)
col_pearson <- colorRamp2(c(0.3, 0.65, 1), c("#0f4259", "white", "#800020"))

cor_pearson_heatmap <- Heatmap(cor_pearson, name = "Pearson correlation",
                               row_names_side = "left",
                               row_dend_side = "left",
                               row_dend_width = unit(30, "mm"),
                               column_names_side = "top", column_names_rot = 45,
                               show_column_dend = FALSE,
                               top_annotation = ha,
                               left_annotation = rowha,
                               rect_gp = gpar(col = "white", lwd = 0.5),
                               col = col_pearson)
cor_pearson_heatmap

# saveHeatmap(heatmap_obj = cor_pearson_heatmap,
#             output_dir = output_dir,
#             output_file = "20190711_corRNA0_12_DEG_FC1p5_pearson.pdf",
#             format = "pdf",
#             width = 17, height = 12)

# Correlation analysis : Spearman method
cor_spearman <- cor(mat, method = "spearman")
min(cor_spearman)
col_spearman <- colorRamp2(c(0.6, 0.8, 1), c("#0f4259", "white", "#800020"))

cor_spearman_heatmap <- Heatmap(cor_spearman, name = "Spearman correlation",
                                row_names_side = "left",
                                row_dend_side = "left",
                                row_dend_width = unit(30, "mm"),
                                column_names_side = "top", column_names_rot = 45,
                                show_column_dend = FALSE,
                                top_annotation = ha,
                                left_annotation = rowha,
                                rect_gp = gpar(col = "white", lwd = 0.5),
                                col = col_spearman)
cor_spearman_heatmap

# saveHeatmap(heatmap_obj = cor_spearman_heatmap,
#             output_dir = output_dir,
#             output_file = "20190711_corRNA0_12_DEG_FC1p5_spearman.pdf",
#             format = "pdf",
#             width = 17, height = 12)

# Continue to format matrix for DPGP
matrix_for_DPGP <- function(matrix, rep, time_point = c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12)) {
  column_order <- paste0(time_point, "h_", rep)
  gene_id_tmp <- matrix$gene_id
  matrep <- matrix %>% dplyr::select(column_order)
  
  matrep <- lapply(matrep, as.character) %>% as.data.frame
  matrep_final <- lapply(matrep, paste0, ".0") %>% as.data.frame
  
  rownames(matrep_final) <- gene_id_tmp
  colnames(matrep_final) <- time_point
  
  return(matrep_final)
}

matrep2 <- matrix_for_DPGP(matfiltered, rep = "rep2", time_point = c(0, 0.5, 1, 2, 3, 4, 5, 6))
matrep3 <- matrix_for_DPGP(matfiltered, rep = "rep3", time_point = c(0, 0.5, 1, 2, 3, 4, 5, 6))
matrep4 <- matrix_for_DPGP(matfiltered, rep = "rep4", time_point = c(0, 0.5, 1, 2, 3, 4, 5, 6))

write.table(matrep2, file = file.path(output_dir, "de_transcripts_A549_0_6h_FC1_rep2.txt"),
            quote = FALSE, sep = "\t")
write.table(matrep3, file = file.path(output_dir, "de_transcripts_A549_0_6h_FC1_rep3.txt"),
            quote = FALSE, sep = "\t")
write.table(matrep4, file = file.path(output_dir, "de_transcripts_A549_0_6h_FC1_rep4.txt"),
            quote = FALSE, sep = "\t")
