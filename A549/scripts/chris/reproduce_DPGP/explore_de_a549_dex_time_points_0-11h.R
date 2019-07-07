setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(tidyverse)
library(knitr)

deg_dir <- "results/a549_dex_time_points_1_11hr"
time_point <- paste0(c(1, 3, 5, 7, 9, 11), "hr")

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
rownames(deg_numbers) <- c("# of transcripts", "# of significant transcripts",
                           "# of upregulated (FC > 0)", "# of downregulated (FC < 0)",
                           "# of upregulated (FC > 0.5)", "# of downregulated (FC < -0.5)",
                           "# of upregulated (FC > 1)", "# of downregulated (FC < -1)",
                           "# of upregulated (FC > 1.5)", "# of downregulated (FC < -1.5)",
                           "# of upregulated (FC > 2)", "# of downregulated (FC < -2)",
                           "# of upregulated (FC > 2.5)", "# of downregulated (FC < -2.5)")
kable(deg_numbers)

### 
for (time in time_point) {
  message("##### ", time)
  deg_table <- deg[[time]]$fdr0p1
  unique_genes <- deg_table$gene_id %>% unique
  message(length(unique_genes), " / ", nrow(deg_table))
  # res <- deg_table %>% dplyr::filter(symbol == "KLF6")
  # print(kable(res))
}

### 
transcript_list <- c()
gene_list <- c()
for (time in time_point) {
  message("##### ", time)
  deg_table <- deg[[time]]$fdr0p1
  print(length(deg_table$transcript_id))
  transcript_list <- c(transcript_list, deg_table$transcript_id)
  gene_list <- c(gene_list, deg_table$gene_id)
}

length(transcript_list)
length(unique(transcript_list))
length(gene_list)
length(unique(gene_list))

### Intended to choose a transcript per gene, but maybe not necessary
# per_transcript <- list()
# for (transcript in unique(transcript_list)[1:10]) {
#   FC_list <- c()
#   pval_list <- c()
#   for (time in time_point) {
#     deg_table <- deg[[time]]$fdr0p1 %>% dplyr::filter(transcript_id == transcript)
#     FC <- deg_table$log2FoldChange
#     pval <- deg_table$padj
#     
#     FC_list <- c(FC_list, FC)
#     pval_list <- c(pval_list, pval)
#   }
#   per_transcript[[transcript]]$FC_list <- FC_list
#   per_transcript[[transcript]]$pval_list <- pval_list
# }
# 
# head(per_transcript)
# head(purrr::map(per_transcript, 2))
# 
# ###
# library(metap)
# pval_list <- per_transcript$ENST00000263707.5$pval_list
# print(pval_list)
# s <- sumlog(pval_list)
# print(s)

# What are the transcripts that are differentially expressed at least at two consecutives time point
transcript_1h <- deg[["1hr"]]$fdr0p1$transcript_id # 1450
transcript_3h <- deg[["3hr"]]$fdr0p1$transcript_id # 3354
inter_1h_3h <- intersect(transcript_1h, transcript_3h) # 1450 | 3354 > 976

transcript_3h <- deg[["3hr"]]$fdr0p1$transcript_id # 3354
transcript_5h <- deg[["5hr"]]$fdr0p1$transcript_id # 3104
inter_3h_5h <- intersect(transcript_3h, transcript_5h) # 3354 | 3104 > 2297

transcript_5h <- deg[["5hr"]]$fdr0p1$transcript_id # 3104
transcript_7h <- deg[["7hr"]]$fdr0p1$transcript_id #  2258
inter_5h_7h <- intersect(transcript_5h, transcript_7h) # 33104 | 2258 > 1589

transcript_7h <- deg[["7hr"]]$fdr0p1$transcript_id # 2258
transcript_9h <- deg[["9hr"]]$fdr0p1$transcript_id # 1496
inter_7h_9h <- intersect(transcript_7h, transcript_9h) # 2258 | 1496 > 1095

transcript_9h <- deg[["9hr"]]$fdr0p1$transcript_id # 1496
transcript_11h <- deg[["11hr"]]$fdr0p1$transcript_id # 125
inter_9h_11h <- intersect(transcript_9h, transcript_11h) # 1496 | 125 > 97

allTranscripts <- c(inter_1h_3h, inter_3h_5h, inter_5h_7h, inter_7h_9h, inter_9h_11h)
uniqueTranscripts <- unique(allTranscripts)
length(uniqueTranscripts)

# Make matrix for DPGP
raw <- read_tsv("input/GSE104714/GSE104714_DEX_EtOH_expression.tsv") %>% dplyr::select(-X2)
matfiltered <- raw %>% dplyr::filter(transcript %in% uniqueTranscripts)

mat <- matfiltered %>% dplyr::select(-transcript) %>% as.matrix
rownames(mat) <- matfiltered$transcript

##### Perform correlation matrix
# mat <- mat[rowSums(mat) > 100,]

# Annotation for heatmap > condition: EtOH or DEX
annot_condition <- strsplit(colnames(mat), split = "_") %>% purrr::map(1) %>% unlist
ha <- HeatmapAnnotation(Condition = annot_condition,
                        col = list(Condition = c("EtOH" = "#16DB93", "DEX" = "#EFEA5A")))
rowha = rowAnnotation(Condition = annot_condition,
                      col = list(Condition = c("EtOH" = "#16DB93", "DEX" = "#EFEA5A")),
                      show_legend = FALSE)

# Correlation analysis : Pearson method
cor_pearson <- cor(mat, method = "pearson")
min(cor_pearson)
col_pearson <- colorRamp2(c(0.8, 0.9, 1), c("#0f4259", "white", "#800020"))

cor_pearson_heatmap <- Heatmap(cor_pearson, name = "Pearson correlation",
                               row_names_side = "right",
                               row_dend_side = "right",
                               row_dend_width = unit(25, "mm"),
                               column_names_side = "top", column_names_rot = 45,
                               show_column_dend = FALSE,
                               top_annotation = ha,
                               right_annotation = rowha,
                               rect_gp = gpar(col = "white", lwd = 0.5),
                               col = col_pearson)
cor_pearson_heatmap

# Correlation analysis : Spearman method
cor_spearman <- cor(mat, method = "spearman")
min(cor_spearman)
col_spearman <- colorRamp2(c(0.9, 0.95, 1), c("#0f4259", "white", "#800020"))

cor_spearman_heatmap <- Heatmap(cor_spearman, name = "Spearman correlation",
                                row_names_side = "right",
                                row_dend_side = "right",
                                row_dend_width = unit(25, "mm"),
                                column_names_side = "top", column_names_rot = 45,
                                show_column_dend = FALSE,
                                top_annotation = ha,
                                right_annotation = rowha,
                                rect_gp = gpar(col = "white", lwd = 0.5),
                                col = col_spearman)
cor_spearman_heatmap

# Continue to format matrix for DPGP
matdex <- matfiltered %>% dplyr::select(transcript, starts_with("DEX"))
matdex <- lapply(matdex, as.character) %>% as.data.frame

transcript_name <- gsub("\\.", "_", matdex$transcript)

matrep1 <- matdex %>% dplyr::select(contains("Rep1"))
matrep1_final <- lapply(matrep1, paste0, ".0") %>% as.data.frame
rownames(matrep1_final) <- transcript_name
colnames(matrep1_final)
colnames(matrep1_final) <- c("1", "3", "5", "7", "9", "11")
head(matrep1_final)

matrep2 <- matdex %>% dplyr::select(contains("Rep2"))
matrep2_final <- lapply(matrep2, paste0, ".0") %>% as.data.frame
rownames(matrep2_final) <- transcript_name
colnames(matrep2_final)
colnames(matrep2_final) <- c("1", "3", "5", "7", "9", "11")
head(matrep2_final)

output_dir <- "output/analyses/DPGP_on_a549_dex_1_11hr"
write.table(matrep1_final, file = file.path(output_dir, "de_transcripts_A549_1_11h_rep1.txt"),
            quote = FALSE, sep = "\t")
write.table(matrep2_final, file = file.path(output_dir, "de_transcripts_A549_1_11h_rep2.txt"),
            quote = FALSE, sep = "\t")

