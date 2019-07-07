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
inter_7h_9h <- intersect(transcript_7h, transcript_9h) # 2258 | 1496 > 

transcript_9h <- deg[["9hr"]]$fdr0p1$transcript_id # 1496
transcript_11h <- deg[["11hr"]]$fdr0p1$transcript_id # 125
inter_9h_11h <- intersect(transcript_9h, transcript_11h) # 1496 | 125 > 

allTranscripts <- c(inter_1h_3h, inter_3h_5h, inter_5h_7h, inter_7h_9h, inter_9h_11h)
uniqueTranscripts <- unique(AllTranscripts)
length(uniqueTranscripts)
