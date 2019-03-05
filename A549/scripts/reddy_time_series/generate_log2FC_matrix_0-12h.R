# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(plyr)

de_0.5h <- read_csv("results/a549_dex_time_points/0.5h") %>% select(gene_id, log2FC_0.5h = log2FoldChange)
de_1h <- read_csv("results/a549_dex_time_points/1h") %>% select(gene_id, log2FC_1h = log2FoldChange)
de_2h <- read_csv("results/a549_dex_time_points/2h") %>% select(gene_id, log2FC_2h = log2FoldChange)
de_3h <- read_csv("results/a549_dex_time_points/3h") %>% select(gene_id, log2FC_3h = log2FoldChange)
de_4h <- read_csv("results/a549_dex_time_points/4h") %>% select(gene_id, log2FC_4h = log2FoldChange)
de_5h <- read_csv("results/a549_dex_time_points/5h") %>% select(gene_id, log2FC_5h = log2FoldChange)
de_6h <- read_csv("results/a549_dex_time_points/6h") %>% select(gene_id, log2FC_6h = log2FoldChange)
de_7h <- read_csv("results/a549_dex_time_points/7h") %>% select(gene_id, log2FC_7h = log2FoldChange)
de_8h <- read_csv("results/a549_dex_time_points/8h") %>% select(gene_id, log2FC_8h = log2FoldChange)
de_10h <- read_csv("results/a549_dex_time_points/10h") %>% select(gene_id, log2FC_10h = log2FoldChange)
de_12h <- read_csv("results/a549_dex_time_points/12h") %>% select(gene_id, log2FC_12h = log2FoldChange)

dfList = list(de_0.5h, de_1h, de_2h, de_3h, de_4h, de_5h, de_6h, de_7h, de_8h, de_10h, de_12h)
FC_mat <- join_all(dfList)

gene_id <- FC_mat$gene_id
FC_mat$gene_id <- 0
rownames(FC_mat) <- gene_id
colnames(FC_mat) <- c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12)

# Multiplication par -1 car inversion
FC_mat <- FC_mat * (-1)

write.table(FC_mat, file="results/a549_dex_time_points/FC_mat.csv", sep=",", quote=FALSE)
