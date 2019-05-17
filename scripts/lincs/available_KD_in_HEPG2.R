# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(knitr)
library(ComplexHeatmap)
library(circlize)
source("scripts/lincs/lincs.utils.R")

all_signatures <- read.csv("input/lincs/LINCS_consensus_gene_KD_Signatures_all_37275.xls", sep = "\t")

### Available KD HEPG2
KD_HEPG2 <- get_target(all_signatures, cell_line = "HEPG2", time = "96 h")

OurFavoritesGenes <-c("ZNF707",
                      "UBP1", "SF3A1", "MARS", "EDC4", "USP7", "RBM6", "SNW1", "WTAP",
                      "ZNP707", "KAP1", "RPN2",
                      "PFN1", "LBR",
                      "SULT2A1", "RPL36", "COPG", "ALDH18A1", "MAP4", "MYO1B", "ACADSB", "TMEM214")

available_KD_HEPG2 <- OurFavoritesGenes %in% KD_HEPG2
names(available_KD_HEPG2) <- OurFavoritesGenes
kable(as.data.frame(available_KD_HEPG2))

####################
#   ZNF707
####################

### ZNF707 KD in 7 cell lines
ZNF707 <- all_signatures %>% filter(TargetGene %in% "ZNF707", Time == "96 h")
signIds_ZNF707 <- ZNF707 %>% pull(SignatureId) %>% as.character
targets_ZNF707 <- paste0(ZNF707 %>% pull(CellLine) %>% as.character, "_ZNF707") %>% as.character
signature_ZNF707 <- downloadSignatureInBatch(signIds_ZNF707, targets_ZNF707)

rownames(signature_ZNF707) <- signature_ZNF707$Name_GeneSymbol
signature_ZNF707 <- signature_ZNF707 %>% select(-Name_GeneSymbol)

mat <- as.matrix(signature_ZNF707)
max(mat)
min(mat)

# Heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("#0f4259", "white", "#800020"))
heatmap_sign_znf707 <- Heatmap(mat, name = "LogDiffExp",
                               row_names_side = "left",
                               row_names_gp = gpar(fontsize = 0),
                               row_dend_side = "right",
                               show_row_dend = FALSE,
                               column_names_side = "top",
                               column_names_gp = gpar(fontsize = 11),
                               column_names_rot = 45,
                               column_dend_side = "bottom",
                               column_dend_height = unit(30, "mm"),
                               col = col_fun,
                               rect_gp = gpar(col = "white", lwd = 0))
pdf(file = "output/analysis/lincs/znf707/heatmap_signature_KD_ZNF707_20190517.pdf", width = 10, height = 12)
print(heatmap_sign_znf707)
dev.off()

# Pearson correlation between these 7 ZNF707
cor.pearson.znf707 <- cor(mat, method = "pearson")
max(cor.pearson.znf707 - diag(nrow(cor.pearson.znf707))) # -diag(n) allow to remove the identity matrix
min(cor.pearson.znf707)

col_fun2 = colorRamp2(c(-1, 0, 1), c("#0f4259", "white", "#800020"))
heatmap_ZNF707 <- Heatmap(cor.pearson.znf707, name = "Pearson correlation",
                          row_names_side = "left",
                          row_dend_side = "right",
                          column_names_side = "top",
                          column_names_rot = 45,
                          column_dend_side = "bottom",
                          row_dend_width = unit(30, "mm"),
                          column_dend_height = unit(30, "mm"),
                          column_dend_reorder = TRUE,
                          col = col_fun2,
                          rect_gp = gpar(col = "white", lwd = 2),
                          cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%.2f", cor.pearson.znf707[i, j]), x, y, gp = gpar(fontsize = 15))
                          })
pdf(file = "output/analysis/lincs/znf707/heatmap_pearson_KD_ZNF707_20190517.pdf", width = 14, height = 11)
print(heatmap_ZNF707)
dev.off()

####################
#   ZNF707 + others
####################

### ZNF707 KD in 7 cell lines + others
ZNF707o <- all_signatures %>% filter(TargetGene %in% OurFavoritesGenes, Time == "96 h")
signIds_ZNF707o <- ZNF707 %>% pull(SignatureId) %>% as.character
targets_ZNF707o <- paste0(ZNF707o %>% pull(CellLine) %>% as.character, "_", ZNF707o %>% pull(TargetGene) %>% as.character)
signature_ZNF707o <- downloadSignatureInBatch(signIds_ZNF707o, targets_ZNF707o)

col_celllines <- ZNF707o %>% pull(CellLine) %>% as.character
colors_celllines <- c("HA1E" = "#A4036F", "A549" = "#048BA8", "ASC" = "#16DB93", "A375" = "#EFEA5A", "HEKTE" = "#F29E4C",
                      "HCC515" = "#0f4259", "HEPG2" = "#800020", "HT29" = "#241023", "MCF7" = "#D5E68D",
                      "NPC" = "#C2E7DA", "PC3" = "#839788", "SKL" = "#BAA898")

ha <- HeatmapAnnotation(CellLine = col_celllines,
                        col = list(CellLine = colors_celllines))
rowha = rowAnnotation(CellLine = col_celllines,
                      col = list(CellLine = colors_celllines),
                      show_legend = FALSE)

rownames(signature_ZNF707o) <- signature_ZNF707o$Name_GeneSymbol
signature_ZNF707o <- signature_ZNF707o %>% select(-Name_GeneSymbol)

mat <- as.matrix(signature_ZNF707o)
max(mat)
min(mat)

# Heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("#0f4259", "white", "#800020"))

heatmap_sign_znf707o <- Heatmap(mat, name = "LogDiffExp",
                                row_names_side = "left",
                                row_names_gp = gpar(fontsize = 0),
                                row_dend_side = "right",
                                show_row_dend = FALSE,
                                column_names_side = "top",
                                column_names_gp = gpar(fontsize = 11),
                                column_names_rot = 45,
                                column_dend_side = "bottom",
                                row_dend_width = unit(50, "mm"),
                                column_dend_height = unit(50, "mm"),
                                col = col_fun,
                                rect_gp = gpar(col = "white", lwd = 0),
                                top_annotation = ha)
pdf(file = "output/analysis/lincs/znf707/heatmap_signature_KD_ZNF707_others_20190517.pdf", width = 25, height = 22)
print(heatmap_sign_znf707o)
dev.off()

# Pearson correlation between these 7 ZNF707
cor.pearson.znf707o <- cor(mat, method = "pearson")
max(cor.pearson.znf707o - diag(nrow(cor.pearson.znf707o))) # -diag(n) allow to remove the identity matrix
min(cor.pearson.znf707o)

col_fun2 = colorRamp2(c(-1, 0, 1), c("#0f4259", "white", "#800020"))
heatmap_ZNF707o <- Heatmap(cor.pearson.znf707o, name = "Pearson correlation",
                          row_names_side = "left",
                          row_dend_side = "right",
                          column_names_side = "top",
                          column_names_rot = 45,
                          column_dend_side = "bottom",
                          row_dend_width = unit(30, "mm"),
                          column_dend_height = unit(30, "mm"),
                          column_dend_reorder = TRUE,
                          col = col_fun2,
                          rect_gp = gpar(col = "white", lwd = 1),
                          top_annotation = ha,
                          left_annotation = rowha)
pdf(file = "output/analysis/lincs/znf707/heatmap_pearson_KD_ZNF707_others_20190517.pdf", width = 25, height = 22)
print(heatmap_ZNF707o)
dev.off()

