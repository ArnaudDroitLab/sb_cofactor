setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(dplyr)
library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(-2, 0, 2), c("#0f4259", "white", "#800020"))

# Manually downloaded : learn how to use API if we continue in that direction
GeneSymbol <- read.csv("input/iLINCS_KD/AFF4_sig_Sat_Apr_13_19_05_48_2019_323144.xls", sep="\t") %>% pull(Name_GeneSymbol) %>% as.character
AFF4 <- read.csv("input/iLINCS_KD/AFF4_sig_Sat_Apr_13_19_05_48_2019_323144.xls", sep="\t") %>% select(Value_LogDiffExp)
BRPF1 <- read.csv("input/iLINCS_KD/BRPF1_sig_Sat_Apr_13_19_39_44_2019_4075752.xls", sep="\t") %>% select(Value_LogDiffExp)
CREBBP <- read.csv("input/iLINCS_KD/CREBBP_sig_Sat_Apr_13_19_38_44_2019_5069733.xls", sep="\t") %>% select(Value_LogDiffExp)
EP300 <- read.csv("input/iLINCS_KD/EP300_sig_Sat_Apr_13_19_39_13_2019_5331544.xls", sep="\t") %>% select(Value_LogDiffExp)
HDAC4 <- read.csv("input/iLINCS_KD/HDAC4_sig_Sat_Apr_13_19_40_14_2019_8177815.xls", sep="\t") %>% select(Value_LogDiffExp)
HDAC8 <- read.csv("input/iLINCS_KD/HDAC8_sig_Sat_Apr_13_19_19_18_2019_8784530.xls", sep="\t") %>% select(Value_LogDiffExp)
KAT6A <- read.csv("input/iLINCS_KD/KAT6A_sig_Sat_Apr_13_19_40_37_2019_2482511.xls", sep="\t") %>% select(Value_LogDiffExp)
KAT6B <- read.csv("input/iLINCS_KD/KAT6B_sig_Sat_Apr_13_19_41_03_2019_6249354.xls", sep="\t") %>% select(Value_LogDiffExp)
NIPBL <- read.csv("input/iLINCS_KD/NIPBL_sig_Sat_Apr_13_19_20_13_2019_5014748.xls", sep="\t") %>% select(Value_LogDiffExp)
SMARCA2 <- read.csv("input/iLINCS_KD/SMARCA2_sig_Sat_Apr_13_19_41_56_2019_6550154.xls", sep="\t") %>% select(Value_LogDiffExp)
SMARCA4 <- read.csv("input/iLINCS_KD/SMARCA4_sig_Sat_Apr_13_19_20_40_2019_2127956.xls", sep="\t") %>% select(Value_LogDiffExp)
SMARCB1 <- read.csv("input/iLINCS_KD/SMARCB1_sig_Sat_Apr_13_19_42_20_2019_4245267.xls", sep="\t") %>% select(Value_LogDiffExp)
SMARCE1 <- read.csv("input/iLINCS_KD/SMARCE1_sig_Sat_Apr_13_19_42_43_2019_8850870.xls", sep="\t") %>% select(Value_LogDiffExp)
SMC3 <- read.csv("input/iLINCS_KD/SMC3_sig_Sat_Apr_13_19_21_41_2019_1818884.xls", sep="\t") %>% select(Value_LogDiffExp)

CTCF <- read.csv("input/iLINCS_KD/CTCF_sig_Sat_Apr_13_19_54_22_2019_8607779.xls", sep="\t") %>% select(Value_LogDiffExp)
POLR2A <- read.csv("input/iLINCS_KD/POLR2A_sig_Sat_Apr_13_19_56_41_2019_7748946.xls", sep="\t") %>% select(Value_LogDiffExp)
CHEK1 <- read.csv("input/iLINCS_KD/CHEK1_sig_Sun_Apr_14_20_24_16_2019_4192854.xls", sep="\t") %>% select(Value_LogDiffExp)
CCL2 <- read.csv("input/iLINCS_KD/CCL2_sig_Sun_Apr_14_20_24_11_2019_2966412.xls", sep="\t") %>% select(Value_LogDiffExp)
B2M <- read.csv("input/iLINCS_KD/B2M_sig_Sun_Apr_14_20_23_38_2019_8911239.xls", sep="\t") %>% select(Value_LogDiffExp)
ADAR <- read.csv("input/iLINCS_KD/ADAR_sig_Sun_Apr_14_20_22_51_2019_8870966.xls", sep="\t") %>% select(Value_LogDiffExp)
ABCA1 <- read.csv("input/iLINCS_KD/ABCA1_sig_Sun_Apr_14_20_21_56_2019_6692401.xls", sep="\t") %>% select(Value_LogDiffExp)
ABAT <- read.csv("input/iLINCS_KD/ABAT_sig_Sun_Apr_14_20_21_51_2019_9531899.xls", sep="\t") %>% select(Value_LogDiffExp)
AATF <- read.csv("input/iLINCS_KD/AATF_sig_Sun_Apr_14_20_21_48_2019_8154118.xls", sep="\t") %>% select(Value_LogDiffExp)
AARS <- read.csv("input/iLINCS_KD/AARS_sig_Sun_Apr_14_20_21_12_2019_3064136.xls", sep="\t") %>% select(Value_LogDiffExp)
A2M <- read.csv("input/iLINCS_KD/A2M_sig_Sun_Apr_14_20_20_19_2019_1039643.xls", sep="\t") %>% select(Value_LogDiffExp)

df <- data.frame(GeneSymbol, AFF4, BRPF1, CREBBP, EP300, HDAC4,
                 HDAC8, KAT6A, KAT6B, NIPBL, SMARCA2, 
                 SMARCA4, SMARCB1, SMARCE1, SMC3, CTCF, POLR2A,
                 CHEK1, CCL2, B2M, ADAR, ABCA1, ABAT, AATF, AARS, A2M)

colnames(df) <- c("Gene", "AFF4", "BRPF1", "CREBBP", "EP300", "HDAC4",
                  "HDAC8", "KAT6A", "KAT6B", "NIPBL", "SMARCA2", 
                  "SMARCA4", "SMARCB1", "SMARCE1", "SMC3", "CTCF", "POLR2A",
                  "CHEK1", "CCL2", "B2M", "ADAR", "ABCA1", "ABAT", "AATF", "AARS", "A2M")

rownames(df) <- GeneSymbol

mat <- as.matrix(df[,2:ncol(df)])

Heatmap(mat,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 0),
        row_dend_side = "right",
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 11),
        column_names_rot = 45,
        column_dend_side = "bottom",
        row_dend_width = unit(50, "mm"),
        column_dend_height = unit(50, "mm"),
        col = col_fun,
        rect_gp = gpar(col = "white", lwd = 0))
