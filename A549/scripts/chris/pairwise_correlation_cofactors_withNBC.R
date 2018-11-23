setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
library(ef.utils)
library(plotly)

# load matrices: ctrl and dex (1 hour) conditions
load("output/analyses/cofactors/mat_cofactors.Rdata")

corr_ctrl <- pairwise_overlap(mat_cofactors_ctrl)
corr_dex <- pairwise_overlap(mat_cofactors_dex)

heatmap(corr_ctrl)
heatmap(corr_dex)

NIPBL_0m <- sort(corr_ctrl["NIPBL",], decreasing = TRUE)
BRD4_0m <- sort(corr_ctrl["BRD4",], decreasing = TRUE)
CDK9_0m <- sort(corr_ctrl["CDK9",], decreasing = TRUE)
NBC_0m <- sort(corr_ctrl["NBC",], decreasing = TRUE)
NBC_speCTRL_0m <- sort(corr_ctrl["NBC_speCTRL",], decreasing = TRUE)
NBC_common_0m <- sort(corr_ctrl["NBC_common",], decreasing = TRUE)
NBC_speDEX_0m <- sort(corr_ctrl["NBC_speDEX",], decreasing = TRUE)

NIPBL_1h <- sort(corr_dex["NIPBL",], decreasing = TRUE)
BRD4_1h <- sort(corr_dex["BRD4",], decreasing = TRUE)
CDK9_1h <- sort(corr_dex["CDK9",], decreasing = TRUE)
NBC_1h <- sort(corr_dex["NBC",], decreasing = TRUE)
NBC_speCTRL_1h <- sort(corr_dex["NBC_speCTRL",], decreasing = TRUE)
NBC_common_1h <- sort(corr_dex["NBC_common",], decreasing = TRUE)
NBC_speDEX_1h <- sort(corr_dex["NBC_speDEX",], decreasing = TRUE)

# check #27 LEVELPLOT WITH LATTICE