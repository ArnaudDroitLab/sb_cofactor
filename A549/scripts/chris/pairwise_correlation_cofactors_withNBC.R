setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
library(ef.utils)
library(plotly)

# load matrices: ctrl and dex (1 hour) conditions
load("output/analyses/cofactors/mat_cofactors.Rdata")

# output_path
output_path <- "output/analyses/cofactors"

corr_ctrl <- pairwise_overlap(mat_cofactors_ctrl)
corr_dex <- pairwise_overlap(mat_cofactors_dex)

heatmap(corr_ctrl)
heatmap(corr_dex)

#### EtOH
NBC_0m <- sort(corr_ctrl["NBC",], decreasing = TRUE)
NIPBL_0m <- sort(corr_ctrl["NIPBL",], decreasing = TRUE)
BRD4_0m <- sort(corr_ctrl["BRD4",], decreasing = TRUE)
CDK9_0m <- sort(corr_ctrl["CDK9",], decreasing = TRUE)


NBC_speCTRL_0m <- sort(corr_ctrl["NBC_speCTRL",], decreasing = TRUE)
NBC_common_0m <- sort(corr_ctrl["NBC_common",], decreasing = TRUE)
NBC_speDEX_0m <- sort(corr_ctrl["NBC_speDEX",], decreasing = TRUE)

overlaps_etoh.1 <- data.frame(names(NBC_0m), NBC_0m,
                            names(NIPBL_0m), NIPBL_0m,
                            names(BRD4_0m), BRD4_0m,
                            names(CDK9_0m), CDK9_0m)

overlaps_etoh.2 <- data.frame(names(NBC_speCTRL_0m), NBC_speCTRL_0m,
                              names(NBC_common_0m), NBC_common_0m,
                              names(NBC_speDEX_0m), NBC_speDEX_0m)

write.table(overlaps_etoh.1, file = file.path(output_path, "overlaps_cofactors_etoh_1.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(overlaps_etoh.2, file = file.path(output_path, "overlaps_cofactors_etoh_2.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

#### DEX 1h
NIPBL_1h <- sort(corr_dex["NIPBL",], decreasing = TRUE)
BRD4_1h <- sort(corr_dex["BRD4",], decreasing = TRUE)
CDK9_1h <- sort(corr_dex["CDK9",], decreasing = TRUE)
NBC_1h <- sort(corr_dex["NBC",], decreasing = TRUE)

NBC_speCTRL_1h <- sort(corr_dex["NBC_speCTRL",], decreasing = TRUE)
NBC_common_1h <- sort(corr_dex["NBC_common",], decreasing = TRUE)
NBC_speDEX_1h <- sort(corr_dex["NBC_speDEX",], decreasing = TRUE)

overlaps_dex.1 <- data.frame(names(NBC_1h), NBC_1h,
                              names(NIPBL_1h), NIPBL_1h,
                              names(BRD4_1h), BRD4_1h,
                              names(CDK9_1h), CDK9_1h)

overlaps_dex.2 <- data.frame(names(NBC_speCTRL_1h), NBC_speCTRL_1h,
                              names(NBC_common_1h), NBC_common_1h,
                              names(NBC_speDEX_1h), NBC_speDEX_1h)

write.table(overlaps_dex.1, file = file.path(output_path, "overlaps_cofactors_dex_1.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(overlaps_dex.2, file = file.path(output_path, "overlaps_cofactors_dex_2.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
#### attempt to data viz
# # check #27 LEVELPLOT WITH LATTICE
# NIPBL_0m
# NIPBL_0m_df <- data.frame(main_cofactor = rep("NIPBL", length(names(NIPBL_0m))), cofactors = names(NIPBL_0m), overlaps = NIPBL_0m)
# levelplot(overlaps ~ main_cofactor*cofactors, data = NIPBL_0m_df,
#           col.regions = heat.colors(100)[length(heat.colors(100)):1])
# 
# #### test lattice
# 
# library("lattice")
# 
# ## Example data
# x <- seq(1,10, length.out=20)
# y <- seq(1,10, length.out=20)
# data <- expand.grid(X=x, Y=y)
# data$Z <- runif(400, 0, 5)
# 
# ## Try it out
# par(mar=c(3,4,2,2))
# levelplot(Z ~ X*Y, data=data, xlab="X" , col.regions = heat.colors(100)[length(heat.colors(100)):1], main="")


