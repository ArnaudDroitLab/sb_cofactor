# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

source("scripts/ckn_utils.R")
source("scripts/load_reddy.R")
library(ef.utils)
library(ChIPseeker)

# Loading peaks
peaks_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NB"
gainNB <- rtracklayer::import(con = file.path(peaks_dir, "NB_specific_DEX.bed")); print(length(gainNB)) # 803
lossNB <- rtracklayer::import(con = file.path(peaks_dir, "NB_specific_CTRL.bed")); print(length(lossNB)) # 5952
commonNB <- rtracklayer::import(con = file.path(peaks_dir, "NB_common.bed")); print(length(commonNB)) # 1680

# TFs

# tf_list <- c("NR3C1", "EP300", "CTCF", "SMC3", "RAD21", "FOSL2", "BCL3", "JUN", "JUNB", "HES2", "CEPBP")

tf_list <- c("NR3C1", "EP300")

# for (tf in tf_list) {
#   tf_regions <- load_reddy_binding_consensus(tf)
#   names(tf_regions)
# }

gr <- load_reddy_binding_consensus("GR")
names(gr)
# 1 hour

ep300 <- load_reddy_binding_consensus("EP300")
names(ep300)
# 0 hour / 0 minute
# 1 hour

ctcf <- load_reddy_binding_consensus("CTCF")
names(ctcf)
# 0 hour
# 1 hour

smc3 <- load_reddy_binding_consensus("SMC3")
names(smc3)
# 0 minute
# 1 hour

rad21 <- load_reddy_binding_consensus("RAD21")
names(rad21)
# 0 minute
# 1 hour

fosl2 <- load_reddy_binding_consensus("FOSL2")
names(fosl2)
# 0 hour
# 1 hour

bcl3 <- load_reddy_binding_consensus("BCL3")
names(bcl3)
# 0 hour
# 1 hour

jun <- load_reddy_binding_consensus("JUN")
names(jun)
# 0 hour / 0 minute
# 1 hour

junb <- load_reddy_binding_consensus("JUNB")
names(junb)
# 0 hour
# 1 hour

hes2 <- load_reddy_binding_consensus("HES2")
names(hes2)
# 0 hour
# 1 hour

cepbp <- load_reddy_binding_consensus("CEBPB")
names(cepbp)
# 0 hour
# 1 hour



