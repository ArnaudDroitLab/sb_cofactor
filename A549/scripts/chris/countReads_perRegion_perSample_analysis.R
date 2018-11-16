setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/ckn_utils.R")


corr(as.matrix(count_total$NIPBL_CTRL, count_total$BRD4_CTRL))
# faire une PCA?
# correlation entre ces deux listes de chiffres
# faire une PCA