source("https://bioconductor.org/biocLite.R")
setwd("/home/chris/Bureau/A549_script")
library(org.Hs.eg.db)

###################
# Infos on org.Hs.eg.db library
# columns(org.Hs.eg.db)
# keytypes(org.Hs.eg.db)

###################
# Retrieve mapIds: get SYMBOL from ENSEMBL IDs


raw_counts <- read.csv("raw_counts.csv", sep=",", header=TRUE, row.names=1)
ensembl_ids <- rownames(raw_counts)
symbol <- mapIds(org.Hs.eg.db, keys=ensembl_ids, column="SYMBOL", keytype="ENSEMBL")
map <- as.data.frame(symbol)

nb_NA <- sum(is.na(map$symbol))
#  35163 NAs ==» 58.137 % of ENSEMBL IDS has no SYMBOLs
# total: 60483 ENSEMBL IDS

###################

write.csv(map, file="map_ensembl_symbol.csv", row.names=TRUE, quote=FALSE)
