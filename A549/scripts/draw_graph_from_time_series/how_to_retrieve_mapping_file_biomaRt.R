source("https://bioconductor.org/biocLite.R")
setwd("/home/chris/Bureau/A549_script")
library(biomaRt)


# How to retrieve mapping file: Ensemble ID - HGNC symbol

listMarts()    # to see which database options are present

ensembl=useMart("ensembl")  # using ensembl database data

listDatasets(ensembl)     # function to see which datasets are present in ensembl

ensembl=useDataset("hsapiens_gene_ensembl", mart=ensembl)   # from ensembl using homosapien gene data

listFilters(ensembl)  # check which filters are available

listAttributes(ensembl) # check attributes are available to select.More information on ensembl data base
# 
# genes.with.id=getBM(attributes=c("ensembl_gene_id", "external_gene_id"),values=gene_names, mart= ensembl) # function to get  gene id's and gene name from data base
#
mapping <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), mart=ensembl)
mapping[mapping == ""] <- NA

nb_NA <- sum(is.na(mapping$hgnc_symbol))

#  22582 NAs ==» 35.26839 % of ENSEMBL IDS has no SYMBOLs
# total: 64029 ENSEMBL IDS

write.csv(mapping, file="map_ensembl_hgnc.csv", row.names=FALSE, quote=FALSE)
