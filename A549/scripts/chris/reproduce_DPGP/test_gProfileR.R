setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(gProfileR)
library(tidyverse)

# Retrieve gene_id in ENSG format
raw <- read_tsv("results/a549_dex_time_points/raw_counts_with_colnames.txt")
gene_id <- raw$gene_id

#
q <- gconvert(query = gene_id[1:100], organism = "hsapiens", filter_na = FALSE)
names(q)

#
gostres <- gprofiler(query = gene_id[1:100], organism = "hsapiens", png_fn = t)