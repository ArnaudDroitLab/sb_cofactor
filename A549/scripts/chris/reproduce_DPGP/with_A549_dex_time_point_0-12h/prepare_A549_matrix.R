setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(ENCODExplorer)
library(tidyverse)
library(knitr)

# Aim: Being able to identify sample (time point + replicate) with the column name in the raw count table

# Retrieve information for RNA samples with ENCODExplorer
all_rnaseq <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_type="tsv", assay="polyA RNA-seq")
a549_dex_samples <- all_rnaseq %>% dplyr::filter(assembly == "GRCh38", lab == "Tim Reddy, Duke", submitted_by == "Ian McDowell")
report_tsv <- a549_dex_samples %>% dplyr::select(accession, file_accession, file_size, submitted_by, target, treatment_duration, treatment_duration_unit, biological_replicates, controls)
report_tsv$treatment_duration[is.na(report_tsv$treatment_duration)] <- 0
report_tsv$treatment_duration[report_tsv$treatment_duration == "30"] <- "0.5"
report_tsv$treatment_duration_unit <- "hour"
report_tsv <- report_tsv %>% arrange(treatment_duration, biological_replicates)
kable(report_tsv)

# Load raw counts and rename columns
raw_counts <- read_csv("results/a549_dex_time_points/raw_counts.csv")

new_colnames <- c()
for (ENCF in colnames(raw_counts)[2:length(colnames(raw_counts))]) {
  line <- report_tsv %>% dplyr::filter(file_accession == ENCF)
  column_name <- paste0(line$treatment_duration, "h_rep", line$biological_replicates)
  message(ENCF, " | ", column_name)
  new_colnames <- c(new_colnames, column_name)
}

colnames(raw_counts) <- c("gene_id", new_colnames)
raw_counts

write.table(raw_counts, file = "results/a549_dex_time_points/raw_counts_with_colnames.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)