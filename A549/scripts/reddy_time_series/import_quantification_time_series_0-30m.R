library(tidyverse)
library(stringr)
library(tools)

files <- dir("input/time_series_0-30m_tsv", pattern = "ENCFF", full.names = TRUE)

import_counts <- function(filename) {
    n <- file_path_sans_ext(filename) %>% basename
    suppressMessages(x <- read_tsv(filename, skip = 1))
    colnames(x) <- c("gene_id", "chr", "start", "end", "strand", "length", "count")
    x <- dplyr::select(x, gene_id, count) %>%
        mutate(gene_id = str_replace(gene_id, "\\.[0-9]+$", ""))
    colnames(x) <- c("gene_id", n)
    x
}
raw_counts_II <- Reduce("left_join", map(files, import_counts))
dir_results <- "results/a549_dex_time_points"

#####
# Import raw counts I

rac_counts_1








write_csv(raw_counts, paste0(dir_results, "/raw_counts_II.csv"))
