library(ENCODExplorer)
library(tidyverse)

encode_df <- get_encode_df()

time_point <- c("0m", "5m", "10m",
                "15m", "20m", "25m")

accession <- c("ENCSR937WIG", "ENCSR482TZY", "ENCSR742VGF",
               "ENCSR964GKZ", "ENCSR831FZM", "ENCSR525HSH")

design <- data.frame(accession = accession, time_point = time_point)

sample_sheet <- dplyr::filter(encode_df, accession %in% design$accession) %>%
    dplyr::filter(file_type == "tsv") %>%
    dplyr::select(accession, file_accession, href) %>%
    mutate(href = paste0("https://www.encodeproject.org", href)) %>%
    left_join(design, by = "accession")

write_csv(sample_sheet, "input/sample_sheet_time_series_0-30m.csv")
