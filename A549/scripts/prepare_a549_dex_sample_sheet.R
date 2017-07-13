library(ENCODExplorer)

accession <- c( "ENCSR070IYV", "ENCSR154TDP", "ENCSR224PTG", "ENCSR255VBV",
               "ENCSR326PTW", "ENCSR385TFN", "ENCSR543QRU", "ENCSR546PPG",
               "ENCSR624RID", "ENCSR632DQP", "ENCSR656FIH", "ENCSR924BHF")

time_point <- c( "0.5h", "12h", "10h", "6h", "4h", "7h", "5h", "8h", "3h",
                "0h", "1h", "2h")

design <- data.frame(accession = accession, time_point = time_point)

sample_sheet <- filter(encode_df, accession %in% design$accession) %>%
    filter(file_type == "tsv") %>%
    dplyr::select(accession, file_accession, href) %>%
    mutate(href = paste0("https://www.encodeproject.org", href)) %>%
    left_join(design, by = "accession")

write_csv(sample_sheet, "input/sample_sheet_ENCSR897XFT.csv")
