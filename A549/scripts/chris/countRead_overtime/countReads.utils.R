library(GenomicRanges)

generate_coordVector <- function(granges) {
  df <- as.data.frame(granges)
  coord_Vector <- paste0(df$seqnames, ":", df$start, "-", df$end)
  return (coord_Vector)
}

get_reduced_peaks_from_list <- function(peak_list) {
  all_peaks <- unlist(GRangesList(peak_list))
  all_peaks_reduced <- reduce(all_peaks)
  return(all_peaks_reduced)
}

make_report_bam <- function(target_name, all_chip_bam) {
  bam <- all_chip_bam %>% filter(target == target_name, assembly == "GRCh38", lab == "Tim Reddy, Duke")
  report_bam <- bam %>% select(accession, file_accession, file_size, submitted_by, target, treatment_duration, treatment_duration_unit, biological_replicates, controls)
  report_bam$treatment_duration[is.na(report_bam$treatment_duration)] <- 0
  report_bam$treatment_duration_unit[is.na(report_bam$treatment_duration_unit)] <- "minute"
  report_bam <- report_bam %>% arrange(desc(treatment_duration_unit), treatment_duration, submitted_by, biological_replicates)
  print(kable(report_bam))
  return(report_bam)
}

generate_bamPath_from_report <- function(report_bam, folder) {
  bamPath <- paste0(folder,
                     paste0(report_bam$target, "_", report_bam$treatment_duration, report_bam$treatment_duration_unit,
                            "_rep", report_bam$biological_replicates,
                            "_", report_bam$file_accession, ".bam"))
  return(bamPath)
}

countRead <- function(regions, peaks_coordVector, bamPath, report_bam) {
  which <- regions
  what <- c("rname", "strand", "pos", "qwidth", "seq")
  param <- ScanBamParam(which=which, what=what)
  
  count_total <- data.frame(peaks_coordVector)
  for (bam in bamPath) {
    message(bam)
    index_bam <- gsub("\\.bam", "", bam)
    count_bySample <- countBam(bam, index = index_bam, param = param)
    records <- count_bySample$records
    count_total <- data.frame(count_tota, records)
  }
  names(count_total) <- c("Coordinates", paste0(report_bam$target, "_", report_bam$treatment_duration, report_bam$treatment_duration_unit,
                                                "_rep", report_bam$biological_replicates,
                                                "_", report_bam$file_accession))
  return(count_total)
}