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
  report_bam <- bam %>% dplyr::select(accession, file_accession, file_size, submitted_by, target, treatment_duration, treatment_duration_unit, biological_replicates, controls)
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

make_report_WCE_bam <- function(report_bam, all_chip_bam){
  accession_WCE_list <- unique(report_bam$controls)
  wce_bam <- all_chip_bam[all_chip_bam$accession %in% accession_WCE_list, ] %>%
    arrange(factor(accession, levels = accession_WCE_list), biological_replicates)
  report_wce_bam <- wce_bam %>% dplyr::select(controls = accession, file_accession, file_size, submitted_by, target, treatment_duration, treatment_duration_unit, biological_replicates)
  report_wce_bam$target <- paste0(unique(report_bam$target), "-WCE")
  for (wce_acc in accession_WCE_list) {
    tmp <- report_bam[report_bam$controls == wce_acc, ]
    tmp_td <- unique(tmp$treatment_duration)
    tmp_tdu <- unique(tmp$treatment_duration_unit)
    report_wce_bam[which(report_wce_bam$controls == wce_acc), ]$treatment_duration <- tmp_td
    report_wce_bam[which(report_wce_bam$controls == wce_acc), ]$treatment_duration_unit <- tmp_tdu 
  }
  print(kable(report_wce_bam))
  return(report_wce_bam)
}

countRead <- function(regions, peaks_coordVector, bamPath, report_bam) {
  which <- regions
  what <- c("rname", "strand", "pos", "qwidth", "seq")
  param <- ScanBamParam(which=which, what=what)
  
  count_total <- data.frame(peaks_coordVector)
  cpt <- 0
  cptmax <- length(bamPath)
  for (bam in bamPath) {
    cpt <- cpt + 1
    message("##### ", cpt, " / ", cptmax, "\n", "# ", bam)
    index_bam <- gsub("\\.bam", "", bam)
    count_bySample <- countBam(bam, index = index_bam, param = param)
    records <- count_bySample$records
    count_total <- data.frame(count_total, records)
  }
  names(count_total) <- c("Coordinates", paste0(report_bam$target, "_", report_bam$treatment_duration, report_bam$treatment_duration_unit,
                                                "_rep", report_bam$biological_replicates,
                                                "_", report_bam$file_accession))
  return(count_total)
}

##### 2019-01-29: Deprecated | use indexBam from Rsamtools instead > do not require samtools to be in path
##### index BAM (samtools must be in path)
# index_bam <- function(bamPath) {
#   for (bam in bamPath) {
#     bai <- gsub("\\.bam", ".bai", bam)
#     cmd_line <- paste("samtools index", bam, bai)
#     cat(cmd_line, "\n")
#     system(cmd_line)
#   }
# }

##### Figures
retrieve_sumcount <- function(countTable_filename) {
  input_path <- "output/analyses/countTable_overtime"
  countTable <- read.table(file.path(input_path, countTable_filename), header = TRUE)
  rownames(countTable) <- countTable$Coordinates
  sumcount <- colSums(countTable[, 2:ncol(countTable)])
  return(sumcount)
}

# Assigned colnames relative to time in minutes (to have proper scaling on x axis)
remove_unit <- function(time_wunit) {
  # remove unit
  # convert minute in hour
  time <- gsub("minute", "", time_wunit)
  if (grepl("hour", time)) {
    time <- gsub("hour", "", time_wunit)
    time <- as.integer(time) * 60
  }
  return(time)
}

assignTime <- function(ratioTable) {
  oldColnames <- names(ratioTable)
  lsplit <- strsplit(names(ratioTable), split = "_")
  time_wunit <- sapply(lsplit, function(x) x[2])
  newColnames <- sapply(time_wunit, remove_unit)
  names(ratioTable) <- newColnames
  return(ratioTable)
}

dfForPlot <- function(scaledRatioTable, protein) {
  proteinVector <- rep(protein, length(scaledRatioTable))
  df <- data.frame(Time = names(scaledRatioTable), Protein = proteinVector, Value = scaledRatioTable)
  df$Time <- as.character(df$Time)
  df$Time <- as.numeric(df$Time)
  return(df)
}

bigdfForPlot <- function(valuesList, protein) {
  bigdf <- data.frame()
  cpt <- 0
  for (values in valuesList) {
    cpt <- cpt + 1
    df <- dfForPlot(values, protein[cpt])
    bigdf <- rbind(bigdf, df)
  }
  return(bigdf)
}

calculate_mean <- function(bigdf) {
  df_mean <- aggregate(bigdf$Value, by=list(bigdf$Time, bigdf$Protein), mean)
  colnames(df_mean) <- c("time", "protein", "mean")
  return(df_mean)
}

plotReadCount <- function(bigdf) {
  df_mean <- calculate_mean(bigdf)
  
  plot <- ggplot(bigdf, aes(x = Time, y = Value)) +
    geom_point(aes(color = Protein), size = 1) +
    geom_line(data=df_mean, aes(x=time, y=mean, group=protein, color=protein)) +
    scale_x_continuous(name="Time",
                       labels = c("0h", "", "", "", "", "", "30m", "1h", "2h", "3h", "4h", "5h", "6h", "7h", "8h", "10h", "12h"),
                       breaks = c(0, 5, 10, 15, 20, 25, 30, 60, 120, 180, 240, 300, 360, 420, 480, 600, 720)) +
    # scale_x_continuous(name="Time", labels=c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12), breaks=c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12)) +
    # scale_y_continuous(name="Value", labels=comma)
    scale_y_continuous(name="Value")
  return(plot)
}