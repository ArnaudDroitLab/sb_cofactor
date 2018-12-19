# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/ckn_utils.R")

####
oldpath <- Sys.getenv("PATH")
newpath <- paste(oldpath, ":/usr/bin/bin", sep='')
Sys.setenv(PATH = newpath)

#####
generate_coordVector <- function(granges) {
  df <- as.data.frame(granges)
  coord_Vector <- paste0(df$seqnames, ":", df$start, "-", df$end)
  return (coord_Vector)
}

generate_bampathVector <- function(cofactors, conditions) {
  bampath <- c()
  bam_dir <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment"
  for (condition in conditions) {
    for (cofactor in cofactors) {
      basename <- paste0("A549", "_", condition, "_", cofactor, "_rep1")
      path <- file.path(bam_dir, basename, paste0(basename, ".sorted.dup.bam"))
      bampath <- c(bampath, path)
    }
  }
  return(bampath)
}

load_specific_cofactor_peaks <- function(cofactor, region_set) {
  peaks_dir <- file.path("output/chip-pipeline-GRCh38/peak_call", paste0("A549_", cofactor))
  filename <- list(CTRL = paste0("A549_", cofactor, "_CTRL_specific.bed"),
                   common = paste0("A549_", cofactor, "_common.bed"),
                   DEX = paste0("A549_", cofactor, "_DEX_specific.bed"))
  peaks_path <- file.path(peaks_dir, filename[[region_set]])
  message(peaks_path)
  peaks <- rtracklayer::import(peaks_path)
  return(peaks)
}

#####
countReads <- function(peaks_set, cofactor) {
  bam_path <- generate_bampathVector(cofactors = cofactor,
                                     conditions = c("CTRL", "DEX"))
  
  peaks_coordVector <- generate_coordVector(peaks_set)
  
  count_total <- data.frame(peaks_coordVector)
  for (bam in bam_path) {
    count_bySample <- c()
    for (region in peaks_coordVector) {
      cmd_line <- paste("samtools view -c", bam, region, sep = " ")
      count <- as.integer(system(cmd_line, intern = TRUE))
      count_bySample <- c(count_bySample, count)
      cat(region, count, "| ")
    }
    cat("\n")
    count_total <- data.frame(count_total, count_bySample)
  }
  names(count_total) <- c("Coordinates", paste0(cofactor, "_CTRL"), paste0(cofactor, "_DEX"))
  return(count_total)
}

getRPM_table <- function(countTable, cofactor) {
  bam_path <- generate_bampathVector(cofactors = cofactor,
                                     conditions = c("CTRL", "DEX"))
  total_read <- c()
  for (bam in bam_path) {
    cmd_line <- paste("samtools view -c", bam, sep = " ")
    cat(cmd_line, "\n")
    total <- as.integer(system(cmd_line, intern = TRUE))
    total_read <- c(total_read, total)
    cat(bam, total, "\n")
  }
  RPM_table <- (countTable[, 2:3] * 1000000) / total_read
  return(RPM_table)
}

sumRPM <- function(RPM_table, cofactor) {
  sumRPM_value <- colSums(RPM_table)
  return(sumRPM_value)
}

###
output_path <- "output/analyses/sumRPM_specific_without_common"
mkdir_cmd_line <- paste("mkdir", output_path, sep = " ")
system(mkdir_cmd_line)

cofactors <- c("NIPBL", "BRD4", "CDK9", "MED1", "SMC1A")
# cofactors <- c("NIPBL", "BRD4")

sumRPM_values <- c()
for (cofactor in cofactors) {
  for (region_set in c("CTRL", "DEX")) {
    message("##########\t", cofactor, " | ", region_set)
    peaks_cond_specific <- load_specific_cofactor_peaks(cofactor, region_set)
    message(length(peaks_cond_specific))
    # peaks_common <- load_specific_cofactor_peaks(cofactor, "common")
    # message(length(peaks_common))
    # peaks <- GenomicRanges::reduce(c(peaks_cond_specific, peaks_common))
    # message(length(peaks))
    
    counts <- countReads(peaks_cond_specific, cofactor)
    write.table(counts, file.path(output_path, paste0("countTable_", cofactor, "_", region_set, "_specific.txt")),
                 sep = "\t", quote = FALSE, row.names = FALSE)

    # counts_reduced <- read.table(file = file.path(output_path, paste0("countTable_", cofactor, "_specific.txt")),
    #                            header = TRUE)

     RPM_table <- getRPM_table(counts, cofactor)
  
   ### boxplot
     boxplot_filename <- file.path(output_path, paste0("RPM_boxplot_", cofactor, "_", region_set, "_specific.png"))
     png(boxplot_filename)
     boxplot(RPM_table, outline = FALSE)
     dev.off()
  
     sumRPM_specific <- sumRPM(RPM_table, cofactor)
  
     sumRPM_values <- c(sumRPM_values, sumRPM_specific)
     message("##################")
     print(sumRPM_values)
     message("##################")
  }
}

save(sumRPM_values, file = file.path(output_path, "sumRPM_specific_values_without_common.RData"))