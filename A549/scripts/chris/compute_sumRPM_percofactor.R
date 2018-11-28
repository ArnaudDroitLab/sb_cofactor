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

load_reduced_cofactor_peaks <- function(cofactor) {
  peaks_dir <- file.path("output/chip-pipeline-GRCh38/peak_call", paste0("A549_", cofactor))
  filename <- paste0("A549_", cofactor, "_reduced.bed")
  peaks_path <- file.path(peaks_dir, filename)
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

sumRPM <- function(countTable, cofactor) {
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
  sumRPM_value <- colSums(RPM_table)
  return(sumRPM_value)
}

###
output_path <- "output/analyses/sumRPM"
mkdir_cmd_line <- paste("mkdir", output_path, sep = " ")
system(mkdir_cmd_line)

cofactors <- c("NIPBL", "BRD4", "CDK9", "MED1", "SMC1A")
sumRPM_values <- c()
for (cofactor in cofactors) {
  peaks_reduced <- load_reduced_cofactor_peaks(cofactor)
  counts_reduced <- countReads(peaks_reduced[1:5], cofactor)
  write.table(counts_reduced, file.path(output_path, paste0("countTable_", cofactor, "_reduced.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)  
  sumRPM_reduced <- sumRPM(counts_reduced, cofactor)
  sumRPM_values <- c(sumRPM_values, sumRPM_reduced)
  message("##################")
  print(sumRPM_values)
  message("##################")
}

save(sumRPM_values, file = file.path(output_path, "sumRPM_values.RData"))

# 
# NIPBL_reduced <- load_reduced_cofactor_peaks("NIPBL")
# counts_NIPBL_reduced <- countReads(NIPBL_reduced[1:10], "NIPBL")
# sumRPM_NIPBL_reduced <- sumRPM(counts_NIPBL_reduced, "NIPBL")
# 
# ###
# output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_perRegion_perSamples"
# 
# speNBC_CTRL <- load_NBC_peaks("ctrl")
# counts_speNBC_CTRL <- countReads_perRegion(speNBC_CTRL)
# save(counts_speNBC_CTRL, file = file.path(output_path, "counts_speNBC_CTRL.RData"))
# write.table(counts_speNBC_CTRL, file = file.path(output_path, "countTable_speNBC_CTRL.txt"),
#             sep = "\t", quote = FALSE, row.names = FALSE)
# rm(counts_speNBC_CTRL)
# 
# NBC_common <- load_NBC_peaks("common")
# counts_NBC_common <- countReads_perRegion(NBC_common)
# save(counts_NBC_common, file = file.path(output_path, "counts_NBC_common.RData"))
# write.table(counts_NBC_common, file = file.path(output_path, "countTable_NBC_common.txt"),
#             sep = "\t", quote = FALSE, row.names = FALSE)
# rm(counts_NBC_common)
# 
# speNBC_DEX <- load_NBC_peaks("dex")
# counts_speNBC_DEX <- countReads_perRegion(speNBC_DEX)
# save(counts_speNBC_DEX, file = file.path(output_path, "counts_speNBC_DEX.RData"))
# write.table(counts_speNBC_DEX, file = file.path(output_path, "countTable_speNBC_DEX.txt"),
#             sep = "\t", quote = FALSE, row.names = FALSE)
# rm(counts_speNBC_DEX)
# ###
