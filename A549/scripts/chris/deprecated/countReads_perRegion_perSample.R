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

load_NBC_peaks <- function(regions_set) {
  peaks_NBC_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NBC"
  filename <- list(ctrl = "A549_NBC_CTRL_specific.bed",
                common = "A549_NBC_common.bed",
                dex = "A549_NBC_DEX_specific.bed")
  NBC_peaks_path <- file.path(peaks_NBC_dir, filename[[regions_set]])
  peaks <- rtracklayer::import(NBC_peaks_path)
  return(peaks)
}
#####
countReads_perRegion <- function(peaks_set) {
  bam_path <- generate_bampathVector(cofactors = c("NIPBL", "BRD4", "CDK9"),
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
    count_total <- data.frame(count_total, count_bySample)
  }
  names(count_total) <- c("Coordinates", "NIPBL_CTRL", "BRD4_CTRL", "CDK9_CTRL",
                          "NIPBL_DEX", "BRD4_DEX", "CDK9_DEX")
  return(count_total)
}

###
output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_perRegion_perSamples"

speNBC_CTRL <- load_NBC_peaks("ctrl")
counts_speNBC_CTRL <- countReads_perRegion(speNBC_CTRL)
save(counts_speNBC_CTRL, file = file.path(output_path, "counts_speNBC_CTRL.RData"))
write.table(counts_speNBC_CTRL, file = file.path(output_path, "countTable_speNBC_CTRL.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
rm(counts_speNBC_CTRL)

NBC_common <- load_NBC_peaks("common")
counts_NBC_common <- countReads_perRegion(NBC_common)
save(counts_NBC_common, file = file.path(output_path, "counts_NBC_common.RData"))
write.table(counts_NBC_common, file = file.path(output_path, "countTable_NBC_common.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
rm(counts_NBC_common)

speNBC_DEX <- load_NBC_peaks("dex")
counts_speNBC_DEX <- countReads_perRegion(speNBC_DEX)
save(counts_speNBC_DEX, file = file.path(output_path, "counts_speNBC_DEX.RData"))
write.table(counts_speNBC_DEX, file = file.path(output_path, "countTable_speNBC_DEX.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
rm(counts_speNBC_DEX)
