# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

options(readr.num_columns = 0)

timepoint <- c(0, 1, 4, 8, 12)

tad_dir <- "input/ENCODE/A549/GRCh38/hic_TAD"

tads <- list()
for (time in timepoint) {
  tad_name <- paste0(time, "h")
  message("###### ", tad_name)
  
  tad_filename <- paste0("hic_dex.t", time, ".tads.bedpe")
  tad_raw <- read_delim(file = file.path(tad_dir, tad_filename), delim = "\t")[, 1:3]
  colnames(tad_raw) <- c("seqnames", "start", "end")
  
  tad <- GRanges(tad_raw)
  print(summary(width(tad)))
  
  tads[[tad_name]] <- tad
  
  message("    # Number of TADs : ", length(tad))
}

all_tads_raw <- unlist(GRangesList(tads)) # 33041 TADs in total
all_tads_reduced <- GenomicRanges::reduce(all_tads_raw) # 4949
all_tads_withoutgap <- GenomicRanges::reduce(all_tads_reduced, min.gapwidth = 25000) # 2916
summary(width(all_tads_withoutgap))





