# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(knitr)
library(DiffBind)
source("scripts/chris/diffbind_GR_0_25m/diffbind.utils.R")

#####
bam_folder <- "input/ENCODE/A549/GRCh38/chip-seq/bam"
bed_folder <- "input/ENCODE/A549/GRCh38/chip-seq/narrow"
sSheet_GR_123 <- build_sSheet("NR3C1", bam_folder, bed_folder, reps = "123")
sSheet_GR_12 <- build_sSheet("NR3C1", bam_folder, bed_folder, reps = "12")
sSheet_EP300 <- build_sSheet("EP300", bam_folder, bed_folder, reps = "123")

##### Perform diffbind for GR and EP300 during 0-25m
# timepoint <- c("0h", "5m", "10m", "15m", "20m", "25m")
timepoint <- c("0m", "5m", "10m", "15m", "20m", "25m")

ltp <- length(timepoint)

# for (i in 1) {
for (i in 1:(ltp-1)) {
  for (j in (i+1):ltp) {
    tp1 <- timepoint[i]
    tp2 <- timepoint[j]
    # perform_diffbind("GR", sSheet_GR, tp1, tp2, reps = "123", output_dir = "output/analyses/GR_diffbind")
    # perform_diffbind("GR", sSheet_GR, tp1, tp2, reps = "12", output_dir = "output/analyses/GR_diffbind")
    # perform_diffbind("EP300", sSheet_EP300, tp1, tp2, reps = "123", output_dir = "output/analyses/EP300_diffbind")
    
    # open_diffBind("GR", tp1, tp2, pval = TF, reps = "123", output_dir = "output/analyses/GR_diffbind")
    # open_diffBind("GR", tp1, tp2, pval = TF, reps = "12", output_dir = "output/analyses/GR_diffbind")
    # open_diffBind("EP300", tp1, tp2, pval = TF, reps = "123", output_dir = "output/analyses/EP300_diffbind")
  }
}

##### Perform diffbind for EP300 0h 1h
# perform_diffbind("EP300", sSheet_EP300, "0h", "1h", reps = "123", output_dir = "output/analyses/EP300_diffbind")
EP300_DOWN <- open_diffBind("EP300", tp1 = "0h", tp2 = "1h", pval = TRUE, reps = "123", output_dir = "output/analyses/EP300_diffbind") %>%
  dplyr::filter(Fold < 0) %>% makeGRangesFromDataFrame

##### Perform diffbind for EP300 0h 0m
# perform_diffbind("EP300", sSheet_EP300, "0h", "0m", reps = "123", output_dir = "output/analyses/EP300_diffbind")
open_diffBind("EP300", tp1 = "0h", tp2 = "0m", pval = TRUE, reps = "123", output_dir = "output/analyses/EP300_diffbind")