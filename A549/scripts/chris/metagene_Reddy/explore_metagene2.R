setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(tidyverse)
library(knitr)
library(metagene2)

##### Get the nth element each vector (from a list of vectors)
# Return a vector
get_nth_element <- function(lst, n) {
  sapply(lst, "[", n)
}

##### Transform "minute" to "m" and "hour" to "h"
format_timepoint <- function(timepoint) {
  if (base::grepl("hour", timepoint)) {
    gsub("hour", "h", timepoint)
  } else {
    gsub("minute", "m", timepoint)
  }
}

##### Make design from list of bam_files
# As all metadata are containing at each bam filename,
# design_metadata can be build from the list of bam_files
make_design_from_bam_list <- function(bam_list) {
  bam_names <- basename(bam_list)
  splitted <- strsplit(bam_names, split = "_")
  
  bam_names_without_ext <- strsplit(bam_names, split = "\\.") %>% get_nth_element(1)
  target <- get_nth_element(splitted, 1)
  target <- gsub("NR3C1", "GR", target)
  timepoint <- get_nth_element(splitted, 2) %>% sapply(format_timepoint, USE.NAMES = FALSE)
  timepoint <- factor(timepoint, levels = c("0m", "5m",  "10m",  "15m", "20m", "25m",
                                            "0h", "30m", "1h", "2h", "3h", "4h", "5h", "6h", "7h", "8h", "10h", "12h"))
  replicate <- get_nth_element(splitted, 3)
  
  design_metadata <- data.frame(design = bam_names_without_ext, target, timepoint, replicate,
                                stringsAsFactors = FALSE)
  
  return(design_metadata)
}

##### Generate a data.frame containing values for plotting
# Output for use as input in plot_metagene_Reddy()
make_df_metagene_Reddy <- function(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks) {
  chip_target <- gsub("GR", "NR3C1", chip_target)
  
  bigdf <- data.frame()
  for (target in chip_target) {
    bam_pattern <- paste0("(", target, ").*\\.bam$")
    bam_files <- list.files(path = bam_folder, pattern = bam_pattern, full.names = TRUE)
    design_meta <- make_design_from_bam_list(bam_files)
    mg <- metagene2$new(regions = peaks,
                        bam_files = bam_files,
                        assay = 'chipseq',
                        cores = 4)
    mg$add_metadata(design_metadata = design_meta)
    df <- mg$get_data_frame()
    
    bigdf <- rbind(bigdf, df)
  }
  bigdf$target <- factor(bigdf$target, levels = c("GR", "EP300", "H3K27ac", "JUN"))
  return(bigdf)
}

# Plot metagene from Reddy time course data
plot_metagene_Reddy <- function(df_metagene, customColors = c("#F5A623", "#4A90E2", "#008000"), title) {
  metagene_plot <- ggplot(df_metagene, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group=replicate)) +
    geom_ribbon(aes(fill = replicate), alpha = 0.3) +
    scale_fill_manual(values = customColors) +
    geom_line(aes(color = replicate), size = 0.5) + 
    scale_color_manual(values = customColors) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    theme(strip.text.y = element_text(angle = 0),
          strip.background = element_rect(colour = NA)) +
    ggtitle(title) +
    xlab("Distance in bins") +
    ylab("Mean coverage (raw)") +
    facet_grid(target ~ timepoint)
  return(metagene_plot)
}

####
saveMetagene <- function(metagene_plot, output_dir, output_file, width_val = 25, height_val = 22, format = "pdf") {
  output_filepath <- file.path(output_dir, paste0(output_file, ".", format))
  pdf(file = output_filepath, width = width_val, height = height_val)
  print(metagene_plot)
  dev.off()
  message(" > Metagene saved in ", output_filepath)
}

# Define one region at first
angptl4_peak <- GRanges("chr19", IRanges(8355327, 8356270))
il11_peak <- GRanges("chr19", IRanges(55372624, 55373186))
gapdh_peak <- GRanges("chr12", IRanges(6532244, 6532808))
mafk_peak <- GRanges("chr7", IRanges(1519343, 1522213))

# BAM
bam_folder <- "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam"

# mafk_metagene only
mafk_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = mafk_peak) 
mafk_plot <- plot_metagene_Reddy(mafk_df_metagene, title = "MAFK")
saveMetagene(metagene_plot = mafk_plot,
            output_dir = "output/analyses/metagene_reddyTimeCourse",
            output_file = "metagene_MAFK",
            format = "pdf",
            width = 23, height = 9)

angptl4_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = angptl4_peak) 
angptl4_plot <- plot_metagene_Reddy(angptl4_df_metagene, title = "ANGPTL4")
saveMetagene(metagene_plot = angptl4_plot,
             output_dir = "output/analyses/metagene_reddyTimeCourse",
             output_file = "metagene_ANGPTL4",
             format = "pdf",
             width = 23, height = 9)

il11_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = il11_peak) 
il11_plot <- plot_metagene_Reddy(il11_df_metagene, title = "IL11")
saveMetagene(metagene_plot = il11_plot,
             output_dir = "output/analyses/metagene_reddyTimeCourse",
             output_file = "metagene_IL11",
             format = "pdf",
             width = 23, height = 9)

gapdh_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = gapdh_peak) 
gapdh_plot <- plot_metagene_Reddy(gapdh_df_metagene, title = "GAPDH")
saveMetagene(metagene_plot = gapdh_plot,
             output_dir = "output/analyses/metagene_reddyTimeCourse",
             output_file = "metagene_GAPDH",
             format = "pdf",
             width = 21, height = 7)
