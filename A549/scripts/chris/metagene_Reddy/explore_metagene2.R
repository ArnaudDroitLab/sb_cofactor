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

##### Make metadata from list of bam_files
# As all metadata are containing at each bam filename,
# design_metadata can be build from the list of bam_files
make_metadata_from_bam_list <- function(bam_list) {
  bam_names <- basename(bam_list)
  splitted <- strsplit(bam_names, split = "_")
  
  bam_names_without_ext <- strsplit(bam_names, split = "\\.") %>% get_nth_element(1)
  target <- get_nth_element(splitted, 1)
  target <- gsub("NR3C1", "GR", target)
  timepoint <- get_nth_element(splitted, 2) %>% sapply(format_timepoint, USE.NAMES = FALSE)
  timepoint <- factor(timepoint, levels = c("0m", "5m",  "10m",  "15m", "20m", "25m",
                                            "0h", "30m", "1h", "2h", "3h", "4h", "5h", "6h", "7h", "8h", "10h", "12h"))
  replicate <- get_nth_element(splitted, 3)
  
  metadata <- data.frame(design = bam_names_without_ext, target, timepoint, replicate,
                                stringsAsFactors = FALSE)
  
  return(metadata)
}

##### Generate a data.frame containing values for plotting
# Output for use as input in plot_metagene_Reddy()
make_df_metagene_Reddy <- function(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks) {
  chip_target <- gsub("GR", "NR3C1", chip_target)
  
  bigdf <- data.frame()
  for (target in chip_target) {
    bam_pattern <- paste0("(", target, ").*\\.bam$")
    bam_files <- list.files(path = bam_folder, pattern = bam_pattern, full.names = TRUE)
    metadata <- make_metadata_from_bam_list(bam_files)
    mg <- metagene2$new(regions = peaks,
                        bam_files = bam_files,
                        assay = 'chipseq',
                        cores = 4)
    mg$add_metadata(design_metadata = metadata)
    df <- mg$get_data_frame()

    bigdf <- rbind(bigdf, df)
  }
  chip_target <- gsub("NR3C1", "GR", chip_target)
  bigdf$target <- factor(bigdf$target, levels = chip_target)
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
    theme(strip.text.x = element_text(size = 15, face = "bold"),
          strip.text.y = element_text(size = 15, face = "bold", angle = 0, hjust = 0),
          strip.background = element_rect(colour = NA)) +
    ggtitle(title) + theme(plot.title = element_text(size = 20, face = "bold")) +
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
# upreg genes
csf3_peak <- GRanges("chr17", IRanges(40021370, 40022011))
gpr1_peak <- GRanges("chr2", IRanges(206213349, 206213709))
socs1_peak <- GRanges("chr16", IRanges(11224153, 11224593))
angptl4_peak <- GRanges("chr19", IRanges(8355327, 8356270))

# downreg genes
dusp8_peak <- GRanges("chr11", IRanges(1559246, 1559846))
il11_peak <- GRanges("chr19", IRanges(55372624, 55373186))
mafk_peak <- GRanges("chr7", IRanges(1519343, 1522213))
rassf10_peak <- GRanges("chr11", IRanges(13007940, 13008299))

# control
gapdh_peak <- GRanges("chr12", IRanges(6532244, 6532808))

# BAM
bam_folder <- "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam"

# At one peak first
# upreg genes
csf3_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = csf3_peak) 
csf3_plot <- plot_metagene_Reddy(csf3_df_metagene, title = "CSF3 at chr17:40021370-40022011")
saveMetagene(metagene_plot = csf3_plot,
             output_dir = "output/analyses/metagene_reddyTimeCourse",
             output_file = "20190624_metagene_CSF3_with_replicates",
             format = "pdf",
             width = 21, height = 7)

gpr1_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = gpr1_peak) 
gpr1_plot <- plot_metagene_Reddy(gpr1_df_metagene, title = "GPR1 at chr2:206212349-206214709")
saveMetagene(metagene_plot = gpr1_plot,
             output_dir = "output/analyses/metagene_reddyTimeCourse",
             output_file = "20190624_metagene_GPR1_with_replicates",
             format = "pdf",
             width = 21, height = 7)

socs1_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = socs1_peak) 
socs1_plot <- plot_metagene_Reddy(socs1_df_metagene, title = "SOCS1 at chr2:206212349-206214709")
saveMetagene(metagene_plot = socs1_plot,
             output_dir = "output/analyses/metagene_reddyTimeCourse",
             output_file = "20190624_metagene_SOCS1_with_replicates",
             format = "pdf",
             width = 21, height = 7)


mafk_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = mafk_peak) 
mafk_plot <- plot_metagene_Reddy(mafk_df_metagene, title = "MAFK at chr7:1519343-1522213")
saveMetagene(metagene_plot = mafk_plot,
            output_dir = "output/analyses/metagene_reddyTimeCourse",
            output_file = "20190624_metagene_MAFK_with_replicates",
            format = "pdf",
            width = 23, height = 9)

angptl4_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = angptl4_peak) 
angptl4_plot <- plot_metagene_Reddy(angptl4_df_metagene, title = "ANGPTL4 at chr19:8355327-8356270")
saveMetagene(metagene_plot = angptl4_plot,
             output_dir = "output/analyses/metagene_reddyTimeCourse",
             output_file = "20190624_metagene_ANGPTL4_with_replicates",
             format = "pdf",
             width = 23, height = 9)

il11_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = il11_peak) 
il11_plot <- plot_metagene_Reddy(il11_df_metagene, title = "IL11 at chr19:55372624-55373186")
saveMetagene(metagene_plot = il11_plot,
             output_dir = "output/analyses/metagene_reddyTimeCourse",
             output_file = "20190624_metagene_IL11_with_replicates",
             format = "pdf",
             width = 23, height = 9)

gapdh_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = gapdh_peak) 
gapdh_plot <- plot_metagene_Reddy(gapdh_df_metagene, title = "GAPDH at chr12:6532244-6532808")
saveMetagene(metagene_plot = gapdh_plot,
             output_dir = "output/analyses/metagene_reddyTimeCourse",
             output_file = "20190624_metagene_GAPDH_with_replicates",
             format = "pdf",
             width = 21, height = 7)

dusp8_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = dusp8_peak) 
dusp8_plot <- plot_metagene_Reddy(dusp8_df_metagene, title = "DUSP8 at chr11:1559246-1559846")
saveMetagene(metagene_plot = dusp8_plot,
             output_dir = "output/analyses/metagene_reddyTimeCourse",
             output_file = "20190624_metagene_DUSP8_with_replicates",
             format = "pdf",
             width = 21, height = 7)

rassf10_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks = rassf10_peak) 
rassf10_plot <- plot_metagene_Reddy(rassf10_df_metagene, title = "RASSF10 at chr11:13007940-13008299")
saveMetagene(metagene_plot = rassf10_plot,
             output_dir = "output/analyses/metagene_reddyTimeCourse",
             output_file = "20190624_metagene_RASSF10_with_replicates",
             format = "pdf",
             width = 21, height = 7)












# TMP
# test at two regions
peaks2 <- list("listpeak1" = c(angptl4_peak, il11_peak),
               "listpeak1" = c(gapdh_peak, mafk_peak))
test_df_metagene <- make_df_metagene_Reddy(chip_target = c("GR", "EP300"), peaks = peaks2)
test_df <- plot_metagene_Reddy(test_df_metagene, title = "test_peaks2")
saveMetagene(metagene_plot = test_df,
             output_dir = "output/analyses/metagene_reddyTimeCourse",
             output_file = "test",
             format = "pdf",
             width = 21, height = 7)

## En cours d'écriture: integrer l'option: merge replicate
make_df_metagene_Reddy <- function(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks, merge_replicates = FALSE) {
  chip_target <- gsub("GR", "NR3C1", chip_target)
  
  bigdf <- data.frame()
  for (target in chip_target) {
    bam_pattern <- paste0("(", target, ").*\\.bam$")
    bam_files <- list.files(path = bam_folder, pattern = bam_pattern, full.names = TRUE)
    if (merge_replicates == TRUE) {
      design <- make_design_from_bam_list(bam_files, target)
    }
    
    metadata <- make_metadata_from_bam_list(bam_files)
    
    mg <- metagene2$new(regions = peaks,
                        bam_files = bam_files,
                        assay = 'chipseq',
                        cores = 4)
    mg$add_metadata(design_metadata = metadata)
    
    mg2 <- mg$group_coverage()
    # 

    # 
    # df <- mg$get_data_frame()
    # 
    # bigdf <- rbind(bigdf, df)
  }
  bigdf$target <- factor(bigdf$target, levels = chip_target)
  return(bigdf)
}


make_design_from_bam_list <- function(bam_list, target) {
  bam_names <- basename(bam_files)
  splitted <- strsplit(bam_names, split = "_")
  
  row_design_names <- paste(get_nth_element(splitted, 1),
                            sapply(get_nth_element(splitted, 2), format_timepoint, USE.NAMES = FALSE),
                            get_nth_element(splitted, 3),
                            sep = "_")
  
  group_bam_names <- paste(get_nth_element(splitted, 1),
                           sapply(get_nth_element(splitted, 2), format_timepoint, USE.NAMES = FALSE),
                           sep = "_") %>% unique
  
  timepoint <- c("0m", "5m",  "10m",  "15m", "20m", "25m",
                 "0h", "30m", "1h", "2h", "3h", "4h", "5h", "6h", "7h", "8h", "10h", "12h")
  ordered_levels <- paste(target, timepoint, sep = "_")
  
  group_bam_names2 <- group_bam_names[order(match(group_bam_names, ordered_levels))]
  
  design <- data.frame(bam_files)
  
  for (bam in group_bam_names2) {
    message(bam)
    TF <- grepl(bam, row_design_names) %>% as.numeric
    design <- cbind(design, TF)
  }
  
  rownames(design) <- row_design_names
  colnames(design) <- c("Samples", group_bam_names2)
  
  return(design)
}