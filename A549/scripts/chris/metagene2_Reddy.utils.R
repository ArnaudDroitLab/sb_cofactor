library(metagene2)

##### Return a vector
get_nth_element <- function(lst, n) {
  sapply(lst, "[", n)
}

#####
format_timepoint <- function(timepoint) {
  if (base::grepl("hour", timepoint)) {
    gsub("hour", "h", timepoint)
  } else {
    gsub("minute", "m", timepoint)
  }
}

#####
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

#####
make_design_from_bam_list <- function(bam_list, target) {
  bam_names <- basename(bam_list)
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
  
  design <- data.frame(bam_list)
  
  for (bam in group_bam_names2) {
    # message(bam)
    TF <- grepl(bam, row_design_names) %>% as.numeric
    design <- cbind(design, TF)
  }
  
  rownames(design) <- row_design_names
  colnames(design) <- c("Samples", group_bam_names2)
  
  return(design)
}

#####
make_metadata_from_design <- function(design) {
  nb_col <- length(colnames(design))
  design <- colnames(design)[2:nb_col]
  
  splitted <- sapply(design, strsplit, split = "_")
  target <- splitted %>% purrr::map(1) %>% unlist
  timepoint <- splitted %>% purrr::map(2) %>% unlist
  timepoint <- factor(timepoint, levels = c("0m", "5m",  "10m",  "15m", "20m", "25m",
                                            "0h", "30m", "1h", "2h", "3h", "4h", "5h", "6h", "7h", "8h", "10h", "12h"))
  replicate <- rep("merged", length(design))
  
  df_design <- data.frame(design, target, timepoint, replicate)
  df_design$target <- gsub("NR3C1", "GR", df_design$target)
  return(df_design)
}

#####
make_df_metagene_Reddy <- function(chip_target = c("GR", "EP300", "H3K27ac", "JUN"), peaks, merge_replicates = FALSE, reps = "123") {
  chip_target <- gsub("GR", "NR3C1", chip_target)
  
  bigdf <- data.frame()
  for (target in chip_target) {
    bam_pattern <- paste0("^", target, "_([0-9]+.*)_rep([", reps, "])_(.*\\.bam$)")
    bam_files <- list.files(path = bam_folder, pattern = bam_pattern, full.names = TRUE)
    # bam_files <- bam_files %>% .[matches("rep[12]", vars = .)] # allow to filter replicate
    metadata <- make_metadata_from_bam_list(bam_files)
    mg <- metagene2$new(regions = peaks,
                        bam_files = bam_files,
                        normalization = "RPM",
                        force_seqlevels = TRUE,
                        assay = 'chipseq',
                        cores = 4)
    mg$add_metadata(design_metadata = metadata)
    df <- mg$get_data_frame()
    
    bigdf <- rbind(bigdf, df)
    
    if (merge_replicates) {
      message("   > Gathering replicates... | ", target)
      design_true <- make_design_from_bam_list(bam_files, target)
      metadata_merge_replicates <- make_metadata_from_design(design_true)
      mg$produce_metagene(design = design_true, design_metadata = metadata_merge_replicates)
      df <- mg$get_data_frame()
      bigdf <- rbind(bigdf, df)
      bigdf$replicate <- factor(bigdf$replicate, levels = c("rep1", "rep2", "rep3", "merged"))
    }
  }
  chip_target <- gsub("NR3C1", "GR", chip_target)
  bigdf$target <- factor(bigdf$target, levels = chip_target)
  return(bigdf)
}

#####
plot_metagene_Reddy <- function(df_metagene, customColors = c("#F5A623", "#4A90E2", "#008000", "#8C001A"), title = "") {
  metagene_plot <- ggplot(df_metagene, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group=replicate)) +
    geom_ribbon(aes(fill = replicate), alpha = 0.3) +
    scale_fill_manual(values = customColors) +
    geom_line(aes(color = replicate, size = replicate, alpha = replicate)) + 
    scale_color_manual(values = customColors) +
    scale_size_manual(values = c(0.3, 0.3, 0.3, 1.3)) +
    scale_alpha_manual(values = c(0.6, 0.6, 0.6, 1)) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    theme(strip.text.x = element_text(size = 15, face = "bold"),
          strip.text.y = element_text(size = 15, face = "bold", angle = 0, hjust = 0),
          strip.background = element_rect(colour = NA)) +
    ggtitle(title) + theme(plot.title = element_text(size = 20, face = "bold")) +
    xlab("Distance in bins") +
    ylab("RPM") +
    facet_grid(target ~ timepoint)
  return(metagene_plot)
}

#####
saveMetagene <- function(metagene_plot, output_dir, output_file, width_val = 25, height_val = 22) {
  output_filepath <- file.path(output_dir, paste0(output_file, ".pdf"))
  pdf(file = output_filepath, width = width_val, height = height_val)
  print(metagene_plot)
  dev.off()
  message(" > Metagene saved in ", output_filepath)
}
