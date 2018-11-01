setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(dplyr)

###### NIPBL Regions
nipbl_regions <- c("output/chip-pipeline-GRCh38/peak_call/A549_NIPBL/A549_NIPBL_CTRL_specific.bed",
                   "output/chip-pipeline-GRCh38/peak_call/A549_NIPBL/A549_NIPBL_common.bed",
                   "output/chip-pipeline-GRCh38/peak_call/A549_NIPBL/A549_NIPBL_DEX_specific.bed")
regions <- paste(nipbl_regions, collapse = " ")
region_labels <- "NIPBL_CTRL NIPBL_common NIPBL_DEX"

###### Samples etoh/dex
targets <- c("BCL3", "CEBPB", "CTCF")
etoh_rep <- c(3, 2, 3)
dex_rep <- c(3, 3, 3)
replicate_nb <- data.frame(targets, etoh_rep, dex_rep)

###### Command line
# 1. computeMatrix
# 2. plotHeatmap

output_dir <- "output/analyses/heatmap_NIPBL_vs_ENCODE_Reddy"
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
bigwig_dir <- "input/ENCODE/A549/GRCh38/chip-seq/bigWig"

generate_sample_path <- function(target, condition, nb_dex) {
  res <- c()
  for (i in 1:nb_dex) {
    filename <- paste0(target, "_", "rep", i, "_", condition, ".bigWig")
    sample_name <- file.path(bigwig_dir, filename)
    res <- c(res, sample_name)
  }
  return(res)
}

compute_matrix <- function(target, replicate_nb) {
  nb <- replicate_nb %>% filter(targets == target)
  nb_etoh <- nb$etoh_rep
  nb_dex <- nb$dex_rep
  
  samples_etoh <- generate_sample_path(target, condition = "etoh", nb_etoh)
  samples_dex <- generate_sample_path(target, condition = "dex_1h", nb_dex)
  sample_scaffold <- c(samples_etoh, samples_dex)
  samples <- paste(sample_scaffold, collapse = " ")
  
  output_path <- matrix_path
  
  cmd_line_scaffold <- c("computeMatrix reference-point --referencePoint center",
                          "--regionsFileName", regions,
                          "--scoreFileName", samples,
                          "--upstream", "1000", "--downstream", "1000", "-p", "8",
                          "--outFileName", output_path)
  
  cmd_line <- paste(cmd_line_scaffold, collapse = " ")
  message(cmd_line)
  system(cmd_line)
}

generate_sample_labels <- function(target, condition, nb_dex) {
  res <- c()
  for (i in 1:nb_dex) {
    label <- paste0(target, "_", "rep", i, "_", condition)
    res <- c(res, label)
  }
  return(res)
}

plot_heatmap <- function(target, replicate_nb) {
  nb <- replicate_nb %>% filter(targets == target)
  nb_etoh <- nb$etoh_rep
  nb_dex <- nb$dex_rep
  
  sample_labels_etoh <- generate_sample_labels(target, condition = "etoh", nb_etoh)
  sample_labels_dex <- generate_sample_labels(target, condition = "dex_1h", nb_dex)
  sample_labels_scaffold <- c(sample_labels_etoh, sample_labels_dex)
  sample_labels <- paste(sample_labels_scaffold, collapse = " ")
  
  output_path <- file.path(output_dir, paste0("nipbl_", target, "_heatmap.png"))
  
  cmd_line_scaffold <- c("plotHeatmap", "--matrixFile", matrix_path,
                         "--colorMap", "rainbow",
                         "--regionsLabel", region_labels,
                         "--samplesLabel", sample_labels,
                         "--outFileName", output_path)
  
  cmd_line <- paste(cmd_line_scaffold, collapse = " ")
  message(cmd_line)
  system(cmd_line)
}

### Main fonction
for (target in targets) {
  message("##########\t", target)
  matrix_path <- file.path(output_dir, paste0("nipbl_", target, "_matrix.gzip"))
  compute_matrix(target, replicate_nb)
  plot_heatmap(target, replicate_nb)
  }