# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(metagene2)
library(GenomicOperations)
library(ComplexHeatmap)
source("scripts/ckn_utils.R")
source("scripts/load_reddy.R")

##### Get bam filepath 
get_bam_filepath <- function(cofactors = c("MED1", "BRD4", "CDK9", "NIPBL", "SMC1A", "GR")) {
  bam_filepath_list <- c()
  
  if ("GR" %in% cofactors) {
    chip_bam_GR_dir <- "input/ENCODE/A549/GRCh38/chip-seq/bam"
    gr_ctrl <- file.path(chip_bam_GR_dir, "NR3C1_0hour_rep2_ENCFF181HLP.bam")
    gr_dex <- file.path(chip_bam_GR_dir, "NR3C1_1hour_rep1_ENCFF331QXR.bam")
    bam_filepath_list <- c(bam_filepath_list, gr_ctrl, gr_dex)
    
    cofactors <- cofactors[!(cofactors == "GR")]
  }

  chip_bam_dir <- "output/chip-pipeline-GRCh38/alignment"
  for (cofactor in cofactors) {
    for (condition in c("CTRL", "DEX")) {
      basename <- paste("A549", condition, cofactor, "rep1", sep = "_")
      bamfilename <- paste0(basename, ".sorted.dup.bam")
      bam_filepath <- file.path(chip_bam_dir, basename, bamfilename)
      bam_filepath_list <- c(bam_filepath_list, bam_filepath)
    }
  }
  return(bam_filepath_list)
}

##### Make metadata from bam filepath list
make_metadata_from_bam_filepath_list <- function(bam_filepath_list) {
  bam_names <- basename(bam_filepath_list)
  splitted <- strsplit(bam_names, split = "_")
  
  basename <- gsub(".bam", "", bam_names)
  
  target <- splitted %>% purrr:::map(3) %>% unlist
  target <- gsub("rep1|rep2", "GR", target)
  
  condition <- splitted %>% purrr:::map(2) %>% unlist
  condition <- gsub("0hour", "CTRL", condition)
  condition <- gsub("1hour", "DEX", condition)
  
  metadata <- data.frame(design = basename, target, condition, stringsAsFactors = FALSE)
  
  return(metadata)
}

#####
plot_metagene_cofactor <- function(df_metagene, customColors = c("#F5A623", "#4A90E2", "#008000", "#8C001A"), title = "") {
  metagene_plot <- ggplot(df_metagene, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group=condition)) +
    geom_ribbon(aes(fill = condition), alpha = 0.3) +
    scale_fill_manual(values = customColors) +
    geom_line(aes(color = condition)) + # , size = condition, alpha = condition)) +
    scale_color_manual(values = customColors) +
    theme_classic() +
    theme(strip.text.x = element_text(size = 11, face = "bold"),
          strip.text.y = element_text(size = 11, face = "bold", angle = 0, hjust = 0),
          axis.text = element_text(size = 7),
          strip.background = element_rect(colour = NA)) +
    scale_x_continuous(breaks = seq(0, 100, 25),
                       labels = c("-1.5", "-0.75", "0", "0.75", "1.5")) +
    facet_grid(target ~ region) + 
    ggtitle(title) + theme(plot.title = element_text(size = 20, face = "bold")) +
    xlab("Distance in kb") +
    ylab("RPM")
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

##### Load cofactors peaks
stchr <- load_cofactor_stdchr_peaks()
sapply(stchr, length)

##### Load GR binding sites
all_gr_regions <- load_reddy_gr_binding_consensus()
gr_regions_list <- GRangesList(c(all_gr_regions[grep("minutes", names(all_gr_regions))], "1 hour" = all_gr_regions[["1 hour"]]))
gr_regions <- GenomicRanges::reduce(unlist(gr_regions_list))

##### Load induced and repressed gene categories
deg <- readRDS(file = "output/analyses/deg.rds")
upreg <- deg$gene_list$FC1$upreg
downreg <- deg$gene_list$FC1$downreg
responsive_genes <- c(upreg, downreg) %>% unique

##### Define bam_files
bam_filepath <- get_bam_filepath(cofactors = c("MED1", "BRD4", "CDK9", "NIPBL", "SMC1A", "GR"))

##### Define cofactors
cofactors <- c("MED1", "BRD4", "CDK9", "NIPBL", "SMC1A")
# cofactors <- c("MED1")

regions <- list()
for (cofactor in cofactors) {
  message("##### ", cofactor)
  regionSets <- stchr[grep(cofactor, names(stchr))]
  print(sapply(regionSets, length))
  
  ctrl_gr <- subsetByOverlaps(regionSets[[1]], gr_regions)
  message("    ", names(regionSets)[1], "\tNumber of regions with GR : ", length(ctrl_gr))
  dex_gr <- subsetByOverlaps(regionSets[[2]], gr_regions)
  message("    ", names(regionSets)[2], "\tNumber of regions with GR : ", length(dex_gr))

  all_peaks_raw <- c(ctrl_gr, dex_gr)
  message("    # Total number of regions with GR : ", length(all_peaks_raw))
  all_peaks <- GenomicRanges::reduce(all_peaks_raw)
  message("    # Total number of non-overlapping regions with GR : ", length(all_peaks))
  
  # linear annotation
  annot <- annotatePeaks(all_peaks, output = "df", tss = 3000, TxDb = most_expressed_TxDb)
  at_promoters <- annot %>% dplyr::filter(Annot == "Promoter")
  message("    Number of regions (at promoters) : ", nrow(at_promoters))
  
  upreg_genes <- at_promoters %>% dplyr::filter(geneId %in% upreg)
  message("    Number of regions (associated with upregulated genes) : ", nrow(upreg_genes))
  message("         > Unique upregulated genes : ", upreg_genes$geneId %>% unique %>% length)
  
  downreg_genes <- at_promoters %>% dplyr::filter(geneId %in% downreg)
  message("    Number of regions (associated with downregulated genes) : ", nrow(downreg_genes))
  message("         > Unique downregulated genes : ", downreg_genes$geneId %>% unique %>% length)
  
  dex_genes <- at_promoters %>% dplyr::filter(geneId %in% responsive_genes)
  message("    Number of regions (associated with responsive genes) : ", nrow(dex_genes))
  message("         > Unique dex responsive genes : ", dex_genes$geneId %>% unique %>% length)
  
  coordinate_regions_raw <- GRangesList("Upregulated genes" = GRanges(upreg_genes),
                                    "Downregulated genes" = GRanges(downreg_genes))
  coordinate_regions <- lapply(coordinate_regions_raw, resize, width = 3000, fix= "center")
  
  # save regions in an object
  upname <- paste(cofactor, "associated_w_upgene", sep = "_")
  downname <- paste(cofactor, "associated_w_downgene", sep = "_")
  regions[[upname]] <- GRanges(upreg_genes)
  regions[[downname]] <- GRanges(downreg_genes)
  
  # # metagene df
  # bam_filepath_list <- get_bam_filepath()
  # metadata <- make_metadata_from_bam_filepath_list(bam_filepath_list)
  # mg <- metagene2$new(regions = coordinate_regions,
  #                     bam_files = bam_filepath,
  #                     normalization = "RPM",
  #                     force_seqlevels = TRUE,
  #                     assay = 'chipseq',
  #                     cores = 4,
  #                     bin_count = 100)
  # mg$add_metadata(design_metadata = metadata)
  # df_cofactor <- mg$get_data_frame()
  # 
  # # metagene_plot
  # title_df <- paste0("Based on ", cofactor, " regions")
  # df_cofactor$target <- factor(df_cofactor$target, levels = c("GR", "MED1", "BRD4", "CDK9", "NIPBL", "SMC1A"))
  # 
  # metagene_plot <- plot_metagene_cofactor(df_cofactor, title = title_df, customColors = c("#008000", "#8C001A"))
  # saveMetagene(metagene_plot,
  #              output_dir = "output/analyses/ecosystem/metagene_cofactors",
  #              output_file = paste("metagene", cofactor, "at_promoters_of_dex_responsive_genes", sep = "_"),
  #              width_val = 6, height_val = 8)
}

# UpSet plot
sapply(regions, length)

comparison <- c(regions, "GR" = gr_regions)
sapply(comparison, length)

inter_cofactors <- GenomicOperations::GenomicOverlaps(GRangesList(comparison))
matrix_cofactors <- intersect_matrix(inter_cofactors)
sum(matrix_cofactors > 1)
matrix_cofactors[matrix_cofactors > 1] <- 1
sum(matrix_cofactors > 1)
colnames(matrix_cofactors)

m4 <- make_comb_mat(matrix_cofactors, remove_empty_comb_set = TRUE)
m4_plot <- displayUpSet(combMat = m4, threshold = 1)
m4_plot


