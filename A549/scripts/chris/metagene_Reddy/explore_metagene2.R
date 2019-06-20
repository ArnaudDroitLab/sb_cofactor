library(tidyverse)
library(knitr)
library(metagene2)
library(wesanderson)

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


### TEST
# regions <- get_demo_regions()
# bam_files <- get_demo_bam_files()
# 
# mg <- metagene2$new(regions = get_demo_regions(), 
#                     bam_files = get_demo_bam_files()[1:4], 
#                     assay = 'chipseq')
# 
# mg$produce_metagene(title = "Demo metagene plot")
# 
# design_meta = data.frame(design=mg$get_design_group_names(),
#                          Align=c("Align1", "Align1", "Align2", "Align2"),
#                          Rep=c(1, 2, 1, 2))
# 
# mg$produce_metagene(design_metadata=design_meta, facet_by=Align~Rep, group_by="region")

# Define one region at first
# angptl4_peak <- GRanges("chr19", IRanges(8355327, 8356270))
# il11_peak <- GRanges("chr19", IRanges(55372624, 55373186))
# gapdh_peak <- GRanges("chr12", IRanges(6532244, 6532808))
mafk_peak <- GRanges("chr7", IRanges(1519343, 1522213))

peak <- mafk_peak

# BAM
bam_folder <- "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam"

# GR
GR_bam_files <- list.files(path = bam_folder, pattern = "(NR3C1).*\\.bam$", full.names = TRUE)
GR_design_meta <- make_design_from_bam_list(GR_bam_files)

GR_mg <- metagene2$new(regions = peak,
                            bam_files = GR_bam_files,
                            assay = 'chipseq')

GR_mg$add_metadata(design_metadata = GR_design_meta)
GR_df <- GR_mg$get_data_frame()

# EP300
EP300_bam_files <- list.files(path = bam_folder, pattern = "(EP300).*\\.bam$", full.names = TRUE)
EP300_design_meta <- make_design_from_bam_list(EP300_bam_files)

EP300_mg <- metagene2$new(regions = peak,
                       bam_files = EP300_bam_files,
                       assay = 'chipseq')

EP300_mg$add_metadata(design_metadata = EP300_design_meta)
EP300_df <- EP300_mg$get_data_frame()

# H3K27ac
H3K27ac_bam_files <- list.files(path = bam_folder, pattern = "(H3K27ac).*\\.bam$", full.names = TRUE)
H3K27ac_design_meta <- make_design_from_bam_list(H3K27ac_bam_files)

H3K27ac_mg <- metagene2$new(regions = peak,
                          bam_files = H3K27ac_bam_files,
                          assay = 'chipseq')

H3K27ac_mg$add_metadata(design_metadata = H3K27ac_design_meta)
H3K27ac_df <- H3K27ac_mg$get_data_frame()

# JUN
JUN_bam_files <- list.files(path = bam_folder, pattern = "(JUN).*\\.bam$", full.names = TRUE)
JUN_design_meta <- make_design_from_bam_list(JUN_bam_files)

JUN_mg <- metagene2$new(regions = peak,
                          bam_files = JUN_bam_files,
                          assay = 'chipseq')

JUN_mg$add_metadata(design_metadata = JUN_design_meta)
JUN_df <- JUN_mg$get_data_frame()

# bigdf
bigdf <- rbind(GR_df, EP300_df, H3K27ac_df, JUN_df)
bigdf$target <- factor(bigdf$target, levels = c("GR", "EP300", "H3K27ac", "JUN") )

CustomColor = wes_palette(n=3, name="Cavalcanti1")
ggplot(bigdf, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group=replicate)) +
  geom_ribbon(aes(fill = replicate), alpha = 0.3) +
  # scale_fill_manual(values=CustomColor) +
  geom_line(aes(color = replicate), size = 0.5) + 
  # scale_color_manual(values=CustomColor) +
  # theme_bw(base_size = 20) +
  # theme(panel.grid.major = element_line(),
  #       panel.grid.minor = element_line(),
  #       panel.background = element_blank()) +
  # theme(legend.position = "bottom",
  #       legend.direction = "vertical") +
  # scale_x_continuous(breaks = seq(0, 100, 25),
  #                    labels = c(c(-1.5, -0.75), "TSS", c(0.75, 1.5))) +
  # theme(axis.text.x = element_text(size=10),
  #       axis.text.y = element_text(size=10)) +
  # ggtitle(title) +
  # xlab("Position from TSS (kb)") +
  # ylab("Mean coverage (RPM)") +
  facet_grid(target ~ timepoint)
