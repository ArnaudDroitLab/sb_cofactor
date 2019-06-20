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
  
  bam_bames_without_ext <- strsplit(bam_names, split = "\\.") %>% get_nth_element(1)
  target <- get_nth_element(splitted, 1)
  gsub("NR3C1", "GR", target)
  timepoint <- get_nth_element(splitted, 2) %>% sapply(format_timepoint, USE.NAMES = FALSE)
  replicate <- get_nth_element(splitted, 3)
  
  design_metadata <- data.frame(design = bam_bames_without_ext, target, timepoint, replicate,
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

#
bam_folder <- "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam"
GR_bam_files <- list.files(path = bam_folder, pattern = "NR3C1.*\\.bam$", full.names = TRUE)
design_meta <- make_design_from_bam_list(GR_bam_files)
angptl4_peak <- GRanges("chr19", IRanges(8355327, 8356270))

mg_angptl4 <- metagene2$new(regions = angptl4_peak,
                            bam_files = GR_bam_files,
                            design_metadata = design_meta,
                            assay = 'chipseq')

# mg_angptl4$produce_metagene(title = "Demo metagene plot")

# design_meta <- data.frame(design = mg_angptl4$get_design_group_names(),
#                           target = rep("GR", 11),
#                           timepoint = c("10h", "10h", "10h", "10m", "10m", "10m",
#                                         "15m", "15m", "15m", "1h", "1h"),
#                           replicate = c("1", "2", "3", "1", "2", "3", "1", "2", "3", "1", "2"))



mg_angptl4$add_metadata(design_metadata = design_meta)

# mg_angptl4$produce_metagene(design_metadata = design_meta)

df <- mg_angptl4$get_data_frame()
# relevel group

CustomColor = wes_palette(n=3, name="Cavalcanti1")
ggplot(df, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group=replicate)) +
  geom_ribbon(aes(fill = replicate), alpha = 0.3) +
  scale_fill_manual(values=CustomColor) +
  geom_line(aes(color = replicate), size = 1) + 
  #scale_color_manual(values=CustomColor) +
  # theme_bw(base_size = 20) +
  # theme(panel.grid.major = element_line(),
  #       panel.grid.minor = element_liCne(),
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
