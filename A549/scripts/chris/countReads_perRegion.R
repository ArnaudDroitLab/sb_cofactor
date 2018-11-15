setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/ckn_utils.R")

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
#####
bam_path <- generate_bampathVector(cofactors = c("NIPBL", "BRD4", "CDK9"),
                                   conditions = c("CTRL", "DEX"))

peaks_NBC_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NBC"

speNBC_CTRL_path <- file.path(peaks_NBC_dir, "A549_NBC_CTRL_specific.bed")
speNBC_CTRL <- rtracklayer::import(speNBC_CTRL_path)

speNBC_CTRL_coordVector <- generate_coordVector(speNBC_CTRL)

# samtools view /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam chr12:104290274-104290284

# function: input: GRranges, output: vector name of region
# system samtools
# samtools view bam region | wc -l
# output: pour chaque set de regions, tableau de counts avec les 6 samples
# faire une PCA?
# correlation entre ces deux listes de chiffres
# faire une PCA
