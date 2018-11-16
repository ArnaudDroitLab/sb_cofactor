setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

# [W::hts_idx_load2] The index file is older than the data file
bam_dir <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment"

cofactors <- c("NIPBL", "BRD4", "CDK9") #, "MED1", "SMC1A")
conditions <- c("CTRL", "DEX")

for (cofactor in cofactors) {
  for (condition in conditions) {
    basename <- paste0("A549", "_", condition, "_", cofactor, "_rep1")
    bam_file <- paste0("A549", "_", condition, "_", cofactor, "_rep1.sorted.dup.bam")
    bai_file <- paste0("A549", "_", condition, "_", cofactor, "_rep1.sorted.dup.bai")
    bam_path <- file.path(bam_dir, basename, bam_file)
    bai_path <- file.path(bam_dir, basename, bai_file)
    
    cmd_line <- paste("samtools index", bam_path, bai_path)
    cat(cmd_line, "\t")
    system(cmd_line)
  }
}