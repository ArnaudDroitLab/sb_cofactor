library(knitr)
library(dplyr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source("metagene_cofactor/function_generate_metagene_object.R")

########################################
# Define directory
########################################
bam_dir <- "output/chip-pipeline-GRCh38/alignment"
output_dir <- "output/chip-pipeline-GRCh38/metagene_cofactor"

###############################################################################
# Define regions over which metagenes will be plotted.
###############################################################################
all_genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Remove all regions not on the main chromosomes.
all_genes = all_genes[!grepl("_", seqnames(all_genes))]

# Remove all genes that are smaller than our defined TSS regions.
all_genes = all_genes[width(all_genes) >= 1500]

# Define TSS regions based on all kept gene locations.
all_TSS = GenomicRanges::promoters(all_genes, upstream=1500, downstream=1500)

region_list = list(AllTSS = all_TSS)

###############################################################################
# Generate the metagene object.
###############################################################################

cofactor_list <- c("BRD4", "CDK9", "NIPBL", "SMC1A", "MED1")

for (cofactor in cofactor_list) {
	generate_metagene_object(cofactor, region_list)
}


