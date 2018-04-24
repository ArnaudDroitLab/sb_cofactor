library(knitr)
library(dplyr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source("scripts/metagene_polII/function_generate_metagene_object.R")

########################################
# Define directory
########################################
bam_dir <- "output/chip-pipeline-PolII-GRCh38/alignment"
output_dir <- "output/chip-pipeline-PolII-GRCh38/metagene_polII"

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

target_list <- c("POL2-ser2", "POL2", "WCE")
sh_list <-  c("shCTRL-1", "shCTRL-2", "shNIPBL-3", "shNIPBL-5")

for (target in target_list) {
	for (sh in sh_list) {
	generate_metagene_object(target, sh, region_list)
	}
}


