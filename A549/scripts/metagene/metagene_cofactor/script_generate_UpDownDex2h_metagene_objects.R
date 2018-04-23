library(knitr)
library(dplyr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source("scripts/metagene/metagene_cofactor/function_generate_metagene_object.R")
source("scripts/load_reddy.R")

########################################
# Define directory
########################################
bam_dir <- "output/chip-pipeline-GRCh38/alignment"
output_dir <- "output/chip-pipeline-GRCh38/metagene/metagene_cofactor"

###############################################################################
# Define regions over which metagenes will be plotted.
###############################################################################

most_expressed = load_most_expressed_transcripts()
most_expressed_TxDb = load_most_expressed_TxDb()
seqlevelsStyle(most_expressed_TxDb) <- "UCSC"

all_genes = most_expressed

# Remove all regions not on the main chromosomes.
all_genes = all_genes[!grepl("_", seqnames(all_genes))]

# Remove all genes that are smaller than our defined TSS regions.
all_genes = all_genes[width(all_genes) >= 200]

# Define TSS regions based on all kept gene locations.
all_TSS = GenomicRanges::promoters(all_genes, upstream=500, downstream=500)

# Determine differentially expressed genes at 2 hour.
de_results = load_reddy_de_list()[["2h"]]$Full
de_results$ENTREZID = mapIds(org.Hs.eg.db, keys=as.character(de_results$gene_id), keytype="ENSEMBL", column="ENTREZID")
de_gene_ids = subset(de_results, abs(log2FoldChange) >= log2(1.5) & padj <= 0.05)$ENTREZID
up_gene_ids = subset(de_results, log2FoldChange <= -log2(1.5) & padj <= 0.05)$ENTREZID
down_gene_ids = subset(de_results, log2FoldChange >= log2(1.5) & padj <= 0.05)$ENTREZID

region_list = list(UpRegulatedTSS=get_tss(all_genes, up_gene_ids),
				   DownRegulatedTSS=get_tss(all_genes, down_gene_ids)

###############################################################################
# Generate the metagene object.
###############################################################################

cofactor_list <- c("BRD4", "CDK9", "NIPBL", "SMC1A", "MED1")

for (cofactor in cofactor_list) {
	generate_metagene_object(cofactor, region_list)
}


