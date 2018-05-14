library(knitr)
library(dplyr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source("scripts/metagene/function_generate_metagene_object.R")
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

# Import GR binding regions.
gr_regions = load_reddy_gr_binding_consensus()

# Determine which genes overlap GR regions.
annotated_gr = ChIPseeker::annotatePeak(gr_regions[["30 minutes"]], TxDb=most_expressed_TxDb)
annotated_gr_df = as.data.frame(annotated_gr)
annotated_gr_df$ENTREZID = mapIds(org.Hs.eg.db, keys=as.character(annotated_gr_df$geneId), keytype="ENSEMBL", column="ENTREZID")
bound_gene_ids = subset(annotated_gr_df, distanceToTSS <= 3000)$ENTREZID
unbound_gene_ids = setdiff(all_genes$entrezgene, bound_gene_ids)

# # Determine differentially expressed genes at 2 hour.
# de_results = load_reddy_de_list()[["2h"]]$Full
# de_results$ENTREZID = mapIds(org.Hs.eg.db, keys=as.character(de_results$gene_id), keytype="ENSEMBL", column="ENTREZID")
# de_gene_ids = subset(de_results, abs(log2FoldChange) >= log2(1.5) & padj <= 0.05)$ENTREZID
# up_gene_ids = subset(de_results, log2FoldChange <= -log2(1.5) & padj <= 0.05)$ENTREZID
# down_gene_ids = subset(de_results, log2FoldChange >= log2(1.5) & padj <= 0.05)$ENTREZID
# 
# # Define group of regions based on DE status AND GR-binding status.
# de_bound_gene_ids = intersect(de_gene_ids, bound_gene_ids)
# up_bound_gene_ids = intersect(up_gene_ids, bound_gene_ids)
# down_bound_gene_ids = intersect(down_gene_ids, bound_gene_ids)
# de_unbound_gene_ids = intersect(de_gene_ids, unbound_gene_ids)
# up_unbound_gene_ids = intersect(up_gene_ids, unbound_gene_ids)
# down_unbound_gene_ids = intersect(down_gene_ids, unbound_gene_ids)


region_list = list(BoundTSS=get_tss(all_genes, bound_gene_ids, flank_size=500),
                   UnboundTSS=get_tss(all_genes, unbound_gene_ids, flank_size=500),
				   AllGRTSS=get_tss(all_genes, annotated_gr_df$ENTREZID, flank_size=500))

###############################################################################
# Generate the metagene object.
###############################################################################

#cofactor_list <- c("BRD4", "CDK9", "NIPBL", "SMC1A", "MED1")

cofactor_list <- c("NIPBL", "SMC1A", "MED1")

for (cofactor in cofactor_list) {
  message("##########     ", cofactor, "     ##########")
	generate_metagene_object(cofactor, region_list, bin=200)
}


