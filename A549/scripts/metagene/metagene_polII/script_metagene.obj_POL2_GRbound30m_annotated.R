library(knitr)
library(dplyr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source("scripts/metagene/metagene_polII/function_generate_metagene_pol2_object.R")
source("scripts/metagene/function_resize_peaks.R")
source("scripts/load_reddy.R")

########################################
# Define directory
########################################
bam_dir <- "output/chip-pipeline-PolII-GRCh38/alignment"
output_dir <- "output/chip-pipeline-PolII-GRCh38/metagene_polII"

###############################################################################
# Define regions over which metagenes will be plotted.
###############################################################################
most_expressed = load_most_expressed_transcripts()
most_expressed_TxDb = load_most_expressed_TxDb()
seqlevelsStyle(most_expressed_TxDb) <- "UCSC"

# Import GR binding regions.
gr_regions = load_reddy_gr_binding_consensus()

#
gr_regions_30 <- gr_regions[["30 minutes"]]
gr_regions_30 <- resize_all_peaks(gr_regions_30, window = 600)

gr_regions_30_annotated <- ChIPseeker::annotatePeak(gr_regions_30, TxDb=most_expressed_TxDb)
gr_regions_30_annotated_df <- as.data.frame(gr_regions_30_annotated)

gr_regions_30_annotated_df$Annot <- gr_regions_30_annotated_df$annotation
gr_regions_30_annotated_df$Annot  <- gsub(" \\(.*\\)", "", gr_regions_30_annotated_df$Annot)
gr_regions_30_annotated_df$Annot  <- as.factor(gr_regions_30_annotated_df$Annot)

###
gr_regions_30_annotated_df_3utr <- GRanges(subset(gr_regions_30_annotated_df, Annot  == "3' UTR")) # 134
gr_regions_30_annotated_df_5utr <- GRanges(subset(gr_regions_30_annotated_df, Annot  == "5' UTR")) # 16
gr_regions_30_annotated_df_distalintergenic <- GRanges(subset(gr_regions_30_annotated_df, Annot  == "Distal Intergenic")) # 6439
gr_regions_30_annotated_df_downstream <- GRanges(subset(gr_regions_30_annotated_df, Annot  == "Downstream")) # 94
gr_regions_30_annotated_df_exon <- GRanges(subset(gr_regions_30_annotated_df, Annot  == "Exon")) # 222
gr_regions_30_annotated_df_intron <- GRanges(subset(gr_regions_30_annotated_df, Annot  == "Intron")) # 3792
gr_regions_30_annotated_df_promoter <- GRanges(subset(gr_regions_30_annotated_df, Annot  == "Promoter")) # 528

# NB: region_list is a list of GRanges object
region_list = list(GR30_3utr = gr_regions_30_annotated_df_3utr,
				   GR30_5utr = gr_regions_30_annotated_df_5utr,
				   GR30_distalintergenic = gr_regions_30_annotated_df_distalintergenic,
				   GR30_downstream = gr_regions_30_annotated_df_downstream,
				   GR30_exon = gr_regions_30_annotated_df_exon,
				   GR30_intron = gr_regions_30_annotated_df_intron,
				   GR30_promoter = gr_regions_30_annotated_df_promoter)

###############################################################################
# Generate the metagene object.
###############################################################################

target_list <- c("POL2-ser2", "POL2", "WCE")
sh_list <-  c("shCTRL-1", "shCTRL-2", "shNIPBL-3", "shNIPBL-5")

for (target in target_list) {
	for (sh in sh_list) {
	generate_metagene_object(target, sh, region_list, bin=200)
	}
}


