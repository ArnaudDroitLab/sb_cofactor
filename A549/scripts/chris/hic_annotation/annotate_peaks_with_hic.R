# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(GenomicInteractions)
library(tidyverse)
library(Gviz)

# load GR_DOWN
GR_diffbind_downreg <- readRDS(file = "output/analyses/annotate_peaks_with_hic/GR_diffbind_downreg.rds")
names(GR_diffbind_downreg)
sapply(GR_diffbind_downreg, length)
GR_25mVS10m <- GR_diffbind_downreg[["25mVS10m_downreg"]] # 4450 regions
names(GR_25mVS10m) <- paste("GR_25mVS10m", paste0("peak", 1:length(GR_25mVS10m)), sep = "_")

# load hg38_promoters
hg38_promoters <- readRDS("output/analyses/annotate_peaks_with_hic/hg38_refseq_promoters.rds")

# 
hic_1h <- readRDS("output/analyses/annotate_peaks_with_hic/hic_1h_GIObject.rds")

# Annotate interactions
annotation.features <- list(promoter = hg38_promoters, enhancer = GR_25mVS10m)
annotateInteractions(hic_1h, annotation.features)

# explore results in hic_1h
regions(hic_1h)

# Node classes
table(regions(hic_1h)$node.class)

# Interactions types
plotInteractionAnnotations(hic_1h, legend = TRUE)
length(hic_1h[isInteractionType(hic_1h, "enhancer", "promoter")])
length(hic_1h[isInteractionType(hic_1h, "enhancer", "distal")])
length(hic_1h[isInteractionType(hic_1h, "enhancer", "enhancer")])

# enhancer-promoter
ep <- hic_1h[isInteractionType(hic_1h, "enhancer", "promoter")]
enh1 <- anchorOne(ep)$enhancer.id %>% unlist %>% na.omit %>% unique
enh2 <- anchorTwo(ep)$enhancer.id %>% unlist %>% na.omit %>% unique
enh12 <- c(enh1, enh2) %>% unique
prom1 <- anchorOne(ep)$promoter.id %>% unlist %>% na.omit %>% unique
prom2 <- anchorTwo(ep)$promoter.id %>% unlist %>% na.omit %>% unique
prom12 <- c(prom1, prom2) %>% unique
sum(grepl("ENST", prom12))
prom12_esng <- prom12[grepl("ENSG", prom12)] %>% as.data.frame # 717 genes
colnames(prom12_esng) <- "gene"
prom12_esng$gene <- as.character(prom12_esng$gene)

# which cluster do they belong?
clusters <- read_tsv("output/analyses/DPGP_on_a549_dex_0_6hr/groenland_FC1/groenland_FC1_optimal_clustering.txt")
per_cluster <- left_join(prom12_esng, clusters, by = "gene")

sum(!is.na(per_cluster$cluster))

# 57 genes belongs to cluster
myClusters <- per_cluster[!is.na(per_cluster$cluster), ]
myClusters$cluster<- paste0("cluster", myClusters$cluster)

recap <- myClusters %>% group_by(cluster) %>% tally


for (geneENSG in myClusters$gene[30:length(myClusters$gene)]) {
  region <- hg38_promoters[geneENSG] %>% as.data.frame
  symbol <- region %>% dplyr::select(geneSymbol) %>% pull(geneSymbol)
  if (symbol == "") {
    symbol = "NA"
  }
  chr <- region %>% dplyr::select(seqnames) %>% pull(seqnames) %>% as.character
  cluster <- myClusters %>% dplyr::filter(gene == geneENSG) %>% pull(cluster)
  message("### ", symbol, " | ", geneENSG, " | ", chr, " | ", cluster)
  title <- paste0(symbol, " | ", geneENSG, " | ", chr, " | ", cluster)
  displayTracks(geneENSG, chr, title)
  Sys.sleep(2)
  saveTracks(tracks = displayTracks(geneENSG, chr, title),
             output_dir = "output/analyses/annotate_peaks_with_hic",
             output_file = paste("20190726", symbol, cluster, geneENSG, sep = "_"),
             format = "pdf",
             width = 20, height = 12)
  Sys.sleep(2)
}

saveTracks <- function(tracks, output_dir, output_file, width_val = 25, height_val = 22, format = "pdf") {
  output_filepath <- file.path(output_dir, paste0(output_file, ".", format))
  pdf(file = output_filepath, width = width_val, height = height_val)
  tracks
  dev.off()
  message(" > Tracks saved in ", output_filepath)
}


displayTracks <- function(GENE_name, chr, title) {
  GENE_region <- resize(hg38_promoters[GENE_name], fix = "center", width = 1000000)
  interaction_track <- InteractionTrack(hic_1h, name = "HiC", chromosome = chr)
  # plotTracks(interaction_track, chromosome="chr9", 
  #            from=start(GENE_region), to=end(GENE_region))
  
  promoterTrack <- AnnotationTrack(hg38_promoters, genome="hg38", name="Promoters",
                                   id=hg38_promoters$geneSymbol,  featureAnnotation="id")
  enhTrack <- AnnotationTrack(GR_25mVS10m, genome="hg38", name="GR_25mVS10m", stacking = "dense")
  
  displayPars(promoterTrack) <- list(fill = "deepskyblue", col = NA, 
                                     fontcolor.feature = "black", fontsize=8,
                                     just.group="below")
  displayPars(enhTrack) <- list(fill = "black", col = NA)
  displayPars(interaction_track) = list(col.interactions="red", 
                                        col.anchors.fill ="blue",
                                        col.anchors.line = "black",
                                        interaction.dimension="height", 
                                        interaction.measure ="counts",
                                        plot.trans=FALSE,
                                        plot.outside = TRUE, 
                                        col.outside="lightblue", 
                                        anchor.height = 0.1)
  gtrack <- GenomeAxisTrack()
  
  plotTracks(list(interaction_track, promoterTrack, enhTrack, gtrack),
             chromosome=chr, from=start(GENE_region), to=end(GENE_region), 
             sizes=c(0.6, 0.2, 0.2, 0.2), main = title)
}



# visualisation 1
anchorOne(ep)
anchorTwo(ep)

ensg = "FRMD3"

# GENE_region <- resize(hg38_refseq_promoters[GENE_name], fix = "center", width = 1000000)
# interaction_track <- InteractionTrack(hic_1h, name = "HiC", chromosome = "chr9")
# plotTracks(interaction_track, chromosome="chr9", 
#            from=start(GENE_region), to=end(GENE_region))
