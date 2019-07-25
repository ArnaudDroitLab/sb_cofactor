# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(GenomicInteractions)
library(GenomicFeatures)
library(Gviz)
source("scripts/ckn_utils.R")

makeGenomicInteractionsFromDF <- function(file_path) {
  raw_hic <- read_tsv(file_path)
  anchor1 <- GRanges(raw_hic[, 1:3] %>% set_names("seqnames", "start", "end"))
  anchor2 <- GRanges(raw_hic[, 4:6] %>% set_names("seqnames", "start", "end"))
  observed <- raw_hic$observed
  metadata <- raw_hic %>% dplyr::select(color, expectedBL, expectedDonut, expectedH, expectedV,
                                 fdrBL, fdrDonut, fdrH, fdrV,
                                 numCollapsed, centroid1, centroid2, radius)
  
  
  hic_interactions <- GenomicInteractions(anchor1, anchor2)
  mcols(hic_interactions) <- cbind(counts = observed, metadata)
  
  return(hic_interactions)
}

hic_1h_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/hic long range interactions/ENCFF385DHX.tsv"
hic_1h <- makeGenomicInteractionsFromDF(hic_1h_path)
mcols(hic_1h)
regions(hic_1h)
anchorOne(hic_1h)
anchorTwo(hic_1h)
summary(width(regions(hic_1h)))
interactionCounts(hic_1h)
summary(interactionCounts(hic_1h))
plot(density(hic_1h$fdrBL))
plot(density(hic_1h$fdrDonut))
plot(density(hic_1h$fdrH))
plot(density(hic_1h$fdrV))

plotCisTrans(hic_1h)
plotCounts(hic_1h, cut = 100)
plotCounts(hic_1h, cut = 500)
plotCounts(hic_1h, cut = 1000)
plotCounts(hic_1h, cut = 1100)

## Annotation

# hg38
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg38.refseq.db <- TxDb.Hsapiens.UCSC.hg38.knownGene
refseq.genes = genes(hg38.refseq.db)
refseq.transcripts = transcriptsBy(hg38.refseq.db, by="gene")
refseq.transcripts = refseq.transcripts[ names(refseq.transcripts) %in% unlist(refseq.genes$gene_id) ] 
hg38_refseq_promoters <- promoters(refseq.transcripts, 2500, 2500)
hg38_refseq_promoters <- unlist(hg38_refseq_promoters)
hg38_refseq_promoters <- unique(hg38_refseq_promoters) # some duplicate promoters from different transcript isoforms
hg38_refseq_promoters$entrez_id <- names(hg38_refseq_promoters)

# QC
(sum(is.na(hg38_refseq_promoters$entrez_id))/length(hg38_refseq_promoters$entrez_id))*100

# get gene symbols
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl", host="www.ensembl.org", verbose = TRUE)
genes <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id"), filter = "entrezgene_id",
               values = hg38_refseq_promoters$entrez_id, mart = mart, verbose = TRUE)
hg38_refseq_promoters$geneSymbol <- genes$hgnc_symbol[match(hg38_refseq_promoters$entrez_id, genes$entrezgene_id)]
hg38_refseq_promoters$ensemblID <- genes$ensembl_gene_id[match(hg38_refseq_promoters$entrez_id, genes$entrezgene_id)]

# quality control
length(hg38_refseq_promoters$geneSymbol)
(sum(is.na(hg38_refseq_promoters$geneSymbol))/length(hg38_refseq_promoters$geneSymbol))*100
(sum(is.na(hg38_refseq_promoters$ensemblID))/length(hg38_refseq_promoters$ensemblID))*100

names(hg38_refseq_promoters) <- hg38_refseq_promoters$geneSymbol
na.symbol <- is.na(names(hg38_refseq_promoters))
names(hg38_refseq_promoters)[na.symbol] <- hg38_refseq_promoters$tx_name[na.symbol]
sum(is.na(names(hg38_refseq_promoters)))

saveRDS(hg38_refseq_promoters, file = "output/hg38_refseq_promoters.rds")

# peaks to annotate using Hi_c
NIPBL_DEX <- load_cofactor_peaks(c("NIPBL"))[["NIPBL_DEX"]]
summary(width(NIPBL_DEX))
names(NIPBL_DEX) <- NIPBL_DEX$name %>% basename
NIPBL_DEX

# Annotate interactions
annotation.features <- list(promoter = hg38_refseq_promoters, enhancer = NIPBL_DEX)
annotateInteractions(hic_1h, annotation.features)

# Explore results
regions(hic_1h)
saveRDS(hic_1h, file = "output/annotated_hic_1h.rds")

# Node classes
table(regions(hic_1h)$node.class)

# Interactions types
plotInteractionAnnotations(hic_1h, legend = TRUE)
length(hic_1h[isInteractionType(hic_1h, "enhancer", "promoter")])
length(hic_1h[isInteractionType(hic_1h, "enhancer", "distal")])
length(hic_1h[isInteractionType(hic_1h, "enhancer", "enhancer")])

#
ep <- hic_1h[isInteractionType(hic_1h, "enhancer", "promoter")]

# visualisation 1
anchorOne(ep[7])
anchorTwo(ep[7])

GENE_name = "FRMD3"

# GENE_region <- resize(hg38_refseq_promoters[GENE_name], fix = "center", width = 1000000)
# interaction_track <- InteractionTrack(hic_1h, name = "HiC", chromosome = "chr9")
# plotTracks(interaction_track, chromosome="chr9", 
#            from=start(GENE_region), to=end(GENE_region))

# visualisation 2
displayTracks <- function(GENE_name) {
  GENE_region <- resize(hg38_refseq_promoters[GENE_name], fix = "center", width = 250000)
  interaction_track <- InteractionTrack(hic_1h, name = "HiC", chromosome = "chr9")
  # plotTracks(interaction_track, chromosome="chr9", 
  #            from=start(GENE_region), to=end(GENE_region))
  
  promoterTrack <- AnnotationTrack(hg38_refseq_promoters, genome="hg38", name="Promoters",
                                   id=names(hg38_refseq_promoters),  featureAnnotation="id")
  enhTrack <- AnnotationTrack(NIPBL_DEX, genome="hg38", name="NIPBL_DEX", stacking = "dense")
  
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
  
  plotTracks(list(interaction_track, promoterTrack, enhTrack),
             chromosome="chr9", from=start(GENE_region), to=end(GENE_region), 
             sizes=c(0.6, 0.2, 0.2))
}

displayTracks("NFX1")
