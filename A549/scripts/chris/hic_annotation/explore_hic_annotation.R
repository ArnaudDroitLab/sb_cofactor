# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(tidyverse)
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

hic_1h_path <- "input/ENCODE/A549/GRCh38/hic long range interactions/ENCFF385DHX.tsv"
hic_1h <- makeGenomicInteractionsFromDF(hic_1h_path)
saveRDS(hic_1h, file = "output/analyses/annotate_peaks_with_hic/hic_1h_GIObject.rds")
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

# # hg38
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
# mart = useMart("ensembl", dataset="hsapiens_gene_ensembl", host="www.ensembl.org", verbose = TRUE)
# genes <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id", "ensembl_transcript_id_version"),
#                filter = "ensembl_transcript_id_version", values = hg38_refseq_promoters$tx_name,
#                mart = mart, verbose = TRUE)
# saveRDS(genes, file = "output/analyses/annotate_peaks_with_hic/gene_nomenclature.rds")
genes <- readRDS("output/analyses/annotate_peaks_with_hic/gene_nomenclature.rds")

hg38_refseq_promoters$geneSymbol <- genes$hgnc_symbol[match(hg38_refseq_promoters$tx_name, genes$ensembl_transcript_id_version)]
hg38_refseq_promoters$ensemblID <- genes$ensembl_gene_id[match(hg38_refseq_promoters$tx_name, genes$ensembl_transcript_id_version)]

(sum(is.na(hg38_refseq_promoters$ensemblID))/length(hg38_refseq_promoters$ensemblID))*100

# quality control
length(hg38_refseq_promoters$geneSymbol)
(sum(is.na(hg38_refseq_promoters$geneSymbol))/length(hg38_refseq_promoters$geneSymbol))*100
(sum(is.na(hg38_refseq_promoters$ensemblID))/length(hg38_refseq_promoters$ensemblID))*100

names(hg38_refseq_promoters) <- hg38_refseq_promoters$ensemblID
na.symbol <- is.na(names(hg38_refseq_promoters))
names(hg38_refseq_promoters)[na.symbol] <- hg38_refseq_promoters$tx_name[na.symbol]
sum(is.na(names(hg38_refseq_promoters)))

saveRDS(hg38_refseq_promoters, file = "output/analyses/annotate_peaks_with_hic/hg38_refseq_promoters.rds")