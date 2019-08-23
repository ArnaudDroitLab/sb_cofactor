# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
source("scripts/ckn_utils.R")

hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
# genes_hg38 <- genes(hg38)
# txdb_genes_hg38 <- makeTxDbFromGRanges(genes_hg38)

#####
`%notin%` <- Negate(`%in%`)

##### Load categories
deg <- readRDS(file = "output/analyses/deg.rds")
upreg <- deg$gene_list$FC1$upreg
downreg <- deg$gene_list$FC1$downreg

##### Count peaks
raw <- load_cofactor_peaks()
stchr <- load_cofactor_stdchr_peaks()
diffbind <- load_diffbind_cofactors_peaks()

sapply(raw, length)
sapply(stchr, length)
sapply(diffbind, length)

##### MED1
MED1 <- diffbind[grep("MED1", names(diffbind))]
sapply(MED1, length)

# linear annotation
linear_annot_MED1 <- lapply(MED1, annotatePeaks, tss = 3000, TxDb = most_expressed_TxDb)

##########################
# analyse MED1_UP
##########################
med1_up_annot <- linear_annot_MED1$MED1_UP
plotAnnotation(med1_up_annot)

med1_up_not_distal <- med1_up_annot %>% dplyr::filter(Annot %notin% c("Distal Intergenic", "Downstream"))
plotAnnotation(med1_up_not_distal)

genes_med1_up_symbol <- med1_up_not_distal %>% dplyr::pull(Symbol) %>% unique %>% sort
genes_med1_up_ensg <- med1_up_not_distal %>% dplyr::select(geneId, Symbol) %>% unique
colnames(genes_med1_up_ensg) <- c("gene", "symbol")

###
# which cluster do they belong?
clusters <- read_tsv("output/analyses/DPGP_on_a549_dex_0_6hr/groenland_FC1/groenland_FC1_optimal_clustering.txt")
per_cluster <- left_join(genes_med1_up_ensg, clusters, by = "gene")
sum(!is.na(per_cluster$cluster))

# 226 genes belongs to cluster out of 1472
myClusters <- per_cluster[!is.na(per_cluster$cluster), ]
myClusters$cluster<- paste0("cluster", myClusters$cluster)

recap <- myClusters %>% group_by(cluster) %>% tally %>% arrange(desc(n))
kable(recap)
###

# which categories do they belong?
res <- cbind(genes_med1_up_ensg, "upreg" = genes_med1_up_ensg$gene %in% upreg, "downreg" = genes_med1_up_ensg$gene %in% downreg)
nb_upreg <- sum(res$upreg)
nb_downreg <- sum(res$downreg)

##########################
# analyse MED1_DOWN
##########################
med1_down_annot <- linear_annot_MED1$MED1_DOWN
plotAnnotation(med1_down_annot)

med1_down_not_distal <- med1_down_annot %>% dplyr::filter(Annot %notin% c("Distal Intergenic", "Downstream"))
plotAnnotation(med1_down_not_distal)

genes_med1_down_symbol <- med1_down_not_distal %>% dplyr::pull(Symbol) %>% unique %>% sort
genes_med1_down_ensg <- med1_down_not_distal %>% dplyr::select(geneId, Symbol) %>% unique
colnames(genes_med1_down_ensg) <- c("gene", "symbol")

# which cluster do they belong?
clusters <- read_tsv("output/analyses/DPGP_on_a549_dex_0_6hr/groenland_FC1/groenland_FC1_optimal_clustering.txt")
per_cluster <- left_join(genes_med1_down_ensg, clusters, by = "gene")
sum(!is.na(per_cluster$cluster))

# 52 genes belongs to cluster out of 703
myClusters <- per_cluster[!is.na(per_cluster$cluster), ]
myClusters$cluster <- paste0("cluster", myClusters$cluster)

recap <- myClusters %>% group_by(cluster) %>% tally %>% arrange(desc(n))
kable(recap)

# which categories do they belong?
res <- cbind(genes_med1_down_ensg, "upreg" = genes_med1_down_ensg$gene %in% upreg, "downreg" = genes_med1_down_ensg$gene %in% downreg)
nb_upreg <- sum(res$upreg)
nb_downreg <- sum(res$downreg)

##########################
# analyse MED1_UNBIASED
##########################
med1_un_annot <- linear_annot_MED1$MED1_UNBIASED
plotAnnotation(med1_un_annot)

med1_un_not_distal <- med1_un_annot %>% dplyr::filter(Annot %notin% c("Distal Intergenic", "Downstream"))
plotAnnotation(med1_un_not_distal)

genes_med1_un_symbol <- med1_un_not_distal %>% dplyr::pull(Symbol) %>% unique %>% sort
genes_med1_un_ensg <- med1_un_not_distal %>% dplyr::select(geneId, Symbol) %>% unique
colnames(genes_med1_un_ensg) <- c("gene", "symbol")

# which cluster do they belong?
clusters <- read_tsv("output/analyses/DPGP_on_a549_dex_0_6hr/groenland_FC1/groenland_FC1_optimal_clustering.txt")
per_cluster <- left_join(genes_med1_un_ensg, clusters, by = "gene")
sum(!is.na(per_cluster$cluster))

# 110 genes belongs to cluster out of 2812
myClusters <- per_cluster[!is.na(per_cluster$cluster), ]
myClusters$cluster <- paste0("cluster", myClusters$cluster)

recap <- myClusters %>% group_by(cluster) %>% tally %>% arrange(desc(n))
kable(recap)

# which categories do they belong?
res <- cbind(genes_med1_un_ensg, "upreg" = genes_med1_un_ensg$gene %in% upreg, "downreg" = genes_med1_un_ensg$gene %in% downreg)
nb_upreg <- sum(res$upreg)
nb_downreg <- sum(res$downreg)
