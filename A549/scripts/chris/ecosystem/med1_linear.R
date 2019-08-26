# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
source("scripts/ckn_utils.R")
source("scripts/load_reddy.R")

hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
# genes_hg38 <- genes(hg38)
# txdb_genes_hg38 <- makeTxDbFromGRanges(genes_hg38)

# Load GR binding sites
all_gr_regions <- load_reddy_gr_binding_consensus()
gr_regions <- GRangesList(c(all_gr_regions[grep("minutes", names(all_gr_regions))], "1 hour" = all_gr_regions[["1 hour"]]))
gr_regions <- unlist(gr_regions)

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

#####
cofactors <- c("MED1", "BRD4", "CDK9", "NIPBL", "SMC1A")
cofactors <- c("MED1")

clusters <- read_tsv("output/analyses/DPGP_on_a549_dex_0_6hr/groenland_FC1/groenland_FC1_optimal_clustering.txt")

for (cofactor in cofactors) {
  message("##### ", cofactor)
  regionSets <- diffbind[grep(cofactor, names(diffbind))]
  print(sapply(regionSets, length))
  
  linear_annot <- lapply(regionSets, annotatePeaks, tss = 3000, TxDb = most_expressed_TxDb)
  
  for(sets in names(regionSets)) {
    message("     # ", sets)
    sets_annot <- linear_annot[[sets]]
    message("     Number of regions : ", nrow(sets_annot))
    message("")
    sets_annot_with_GR <- subsetByOverlaps(GRanges(sets_annot), gr_regions)
    message("     Number of regions with GR : ", length(sets_annot_with_GR))
    sets_annot_not_distal <- as.data.frame(sets_annot_with_GR) %>% dplyr::filter(Annot %in% c("Promoter"))
    message("     Number of regions (at promoters) : ", nrow(sets_annot_not_distal))
    genes_sets_ID <- sets_annot_not_distal %>% dplyr::select(geneId, Symbol) %>% unique
    colnames(genes_sets_ID) <- c("gene", "symbol")
    message("     Number of unique genes : ", nrow(genes_sets_ID))
    
    # # cluster
    # per_cluster <- left_join(genes_sets_ID, clusters, by = "gene")
    # sum(!is.na(per_cluster$cluster))
    # 
    # myClusters <- per_cluster[!is.na(per_cluster$cluster), ]
    # myClusters$cluster<- paste0("cluster", myClusters$cluster)
    # # message("     Number of genes belonging to a cluster : ", nrow(myClusters))
    # 
    # recap <- myClusters %>% group_by(cluster) %>% tally %>% arrange(desc(n))
    # # print(kable(recap))
    
    # activated/repressed over 12h
    res <- cbind(genes_sets_ID, "upreg" = genes_sets_ID$gene %in% upreg, "downreg" = genes_sets_ID$gene %in% downreg)
    nb_upreg <- sum(res$upreg)
    nb_downreg <- sum(res$downreg)
    message("     Associated with ",  nb_upreg, " induced genes > ", round(nb_upreg/nrow(res)*100, 2), " %")
    message("     Associated with ",  nb_downreg, " repressed genes > ", round(nb_downreg/nrow(res)*100, 2), " %")
    message(" ")
  }
}
