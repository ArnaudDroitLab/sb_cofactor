# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(tidyverse)
source("scripts/ckn_utils.R")
source("scripts/load_reddy.R")

# Load GR binding sites
all_gr_regions <- load_reddy_gr_binding_consensus()
gr_regions_list <- GRangesList(c(all_gr_regions[grep("minutes", names(all_gr_regions))], "1 hour" = all_gr_regions[["1 hour"]]))
gr_regions <- GenomicRanges::reduce(unlist(gr_regions_list))

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
# cofactors <- c("MED1", "BRD4")

clusters <- read_tsv("output/analyses/DPGP_on_a549_dex_0_6hr/groenland_FC1/groenland_FC1_optimal_clustering.txt")

genes_reg_by_cofactors_viaLinear <- list()
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
    
    # save genes_list
    genes_reg_by_cofactors_viaLinear[[sets]] <- genes_sets_ID %>% pull(gene) %>% as.character
  }
}

sapply(genes_reg_by_cofactors_viaLinear, length)
# saveRDS(genes_reg_by_cofactors_viaLinear, file = "output/analyses/ecosystem/genes_reg_by_cofactors_viaLinear.rds")

# HCR
hcr <- GenomicRanges::reduce(gr_regions, min.gapwidth = 25000) # 15102
summary(width(hcr))
hcr <- hcr[width(hcr) > 0]
summary(width(hcr))

raw_promoters_A549 <- promoters(most_expressed_TxDb, columns = c("gene_id"))
stdchr_promoters_A549 <- keepStdChr(raw_promoters_A549)
stdchr_promoters_A549$tx_id <- names(stdchr_promoters_A549)
names(stdchr_promoters_A549) <- stdchr_promoters_A549$gene_id

hcr_vs_genes <- list()
for (cofactor in c("BRD4")) {
  geneup <- genes_reg_by_cofactors_viaLinear[[paste(cofactor, "UP", sep = "_")]]
  tss_up <- stdchr_promoters_A549[geneup]
  tss_up$direction <- "up"
  genedown <- genes_reg_by_cofactors_viaLinear[[paste(cofactor, "DOWN", sep = "_")]]
  tss_down <- stdchr_promoters_A549[genedown]
  tss_down$direction <- "down"

  hcr_cofactor <- hcr_vs_genes[[cofactor]]

  for (i in 1:length(hcr)) {
    hcr_name <- paste("HCR", i, sep = "_")
    hcr_x <- hcr[i] + 500000

    ov <- GRanges()
    ov <- c(ov,
            subsetByOverlaps(tss_up, hcr_x),
            subsetByOverlaps(tss_down, hcr_x))

    if (length(ov) > 1) {
      names(ov) <- c()
      hcr_vs_genes[[cofactor]][[hcr_name]] <- as.data.frame(ov)
      message("##### ", cofactor, " | ", hcr_name, " | ", width(hcr_x) - 1000000, " | ", length(ov), " !!!!!!!!!!!!! ")
    } else if (length(ov) == 1) {
      message("##### ", cofactor, " | ", hcr_name, " | ", width(hcr_x) - 1000000, " | ", length(ov))
    }
  }
}

sapply(hcr_vs_genes, length)

hcr_brd4 <- hcr_vs_genes$BRD4
sapply(hcr_brd4, nrow) %>% table

for (elt in names(hcr_brd4)) {
  df <- hcr_brd4[[elt]]
  if (nrow(df) > 2) {
    print(knitr::kable(df))
  }
}

# 
# saveRDS(hcr_vs_genes, file = "output/analyses/hcr/hcr_vs_genes.rds")

### Work with TAD
# tad_1h_raw <- read_delim(file = "input/ENCODE/A549/GRCh38/hic_TAD/hic_dex.t1.tads.bedpe", delim = "\t")[, 1:3]
# colnames(tad_1h_raw) <- c("seqnames", "start", "end")
# tad_1h <- GRanges(tad_1h_raw)
# 
# tad_vs_genes <- list()
# for (cofactor in cofactors) {
#   geneup <- genes_reg_by_cofactors_viaLinear[[paste(cofactor, "UP", sep = "_")]]
#   tss_up <- stdchr_promoters_A549[geneup]
#   tss_up$direction <- "up"
#   genedown <- genes_reg_by_cofactors_viaLinear[[paste(cofactor, "DOWN", sep = "_")]]
#   tss_down <- stdchr_promoters_A549[genedown]
#   tss_down$direction <- "down"
# 
#   tad_cofactor <- tad_vs_genes[[cofactor]]
# 
#   for (i in 1:length(tad_1h)) {
#     tad_name <- paste("TAD", i, sep = "_")
#     tad_x <- tad_1h[i]
# 
#     ov <- GRanges()
#     ov <- c(ov,
#             subsetByOverlaps(tss_up, tad_x),
#             subsetByOverlaps(tss_down, tad_x))
# 
#     if (length(ov) > 1) {
#       names(ov) <- c()
#       tad_vs_genes[[cofactor]][[tad_name]] <- as.data.frame(ov)
#       message("##### ", cofactor, " | ", tad_name, " | ", width(tad_x), " | ", length(ov), " !!!!!!!!!!!!! ")
#     } else {
#       message("##### ", cofactor, " | ", tad_name, " | ", width(tad_x), " | ", length(ov))
#     }
#   }
# }
# 
# sapply(tad_vs_genes, length)
# saveRDS(tad_vs_genes, file = "output/analyses/hcr/tad_vs_genes_v1.rds")

# reperform TAD identification?
