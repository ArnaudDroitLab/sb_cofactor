# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(GenomicInteractions)
source("scripts/ckn_utils.R")
source("scripts/load_reddy.R")

### Load cofactors peaks
diffbind <- load_diffbind_cofactors_peaks()
sapply(diffbind, length)

### Load GR binding sites
all_gr_regions <- load_reddy_gr_binding_consensus()
gr_regions <- GRangesList(c(all_gr_regions[grep("minutes", names(all_gr_regions))], "1 hour" = all_gr_regions[["1 hour"]]))
gr_regions <- unlist(gr_regions)

### Load reference genome
raw_promoters_A549 <- promoters(most_expressed_TxDb, columns = "gene_id")
stdchr_promoters_A549 <- keepStdChr(raw_promoters_A549)
names(stdchr_promoters_A549) <- stdchr_promoters_A549$gene_id

### Load Hi-C data 1 hour
hic_1h <- readRDS("output/analyses/annotate_peaks_with_hic/hic_1h_GIObject.rds")

#####
# `%notin%` <- Negate(`%in%`)

##### Load induced and repressed gene categories
deg <- readRDS(file = "output/analyses/deg.rds")
upreg <- deg$gene_list$FC1$upreg
downreg <- deg$gene_list$FC1$downreg

#### Define cofactors
cofactors <- c("MED1", "BRD4", "CDK9", "NIPBL", "SMC1A")
cofactors <- c("MED1")

for (cofactor in cofactors) {
  message("##### ", cofactor)
  regionSets <- diffbind[grep(cofactor, names(diffbind))]
  print(sapply(regionSets, length))
  
  for(sets in names(regionSets)) {
    message("     ####################")
    message("     ### ", sets)
    message("     ####################")
    
    sets <- regionSets[[sets]]
    message("     Number of regions : ", length(sets))
    
    sets_with_GR <- subsetByOverlaps(sets, gr_regions)
    message("     Number of regions with GR : ", length(sets_with_GR))
    names(sets_with_GR) <- sets_with_GR$name
    message("")
    
    # Annotate with Hi-C
    annotation.features <- list(promoter = stdchr_promoters_A549, enhancer = sets_with_GR)
    annotateInteractions(hic_1h, annotation.features)
    
    ### Explore interactions
    message("Number of annotations : ")
    print(table(regions(hic_1h)$node.class))
    
    # Interactions types
    plotInteractionAnnotations(hic_1h, legend = TRUE, viewpoints = "enhancer")
    nb_ep <- length(hic_1h[isInteractionType(hic_1h, "enhancer", "promoter")])
    nb_ed <- length(hic_1h[isInteractionType(hic_1h, "enhancer", "distal")])
    nb_ee <- length(hic_1h[isInteractionType(hic_1h, "enhancer", "enhancer")])
    message("  Number of enhancer-promoter interactions : ", nb_ep)
    message("  Number of enhancer-distal interactions : ", nb_ed)
    message("  Number of enhancer-enhancer interactions : ", nb_ee)
    message("")
    
    # Focus on enhancer-promoter
    ep <- hic_1h[isInteractionType(hic_1h, "enhancer", "promoter")]
    
    # Retrieve gene name from promoters
    prom1 <- anchorOne(ep)$promoter.id %>% unlist %>% na.omit %>% unique
    prom2 <- anchorTwo(ep)$promoter.id %>% unlist %>% na.omit %>% unique
    prom12 <- c(prom1, prom2) %>% unique %>% as.data.frame
    colnames(prom12) <- c("gene")
    message("Number of unique genes : ", nrow(prom12))
    
    # activated/repressed over 12h
    res <- cbind(prom12, "upreg" = prom12$gene %in% upreg, "downreg" = prom12$gene %in% downreg)
    nb_upreg <- sum(res$upreg)
    nb_downreg <- sum(res$downreg)
    message("     Associated with ",  nb_upreg, " induced genes > ", round(nb_upreg/nrow(res)*100, 2), " %")
    message("     Associated with ",  nb_downreg, " repressed genes > ", round(nb_downreg/nrow(res)*100, 2), " %")
    message(" ")
  }
}
