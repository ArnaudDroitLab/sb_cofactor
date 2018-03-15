#setwd("C:/Dev/Projects/sb_cofactor/A549")

# Load libraries
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# Load DB regions
csaw_regions = GRanges(read.table("output/analyses/ENCODE-csaw-results.txt", sep="\t", header=TRUE))

# Load and annotate HiC contacts
hic_contacts = read.table("input/ENCFF385DHX.tsv", sep="\t", header=TRUE)

hic_left = GRanges(data.frame(seqnames=hic_contacts$chr1, start=hic_contacts$x1, end=hic_contacts$x2))
hic_right = GRanges(data.frame(seqnames=hic_contacts$chr2, start=hic_contacts$y1, end=hic_contacts$y2)) 

hic_left_annot = as.data.frame(annotatePeak(hic_left, tssRegion = c(-1000, 1000), 
                                            TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, level="gene",
                                            annoDb="org.Hs.eg.db", addFlankGeneInfo = TRUE))
hic_right_annot = as.data.frame(annotatePeak(hic_right, tssRegion = c(-1000, 1000), 
                                             TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, level="gene",
                                             annoDb="org.Hs.eg.db", addFlankGeneInfo = TRUE))

# Match DB regions and contacts
overlap_left = findOverlaps(csaw_regions, hic_left)
overlap_right = findOverlaps(csaw_regions, hic_right)
                                             

# Function to associate a DB region with the gene it's in contact with.
get_contact_gene <- function(overlap_results, hic_main, hic_other, query_length) {
    # By default, no association.
    results = rep(NA, query_length)
    
    # Loop over all DB-region -> HiC contact overlaps.
    for(i in 1:length(overlap_results)) {
        # Get indices of overlap.
        query_index = queryHits(overlap_results)[i]
        hic_index = subjectHits(overlap_results)[i]
        
        # See if the HiC regions overlap a gene's TSS.
        main_symbol = ifelse(hic_main$annotation[hic_index]=="Promoter", hic_main$SYMBOL[hic_index], NA)
        other_symbol = ifelse(hic_other$annotation[hic_index]=="Promoter", hic_other$SYMBOL[hic_index], NA)
        
        # If the DB-region didn't already have a gene assignment, assign the one
        # we found.
        if(is.na(results[query_index])) {
            results[query_index] = ifelse(!is.na(other_symbol), other_symbol, main_symbol)
        }
    }

    return(results)
}
              
get_contact_gene_list <- function(overlap_results, hic_main, hic_other, query_length) {
    # By default, no association.
    results = rep(list(c()), query_length)
    
    # Loop over all DB-region -> HiC contact overlaps.
    for(i in 1:length(overlap_results)) {
        # Get indices of overlap.
        query_index = queryHits(overlap_results)[i]
        hic_index = subjectHits(overlap_results)[i]
        
        # See if the HiC regions overlap a gene's TSS.
        main_symbol = ifelse(hic_main$annotation[hic_index]=="Promoter", hic_main$SYMBOL[hic_index], NA)
        other_symbol = ifelse(hic_other$annotation[hic_index]=="Promoter", hic_other$SYMBOL[hic_index], NA)
        
        # If the DB-region didn't already have a gene assignment, assign the one
        # we found.
        if(!is.na(main_symbol)) {
            results[[query_index]] = c(results[[query_index]], main_symbol)
        }
        
         if(!is.na(other_symbol)) {
            results[[query_index]] = c(results[[query_index]], other_symbol)
         }
    }

    return(results)
}
              
# Associate DB regions to genes:                               
direct_contact = rep(NA, length(csaw_regions))
direct_contact[csaw_regions$annotation == "Promoter"] = as.character(csaw_regions$SYMBOL[csaw_regions$annotation == "Promoter"])

left_contact = get_contact_gene_list(overlap_left, hic_left_annot, hic_right_annot, length(csaw_regions))
right_contact = get_contact_gene_list(overlap_right, hic_right_annot, hic_left_annot, length(csaw_regions))

# Summarize DB -> gene associations into a single one by priority.
csaw_regions$DirectContact = direct_contact
csaw_regions$LeftContact = left_contact
csaw_regions$RightContact = right_contact
csaw_regions$GeneAssignment = ifelse(!is.na(direct_contact), direct_contact,
                              ifelse(!is.na(left_contact), left_contact,
                              ifelse(!is.na(right_contact), right_contact, NA)))
                             
de_results = read.table("results/a549_dex_time_points/6h", header=TRUE, sep=",")

db_genes = na.omit(unique(csaw_regions$GeneAssignment))
de_genes = na.omit(unique(de_results$symbol[de_results$padj <= 0.05 & abs(de_results$log2FoldChange) >= 1]))
#de_genes = na.omit(unique(de_results$symbol[de_results$padj <= 0.05]))
grid.draw(venn.diagram(list(DB=db_genes, DE=de_genes), filename=NULL))
     

sorted_regions = csaw_regions[order(abs(csaw_regions$distanceToTSS)),]     
sorted_regions = csaw_regions[order(abs(csaw_regions$distanceToTSS)),]
de_results$ClosestPeak = sorted_regions$distanceToTSS[match(de_results$symbol, sorted_regions$SYMBOL)]

mcols(csaw_regions) = NULL
csaw_regions = as.data.frame(annotatePeak(csaw_regions, tssRegion = c(-3000, 3000), 
                                          TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, level="gene",
                                          annoDb="org.Hs.eg.db", addFlankGeneInfo = TRUE))
csaw_regions = GRanges(csaw_regions)