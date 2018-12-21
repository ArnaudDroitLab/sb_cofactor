library(ef.utils)
library(dplyr)

source("scripts/load_reddy.R")

USE_MACS2_BGDIFF=TRUE
if(!USE_MACS2_BGDIFF) {
    # Load peaks.
    cofactor_peaks = load_cofactor_binding(consensus_method="replicate_1")

    # Determine losses/gains in individual cofactors.
    peak_diff = list(Gain=GRangesList(), Loss=GRangesList(), Common=GRangesList())
    for(cofactor in c("NIPBL", "BRD4", "CDK9")) {
        intersect_obj = build_intersect(GRangesList(DEX=cofactor_peaks[["DEX"]][[cofactor]], 
                                                    CTRL=cofactor_peaks[["CTRL"]][[cofactor]]))
        gains = intersect_obj$Matrix[,"DEX"] > 0 & intersect_obj$Matrix[,"CTRL"] == 0
        losses = intersect_obj$Matrix[,"DEX"] == 0 & intersect_obj$Matrix[,"CTRL"] > 0
        common = intersect_obj$Matrix[,"DEX"] > 0 & intersect_obj$Matrix[,"CTRL"] > 0
        peak_diff[["Gain"]][[cofactor]] = intersect_obj$Regions[gains]
        peak_diff[["Loss"]][[cofactor]] = intersect_obj$Regions[losses]
        peak_diff[["Common"]][[cofactor]] = intersect_obj$Regions[common]
    }
} else {
    peak_diff = load_macs2_diffpeaks(threshold="0.7")
}

# Identify losses/gains in all three cofactors.
all_diff = GRangesList()
for(type in c("Gain", "Loss")) {
    intersect_obj = build_intersect(peak_diff[[type]][c("NIPBL", "BRD4", "CDK9")])
    all_diff[[type]] = intersect_overlap(intersect_obj)
}




# Annotate peaks.
annotation_set = select_annotations("hg38")
most_expressed_TxDb = load_most_expressed_TxDb()
seqlevelsStyle(most_expressed_TxDb) <- "UCSC"
#annotated_regions = lapply(all_diff, annotate_region, annotations.list=annotation_set)
annotated_regions = lapply(all_diff, ChIPseeker::annotatePeak, 
                                         tssRegion=c(-3000, 3000),
                                         TxDb=most_expressed_TxDb,
                                         annoDb=annotation_set$OrgDbStr)

# Fetch name of genes where peaks are found.
bound_genes = lapply(annotated_regions, function(x) {
    as.data.frame(x) %>% 
    filter(grepl("Promoter", as.character(annotation))) %>% 
    pull(transcriptId) %>% 
    unique 
})

# Fetch transcript coordinates for chosen genes.
bound_transcripts = lapply(bound_genes, function(x) {
    transcripts(most_expressed_TxDb, filter=list(tx_name=x))
})

# Get TSSes and TESes.
bound_TSS = lapply(bound_transcripts, flank, width=500, start=TRUE, both=TRUE)
bound_TES = lapply(bound_transcripts, flank, width=500, start=FALSE, both=TRUE)

# Put both TSSes and TESes in a single list.
bound_grl = GRangesList(c(bound_TSS, bound_TES))
grl_lengths = unlist(lapply(bound_grl, length))
single_gr = unlist(bound_grl)
single_gr$GeneExtremity = rep(c("TSS", "TSS", "TES", "TES"), times=grl_lengths)
single_gr$CofactorType = rep(c("Gain", "Loss", "Gain", "Loss"), times=grl_lengths)
seqlevels(single_gr) = seqlevels(single_gr)[grepl("^chr", seqlevels(single_gr))]


# Retrieve bam file names and their sample sheet.
sample_sheet = read.table("raw/PolII/readset.txt", sep="\t", header=TRUE)
bamnames = paste0(sample_sheet$Sample, ".sorted.dup.bam")
sample_sheet$Bam = file.path("output/chip-pipeline-PolII-GRCh38/alignment", sample_sheet$Sample, bamnames)
sample_sheet$design = basename(gsub(".bam", "", sample_sheet$Bam))


# Generate the metagenes.
library(metagene2)
mg = metagene2$new(bam_files=sample_sheet$Bam, regions=single_gr, cores=4, normalization="RPM",
                   design_metadata=sample_sheet, 
                   split_by=c("GeneExtremity", "CofactorType"))
                   
dir.create("output/analyses/PolII/metagenes/", recursive=TRUE, showWarnings=FALSE)
pdf("output/analyses/PolII/metagenes/NBC Gains.pdf", width=14, height=7)
mg$produce_metagene(facet_by=Target~sh*GeneExtremity, group_by="Condition", region_filter=quo(CofactorType=="Gain"))
dev.off()

pdf("output/analyses/PolII/metagenes/NBC Losses.pdf", width=14, height=7)
mg$produce_metagene(facet_by=Target~sh*GeneExtremity, group_by="Condition", region_filter=quo(CofactorType=="Loss"))
dev.off()
