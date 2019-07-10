library(InteractionSet)
library(purrr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(FantomEnhancers.hg19)

source("scripts/load_reddy.R")

####### Step 1: Load and annotate HiC data. 

# Load HiC data.
raw_hic_0h = read.table("input/ENCFF803ZOW.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)

left = GRanges(raw_hic_0h[,1:3] %>% set_names("seqnames", "start", "end"))
right = GRanges(raw_hic_0h[,4:6] %>% set_names("seqnames", "start", "end"))

hic_0h = GInteractions(left, right, mode="strict")

# Add enhancer annotation to HiC interactions.
hg19_enhancers = get_fantom_enhancers_tpm(cell_lines = "A549")
chain = rtracklayer::import.chain("input/hg19ToHg38.over.chain(2)")
hg38_enhancers = rtracklayer::liftOver(hg19_enhancers, chain)

enh_overlap = countOverlaps(regions(hic_0h), hg38_enhancers)
regions(hic_0h)$Enhancer = enh_overlap

# Add GR annotation to HiC interactions.
gr_regions = load_reddy_gr_binding_consensus()

identify_early_binding = function(gr_regions, target_regions) {
    for(i in names(gr_regions)) {
        gr_overlap = countOverlaps(target_regions, gr_regions[[i]])
        mcols(target_regions)[[i]] = gr_overlap
    }
    
    early_columns = names(gr_regions)[grepl("minutes", names(gr_regions))]
    early_matrix = as.matrix(mcols(target_regions)[, early_columns])
    early_binding = apply(early_matrix > 0, 1, sum)
    
    return(early_binding)
}

regions(hic_0h)$EarlyBinding = identify_early_binding(gr_regions, regions(hic_0h))



####### Step 2: Load and annotate promoter regions. 

# Load promoter regions and annotate with GR binding.
promoter_regions = load_annotated_most_expressed_promoters(fix_chr=TRUE)

# Identify genes with GR binding in their promoter region.
promoter_regions$EarlyBinding = identify_early_binding(gr_regions, promoter_regions)



####### Step 3: Identify gene categories based on interactions.

# 3a. Identify genes with GR binding in their 5/10K window.
gene_overlaps = findOverlaps(promoter_regions, regions(hic_0h))

promoter_regions$WindowBound = FALSE
hit_binding = regions(hic_0h)[subjectHits(gene_overlaps)]$EarlyBinding > 0
    
hit_df = data.frame(query=queryHits(gene_overlaps),
               subject=subjectHits(gene_overlaps),
               binding=hit_binding)
               
bind_status = hit_df %>% group_by(query) %>% summarize(Bound=any(binding))
promoter_regions$WindowBound[bind_status$query] = bind_status$Bound

# 3b. Identify genes with GR binding in contacts.
promoter_hic_overlap = findOverlaps(hic_0h, promoter_regions)

f_a = anchors(hic_0h, type="first")
s_a = anchors(hic_0h, type="second")
GR_in_interaction = f_a$EarlyBinding | s_a$EarlyBinding

promoter_regions$ContactBound = FALSE
promoter_regions$ContactBound[subjectHits(promoter_hic_overlap)] = GR_in_interaction[queryHits(promoter_hic_overlap)]

# 3c. Identify genes with GR binding in TADs
tad_raw = read.table("input/hic_dex.t0.tads.bedpe", sep="\t", header=TRUE)
tad_gr = GRanges(tad_raw[,1:3] %>% set_names("seqnames", "start", "end"))
tad_gr$EarlyBinding = identify_early_binding(gr_regions, tad_gr)

tad_overlap = findOverlaps(promoter_regions, tad_gr)
promoter_regions$TADBound = FALSE
promoter_regions$TADBound[queryHits(tad_overlap)] = tad_gr$EarlyBinding[subjectHits(tad_overlap)] > 0

# Assess intra-TAD vs inter-TAD contacts.
f_tad = findOverlaps(f_a, tad_gr)
s_tad = findOverlaps(s_a, tad_gr)

tad_df = data.frame(f=rep(NA, length(f_a)), s=rep(NA, length(s_a)))
tad_df$f[queryHits(f_tad)] = subjectHits(f_tad)
tad_df$s[queryHits(s_tad)] = subjectHits(s_tad)

# Test plot
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

KLF6_region = promoter_regions[promoter_regions$external_gene_name=="KLF6"]
start(KLF6_region) = start(KLF6_region) - 200000
end(KLF6_region) = end(KLF6_region) + 200000

gene_models = autoplot(TxDb.Hsapiens.UCSC.hg38.knownGene, which=KLF6_region)
TAD_track = autoplot(subsetByOverlaps(tad_gr, KLF6_region), which=KLF6_region)
gr_subset = lapply(gr_regions, subsetByOverlaps, ranges=KLF6_region)
GR_track = autoplot(GRangesList(gr_subset), which=KLF6_region)
relevant_interactions = subsetByOverlaps(hic_0h, KLF6_region)
relevant_interactions = relevant_interactions[countOverlaps(anchors(relevant_interactions, type='first'), KLF6_region) > 0 &
                                              countOverlaps(anchors(relevant_interactions, type='second'), KLF6_region) > 0]
rf_a = anchors(relevant_interactions, type='first')
rs_a = anchors(relevant_interactions, type='second')
interaction_gr = data.frame(seqnames=as.character(seqnames(rf_a)),
                         start=start(rf_a) + floor(width(rf_a) / 2),
                         end=start(rs_a) + floor(width(rs_a) / 2),
                         strand='+')
                         
interaction_plot = autoplot(GRanges(interaction_gr), geom="arch")

tracks(gene_models, TAD_track, GR_track, interaction_plot)
