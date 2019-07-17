library(InteractionSet)
library(purrr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(FantomEnhancers.hg19)
library(ENCODExplorer)
source("scripts/load_reddy.R")

####### Step 1: Load and annotate HiC data. 

load_long_range_interactions = function(input_file) {
    raw_hic_0h = read.table(input_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    left = GRanges(raw_hic_0h[,1:3] %>% set_names("seqnames", "start", "end"))
    right = GRanges(raw_hic_0h[,4:6] %>% set_names("seqnames", "start", "end"))
    
    return(GInteractions(left, right, mode="strict"))
}

load_bedpe_tad = function(input_file) {
    tad_raw = read.table(input_file, sep="\t", header=TRUE)
    tad_gr = GRanges(tad_raw[,1:3] %>% set_names("seqnames", "start", "end"))
    
    return(tad_gr)
}

load_encode_file_internal = function(filetype, load_fun, dir=".") {
    # Retrieve all of Tim Reddy's HiC results from ENCODE
    hic_df = ENCODExplorer::queryEncodeGeneric(biosample_name="A549",
                                               assay="Hi-C", 
                                               assembly="GRCh38", 
                                               lab="Tim Reddy, Duke",
                                               output_type=filetype)
                                                  
    # Download the long-range interaction files.                                                  
    dir.create(dir, recursive=TRUE, showWarnings=FALSE)
    encode_files = ENCODExplorer::downloadEncode(hic_df$file_accession,
                                                 dir=dir, force=FALSE)
    loaded_files = lapply(encode_files, load_fun)
    names(loaded_files) = paste(hic_df$treatment_duration,
                                hic_df$treatment_duration_unit)
    names(loaded_files)[is.na(hic_df$treatment_duration)] = "Ctrl"
    
    return(loaded_files)
}

load_reddy_hic = function() {
    # Download the long-range interaction files.  
    lri_gi = load_encode_file_internal(filetype="long range chromatin interactions",
                                       load_fun=load_long_range_interactions,
                                       dir="input/ENCODE/A549/hg38/hic long range interactions")
    
    # Download TAD files                      
    tad_gr = load_encode_file_internal(filetype="topologically associated domains",
                                       load_fun=load_bedpe_tad,
                                       dir="input/ENCODE/A549/hg38/hic TAD")                      
    # Group by time point                                              
    res=list()
    for(time_point in unique(names(lri_gi))) {
        res[[time_point]]=list(LRI=lri_gi[[time_point]],
                               TAD=tad_gr[[time_point]])
    }
    
    return(res)
}

overlap_early_binding = function(gr_regions, target_regions) {
    for(i in names(gr_regions)) {
        gr_overlap = countOverlaps(target_regions, gr_regions[[i]])
        mcols(target_regions)[[i]] = gr_overlap
    }
    
    return(target_regions)
}

identify_early_binding = function(gr_regions, target_regions) {
    target_regions = overlap_early_binding(gr_regions, target_regions)
    
    early_columns = names(gr_regions)[grepl("minutes", names(gr_regions))]
    early_matrix = as.matrix(mcols(target_regions)[, early_columns])
    early_binding = apply(early_matrix > 0, 1, sum)
    
    return(early_binding)
}

load_hg38_fantom_enhancers = function() {
    hg19_enhancers = get_fantom_enhancers_tpm(cell_lines = "A549")
    chain = rtracklayer::import.chain("input/hg19ToHg38.over.chain")
    hg38_enhancers = rtracklayer::liftOver(hg19_enhancers, chain)
    
    return(hg38_enhancers)
}             

annotate_regions_with_enh_gr = function(query_regions, 
                            ehn_regions=load_hg38_fantom_enhancers(),
                            gr_regions=load_reddy_gr_binding_consensus()) {
    enh_overlap = countOverlaps(query_regions, ehn_regions)
    query_regions$EarlyBinding = identify_early_binding(gr_regions, query_regions)
    query_regions = overlap_early_binding(gr_regions, query_regions) 

    return(query_regions)
}

annotate_genes_with_reddy_de = function(query_tss) {
    # Add DE-status annotation
    time_points = Sys.glob("results/a549_dex_time_points/*")
    regex_match = grepl("^[0-9\\.]+[hm]$", basename(time_points))
    time_points = time_points[regex_match]
    names(time_points) = basename(time_points)
    
    de_list = lapply(time_points, function(x) { 
        raw_data = read.table(x, sep=",", header=TRUE)
        de_cat = ifelse(raw_data$padj > 0.05 | is.na(raw_data$padj), "Stable", 
                    ifelse(raw_data$log2FoldChange <= -log2(1.5), "Down",
                        ifelse(raw_data$log2FoldChange >= log2(1.5), "Up", "Stable")))
        raw_data$DE = de_cat
        
        raw_data
    })
    
    de_list = lapply(de_list, function(x) {
        x[match(query_tss$ensembl_gene_id, x$gene_id),]
    })
    
    de_df = as.data.frame(lapply(de_list, '[[', "DE"))
    de_df[is.na(de_df)] <- "Stable"
    names(de_df) = gsub("X", "DE-", names(de_df))
    
    mcols(query_tss) = cbind(mcols(query_tss), de_df)
    
    return(query_tss)
}
             
annotate_distant_gr_binding = function(promoter_regions, input_hic, input_TAD) {
    # 3a. Identify genes with GR binding in their 5/10K window.
    gene_overlaps = findOverlaps(promoter_regions, regions(input_hic))
    
    promoter_regions$WindowBound = FALSE
    hit_binding = regions(input_hic)[subjectHits(gene_overlaps)]$EarlyBinding > 0
        
    hit_df = data.frame(query=queryHits(gene_overlaps),
                   subject=subjectHits(gene_overlaps),
                   binding=hit_binding)
                   
    bind_status = hit_df %>% group_by(query) %>% summarize(Bound=any(binding))
    promoter_regions$WindowBound[bind_status$query] = bind_status$Bound
    
    # 3b. Identify genes with GR binding in contacts.
    promoter_hic_overlap = findOverlaps(input_hic, promoter_regions)
    
    f_a = anchors(input_hic, type="first")
    s_a = anchors(input_hic, type="second")
    GR_in_interaction = f_a$EarlyBinding | s_a$EarlyBinding
    
    promoter_regions$ContactBound = FALSE
    promoter_regions$ContactBound[subjectHits(promoter_hic_overlap)] = GR_in_interaction[queryHits(promoter_hic_overlap)]
    
    # 3c. Identify genes with GR binding in TADs
    input_TAD$EarlyBinding = identify_early_binding(gr_regions, input_TAD)
    
    tad_overlap = findOverlaps(promoter_regions, input_TAD)
    promoter_regions$TADBound = FALSE
    promoter_regions$TADBound[queryHits(tad_overlap)] = input_TAD$EarlyBinding[subjectHits(tad_overlap)] > 0
    
    return(promoter_regions)
}
get_consensual_early_gr_binding = function(gr_regions) {
    # Build list of GR sites which are "consensual" in first 30 minutes 
    early_gr = GRangesList(gr_regions[grepl("minutes", names(gr_regions))])
    gr_intersect = build_intersect(early_gr)
    consensual_binding = apply(gr_intersect$Matrix > 0, 1, sum) >= 4
    consensual_gr_regions = gr_intersect$Regions[consensual_binding]

    return(consensual_gr_regions)
}

# Define a function to go from promoter/enhancer indices to a GRangesList
# of GR sites.
build_list_of_gr_sites = function(promoter_regions, 
                                  promoter_ids, 
                                  gr_regions,
                                  gr_ids) {
    unique_query = unique(promoter_ids)
    gr_by_gene = lapply(unique_query, function(x) {
        gr_regions[gr_ids[promoter_ids==x] ]
    })
    names(gr_by_gene) = promoter_regions$ensembl_gene_id[unique_query]
    gr_by_gene = GRangesList(gr_by_gene)

    gr_by_gene
}

find_per_gene_gr_sites = function(promoter_regions, gr_regions, input_hic) {
    consensual_gr_regions = get_consensual_early_gr_binding(gr_regions)
        
    # Direct binding
    promoter_gr_overlap = findOverlaps(promoter_regions, consensual_gr_regions)
    gr_by_gene_direct = build_list_of_gr_sites(promoter_regions,
                                               queryHits(promoter_gr_overlap),
                                               consensual_gr_regions,
                                               subjectHits(promoter_gr_overlap))
    
    # In-window binding
    gene_overlaps = findOverlaps(promoter_regions, regions(input_hic))
    window_gr_overlap = findOverlaps(regions(input_hic), consensual_gr_regions)
    overlap_df = inner_join(as.data.frame(gene_overlaps), 
                            as.data.frame(window_gr_overlap),
                            c(subjectHits="queryHits"))
    gr_by_gene_window = build_list_of_gr_sites(promoter_regions,
                                               overlap_df$queryHits,
                                               consensual_gr_regions,
                                               overlap_df$subjectHits.y)
    
    # Distant binding
    promoter_gr_distant = linkOverlaps(input_hic, promoter_regions, consensual_gr_regions)
    gr_by_gene_indirect = build_list_of_gr_sites(promoter_regions,
                                                 promoter_gr_distant$subject1,
                                                 consensual_gr_regions,
                                                 promoter_gr_distant$subject2)
    
    unique_genes = unique(c(names(gr_by_gene_direct), 
                            names(gr_by_gene_window), 
                            names(gr_by_gene_indirect)))
    
    # All binding
    gr_by_gene_any = list()
    for(gene in unique_genes) {
        grl = unlist(GRangesList(c(gr_by_gene_direct[[gene]],
                                   gr_by_gene_window[[gene]],
                                   gr_by_gene_indirect[[gene]])))
        gr_by_gene_any[[gene]] = grl
    }

    return(gr_by_gene_any)
}

load_connection_data = function(hic_timepoint) {
     all_hic = load_reddy_hic()
     hic_res = all_hic[[hic_timepoint]]$LRI
     
     # Add enhancer annotation to HiC interactions.
     hg38_enhancers = load_hg38_fantom_enhancers()
     gr_regions = load_reddy_gr_binding_consensus()
     
     regions(hic_res) = annotate_regions_with_enh_gr(regions(hic_res),
                                                    hg38_enhancers,
                                                    gr_regions)
     
     # Load promoter regions and annotate with GR binding, DE status.
     promoter_regions = load_annotated_most_expressed_promoters(fix_chr=TRUE)
     
     promoter_regions = annotate_regions_with_enh_gr(promoter_regions,
                                                    hg38_enhancers,
                                                    gr_regions)
     promoter_regions = annotate_genes_with_reddy_de(promoter_regions)
     
     # Identify gene categories based on interactions.
     promoter_regions = annotate_distant_gr_binding(promoter_regions, 
                                                    hic_res,
                                                    all_hic[[hic_timepoint]]$TAD)

     # Identify per-gene GR sites.
     gr_by_gene_any = find_per_gene_gr_sites(promoter_regions, gr_regions, hic_res)
 
    return(list(HiC=all_hic,
                Promoters=promoter_regions,
                Enhancers=hg38_enhancers,
                GR_all=gr_regions,
                GR_by_gene))
}
     
motif_by_de = function(promoter_regions, gr_by_gene, de_time, de_class) {
    de_genes = promoter_regions$ensembl_gene_id[promoter_regions[[de_time]] %in% de_class]
    de_gr = gr_by_gene[names(gr_by_gene) %in% de_genes]
    de_gr_gr = unlist(GRangesList(unlist(de_gr)))

    res = motif_enrichment(de_gr_gr, select_annotations("hg38"))

    return(res)
}

# rmarkdown::render("scripts/hic_analysis.Rmd", knit_root_dir=getwd(), output_format="html_document", output_file="output/hic_analysis.html", output_dir=getwd())

# res_down = motif_by_de("2h", "Down")
# res_up = motif_by_de("2h", "Up")

# Test plot
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggbio)

visualize_region <- function(input_region, promoter_regions, gr_regions, tad_gr,
                             hic_interactions, output_file) {
    input_region = promoter_regions[promoter_regions$external_gene_name=="KLF6"]
    #start(input_region) = start(input_region) - 200000
    #end(input_region) = end(input_region) + 200000

    start(input_region) = 3725000
    end(input_region) = 4780000


    gene_models = autoplot(TxDb.Hsapiens.UCSC.hg38.knownGene, which=input_region)
    TAD_track = autoplot(subsetByOverlaps(tad_gr, input_region), which=input_region)
    gr_subset = lapply(gr_regions, subsetByOverlaps, ranges=input_region)
    gr_tracks = lapply(gr_subset, function(x) { autoplot(x, geom="rect", which=input_region)})
    names(gr_tracks) = paste0("GR-", names(gr_tracks))
    names(gr_tracks) = gsub(" minutes", "m", names(gr_tracks))
    names(gr_tracks) = gsub(" hours*", "h", names(gr_tracks))

    GR_track = autoplot(GRangesList(gr_subset), geom="rect", which=input_region)
    relevant_interactions = subsetByOverlaps(hic_interactions, input_region)
    relevant_interactions = relevant_interactions[countOverlaps(anchors(relevant_interactions, type='first'), input_region) > 0 &
                                                  countOverlaps(anchors(relevant_interactions, type='second'), input_region) > 0]
    rf_a = anchors(relevant_interactions, type='first')
    rs_a = anchors(relevant_interactions, type='second')
    interaction_gr = data.frame(seqnames=as.character(seqnames(rf_a)),
                             start=start(rf_a) + floor(width(rf_a) / 2),
                             end=start(rs_a) + floor(width(rs_a) / 2),
                             strand='+')

    interaction_plot = autoplot(GRanges(interaction_gr), geom="arch")

    relevant_ids = unlist(anchorIds(relevant_interactions))
    relevant_regions = unique(regions(relevant_interactions)[relevant_ids])
    interaction_regions_plot = autoplot(relevant_regions, which=input_region)

    relevant_promoters = subsetByOverlaps(promoter_regions, input_region)
    promoter_track = autoplot(relevant_promoters, which=input_region)

    track_list = list(HiC=interaction_plot,
                      HiCR=interaction_regions_plot,
                      Genes=gene_models, 
                      Prom=promoter_track,
                      TADs=TAD_track)
    track_list = c(track_list, gr_tracks)
    track_args = c(track_list, 
                   list(heights=c(1, 1, 3, 1, 1, rep(0.5, length(gr_tracks))),
                        xlim=input_region),
                        label.text.angle = 0)
    pdf(output_file, height=14, width=7)
    do.call(tracks, track_args)
    dev.off()
}

# KLF6_region = promoter_regions[promoter_regions$external_gene_name=="KLF6"]
# #start(KLF6_region) = start(KLF6_region) - 200000
# #end(KLF6_region) = end(KLF6_region) + 200000
# 
# start(KLF6_region) = 3725000
# end(KLF6_region) = 4780000
