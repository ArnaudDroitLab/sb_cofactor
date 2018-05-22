###############################################################################
# Sets of function for loading and interrogating the time-course data produced
# by Tim Reddy's lab.
###############################################################################

library(ENCODExplorer)
library(ef.utils)
library(GenomicFeatures)
library(org.Hs.eg.db)

###############################################################################
# Functions to load and query GR-binding data.
###############################################################################

# Download GR-binding peak calls at 16 time points. Each time-point is in 
# triplicate. Only regions identified by two replicates are kept.
# Results are provided as a named-list of GRanges object, in chronological
# order.
load_reddy_gr_binding_consensus <- function(diagnostic_dir=NULL) {
    # Get ENCODE accessions for GR binding time series
    gr_accession = read.table("input/ENCODE_Reddy_GR_ChIP_experiments.txt", header=TRUE, sep="\t")
    gr_accession = gr_accession[order(gr_accession$Order),]
    gr_accession$Time = factor(gr_accession$Time, levels=gr_accession$Time)
    
    chip_dir = "input/ENCODE/A549/GRCh38/chip-seq/narrow"
    dir.create(chip_dir, recursive=TRUE, showWarnings=FALSE)
    
    # Loop over all ENCODE accessions, downloading all replicates and 
    # extracting consensus regions.
    gr_regions = list()
    for(i in 1:nrow(gr_accession)) {
        encode_accession = gr_accession$Experiment[i]
        time_point = as.character(gr_accession$Time[i])
        
        # Download and import peak calls for the time point.
        encodeResults = ENCODExplorer::queryEncodeGeneric(accession=encode_accession, file_type="bed narrowPeak", assembly="GRCh38")
        downloaded_files = ENCODExplorer::downloadEncode(encodeResults, dir=chip_dir, force=FALSE)
        names(downloaded_files) = gsub(".*\\/(.*).bed.gz", "\\1", downloaded_files)
        
        extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
        binding_sites = GRangesList(lapply(downloaded_files, rtracklayer::import, extraCols=extraCols))
        
        # Only keep peak calls found in at least two replicates.
        intersect_object = build_intersect(binding_sites)
        two_replicates = rowSums(intersect_object$Matrix) >= 2
        gr_regions[[time_point]] = intersect_object$Regions[two_replicates]
        if(!is.null(diagnostic_dir)) {
            # Turn off VennDiagram logging.
            futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
        
            dir.create(file.path(diagnostic_dir, "Venn diagrams of peak calls"), recursive=TRUE, showWarnings=FALSE)
            intersect_venn_plot(intersect_object, file.path(diagnostic_dir, "Venn diagrams of peak calls", paste0(time_point, ".tiff")))
        }
    
        # Get rid of incomplete/alternate scaffolds, since they interfere with annotation later.
        gr_regions[[time_point]] = gr_regions[[time_point]][!grepl("_", seqnames(gr_regions[[time_point]]))]
    }
    
    return(gr_regions)
}

# Loads GR-binding data the same way load_reddy_gr_binding_consensus does,
# but create an intersection object out of all binding regions and annotate
# them.
load_reddy_gr_binding_intersect <- function(diagnostic_dir=NULL, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene) {
    gr_regions = load_reddy_gr_binding_consensus(diagnostic_dir)

    # Get ENCODE accessions for GR binding time series
    intersect_all = build_intersect(GRangesList(gr_regions))
    all_regions_annotation = ChIPseeker::annotatePeak(intersect_all$Regions, TxDb=TxDb)
    intersect_all$Regions = GRanges(as.data.frame(all_regions_annotation))

    # Add gene id columns for unified post-processing.
    if(grepl("^ENSG", intersect_all$Regions$geneId[1])) {
        # Ids are ENSEMBL. Add ENTREZ ids.
        intersect_all$Regions$ENSEMBLGENE = intersect_all$Regions$geneId
        intersect_all$Regions$ENTREZID = mapIds(org.Hs.eg.db, intersect_all$Regions$geneId, column="ENTREZID", keytype="ENSEMBL")
    } else {
        # Ids are ENTREZ. Add ENSEMBL ids.
        intersect_all$Regions$ENSEMBLGENE = mapIds(org.Hs.eg.db, intersect_all$Regions$geneId, column="ENSEMBL", keytype="ENTREZID")
        intersect_all$Regions$ENTREZID = intersect_all$Regions$geneId
    }
    
    return(intersect_all)
}

# Returns a vector of ENTREZ gene ids of genes whose promoter is 
# GR-bound (1kb from start site) at a given time point.
get_gr_bound_genes_at_time_point <- function(gr_intersect, gr_time, id="geneId", distance=1000) {
    gr_ranges = gr_intersect$Regions[gr_intersect$List[[as.character(gr_time)]]]
    gr_promoters = gr_ranges[gr_ranges$distanceFromTSS <= distance]
    
    # Get the total number of GR-bound genes
    gr_genes = mcols(gr_promoters)[[id]]
    
    return(gr_genes)
}

# Returns a vector of ENTREZ gene ids of genes whose promoter is 
# GR-bound at or before a given time point.
get_gr_bound_genes_at_or_before_time_point <- function(gr_intersect, gr_time, id="geneId", distance=1000) {
    # Find all applicable time points
    before_levels = 1:which(as.character(gr_time)==gr_intersect$Name)
    before_times = gr_intersect$Name[before_levels]
    
    before_regions = c()
    for(i in before_times) {
        before_regions = c(before_regions, gr_intersect$List[[as.character(i)]])
    }
    
    gr_ranges = gr_intersect$Regions[unique(before_regions)]
    gr_promoters = gr_ranges[gr_ranges$distanceFromTSS <= distance]
    
    # Get the total number of GR-bound genes
    gr_genes = mcols(gr_promoters)[[id]]
    
    return(gr_genes)
}

# Returns a vector of ENTREZ gene ids of genes whose promoter is 
# GR-bound at any time point.
get_gr_bound_genes_any_time_point <- function(gr_intersect, gr_time, id="geneId", distance=1000) {
    gr_promoters = gr_intersect$Regions[gr_intersect$Regions$distanceFromTSS <= distance]
    
    # Get the total number of GR-bound genes
    gr_genes = mcols(gr_promoters)[[id]]
    
    return(gr_genes)
}

###############################################################################
# Functions to load and query differentiall expressed genes data.
# DE results are a reanalysis of the counts provided in ENCODE for the
# Reddy experiments.
###############################################################################

# Loads all DE results and return them as a named list.
# Each element is named after its time-point, and is itself a list with
# four elements: Full, DE, Up and Down.
load_reddy_de_list <- function() {
     # Loop over all result files.
    de_results = list()
    for(de_result_file in Sys.glob("results/a549_dex_time_points/*h")) {
        time_point = gsub(".*points\\/(.*h)$", "\\1", de_result_file)
        de_results[[time_point]] = list()
        de_results[[time_point]]$Full = read.csv(de_result_file)
        de_results[[time_point]]$Full$ENTREZID = mapIds(org.Hs.eg.db, keys=as.character(de_results[[time_point]]$Full$gene_id), keytype="ENSEMBL", column="ENTREZID")
        de_results[[time_point]]$DE_indices = with(de_results[[time_point]]$Full, abs(log2FoldChange) > log2(1.5) & padj < 0.05)
        de_results[[time_point]]$DE = de_results[[time_point]]$Full[de_results[[time_point]]$DE_indices,]
        
        de_results[[time_point]]$Up_indices = with(de_results[[time_point]]$Full, log2FoldChange < -log2(1.5) & padj < 0.05)
        de_results[[time_point]]$Up = de_results[[time_point]]$Full[de_results[[time_point]]$Up_indices,]
        
        de_results[[time_point]]$Down_indices = with(de_results[[time_point]]$Full, log2FoldChange > log2(1.5) & padj < 0.05)
        de_results[[time_point]]$Down = de_results[[time_point]]$Full[de_results[[time_point]]$Down_indices,]    
    }    
    
    return(de_results)
}

# Return a list of DE genes for a given time point and direction (Up, Down,
# or "DE" for all). Can return eiher ENTREZID (id_type=ENTREZID)or ENSEMBL ids
# (id_type=ENSEMBLID)
get_de_genes <- function(de_results, time_point, direction="DE", id_type="ENTREZID") {
    if(id_type=="ENSEMBLID") {
        id_type="gene_id"
    }
    return(de_results[[time_point]][[direction]][[id_type]])
}

# Return a data-frame including the fold-changes of all DE genes.
get_reddy_fc_dataframe <- function(de_results) {
    all_de = unique(unlist(lapply(de_results, function(x) { c(as.character(x[["Up"]]$ENTREZID), as.character(x[["Down"]]$ENTREZID)) })))
    all_fc = data.frame(ENTREZID=all_de)
    for(de_item in names(de_results)) {
        cur_fc = de_results[[de_item]]$Full
        all_fc[[de_item]] = cur_fc$log2FoldChange[match(all_fc$ENTREZID, cur_fc$ENTREZID)]
    }   
    return(all_fc)
}


load_most_expressed_transcripts <- function() {
    return(GRanges(read.table("output/analyses/tss_gene_coordinates.txt", sep="\t", header=TRUE)))
}

load_most_expressed_TxDb <- function() {
    if(!exists("most_expressed")) {
        most_expressed = load_most_expressed_transcripts()
    }

    cache_path = "output/analyses/most_expressed_TxDb.RData"
    if(!file.exists(cache_path)) {
        most_expressed_TxDb = makeTxDbFromBiomart(transcript_ids=as.character(most_expressed$ensembl_transcript_id), host="useast.ensembl.org")
        AnnotationDbi::saveDb(most_expressed_TxDb, file=cache_path)
    } else {
        most_expressed_TxDb = AnnotationDbi::loadDb(cache_path)
    }
    
    return(most_expressed_TxDb)
}

get_gene_bodies <- function(all_genes, ..., gene_id="entrezgene") {
    if(length(list(...)) == 0) {
        return(all_genes)
    } else {
        return(all_genes[mcols(all_genes)[[gene_id]] %in% Reduce(intersect, list(...))])
    }
}

get_tss <- function(all_genes, ..., flank_size=200) {
    return(GenomicRanges::promoters(get_gene_bodies(all_genes, ...), upstream=flank_size, downstream=flank_size))
}

# Get GR peaks close to gene list
get_peak_near_gene <- function() {

}

load_cofactor_binding <- function(consensus_method="union", diagnostic_dir=NULL) {
    # Import all peaks
    all_peaks = import_into_grl(input.dir="output/chip-pipeline-GRCh38/", file.format="narrow", dir.type="mugqic")

    results = list()
    for(condition in c("DEX", "CTRL")) {
        condition_prefix = paste0("A549_", condition, "_")
        if(consensus_method=="pool") {
            # Remove single replicates
            results[[condition]] = all_peaks[!grepl("rep", names(all_peaks)) & grepl(condition, names(all_peaks))]
            names(results[[condition]]) = gsub(condition_prefix, "", names(results[[condition]]))
        } else {
            # Remove pooled peak calls
            relevant_peaks = all_peaks[grepl("rep", names(all_peaks)) & grepl(condition, names(all_peaks))]
            cofactors = gsub("_rep.*", "", gsub(condition_prefix, "", names(relevant_peaks)))
            results[[condition]] = GRangesList()

            # Loop over cofactors, and either intersect/union both replicates.
            for(cofactor in unique(cofactors)) {
                if(consensus_method=="intersect") {
                    # Remove the prefix within the name so that the venn diagrams
                    # will be nicer.
                    cofactor_peaks = relevant_peaks[cofactors==cofactor]
                    names(cofactor_peaks) = gsub(".*_(rep.)", "\\1", names(cofactor_peaks))

                    # Build an intersect object out of both replicates, and plot a venn
                    # diagram if requested.
                    intersect_obj = build_intersect(cofactor_peaks)
                    if(!is.null(diagnostic_dir)) {
                       venn_label = paste(cofactor, condition)
                       intersect_venn_plot(intersect_obj,
                                            file=file.path(diagnostic_dir, paste0(venn_label, ".tiff")),
                                            title=venn_label)
                    }

                    # Return the overlap between replicates.
                    results[[condition]][[cofactor]] = intersect_overlap(intersect_obj)
                } else {    # Mode is presumed to be "union".
                    results[[condition]][[cofactor]] = reduce(unlist(relevant_peaks[cofactors==cofactor]))
                }
            }
        }
    }
    
    return(results)
}

