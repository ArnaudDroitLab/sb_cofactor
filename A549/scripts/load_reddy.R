###############################################################################
# Sets of function for loading and interrogating the time-course data produced
# by Tim Reddy's lab.
###############################################################################

library(ENCODExplorer)
library(ef.utils)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(knitr)
library(dplyr)
library("biomaRt")
    
###############################################################################
# Function to make ENCODE_Reddy_ChIP_experiments.txt file
###############################################################################
make_ENCODE_Reddy_ChIP_experiments_file <- function(target_name) {
    all_chip_bed <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bed", assay="ChIP-seq")
    target_bed <- all_chip_bed %>% filter(target == target_name, assembly == "GRCh38", lab == "Tim Reddy, Duke")
    report_target_bed <- target_bed %>% dplyr::select(accession, file_accession, submitted_by, file_format, file_type, target, treatment, treatment_duration, treatment_duration_unit, biological_replicates, controls)
    report_target_bed$treatment_duration[is.na(report_target_bed$treatment_duration)] <- 0
    report_target_bed[which(report_target_bed$submitted_by == "Alejandro Barrera"), ]$treatment_duration_unit[is.na(report_target_bed[which(report_target_bed$submitted_by == "Alejandro Barrera"), ]$treatment_duration_unit)] <- "minute"
    report_target_bed[which(report_target_bed$submitted_by == "Ian McDowell"), ]$treatment_duration_unit[is.na(report_target_bed[which(report_target_bed$submitted_by == "Ian McDowell"), ]$treatment_duration_unit)] <- "hour"
    report_target_bed <- report_target_bed %>% arrange(desc(treatment_duration_unit), treatment_duration, submitted_by, biological_replicates)
    print(kable(report_target_bed))
    
    exp_tmp <- report_target_bed %>% dplyr::select(accession, treatment_duration, treatment_duration_unit) %>% unique
    exp_target <- data.frame(Experiment = exp_tmp$accession, Time = paste0(exp_tmp$treatment_duration, " ", exp_tmp$treatment_duration_unit), Order = 1:nrow(exp_tmp))
    print(kable(exp_target))
    
    filename <- paste0("ENCODE_Reddy_", target_name, "_ChIP_experiments.txt")
    write.table(exp_target, file = file.path("input", filename), row.names = FALSE, quote = FALSE, sep = "\t")
    message(filename, " saved in ", file.path(getwd(), "input"))
}

###############################################################################
# Function to load and query DNA binding protein data
###############################################################################
load_reddy_binding_consensus <- function(target_name, diagnostic_dir=NULL) {
    # Get ENCODE accessions for target binding time series
    filename <- paste0("ENCODE_Reddy_", target_name, "_ChIP_experiments.txt")
    target_accession = read.table(file.path("input", filename), header=TRUE, sep="\t")
    target_accession = target_accession[order(target_accession$Order),]
    target_accession$Time = factor(target_accession$Time, levels=target_accession$Time)
    
    chip_dir = "input/ENCODE/A549/GRCh38/chip-seq/narrow"
    dir.create(chip_dir, recursive=TRUE, showWarnings=FALSE)
    
    # Loop over all ENCODE accessions, downloading all replicates and 
    # extracting consensus regions.
    target_regions = list()
    for(i in 1:nrow(target_accession)) {
        encode_accession = target_accession$Experiment[i]
        time_point = as.character(target_accession$Time[i])
        
        # Download and import peak calls for the time point.
        encodeResults = ENCODExplorer::queryEncodeGeneric(accession=encode_accession, file_type="bed narrowPeak", assembly="GRCh38")
        downloaded_files = ENCODExplorer::downloadEncode(encodeResults, dir=chip_dir, force=FALSE)
        names(downloaded_files) = gsub(".*\\/(.*).bed.gz", "\\1", downloaded_files)
        
        extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
        binding_sites = GRangesList(lapply(downloaded_files, rtracklayer::import, extraCols=extraCols))
        
        # Only keep peak calls found in at least two replicates.
        intersect_object = build_intersect(binding_sites)
        two_replicates = rowSums(intersect_object$Matrix) >= 2
        target_regions[[time_point]] = intersect_object$Regions[two_replicates]
        if(!is.null(diagnostic_dir)) {
            # Turn off VennDiagram logging.
            futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
            
            dir.create(file.path(diagnostic_dir, "Venn diagrams of peak calls"), recursive=TRUE, showWarnings=FALSE)
            intersect_venn_plot(intersect_object, file.path(diagnostic_dir, "Venn diagrams of peak calls", paste0(time_point, ".tiff")))
        }
        
        # Get rid of incomplete/alternate scaffolds, since they interfere with annotation later.
        target_regions[[time_point]] = target_regions[[time_point]][!grepl("_", seqnames(target_regions[[time_point]]))]
    }
    
    return(target_regions)
}


###############################################################################
# Functions to load and query GR-binding data.
###############################################################################

# Download GR-binding peak calls at 16 time points. Each time-point is in 
# triplicate. Only regions identified by two replicates are kept.
# Results are provided as a named-list of GRanges object, in chronological
# order.
load_reddy_gr_binding_consensus <- function(diagnostic_dir=NULL) {
    gr_regions <- load_reddy_binding_consensus("GR", diagnostic_dir)
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
        # Filter out deprecated transcript ids
        all_ids = getBM(attributes="ensembl_transcript_id", mart=useMart("ensembl", dataset="hsapiens_gene_ensembl"))
        existing_ids = intersect(most_expressed$ensembl_transcript_id, all_ids$ensembl_transcript_id)
        most_expressed_TxDb = makeTxDbFromBiomart(transcript_ids=existing_ids, host="useast.ensembl.org")
        AnnotationDbi::saveDb(most_expressed_TxDb, file=cache_path)
    } else {
        most_expressed_TxDb = AnnotationDbi::loadDb(cache_path)
    }
    
    return(most_expressed_TxDb)
}

load_most_expressed_promoters <- function(fix_chr=FALSE, upstream=1000, downstream=1000) {
    # Load most expressed TxDb and generate a promoter GRanges object.
    most_expressed_TxDb = load_most_expressed_TxDb()
    promoter_regions = promoters(most_expressed_TxDb, upstream=upstream, downstream=downstream)
    if(fix_chr) {
        promoter_regions = promoter_regions[grepl("^(\\d+|XY)$", seqnames(promoter_regions))]
        seqlevels(promoter_regions) = gsub("^", "chr", seqlevels(promoter_regions))
    }
    
    return(promoter_regions)
}

load_annotated_most_expressed_promoters <- function(fix_chr=FALSE, upstream=1000, downstream=1000) {
    promoter_regions = load_most_expressed_promoters(fix_chr, upstream, downstream)
    
    # Get annotations through biomart.
    ensembl = biomaRt::useMart("ensembl", host="uswest.ensembl.org")
    ensembl = biomaRt::useDataset("hsapiens_gene_ensembl",mart=ensembl)
    biomart_attr = c("ensembl_transcript_id", "external_gene_name",
                     "description", "hgnc_symbol", "gene_biotype",
                     "transcript_biotype")
    annotations = biomaRt::getBM(attributes=biomart_attr,
                                 filters="ensembl_transcript_id",
                                 values=names(promoter_regions),
                                 mart=ensembl)
    
    # Some transcript ids may appear more than once due to changing annotations
    annotations = annotations[!duplicated(annotations$ensembl_transcript_id),]
                        
    concat_annot = dplyr::left_join(as.data.frame(mcols(promoter_regions)),
                                    annotations,
                                    c("tx_name"="ensembl_transcript_id"))                    
    mcols(promoter_regions) = concat_annot

    return(promoter_regions)
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
            rep_number = gsub(".*_rep(.)", "\\1", names(relevant_peaks))
            
            if(consensus_method != "separate") {
                results[[condition]] = GRangesList()
            } else {
                results[[condition]] = list()
            }

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
                } else if(consensus_method=="union") { 
                    results[[condition]][[cofactor]] = reduce(unlist(relevant_peaks[cofactors==cofactor]))
                } else if(consensus_method=="separate") {
                    results[[condition]][[cofactor]] = relevant_peaks[cofactors==cofactor]
                    names(results[[condition]][[cofactor]]) = gsub(".*(rep.)", "\\1", names(results[[condition]][[cofactor]]))
                } else if(consensus_method=="replicate_1") {
                     results[[condition]][[cofactor]] = relevant_peaks[cofactors==cofactor & rep_number=="1"]
                } else {
                    stop("Invalid consensus_method")
                }
            }
        }
    }
    
    return(results)
}

load_macs2_diffpeaks <- function(threshold="1.0") {
    output_bed = Sys.glob(paste0("output/chip-pipeline-GRCh38/binding_diff/*/output_filters/*", threshold,"*.bed"))
    cofactor = gsub(".*filters/A549_(.*)_(.*)_rep1_.*Peak.?_M_(.*)_biased.*", "\\2", output_bed)
    binding_status = gsub(".*filters/A549_(.*)_(.*)_rep1_.*Peak.?_M_(.*)_biased.*", "\\3", output_bed)
    binding_status = ifelse(grepl("below", binding_status), "Loss", "Gain")
    
    imported_bed = lapply(output_bed, rtracklayer::import, format="bed")
    names(imported_bed) <- cofactor
    
    list(Gain=GRangesList(imported_bed[binding_status=="Gain"]), 
         Loss=GRangesList(imported_bed[binding_status=="Loss"]))
}


intersect_ranges <- function(range_list, min_replicate=2, diagnostic_dir=NULL, label=NULL) {
    intersect_object = ef.utils::build_intersect(range_list)
    enough_replicates = rowSums(intersect_object$Matrix) >= min_replicate
    
    if(!is.null(diagnostic_dir)) {
        # Turn off VennDiagram logging.
        futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    
        dir.create(file.path(diagnostic_dir, "Venn diagrams of peak calls"), recursive=TRUE, showWarnings=FALSE)
        ef.utils::intersect_venn_plot(intersect_object, file.path(diagnostic_dir, "Venn diagrams of peak calls", paste0(label, ".tiff")))
    }

    return(intersect_object$Regions[enough_replicates])
}

# Loads all peaks calls for a given target and build consensus regions
# for each time point.
summarize_GGR_chip <- function(encode_results, diagnostic_dir=NULL) {
    # Make sure the passed in ENCODE subset ahs only one target.
    target = unique(encode_results$target)
    stopifnot(length(target)==1)
    
    # Determine if the target has broad or narrow peaks so we can
    # load the peaks in an appropriate manner.
    if(all(unique(encode_results$file_type)=="bed narrowPeak")) {
        bed_type="narrow"
    } else if(all(unique(encode_results$file_type)=="bed broadPeak")) {
        bed_type="broad"
    } else {
        stop("Not all downloaded files are of teh same type!")
    }
    
    # Download the peak files.
    download_dir = file.path("input/ENCODE/A549/GRCh38/chip-seq", bed_type)
    dir.create(download_dir, recursive=TRUE, showWarnings=FALSE)
    downloaded_files = ENCODExplorer::downloadEncode(encode_results, dir=download_dir, force=FALSE)
    
    # Import the peak files into a GRangesList object.
    names(downloaded_files) <- gsub(".bed.gz", "", basename(downloaded_files))
    all_regions = ef.utils::import_files_into_grl(downloaded_files, file.format=bed_type, 
                                                  file.ext=".bed.gz", discard.metadata=TRUE)
    
    # Convert to GRCh38 if necessary
    for(i in which(encode_results$assembly=="hg19")) {
        hg19_to_hg38_chain = rtracklayer::import.chain("input/hg19ToHg38.over.chain")
        all_regions[[i]] = reduce(unlist(rtracklayer::liftOver(all_regions[[i]], hg19_to_hg38_chain)))
    }
    
    # # Determine the time_point for all downloaded files.
    is_ctrl = is.na(encode_results$treatment) | (encode_results$treatment=="ethanol")
    time_label = paste(encode_results$treatment_duration, encode_results$treatment_duration_unit)
    time_points = ifelse(is_ctrl, "0 minute", time_label)
    names(time_points) = encode_results$file_accession
    
    # Group files by time point, then intersect them.
    results = list()
    for(time_point in unique(time_points)) {
        time_point_regions = all_regions[names(time_points[time_points==time_point])]
        if(length(time_point_regions) > 1) {
            results[[time_point]] = intersect_ranges(time_point_regions, min_replicate=2,
                                                     diagnostic_dir=diagnostic_dir, label=paste0(target, "-", time_point))
        } else {
            warning("Not enough replicate for ", target, "-", time_point, ". Using all regions.")
            results[[time_point]] = time_point_regions[[1]]
        }
    }
    
    return(GenomicRanges::GRangesList(results))
}
