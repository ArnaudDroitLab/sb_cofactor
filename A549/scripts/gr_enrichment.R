library(ENCODExplorer)
library(ef.utils)
library(GenomicRanges)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(org.Hs.eg.db)
library(reshape2)
library(dplyr)
library(gplots)

# Turn off VennDiagram logging.
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

output_dir="output/analyses/gr_enrichment"
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

###############################################################################
# Load and annotate GR binding data.
###############################################################################

# Get GR binding time series
gr_accession = read.table("input/ENCODE_Reddy_GR_ChIP_experiments.txt", header=TRUE, sep="\t")
gr_accession = gr_accession[order(gr_accession$Order),]
gr_accession$Time = factor(gr_accession$Time, levels=gr_accession$Time)

chip_dir = "input/ENCODE/A549/GRCh38/chip-seq/narrow"
dir.create(chip_dir, recursive=TRUE, showWarnings=FALSE)
gr_regions = list()
for(i in 1:nrow(gr_accession)) {
    encode_accession = gr_accession$Experiment[i]
    time_point = as.character(gr_accession$Time[i])
    
    # Download and import peak calls for the time point.
    encodeResults = queryEncodeGeneric(accession=encode_accession, file_type="bed narrowPeak", assembly="GRCh38")
    downloaded_files = downloadEncode(encodeResults, dir=chip_dir, force=FALSE)
    names(downloaded_files) = gsub(".*\\/(.*).bed.gz", "\\1", downloaded_files)
    
    extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
    binding_sites = GRangesList(lapply(downloaded_files, rtracklayer::import, extraCols=extraCols))
    
    # Only keep peak calls found in at least two replicates.
    intersect_object = build_intersect(binding_sites)
    two_replicates = rowSums(intersect_object$Matrix) >= 2
    gr_regions[[time_point]] = intersect_object$Regions[two_replicates]
    dir.create(file.path(output_dir, "Venn diagrams of peak calls"))
    intersect_venn_plot(intersect_object, file.path(output_dir, "Venn diagrams of peak calls", paste0(time_point, ".tiff")))

    # Get rid of incomplete/alternate scaffolds, since they interfere with annotation later.
    gr_regions[[time_point]] = gr_regions[[time_point]][!grepl("_", seqnames(gr_regions[[time_point]]))]
}
    

# Are binding regions conserved across time points?
intersect_all = build_intersect(GRangesList(gr_regions))
all_regions_annotation = ChIPseeker::annotatePeak(intersect_all$Regions, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
all_gr_regions = GRanges(as.data.frame(all_regions_annotation))

###############################################################################
# Plot some basic metrics of the GR binding data.
###############################################################################

pdf(file.path(output_dir, "Heatmap of GR-binding sites.pdf"))
heatmap.2(intersect_all$Matrix, dendrogram="column", scale="none", trace="none")
dev.off()    
    
# Do the type of GR binding change over time?
gr_binding_type_df = data.frame()
for(i in names(intersect_all$List)) {
    which_regions = all_gr_regions[intersect_all$List[[i]]]
    row_df = data.frame(Time=i,
                        Annotation=gsub(" \\(.*\\)", "", which_regions$annotation))
    gr_binding_type_df = rbind(gr_binding_type_df, row_df)                    
                        
}
gr_binding_type_df$Time = factor(gr_binding_type_df$Time, levels=gr_accession$Time)

# Plot GR binding type in absolute number of regions.
ggplot(gr_binding_type_df, aes(x=Time, fill=Annotation)) +
    geom_bar(color="black") +
    labs(x="Time", y="Number of identified GR-binding sites", title="GR-binding over time in Reddy dataset") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1), plot.title = element_text(hjust=0.5))
ggsave(file.path(output_dir, "GR-bound regions over time in Reddy dataset (absolute).pdf"))
    
# Plot GR binding type as proportion of current GR-bindings.
gr_proportions = gr_binding_type_df %>% 
                    group_by(Time, Annotation) %>% 
                    summarise(n=n()) %>% 
                    mutate(freq = n / sum(n))
ggplot(gr_proportions, aes(x=Time, fill=Annotation, y=freq*100)) +
    geom_bar(color="black", stat="identity") +
    labs(x="Time", y="Percentage of identified GR-binding sites", title="GR-binding over time in Reddy dataset") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1), plot.title = element_text(hjust=0.5))
ggsave(file.path(output_dir, "GR-bound regions over time in Reddy dataset (proportions).pdf"))

# Plot the number of genes bound by GR
gr_bound_genes = data.frame(Gene=unique(all_gr_regions[all_gr_regions$annotation=="Promoter (<=1kb)"]$geneId))
for(i in names(intersect_all$List)) {
    which_regions = all_gr_regions[intersect_all$List[[i]]]
    gr_bound_genes[[i]] = gr_bound_genes$Gene %in% which_regions[which_regions$annotation=="Promoter (<=1kb)"]$geneId
}

# Plot a heatmap of GR-bound genes.
pdf(file.path(output_dir, "Heatmap of GR-bound genes.pdf"))
heatmap.2(as.matrix(gr_bound_genes[,-1])+0, dendrogram="column", scale="none", trace="none")
dev.off() 

# Plot a bar graph of the number of GR-bound genes.
gr_bound_number = data.frame(Time=colnames(gr_bound_genes)[-1], 
                             Genes=apply(gr_bound_genes[-1], 2, sum))
gr_bound_number$Time = factor(gr_bound_number$Time, levels=gr_accession$Time)
ggplot(gr_bound_number, aes(x=Time, y=Genes)) +
    geom_bar(color="black", stat="identity", fill="#4FC3F7") +
    labs(x="Time", y="Number of genes with GR-binding at their TSS", title="GR-binding over time in Reddy dataset") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1), plot.title = element_text(hjust=0.5))    
ggsave(file.path(output_dir, "GR-bound genes over time in Reddy dataset.pdf"))   
    
###############################################################################
# Load DE results
###############################################################################

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

###############################################################################
# Plot simple DE metrics.
###############################################################################

# Plot the number of DE genes per time point.
de_over_time = data.frame(Time=names(de_results))
de_over_time$Time = factor(de_over_time$Time, levels=c("0.5h", "1h", "2h", "3h", "4h", "5h", "6h", "7h", "8h", "10h", "12h"))
de_over_time$Up = unlist(lapply(de_results, function(x) { sum(x$Up_indices, na.rm=TRUE)}))
de_over_time$Down = unlist(lapply(de_results, function(x) { sum(x$Down_indices, na.rm=TRUE)}))
de_over_time = melt(de_over_time, variable.name="Direction", value.name="NumberOfGenes")

ggplot(de_over_time, aes(x=Time, y=NumberOfGenes, fill=Direction)) + 
    geom_bar(stat="identity", color="black") +
    labs(x="Time", y="Number of identified DE genes", title="Differential expression over time in Reddy dataset") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1), plot.title = element_text(hjust=0.5))
ggsave(file.path(output_dir, "DE over time in Reddy dataset.pdf"))

# Make a matrix of fold-changes so we can plot a heatmap.
all_de = unique(unlist(lapply(de_results, function(x) { c(as.character(x[["Up"]]$ENTREZID), as.character(x[["Down"]]$ENTREZID)) })))
all_fc = data.frame(ENTREZID=all_de)
for(de_item in names(de_results)) {
    cur_fc = de_results[[de_item]]$Full
    all_fc[[de_item]] = cur_fc$log2FoldChange[match(all_fc$ENTREZID, cur_fc$ENTREZID)]
}   

# Remove NAs, and cap FC values at +/-5 to keep the color scale more useful.
heatmap_matrix = as.matrix(all_fc[,-1])
heatmap_matrix[is.na(heatmap_matrix)] = 0
heatmap_matrix[heatmap_matrix < -5] = -5
heatmap_matrix[heatmap_matrix > 5] = 5

color_breaks = seq(-5,5,length=40)
color_palette = colorRampPalette(c("red", "white", "green"))(n = 39)

pdf(file.path(output_dir, "DE heatmap.pdf"))
heatmap.2(heatmap_matrix, dendrogram="column", scale="none", trace="none", col=color_palette, breaks=color_breaks)
dev.off()
    
###############################################################################
# Plot a combined GR-binding/DE heatmap.
###############################################################################
    
pdf(file.path(output_dir, "Heatmap of GR-bound genes (2).pdf"))
# Get an overall DE status for all genes.
directions = ifelse(apply(all_fc[,-1], 1, mean, na.rm=TRUE)<0, "Green", "Red")
gr_genes_directions = directions[match(gr_bound_genes$Gene, all_fc$ENTREZID)]
heatmap.2(as.matrix(gr_bound_genes[,-1])+0, dendrogram="column", scale="none", trace="none", RowSideColors=gr_genes_directions)
dev.off()     
###############################################################################
# Perform enrichment of GR binding.
############################################################################    

get_gr_genes_at_time_point <- function(gr_time) {
    gr_ranges = all_gr_regions[intersect_all$List[[as.character(gr_time)]]]
    gr_promoters = gr_ranges[gr_ranges$annotation=="Promoter (<=1kb)"]
    
    # Get the total number of GR-bound genes
    gr_genes = gr_promoters$geneId
    
    return(gr_genes)
}

get_gr_genes_at_or_before_time_point <- function(gr_time) {
    # Find all applicable time points
    before_levels = 1:which(as.character(gr_time)==levels(gr_accession$Time))
    before_times = levels(gr_accession$Time)[before_levels]
    
    before_regions = c()
    for(i in before_times) {
        before_regions = c(before_regions, intersect_all$List[[as.character(i)]])
    }
    
    gr_ranges = all_gr_regions[unique(before_regions)]
    gr_promoters = gr_ranges[gr_ranges$annotation=="Promoter (<=1kb)"]
    
    # Get the total number of GR-bound genes
    gr_genes = gr_promoters$geneId
    
    return(gr_genes)
}

get_gr_genes_any_time_point <- function(gr_time) {
    gr_promoters = all_gr_regions[all_gr_regions$annotation=="Promoter (<=1kb)"]
    
    # Get the total number of GR-bound genes
    gr_genes = gr_promoters$geneId
    
    return(gr_genes)
}
    
perform_gr_enrichment <- function(gr_function, time_match_column, label) {
    time_matches = read.table("input/TimeMatches.txt", sep="\t", header=TRUE)
    enrich_df = data.frame()
    for(de_index in 1:nrow(time_matches)) {
        for(direction in c("Up", "Down")) {
            de_time = time_matches$DETime[de_index]
            gr_time = time_matches[[time_match_column]][de_index]
            
            # Get total number of genes.
            n_total_genes = length(unique(de_results[[de_time]]$Full$ENTREZID))
            
            # Get the number of DEregulated genes.
            de_genes = de_results[[de_time]][[direction]]$ENTREZID
            n_de_genes = length(unique(de_genes))
            
            # Get the list of GR-bound genes.
            gr_genes = gr_function(gr_time)
            n_gr_genes = length(unique(gr_genes))
            
            # Get the total number of GR-bound de-regulated genes.
            n_gr_de_genes = length(unique(intersect(de_genes, gr_genes)))
            
            # Do hyper geometric enrichment.
            enrichment_pval = phyper(n_gr_de_genes, n_gr_genes, n_total_genes - n_gr_genes, n_de_genes, lower.tail=FALSE)
            
            row_df = data.frame(Direction=direction,
                                DETime=de_time,
                                GRTime=gr_time,
                                DEGenes = n_de_genes,
                                GRGenes = n_gr_genes,
                                GRDEGenes = n_gr_de_genes,
                                TotalGenes = n_total_genes, 
                                DEProportion = n_de_genes / n_total_genes,
                                GRProportionAll = n_gr_genes / n_total_genes,
                                GRProportionDE = n_gr_de_genes / n_de_genes,
                                PVal=enrichment_pval)
                                
            enrich_df = rbind(enrich_df, row_df)
        }
    }
        
    enrich_df$Enrichment = log2(enrich_df$GRProportionDE / enrich_df$GRProportionAll)
    write.table(enrich_df, row.names=FALSE, col.names=TRUE, sep="\t", file=file.path(output_dir, paste0("Enrichment ", label, ".txt")))
}

perform_gr_enrichment(get_gr_genes_at_time_point, "HourBefore", "GR-bound 1 hour before")
perform_gr_enrichment(get_gr_genes_any_time_point, "HourBefore", "any GR binding")
perform_gr_enrichment(get_gr_genes_at_or_before_time_point, "SameHour", "GR-bound any time before")
