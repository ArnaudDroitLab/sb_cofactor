library(ENCODExplorer)
library(ef.utils)
library(GenomicRanges)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(org.Hs.eg.db)
library(reshape2)

output_dir="output/analyses/gr_enrichment"
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

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
    downloaded_files = downloadEncode(encodeResults, dir=chip_dir)
    names(downloaded_files) = gsub(".*\\/(.*).bed.gz", "\\1", downloaded_files)
    
    extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
    binding_sites = GRangesList(lapply(downloaded_files, rtracklayer::import, extraCols=extraCols))
    
    # Only keep peak calls found in at least two replicates.
    intersect_object = build_intersect(binding_sites)
    two_replicates = rowSums(intersect_object$Matrix) >= 2
    gr_regions[[time_point]] = intersect_object$Regions[two_replicates]
    intersect_venn_plot(intersect_object, file.path(output_dir, paste0(time_point, ".tiff")))

    # Get rid of incomplete/alternate scaffolds, since they interfere with annotation later.
    gr_regions[[time_point]] = gr_regions[[time_point]][!grepl("_", seqnames(gr_regions[[time_point]]))]
}
    
# Plot the number of binding sites over time.
binding_over_time = data.frame(NumberOfRegions=unlist(lapply(gr_regions, length)),
                               Time=factor(names(gr_regions), levels=gr_accession$Time))
ggplot(binding_over_time, aes(x=Time, y=NumberOfRegions)) + 
    geom_bar(stat="identity", color="black", fill="#4FC3F7") +
    labs(x="Time", y="Number of identified GR-binding sites", title="GR-binding over time in Reddy dataset") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1), plot.title = element_text(hjust=0.5))
ggsave(file.path(output_dir, "GR-binding over time in Reddy dataset.pdf"))

# Are binding regions conserved across time points?
intersect_all = build_intersect(GRangesList(gr_regions))
all_regions_annotation = ChIPseeker::annotatePeak(intersect_all$Regions, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
all_gr_regions = GRanges(as.data.frame(all_regions_annotation))

pdf(file.path(output_dir, "Heatmap of GR-binding sites.pdf"))
heatmap(intersect_all$Matrix, Rowv=NA, Colv=NULL)
dev.off()    
    
# Load DE results
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

# Plot the number of DE genes per time point.
de_over_time = data.frame(Time=names(de_results))
de_over_time$Time = factor(de_over_time$Time, levels=c("0.5h", "1h", "2h", "3h", "4h", "5h", "6h", "7h", "8h", "10h", "12h"))
#de_over_time$DE = unlist(lapply(de_results, function(x) { sum(x$DE_indices, na.rm=TRUE)}))
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
all_de = c()
for(de_item in de_results) {
    all_de = c(all_de, as.character(de_item$DE$gene_id))
}
all_de = unique(all_de)
   
all_fc = data.frame(gene_id=all_de)
for(de_item in names(de_results)) {
    all_fc[[de_item]] = de_results[[de_item]]$Full$log2FoldChange[match(all_fc$gene_id, de_results[[de_item]]$Full$gene_id)]
}   

pdf(file.path(output_dir, "DE heatmap.pdf"))
heatmap(as.matrix(all_fc[,-1]), Rowv=NA, Colv=NULL)
dev.off()
    
# Perform enrichment of GR binding
time_matches = read.table("input/TimeMatches.txt", sep="\t", header=TRUE)
enrich_df = data.frame()
for(de_index in 1:nrow(time_matches)) {
    for(direction in c("Up", "Down")) {
        de_time = time_matches$DETime[de_index]
        gr_time = time_matches$HourBefore[de_index]
        
        # Get total number of genes.
        n_total_genes = length(unique(de_results[[de_time]]$Full$ENTREZID))
        
        # Get the number of DEregulated genes.
        de_genes = de_results[[de_time]][[direction]]$ENTREZID
        n_de_genes = length(unique(de_genes))
        
        # Get the list of GR-bound genes.
        gr_ranges = all_gr_regions[intersect_all$List[[gr_time]]]
        gr_promoters = gr_ranges[gr_ranges$annotation=="Promoter (<=1kb)"]
        
        # Get the total number of GR-bound genes
        gr_genes = gr_promoters$geneId
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
write.table(enrich_df, row.names=FALSE, col.names=TRUE, sep="\t", file=file.path(output_dir, "Enrichment.txt"))