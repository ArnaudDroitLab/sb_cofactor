library(metagene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(ggplot2)
library(org.Hs.eg.db)

###############################################################################
# Define functions for alter analysis.
###############################################################################

# Function for subsetting the metagene data-frame and associating
# sample conditions to profiles.
subset_metagene_df <- function(region_name) {
    metagene.obj = metagenes[[region_name]]
    
    sh_subset = subset(metagene.obj$get_data_frame(), grepl("sh(C..L|NIPBL).._", group))
    sh_subset$group = gsub("_regions", "", sh_subset$group)
    sh_subset$group = gsub("CRTL", "CTRL", sh_subset$group)
    sh_subset$Condition = ifelse(grepl("Dex", sh_subset$group), "Dex", "EtOH")
    sh_subset$Antibody = ifelse(grepl("POL2.ser2", sh_subset$group), "PolII-ser2", "PolII")
    sh_subset$sh_name = gsub("POL.*(sh.*)_.*", "\\1", sh_subset$group)

    return(sh_subset)
}

# Function for plotting the sh effect on a group of regions.
plot_sh_effect <- function(region_name) {
    color_palette = c(shCTRL.1="#40C4FF", shCTRL.2="#00B0FF", shNIPBL.3="#FF8A80", shNIPBL.5="#FF5252")
    do_meta_ggplot_single(region_name, "sh_name", "Condition", color_palette, "sh")
}

# Function for plotting the dex effect on a group of region.
plot_dex_effect <- function(region_name) {
    do_meta_ggplot_single(region_name, "Condition", "sh_name", c(EtOH="#40C4FF", Dex="#FF8A80"), "Dex")
}

do_meta_ggplot_single <- function(region_name, group_var, facet_var, color_palette, file_label) {
    sh_subset = subset_metagene_df(region_name)
    do_meta_ggplot_generic(sh_subset, region_name, group_var, facet_var, color_palette, file_label, FacetVar~Antibody)
}

do_meta_ggplot_double <- function(region_name1, region_name2, label1, label2, facet_var, color_palette, file_label, facet_formula, group_label) {
    subset_1 = subset_metagene_df(region_name1)
    subset_2 = subset_metagene_df(region_name2)
    subset_1$RegionName = label1
    subset_2$RegionName = label2
    
    sh_subset = rbind(subset_1, subset_2)
    do_meta_ggplot_generic(sh_subset, group_label, "RegionName", facet_var, color_palette, file_label, facet_formula)
}

do_meta_ggplot_generic <- function(sh_subset, region_label, group_var, facet_var, color_palette, file_label, facet_formula) {
    sh_subset$GroupVar = sh_subset[[group_var]]
    sh_subset$FacetVar = sh_subset[[facet_var]]
    
    if(grepl("TSS", region_label)) {
        sh_subset$bin = (sh_subset$bin - 50) * 4
        xlabel = "Distance from TSS (bp)"
    } else {
        xlabel = "Relative distance from TSS (0) to TES (100)"
    }
    
    title_str = paste0("Effect of ", file_label, " on ", region_label)
    ggplot(sh_subset, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group=GroupVar, color=GroupVar, fill=GroupVar)) +
        geom_line() +
        geom_ribbon(alpha = 0.1) +
        facet_grid(facet_formula) +
        scale_color_manual(name=group_var, values=color_palette) +
        scale_fill_manual(name=group_var, values=color_palette) +
        ylab("Mean coverage (RPM)") +
        xlab(xlabel) +
        ggtitle(title_str) +
        theme(plot.title = element_text(hjust = 0.5))
    
    file_name = paste0("Metagene - ", title_str, ".pdf")
    ggsave(file.path(output_dir, file_name))
}

###############################################################################
# Define input/output files and directories.
###############################################################################

# Define input/output directories.
mugqic_dir = "output/chip-pipeline-PolII-GRCh38"
output_dir = "output/analyses/PolII/metagenes"

# Read the design file, which contains the names of the BAM to import,
design = read.table("raw/PolII/design.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Remove narrow columns, suffixes.
design = design[,!grepl("N.N", colnames(design))]
colnames(design) = gsub("\\.B$", "", colnames(design))

rownames(design) = design$Sample
design$Sample = file.path(mugqic_dir, "alignment", design$Sample, design$Sample)
design$Sample = paste0(design$Sample, ".sorted.dup.bam")

# MUGQIC pipeline uses 1 for control, 2 for treatment. Metagene uses the opposite.
design[design==1] = 3
design[design==2] = 1
design[design==3] = 2



###############################################################################
# Define regions over which metagenes will be plotted.
###############################################################################

all_genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Remove all regions not on the main chromosomes.
all_genes = all_genes[!grepl("_", seqnames(all_genes))]

# Remove all genes that are smaller than our defined TSS regions.
all_genes = all_genes[width(all_genes) >= 200]

# Define TSS regions based on all kept gene locations.
all_TSS = GenomicRanges::promoters(all_genes, upstream=200, downstream=200)

# Import GR binding regions.
extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
gr_regions = rtracklayer::import("output/ENCODE-chip/peak_call/GR_1hr/GR_1hr_peaks.narrowPeak", format="bed", extraCols=extraCols)

# Determine which genes overlap GR regions.
annotated_gr = ChIPseeker::annotatePeak(gr_regions, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
bound_gene_ids = subset(as.data.frame(annotated_gr), annotation=="Promoter (<=1kb)")$geneId

# Determine which genes are bound by GR.
bound_genes = all_genes[all_genes$gene_id %in% bound_gene_ids]
bound_TSS = GenomicRanges::promoters(bound_genes, upstream=200, downstream=200)
unbound_genes = all_genes[!(all_genes$gene_id %in% bound_gene_ids)]
unbound_TSS = GenomicRanges::promoters(unbound_genes, upstream=200, downstream=200)

# Determine differentially expressed genes at 1 hour.
de_results = read.csv("results/a549_dex_time_points/1h", header=TRUE)
de_results$ENTREZID = mapIds(org.Hs.eg.db, keys=as.character(de_results$gene_id), keytype="ENSEMBL", column="ENTREZID")
de_gene_ids = subset(de_results, abs(log2FoldChange) >= log2(1.5) & padj <= 0.05)$ENTREZID
up_gene_ids = subset(de_results, log2FoldChange <= -log2(1.5) & padj <= 0.05)$ENTREZID
down_gene_ids = subset(de_results, log2FoldChange >= log2(1.5) & padj <= 0.05)$ENTREZID

# Define group of regions based on DE status.
de_genes = all_genes[all_genes$gene_id %in% de_gene_ids]
de_TSS = GenomicRanges::promoters(de_genes, upstream=200, downstream=200)
up_genes = all_genes[all_genes$gene_id %in% up_gene_ids]
up_TSS = GenomicRanges::promoters(up_genes, upstream=200, downstream=200)
down_genes = all_genes[all_genes$gene_id %in% down_gene_ids]
down_TSS = GenomicRanges::promoters(down_genes, upstream=200, downstream=200)

region_list = list(AllGeneBodies=all_genes, AllTSS=all_TSS, 
                   BoundGeneBodies=bound_genes, BoundTSS=bound_TSS, 
                   UnboundGeneBodies=unbound_genes, UnboundTSS=unbound_TSS,
                   DEGeneBodies=de_genes, DETSS=de_TSS,
                   UpRegulatedGeneBodies=up_genes, UpRegulatedTSS=up_TSS,
                   DownRegulatedGeneBodies=down_genes, DownRegulatedTSS=down_TSS)
         

###############################################################################
# Generate the metagene objects and plots for all regions.
###############################################################################
         
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
metagenes = list()
for(region_name in names(region_list)) {
    loaded.cache.filename = file.path(output_dir, paste0(region_name, " Loaded.RData"))
    if(!file.exists(loaded.cache.filename)) {
        metagene.obj = metagene$new(regions=region_list[[region_name]], bam_files=design$Sample, assay='chipseq', force_seqlevels=TRUE)
        save(metagene.obj, file=loaded.cache.filename)
    } else {
        load(loaded.cache.filename)
    }
    
    matrix.cache.filename = file.path(output_dir, paste0(region_name, " matrix.RData"))
    if(!file.exists(matrix.cache.filename)) {
        bin_count = min(c(unlist(lapply(region_list[[region_name]], width)), 100))
        metagene.obj$produce_table(design = design, normalization="RPM", flip_regions=TRUE, bin_count=bin_count)
        save(metagene.obj, file=matrix.cache.filename)
    } else {
        load(matrix.cache.filename)
    }
    
    df.cache.filename = file.path(output_dir, paste0(region_name, " data-frame.RData"))
    if(!file.exists(df.cache.filename)) {
        metagene.obj$produce_data_frame()
        save(metagene.obj, file=df.cache.filename)
    } else {
        load(df.cache.filename)
    }
    
    metagenes[[region_name]] = metagene.obj
    plot_sh_effect(region_name)
    plot_dex_effect(region_name)
}

# Make combined plots comparing Upregulated to Downregulated genes.
for(region_type in c("GeneBodies", "TSS")) {
    color_palette = c(Downregulated="#40C4FF", Upregulated="#FF8A80")
    do_meta_ggplot_double(region_name1=paste0("UpRegulated", region_type), region_name2=paste0("DownRegulated", region_type),
                          label1="Upregulated", label2="Downregulated", 
                          facet_var="sh_name", color_palette=color_palette, 
                          file_label="DE direction", facet_formula=sh_name~Antibody*Condition,
                          group_label=region_type)
}

# Make combined plots comparing GR-bound genes to those not bound by GR.
for(region_type in c("GeneBodies", "TSS")) {
    color_palette = c("GR-Unbound"="#40C4FF", "GR-Bound"="#FF8A80")
    do_meta_ggplot_double(region_name1=paste0("Bound", region_type), region_name2=paste0("Unbound", region_type),
                          label1="GR-Bound", label2="GR-Unbound", 
                          facet_var="sh_name", color_palette=color_palette, 
                          file_label="GR binding", facet_formula=sh_name~Antibody*Condition,
                          group_label=region_type)
}
