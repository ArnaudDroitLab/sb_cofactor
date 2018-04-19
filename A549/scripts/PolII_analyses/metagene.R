library(metagene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(ggplot2)
library(org.Hs.eg.db)
library(dplyr)

source("scripts/load_reddy.R")
###############################################################################
# Define functions for alter analysis.
###############################################################################

# Function for subsetting the metagene data-frame and associating
# sample conditions to profiles.
subset_metagene_df <- function(region_name, antibody=NA) {
    metagene.obj = metagenes[[region_name]]
    
    sh_subset = subset(metagene.obj$get_data_frame(), grepl("sh(C..L|NIPBL).._", group))
    sh_subset$group = gsub("_regions", "", sh_subset$group)
    sh_subset$group = gsub("CRTL", "CTRL", sh_subset$group)
    sh_subset$Condition = ifelse(grepl("Dex", sh_subset$group), "Dex", "EtOH")
    sh_subset$Antibody = ifelse(grepl("POL2.ser2", sh_subset$group), "PolII-ser2", "PolII")
    sh_subset$sh_name = gsub("POL.*(sh.*)_.*", "\\1", sh_subset$group)
    sh_subset$Direction = ifelse(grepl("Regulated", region_name), gsub("(.*Regulated).*", "\\1", region_name), NA)
    sh_subset$Bound = ifelse(grepl("ound", region_name), gsub(".*((Unb|B)ound).*", "\\1", region_name), NA)
    
    if(!is.na(antibody)) {
        sh_subset = sh_subset[sh_subset$Antibody==antibody,]
    }
    
    return(sh_subset)
}

determine_bp_per_bin <- function(region_name, bin_number=100) {
    region_widths = width(metagenes[[region_name]]$get_regions[[1]])
    if(all(region_widths == mean(region_widths))) {
        # All regions are the same length.
        return(region_widths[1] / bin_number)
    } else {
        return(NA)
    }
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

do_meta_ggplot_single <- function(region_name, group_var, facet_var, color_palette, file_label, antibody=NA) {
    sh_subset = subset_metagene_df(region_name, antibody)
    bp_per_bin = determine_bp_per_bin(region_name) 
    do_meta_ggplot_generic(sh_subset, region_name, group_var, facet_var, color_palette, file_label, FacetVar~Antibody, bp_per_bin=bp_per_bin)
}

do_meta_ggplot_double <- function(region_name1, region_name2, label1, label2, facet_var, color_palette, file_label, facet_formula, group_label, antibody=NA) {
    subset_1 = subset_metagene_df(region_name1, antibody)
    subset_2 = subset_metagene_df(region_name2, antibody)
    subset_1$RegionName = label1
    subset_2$RegionName = label2
    
    sh_subset = rbind(subset_1, subset_2)
    bp_per_bin = determine_bp_per_bin(region_name1)    
    do_meta_ggplot_generic(sh_subset, group_label, "RegionName", facet_var, color_palette, file_label, facet_formula, bp_per_bin=bp_per_bin)
}

do_meta_ggplot_any <- function(name_list, facet_var, color_palette, file_label, facet_formula, group_label, antibody="both") {
    all_subset = data.frame()
    for(name_index in 1:length(name_list)) {
        this_subset = subset_metagene_df(name_list[[name_index]], antibody)
        this_subset$RegionName = names(name_list)[name_index]
        all_subset = rbind(all_subset, this_subset)
    }

    bp_per_bin = determine_bp_per_bin(name_list[[1]])
    do_meta_ggplot_generic(all_subset, group_label, "RegionName", facet_var, color_palette, file_label, facet_formula, bp_per_bin=bp_per_bin)
}


do_meta_ggplot_generic <- function(sh_subset, region_label, group_var, facet_var, color_palette, file_label, facet_formula, bp_per_bin=4) {
    sh_subset$GroupVar = sh_subset[[group_var]]
    sh_subset$FacetVar = sh_subset[[facet_var]]
    
    if(grepl("TSS", region_label)) {
        sh_subset$bin = (sh_subset$bin - 50) * bp_per_bin
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

most_expressed = load_most_expressed_transcripts()
most_expressed_TxDb = load_most_expressed_TxDb()
seqlevelsStyle(most_expressed_TxDb) <- "UCSC"

all_genes = most_expressed

# Remove all regions not on the main chromosomes.
all_genes = all_genes[!grepl("_", seqnames(all_genes))]

# Remove all genes that are smaller than our defined TSS regions.
all_genes = all_genes[width(all_genes) >= 200]

# Define TSS regions based on all kept gene locations.
all_TSS = GenomicRanges::promoters(all_genes, upstream=1000, downstream=1000)

# Import GR binding regions.
gr_regions = load_reddy_gr_binding_consensus()

# Determine which genes overlap GR regions.
annotated_gr = ChIPseeker::annotatePeak(gr_regions[["30 minutes"]], TxDb=most_expressed_TxDb)
annotated_gr_df = as.data.frame(annotated_gr)
annotated_gr_df$ENTREZID = mapIds(org.Hs.eg.db, keys=as.character(annotated_gr_df$geneId), keytype="ENSEMBL", column="ENTREZID")
bound_gene_ids = subset(annotated_gr_df, distanceToTSS <= 3000)$ENTREZID
unbound_gene_ids = setdiff(all_genes$entrezgene, bound_gene_ids)

# Determine differentially expressed genes at 1 hour.
de_results = load_reddy_de_list()[["2h"]]$Full
de_results$ENTREZID = mapIds(org.Hs.eg.db, keys=as.character(de_results$gene_id), keytype="ENSEMBL", column="ENTREZID")
de_gene_ids = subset(de_results, abs(log2FoldChange) >= log2(1.5) & padj <= 0.05)$ENTREZID
up_gene_ids = subset(de_results, log2FoldChange <= -log2(1.5) & padj <= 0.05)$ENTREZID
down_gene_ids = subset(de_results, log2FoldChange >= log2(1.5) & padj <= 0.05)$ENTREZID

# Define group of regions based on DE status AND GR-binding status.
de_bound_gene_ids = intersect(de_gene_ids, bound_gene_ids)
up_bound_gene_ids = intersect(up_gene_ids, bound_gene_ids)
down_bound_gene_ids = intersect(down_gene_ids, bound_gene_ids)
de_unbound_gene_ids = intersect(de_gene_ids, unbound_gene_ids)
up_unbound_gene_ids = intersect(up_gene_ids, unbound_gene_ids)
down_unbound_gene_ids = intersect(down_gene_ids, unbound_gene_ids)

# Import regions based on expression level
raw_counts = read.csv("results/a549_dex_time_points/raw_counts.csv")
filter_df = data.frame(gene_id=raw_counts$gene_id, mean=rowMeans(raw_counts[,-1]))
filter_df = filter_df[filter_df$mean >= 10,]
filter_df$ENTREZID = mapIds(org.Hs.eg.db, keys=as.character(filter_df$gene_id), keytype="ENSEMBL", column="ENTREZID")
filter_df$quartile = cut(filter_df$mean, breaks=quantile(filter_df$mean, probs=seq(0,1, 1/15)), labels=paste0("Q", 1:15))


q1_gene_ids = filter_df$ENTREZID[filter_df$quartile=="Q1"]
q8_gene_ids = filter_df$ENTREZID[filter_df$quartile=="Q8"]
q15_gene_ids = filter_df$ENTREZID[filter_df$quartile=="Q15"]


region_list = list(#AllGeneBodies=all_genes,
                   #AllTSS=all_TSS, 
                   #BoundGeneBodies=get_gene_bodies(all_genes, bound_gene_ids),
                   #BoundTSS=get_tss(all_genes, bound_gene_ids), 
                   #UnboundGeneBodies=get_gene_bodies(all_genes, unbound_gene_ids),
                   #UnboundTSS=get_tss(all_genes, unbound_gene_ids),
                   #DEGeneBodies=get_gene_bodies(all_genes, de_gene_ids),
                   #DETSS=get_tss(all_genes, de_gene_ids),
                   #UpRegulatedGeneBodies=get_gene_bodies(all_genes, up_gene_ids),
                   #UpRegulatedTSS=get_tss(all_genes, up_gene_ids),
                   #DownRegulatedGeneBodies=get_gene_bodies(all_genes, down_gene_ids),
                   #DownRegulatedTSS=get_tss(all_genes, down_gene_ids),
                   UpRegulatedBoundGeneBodies=get_gene_bodies(all_genes, up_gene_ids, bound_gene_ids),
                   UpRegulatedBoundTSS=get_tss(all_genes, up_gene_ids, bound_gene_ids, flank_size=1000),
                   UpRegulatedUnboundGeneBodies=get_gene_bodies(all_genes, up_gene_ids, unbound_gene_ids),
                   UpRegulatedUnboundTSS=get_tss(all_genes, up_gene_ids, unbound_gene_ids, flank_size=1000),                   
                   DownRegulatedBoundGeneBodies=get_gene_bodies(all_genes, down_gene_ids, bound_gene_ids),
                   DownRegulatedBoundTSS=get_tss(all_genes, down_gene_ids, bound_gene_ids, flank_size=1000),
                   DownRegulatedUnboundGeneBodies=get_gene_bodies(all_genes, down_gene_ids, unbound_gene_ids),
                   DownRegulatedUnboundTSS=get_tss(all_genes, down_gene_ids, unbound_gene_ids, flank_size=1000),
                   Q1TSS=get_tss(all_genes, q1_gene_ids, flank_size=1000),
                   Q8TSS=get_tss(all_genes, q8_gene_ids, flank_size=1000),
                   Q15TSS=get_tss(all_genes, q15_gene_ids, flank_size=1000),
                   ACTBTSS=get_tss(all_genes, 60, flank_size=1000))
         

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
#for(region_type in c("GeneBodies", "TSS")) {
#    color_palette = c(Downregulated="#40C4FF", Upregulated="#FF8A80")
#    do_meta_ggplot_double(region_name1=paste0("UpRegulated", region_type), region_name2=paste0("DownRegulated", region_type),
#                          label1="Upregulated", label2="Downregulated", 
#                          facet_var="sh_name", color_palette=color_palette, 
#                          file_label="DE direction", facet_formula=sh_name~Antibody*Condition,
#                          group_label=region_type)
#}
#
## Make combined plots comparing GR-bound genes to those not bound by GR.
#for(region_type in c("GeneBodies", "TSS")) {
#    color_palette = c("GR-Unbound"="#40C4FF", "GR-Bound"="#FF8A80")
#    do_meta_ggplot_double(region_name1=paste0("Bound", region_type), region_name2=paste0("Unbound", region_type),
#                          label1="GR-Bound", label2="GR-Unbound", 
#                          facet_var="sh_name", color_palette=color_palette, 
#                          file_label="GR binding", facet_formula=sh_name~Antibody*Condition,
#                          group_label=region_type)
#}

# Compare GR-Bound*Direction interaction.
for(antibody in c("PolII", "PolII-ser2")) {
    for(region_type in c("TSS", "GeneBodies")) {
        region_list = list("GR-Bound"   = paste0("UpRegulatedBound", region_type),
                           "GR-Unbound" = paste0("UpRegulatedUnbound", region_type),
                           "GR-Bound"   = paste0("DownRegulatedBound", region_type),
                           "GR-Unbound" = paste0("DownRegulatedUnbound", region_type))
        file_label = paste0("on ", antibody, " of GR-Binding")
        do_meta_ggplot_any(region_list,
                           facet_var="Direction", color_palette=c("GR-Unbound"="#40C4FF", "GR-Bound"="#FF8A80"),
                           file_label=file_label, facet_formula=sh_name~Direction*Condition,
                           group_label=region_type, antibody=antibody)
    }
}

# Plot Q1, Q8 Q15 against each other.
q_list = list("Q1"="Q1TSS", "Q8"="Q8TSS", "Q15"="Q15TSS")
q_palette = c(Q1="#4FC3F7", Q8="#039BE5", Q15="#01579B")
do_meta_ggplot_any(q_list, facet_var="sh_name", color_palette=q_palette,
                   file_label="gene expression level", facet_formula=sh_name~Condition,
                   group_label="Quantile (Out of 15)", antibody="PolII")

do_meta_ggplot_any(list("ACTB"="ACTBTSS"), facet_var="sh_name", color_palette=c(ATCB="#039BE5"),
                   file_label="ACTB", facet_formula=sh_name~Condition,
                   group_label="ACTB", antibody="PolII")        
                   