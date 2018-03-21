library(metagene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(ggplot2)

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
    sh_subset = subset_metagene_df(region_name)
    
    color_palette = c(shCTRL.1="#40C4FF", shCTRL.2="#00B0FF", shNIPBL.3="#FF8A80", shNIPBL.5="#FF5252")
    ggplot(sh_subset, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group=sh_name, color=sh_name, fill=sh_name)) +
        geom_line() +
        geom_ribbon(alpha = 0.1) +
        facet_grid(Condition~Antibody) +
        scale_color_manual(values=color_palette) +
        scale_fill_manual(values=color_palette)
    file_name = paste0("Metagene sh effect ", region_name, ".pdf")
    ggsave(file.path(output_dir, file_name))
}

# Function for plotting the dex effect on a group of region.
plot_dex_effect <- function(region_name) {
    sh_subset = subset_metagene_df(region_name)
    
    color_palette = c(EtOH="#40C4FF", Dex="#FF8A80")
    ggplot(sh_subset, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group=Condition, color=Condition, fill=Condition)) +
        geom_line() +
        geom_ribbon(alpha = 0.1) +
        facet_grid(sh_name~Antibody) +
        scale_color_manual(values=color_palette) +
        scale_fill_manual(values=color_palette)
    file_name = paste0("Metagene dex effect ", region_name, ".pdf")
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

region_list = list(AllGeneBodies=all_genes, AllTSS=all_TSS, BoundGeneBodies=bound_genes, BoundTSS=bound_TSS, UnboundGeneBodies=unbound_genes, UnboundTSS=unbound_TSS)
         

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


# Plot bound/unbound differences.
for(region_type in c("GeneBodies", "TSS")) {
    bound_subset = subset_metagene_df(paste0("Bound", region_type))
    unbound_subset = subset_metagene_df(paste0("Unbound", region_type))
    bound_subset$Bound = "GR-Bound"
    unbound_subset$Bound = "GR-Unbound"
    
    sh_subset = rbind(bound_subset, unbound_subset)
    
    # Plot sh effect
    color_palette = c(shCTRL.1="#40C4FF", shCTRL.2="#00B0FF", shNIPBL.3="#FF8A80", shNIPBL.5="#FF5252")
    ggplot(sh_subset, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group=sh_name, color=sh_name, fill=sh_name, linetype=Bound)) +
        geom_line() +
        geom_ribbon(alpha = 0.1) +
        facet_grid(Condition~Antibody) +
        scale_color_manual(values=color_palette) +
        scale_fill_manual(values=color_palette)
    file_name = paste0("Metagene sh effect GR-bound vs GR-unbound ", region_type, ".pdf")
    ggsave(file.path(output_dir, file_name))
    
    # Plot dex effect
    color_palette = c(EtOH="#40C4FF", Dex="#FF8A80")
    ggplot(sh_subset, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group=Condition, color=Condition, fill=Condition, linetype=Bound)) +
        geom_line() +
        geom_ribbon(alpha = 0.1) +
        facet_grid(sh_name~Antibody) +
        scale_color_manual(values=color_palette) +
        scale_fill_manual(values=color_palette)
    file_name = paste0("Metagene dex effect GR-bound vs GR-unbound ", region_name, ".pdf")
    ggsave(file.path(output_dir, file_name))
}
