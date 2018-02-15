library(metagene)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)

run.options <- commandArgs(trailingOnly=TRUE)
if(length(run.options) < 2) {
    # Not enough arguments, use default
    nrs = c("GR", "AR", "LCC9")
    regions = c("Promoter", "First", "Other", "Intergenic")
} else {
    nrs = run.options[1]
    regions = run.options[2]
}

annoSplit <- function(csAnno, DE_genes) {
    anno.gr = as.GRanges(csAnno)
    
    promoter = grepl("Promoter", anno.gr$annotation)
    promoter_DE = grepl("Promoter", anno.gr$annotation) & (anno.gr$SYMBOL %in% names(DE_genes))
    promoter_DE_up = grepl("Promoter", anno.gr$annotation) & (anno.gr$SYMBOL %in% names(DE_genes)) & !is.na(DE_genes[anno.gr$SYMBOL]) & (DE_genes[anno.gr$SYMBOL] > 0)
    promoter_DE_down = grepl("Promoter", anno.gr$annotation) & (anno.gr$SYMBOL %in% names(DE_genes)) & !is.na(DE_genes[anno.gr$SYMBOL]) & (DE_genes[anno.gr$SYMBOL] < 0)
    promoter_NotDE = grepl("Promoter", anno.gr$annotation) & !(anno.gr$SYMBOL %in% names(DE_genes))
    first.exon.intron = grepl("exon 1 of", anno.gr$annotation) | grepl("intron 1 of", anno.gr$annotation)
    intergenic = grepl("intergenic", anno.gr$annotation) | grepl("Downstream", anno.gr$annotation)
    other.exon = !promoter & !first.exon.intron & !intergenic
    
    return(GRangesList(list(Promoter=anno.gr[promoter],
                            Promoter_DE=anno.gr[promoter_DE],
                            Promoter_DE_Up = anno.gr[promoter_DE_up],
                            Promoter_DE_Down = anno.gr[promoter_DE_down],
                            Promoter_NotDE=anno.gr[promoter_NotDE],
                            First=anno.gr[first.exon.intron],
                            Other=anno.gr[intergenic],
                            Intergenic=anno.gr[other.exon])))
}

# Split the design in two: MCF7
perform_metagene <- function (design, design_regex, regions, out.path, out.label) {
    design_subset = design[grepl(design_regex, design$sample),]
    design_subset = design_subset[,!apply(design_subset, 2, function(x) { all(x==0) })]

    cache.path = file.path(out.path, "cache")
    dir.create(cache.path, recursive=TRUE, showWarnings=FALSE)
    loaded.cache.filename = file.path(cache.path, paste0(out.label, " Loaded.RData"))
    if(!file.exists(loaded.cache.filename)) {
        metagene.obj = metagene$new(regions=regions, bam_files=design_subset$sample)
        save(metagene.obj, file=loaded.cache.filename)
    } else {
        load(loaded.cache.filename)
    }
    
    matrix.cache.filename = file.path(cache.path, paste0(out.label, " matrix.RData"))
    if(!file.exists(matrix.cache.filename)) {
        if(class(GRangesList(regions, regions))=="GRangesList") {
            bin_count = min(c(unlist(lapply(regions, width)), 200))
        } else {
            bin_count = min(c(width(regions), 200))
        }
        metagene.obj$produce_matrices(design = design_subset, normalization="RPM", flip_regions=FALSE, bin_count=bin_count)
        save(metagene.obj, file=matrix.cache.filename)
    } else {
        load(matrix.cache.filename)
    }
    
    df.cache.filename = file.path(cache.path, paste0(out.label, " data-frame.RData"))
    if(!file.exists(df.cache.filename)) {
        metagene.obj$produce_data_frame(stat = "bootstrap")
        save(metagene.obj, file=df.cache.filename)
    } else {
        load(df.cache.filename)
    }
    
    return(metagene.obj)
}

plot_metagene <- function(meta_obj, label, outdir) {
    meta_df = meta_obj$get_data_frame()
    meta_df$Cell = gsub("(.*)_(.*)_(.*)_regions", "\\1", meta_df$group)
    meta_df$Condition = gsub("(.*)_(.*)_(.*)_regions", "\\2", meta_df$group)
    meta_df$Target = gsub("(.*)_(.*)_(.*)_regions", "\\3", meta_df$group)
    
    ggplot(meta_df, aes(x=position, y=value, ymin=qinf, ymax=qsup, color=Condition, fill=Condition)) +
        geom_line(mapping=aes(group=Condition)) +
        geom_ribbon(alpha=0.5) +
        facet_grid(Target~., scales="free")
    ggsave(file.path(outdir, paste0(label, ".pdf")))
}

# Load GR- and ER-bound regions in treatment condition.
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
nr.regions = list(GR=import("A549/output/ENCODE-chip/peak_call/GR_1hr/GR_1hr_peaks.narrowPeak", format = "BED", extraCols = extraCols_narrowPeak),
                  ER=import("MCF7/output/chip-pipeline-GRCh38/peak_call/MCF7_E2_ERA/MCF7_E2_ERA_peaks.narrowPeak", format = "BED", extraCols = extraCols_narrowPeak),
                  LCC9=import("LCC9/output/chip-pipeline-GRCh38/peak_call/LCC9_E2_ERA/LCC9_E2_ERA_peaks.narrowPeak", format = "BED", extraCols = extraCols_narrowPeak))
nr.annotated = lapply(nr.regions, annotatePeak, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")

#### Fix LCC9 because data is actually hg19 for now
nr.annotated[["LCC9"]] = annotatePeak(nr.regions[["LCC9"]], TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
#### End fix/hack

# Read DE genes
a549_table = read.table("input/A549_DE.txt", sep=",", header=TRUE, stringsAsFactors=FALSE)
a549_table = a549_table[a549_table$padj <= 0.05 & abs(a549_table$log2FoldChange) >= log2(1.5) & !is.na(a549_table$padj),]
a549_DE = a549_table$log2FoldChange
names(a549_DE) = a549_table$symbol

mcf7_table = read.table("input/MCF7_E2_32h.txt", header=TRUE, stringsAsFactors=FALSE)
mcf7_DE = mcf7_table$Fold_change_32h
names(mcf7_DE) = mcf7_table$Gene

nr.DE = list(GR=a549_DE, ER=mcf7_DE, LCC9=mcf7_DE)

nr.split = list()
for(nr in names(nr.annotated)) {
    nr.split[[nr]] = annoSplit(nr.annotated[[nr]], nr.DE[[nr]])
}
                          
# Read the design file, which contains the names of the BAM to import,
design = read.table("input/design.txt", stringsAsFactors=FALSE, header=TRUE, sep="\t")

nr.regex = list(GR="A549", ER="MCF7", LCC9="LCC9")

metagenes = list()
for(nr in nrs) {
    for(region in regions) {
        label = paste(nr, region, sep="-")
        metagenes[[label]] = perform_metagene(design, nr.regex[[nr]], nr.split[[nr]][[region]], "output/metagenes", label)
        plot_metagene(metagenes[[label]], label, "output/metagenes")
    }
}
