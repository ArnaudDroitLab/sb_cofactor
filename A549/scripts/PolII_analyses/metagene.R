library(metagene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)

output_dir = "output/chip-pipeline-GRCh38"

# Read the design file, which contains the names of the BAM to import,
design = read.table("raw/design.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Remove narrow columns, suffixes.
design = design[,!grepl("N.N", colnames(design))]
colnames(design) = gsub("\\.B$", "", colnames(design))

rownames(design) = design$Sample
design$Sample = file.path(output_dir, "alignment", design$Sample, design$Sample)
design$Sample = paste0(design$Sample, ".sorted.dup.bam")

# MUGQIC pipeline uses 1 for control, 2 for treatment. Metagene uses the opposite.
design[design==1] = 3
design[design==2] = 1
design[design==3] = 2


all_genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Remove all regions not on the main chromosomes.
all_genes = all_genes[!grepl("_", seqnames(all_genes))]

# Remove all genes that are smaller than our defined TSS regions.
all_genes = all_genes[width(all_genes) >= 200]

# Define TSS regions based on all kept gene locations.
all_TSS = GenomicRanges::promoters(all_genes, upstream=200, downstream=200)

# Import GR binding regions.
extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
gr_regions = rtracklayer::import("../CofactorHR/A549/output/ENCODE-chip/peak_call/GR_1hr/GR_1hr_peaks.narrowPeak", format="bed", extraCols=extraCols)

# Determine which genes overlap GR regions.
annotated_gr = ChIPseeker::annotatePeak(gr_regions, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
bound_gene_ids = subset(as.data.frame(annotated_gr), annotation=="Promoter (<=1kb)")$geneId

# Determine which genes are bound by GR.
bound_genes = all_genes[all_genes$gene_id %in% bound_gene_ids]
bound_TSS = GenomicRanges::promoters(bound_genes, upstream=200, downstream=200)

region_list = list(AllGeneBodies=all_genes, AllTSS=all_TSS, BoundGeneBodies=bound_genes, BoundTSS=bound_TSS)
                 
output_path = "output/metagenes"
dir.create(output_path, recursive=TRUE, showWarnings=FALSE)
metagenes = list()
for(region_name in names(region_list)) {
    loaded.cache.filename = file.path(output_path, paste0(region_name, " Loaded.RData"))
    if(!file.exists(loaded.cache.filename)) {
        metagene.obj = metagene$new(regions=region_list[[region_name]], bam_files=design$Sample, assay='chipseq', force_seqlevels=TRUE)
        save(metagene.obj, file=loaded.cache.filename)
    } else {
        load(loaded.cache.filename)
    }
    
    matrix.cache.filename = file.path(output_path, paste0(region_name, " matrix.RData"))
    if(!file.exists(matrix.cache.filename)) {
        bin_count = min(c(unlist(lapply(region_list[[region_name]], width)), 100))
        metagene.obj$produce_table(design = design, normalization="RPM", flip_regions=TRUE, bin_count=bin_count)
        save(metagene.obj, file=matrix.cache.filename)
    } else {
        load(matrix.cache.filename)
    }
    
    df.cache.filename = file.path(output_path, paste0(region_name, " data-frame.RData"))
    if(!file.exists(df.cache.filename)) {
        metagene.obj$produce_data_frame()
        save(metagene.obj, file=df.cache.filename)
    } else {
        load(df.cache.filename)
    }
    
    metagenes[[region_name]] = metagene.obj
}

test = subset(metagene.obj$get_data_frame(), group %in% c("POL2_shCRTL.1_Dex_regions", "POL2_shCTRL.2_Dex_regions", "POL2_shNIPBL.3_Dex_regions", "POL2_shNIPBL.5_Dex_regions"))
pdf("Pol2_Dex.pdf")
ggplot(test, aes(x=(bin-50)*4, y=value, ymin=qinf, ymax=qsup, group=group, color=group, fill=group)) + geom_line()
dev.off()

test = subset(metagene.obj$get_data_frame(), group %in% c("POL2_shCRTL.1_EtOH_regions", "POL2_shCTRL.2_EtOH_regions", "POL2_shNIPBL.3_EtOH_regions", "POL2_shNIPBL.5_EtOH_regions"))
pdf("Pol2_EtOH.pdf")
ggplot(test, aes(x=(bin-50)*4, y=value, ymin=qinf, ymax=qsup, group=group, color=group, fill=group)) + geom_line()
dev.off()

test = subset(metagene.obj$get_data_frame(), group %in% c("POL2.ser2_shCRTL.1_Dex_regions", "POL2.ser2_shCTRL.2_Dex_regions", "POL2.ser2_shNIPBL.3_Dex_regions", "POL2.ser2_shNIPBL.5_Dex_regions"))
pdf("Pol2-ser2_Dex.pdf")
ggplot(test, aes(x=(bin-50)*4, y=value, ymin=qinf, ymax=qsup, group=group, color=group, fill=group)) + geom_line()
dev.off()

test = subset(metagene.obj$get_data_frame(), group %in% c("POL2.ser2_shCRTL.1_EtOH_regions", "POL2.ser2_shCTRL.2_EtOH_regions", "POL2.ser2_shNIPBL.3_EtOH_regions", "POL2.ser2_shNIPBL.5_EtOH_regions"))
pdf("Pol2-ser2_EtOH.pdf")
ggplot(test, aes(x=(bin-50)*4, y=value, ymin=qinf, ymax=qsup, group=group, color=group, fill=group)) + geom_line()
dev.off()
