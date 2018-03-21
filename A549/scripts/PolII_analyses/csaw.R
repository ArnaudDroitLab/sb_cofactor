require(csaw)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
require(edgeR)

output_dir = "output/chip-pipeline-GRCh38"

# Determine bam file paths and conditions.
bam_files = Sys.glob(file.path(output_dir, "alignment/*/*.sorted.dup.bam"))
bam_files = bam_files[grepl("POL2", bam_files)]

sample_name = gsub(".*alignment\\/(.*)\\/.*", "\\1", bam_files)
# Fix spelling error
sample_name = gsub("shCRTL", "shCTRL", sample_name)

antibody = gsub("_sh.*", "", sample_name)
sh_name = gsub(".*_(sh.*)_.*_Rep.*", "\\1", sample_name)
sh_target = gsub(".*_(sh.*)-._.*_Rep.*", "\\1", sample_name)
sh_target = factor(sh_target, levels=c("shCTRL", "shNIPBL"))
condition = gsub(".*_(EtOH|Dex)_.*", "\\1", sample_name)
condition = factor(condition, levels=c("EtOH", "Dex"))

# Determine bed file paths and grab all regions.
bed_files = Sys.glob(file.path(output_dir, "peak_call/*/*.broadPeak"))
extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric")
bed_regions_list= lapply(bed_files, rtracklayer::import, format="BED", extraCols=extraCols)    
all_peak_regions = reduce(unlist(GRangesList(bed_regions_list)))


param <- readParam(minq=50)
count_data <- windowCounts(bam_files, ext=200, width=10, param=param)

# 2. Filtering out uninteresting regions.

#keep <- aveLogCPM(asDGEList(count_data)) >= -1
keep <- overlapsAny(rowRanges(count_data), all_peak_regions)
count_data <- count_data[keep,]

# 3. Calculating normalization factors.
binned <- windowCounts(bam_files, bin=TRUE, width=10000, param=param)
normfacs <- normOffsets(binned)

identify_db_windows <- function(sample_subset, design_variable, label) {
    count_data_subset = count_data[,sample_subset]
    norm_factors_subset = normfacs[sample_subset]
    
    sh=sh_target[sample_subset]
    cond=condition[sample_subset]
    if(design_variable=="sh") {
        design=model.matrix(~sh)
    } else {
        design=model.matrix(~cond)
    }
    
    y <- asDGEList(count_data_subset, norm.factors=norm_factors_subset)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    results <- glmQLFTest(fit)
    
    # 5. Correcting for multiple testing.
    merged <- mergeWindows(rowRanges(count_data_subset), tol=1000L)
    tabcom <- combineTests(merged$id, results$table)
    
    gr.results = merged$region
    mcols(gr.results) = tabcom
    
    annotated.results = annotatePeak(gr.results, tssRegion = c(-1000, 1000), 
                                     TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, level="gene",
                                     annoDb="org.Hs.eg.db", addFlankGeneInfo = TRUE)
    final.results = as.data.frame(annotated.results)
    if(sum(final.results$FDR <= 0.05) != 0) {
        final.results = final.results[final.results$FDR <= 0.05,]
        write.table(final.results, file=paste0("csaw results ", label, ".txt"), sep="\t", col.names=TRUE, row.names=FALSE)
        export(GRanges(final.results), paste0("csaw results ", label, ".bed"), format="bed")
    } else {
        cat("No results found for ", label)
    }
}

# 4. Identifying DB windows for sh effects.
for(antibody_name in unique(antibody)) {
    antibody_subset = antibody==antibody_name
    for(condition_name in unique(condition)) {
        condition_subset = condition==condition_name
        sample_subset = antibody_subset & condition_subset

        identify_db_windows(sample_subset, "sh", paste0(antibody_name, " DB under the effect of shNIPBL in ", condition_name))
    }
}

# 5. Identifying DB windows for condition effects
for(antibody_name in unique(antibody)) {
    antibody_subset = antibody==antibody_name
    for(sh_target_name in unique(sh_target)) {
        sh_subset = sh_target==sh_target_name
        sample_subset = antibody_subset & sh_subset

        identify_db_windows(sample_subset, "cond", paste0(antibody_name, " DB under the effect of Dex in ", sh_target_name))
    }
}
