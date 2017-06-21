require(csaw)

bam.files = Sys.glob("output/ENCODE-chip/alignment/*/*.sorted.dup.bam")
bam.files = bam.files[grepl("GR_.hr", bam.files)]

condition = gsub(".*\\/(.*)_..sorted.dup.bam", "\\1", bam.files)

param <- readParam(minq=50)
count.data <- windowCounts(bam.files, ext=200, width=10, param=param)

# 2. Filtering out uninteresting regions.
require(edgeR)
keep <- aveLogCPM(asDGEList(count.data)) >= -1
count.data <- count.data[keep,]

# 3. Calculating normalization factors.
binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
normfacs <- normOffsets(binned)

# 4. Identifying DB windows.
design <- model.matrix(~factor(condition))
y <- asDGEList(count.data, norm.factors=normfacs)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit)

# 5. Correcting for multiple testing.
merged <- mergeWindows(rowRanges(count.data), tol=1000L)
tabcom <- combineTests(merged$id, results$table)

gr.results = merged$region
mcols(gr.results) = tabcom

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)

annotated.results = annotatePeak(gr.results, tssRegion = c(-1000, 1000), 
                                 TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, level="gene",
                                 annoDb="org.Hs.eg.db", addFlankGeneInfo = TRUE)
final.results = as.data.frame(annotated.results)
final.results = final.results[final.results$FDR <= 0.05,]
write.table(final.results, file="output/analyses/ENCODE-csaw-results.txt", sep="\t", col.names=TRUE, row.names=FALSE)
export(GRanges(final.results), "output/analyses/ENCODE-csaw-results.bed", format="bed")
