# This script generates a simple correlation/clustering
# of our RNA-Seq results to help diagnose issues
# related to the lack of DE genes found in DEX condition shNIPBL.
library(reshape2)
library(GGally)

rpkm = read.table("output/rna-pipeline-GRCh38/raw_counts/matrixRPKM.tsv")

png("Correlation matrix.png", width=1000, height=1000)
ggpairs(log2(rpkm+0.01))
dev.off()

pdf("Hclust.pdf")
plot(hclust(dist(log2(t(rpkm+0.01)))))
dev.off()
