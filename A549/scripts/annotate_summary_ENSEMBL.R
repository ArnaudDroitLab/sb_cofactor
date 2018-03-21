source("/gs/home/efournier/mugqic_utils/summarize_dge.R")

dge = summarize.dge("output/rna-pipeline", "edgeR", subdir="DGE_ENSEMBL")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
dge$EntrezID = mapIds(TxDb.Hsapiens.UCSC.hg19.knownGene, dge$id, "GENEID", "TXNAME")
dge$SYMBOL = mapIds(org.Hs.eg.db, dge$EntrezID, "SYMBOL", "ENTREZID")

test$SYMBOL= symbol.list[test$EntrezID]
write.table(as.data.frame(test), "/gs/scratch/efournier/CofactorHR/A549/output/rna-pipeline/DGE/summary.tsv",
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)