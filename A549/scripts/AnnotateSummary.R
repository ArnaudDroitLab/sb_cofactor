source("/gs/home/efournier/mugqic_utils/summarize_dge.R/)

test = summarize.dge("/gs/scratch/efournier/CofactorHR/A549/output/rna-pipeline", "edgeR")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
test$EntrezID = mapIds(TxDb.Hsapiens.UCSC.hg19.knownGene, test$id, 
                       "GENEID", "TXNAME") symbol.list = unlist(mapIds(org.Hs.eg.db, test$EntrezID, "SYMBOL", "ENTREZID"))
test$SYMBOL= symbol.list[test$EntrezID]
write.table(as.data.frame(test), "/gs/scratch/efournier/CofactorHR/A549/output/rna-pipeline/DGE/summary.tsv",
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
