source("../mugqic_utils/summarize_dge.R")
all_dge = summarize.dge("A549/output/rna-pipeline-GRCh38", "edgeR")

get_sig <- function(contrast, direction="both", fc_threshold=log2(1.5), p_threshold=0.05) {
    fc_column = paste0("log_FC.", contrast)
    p_column = paste0("edger.adj.p.value.", contrast)
    
    dir_logical = TRUE
    if(direction=="up") {
        dir_logical = all_dge[[fc_column]] > 0
    } else if(direction=="down"){
        dir_logical = all_dge[[fc_column]] < 0
    }
    
    return (all_dge[[p_column]] <= p_threshold & !is.na(all_dge[[p_column]]) &
            abs(all_dge[[fc_column]]) >= fc_threshold & !is.na(all_dge[[fc_column]]) &
            dir_logical)
}

get_sig_entrez <- functon(contrast, direction="both", fc_threshold=log2(1.5), p_threshold=0.05) {
    
}

sum(get_sig("shNIPBL_vs_shMamm-EtOH", "up"), na.rm=TRUE)
sum(get_sig("shNIPBL_vs_shMamm-EtOH", "down"), na.rm=TRUE)
sum(get_sig("shNIPBL_vs_shMamm-EtOH"), na.rm=TRUE)

library(ef.utils)
library(ChIAnalysis)

chia.params = build_chia_params(biosample = "A549",
                                genome.build = "hg19",
                                tssRegion = c(-1000, 1000))
cache(chia.params <- add_encode_data(chia.params))

# Get ENCODE TF bindings.

library(org.Hs.eg.db)
all_dge$entrez_ids = mapIds(org.Hs.eg.db, keys=all_dge$id, keytype="ENSEMBL", column="ENTREZID")

annot_list = select_annotations("hg19")

# Get promoters of all genes
# Get the transcription regions from the database.
tx_regions = AnnotationDbi::select(annot_list$TxDb, all_dge$entrez_ids, c("TXCHROM", "TXSTART", "TXEND", "TXSTRAND"), "GENEID")

# Keep only the first record for each gene.
tx_regions = tx_regions[match(all_dge$entrez_ids, tx_regions$GENEID),]

# Remove NAs
tx_regions = tx_regions[!is.na(tx_regions$TXCHROM),]
tx_regions$Ensembl = mapIds(org.Hs.eg.db, keys=tx_regions$GENEID, keytype="ENTREZID", column="ENSEMBL")

with.tf = associate_tf(GRanges(tx_regions), chia.params$tf.regions)

tf_names = gsub("TF.overlap.", "", grep("TF.overlap", names(mcols(with.tf)), value=TRUE))
diff_ids = 

res_df = NULL
for(tf in tf_names) {
    # Determien which genes have TF binding sites in their promoter.
    tf_presence = mcols(with.tf)[[paste0("TF.overlap.", tf)]] > 0

    # Determine which genes are DE.
    diff_genes = all_dge$id[get_sig("shNIPBL_vs_shMamm-EtOH", "down")]
    diff_genes_indices = mcols(with.tf)$Ensembl %in% diff_genes

    # Get metrics.
    total_with_TF = sum(tf_presence)
    total_without_TF = sum(!tf_presence)
    total_diff = sum(diff_genes_indices)
    total_not_diff = sum(!diff_genes_indices)
    diff_with_tf = sum(tf_presence & diff_genes_indices)
    diff_without_tf = sum(!tf_presence & diff_genes_indices)
    total_proportion = total_with_TF/(total_with_TF+total_without_TF)
    diff_proportion = diff_with_tf/(diff_with_tf + diff_without_tf)
    if(total_proportion != 0) {
        enrichment = log2(diff_proportion / total_proportion)    
    } else {
        enrichment = NA
    }

    p_val = phyper(diff_with_tf, total_with_TF, total_without_TF, total_diff, lower.tail=FALSE)
    partial_df = data.frame(TF=tf, GenesWithBinding=total_with_TF, GenesWithoutBinding=total_without_TF,
                            DiffGeneWithBinding=diff_with_tf, DiffGeneWithoutBinding=diff_without_tf,
                            TotalProportion=total_proportion, DiffProportion=diff_proportion,
                            Enrichment = enrichment, PVal=p_val)
                            
    if(is.null(res_df)) {
        res_df = partial_df
    } else {
        res_df = rbind(res_df, partial_df)
    }
}
