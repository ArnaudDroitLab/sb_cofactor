# The start coordinates provided by the TxDb package are far sometimes far-removed from
# the actual TSSes from the gene. A particularly egregious offender is the ACTB gene,
# whose most used TSS is more than 20 kb away from the most used start site, or more than
# 7 time the length of the whole gene.
#
# To remedy this, we obtain per-transcript expression levels obtained from ENCODE
# (polyA mRNA data in A549 cells treated with 500pM dexamethasone). The dataset is 
# deprecated and based on the hg19 assembly of the human genome, but should be
# good enough for choosing between isoforms in the vast majority of cases.
#
# We annotate the data using bioMart, and filter out transcripts based on evidence
# (keeping only those with the strongest evidence) and removing those without
# an ENTREZ ID (to remove some non-coding/non-descript transcripts that we
# probably wouldn'T care about anyway).
# 
# This relies on hitting the internet twice, and thus cannot be executed from
# a node on cedar or graham.

library(biomaRt)
library(dplyr)
require(ENCODExplorer)
require(GenomicRanges)


# Read in per-transcript expression
all_res = ENCODExplorer::queryEncodeGeneric(biosample_name="A549",
                                            output_type="transcript quantifications",
                                            assembly="GRCh38")
downloaded_files = ENCODExplorer::downloadEncode(file_acc=all_res$file_accession)
rna_results = lapply(downloaded_files, read.table, header=TRUE, stringsAsFactors=FALSE)
all_results = do.call(cbind, rna_results)

# Make sure all row ids match.
id_columns = colnames(all_results) == "transcript_id"
single_value = function(x) { length(unique(x))==1 }
stopifnot(all(apply(all_results[,id_columns], 1, single_value)))

# Average TPMs
tpm_columns = colnames(all_results) == "TPM"
tpm_mean = apply(all_results[,tpm_columns], 1, mean)

transcript_expression = data.frame(ensembl_transcript_id=all_results$transcript_id,
                                   TPM=tpm_mean,
                                   stringsAsFactors=FALSE)

#ENCODExplorer::downloadEncode(file_acc="ENCFF289BZX")
#transcript_expression = read.table("ENCFF289BZX.tsv", sep="\t", header=TRUE)
#transcript_expression$ensembl_transcript_id = gsub("\\..*$", "", transcript_expression$transcript_id)

# Annotate individual transcripts.
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
attribute_set = c('ensembl_transcript_id', "entrezgene_id", 'external_gene_name',
                  'chromosome_name', 'start_position', 'end_position', 'strand',
                  'transcript_start', 'transcript_end', 'transcription_start_site', 
                  'transcript_tsl', "entrezgene_trans_name", "ensembl_gene_id")
# Ensembl is very buggy lately and keeps ond isconnecting before the 
# requests are done. We'll split up the request in ~10 so all of them complete under 30 seconds.

transcript_expression$ensembl_transcript_id_no_version = gsub("\\..*$", "", transcript_expression$ensembl_transcript_id)
cut_intervals = split(transcript_expression$ensembl_transcript_id_no_version,
                      cut(seq_along(transcript_expression$ensembl_transcript_id_no_version), breaks=40))
all_bm = lapply(cut_intervals, function(x) {
    res <- NULL
    attempt <- 0
    while( is.null(res) && attempt <= 5 ) {
      attempt <- attempt + 1
      message("attempt ", attempt)
      try({
        res = getBM(attributes=attribute_set, 
              filters = 'ensembl_transcript_id', 
              values = x, 
              mart = ensembl)
      })
    } 
    
    res
})

all_annot = do.call(rbind, all_bm)

transcript_expression = dplyr::left_join(transcript_expression, all_annot, by=c("ensembl_transcript_id_no_version"="ensembl_transcript_id"))
                         
# Filter transcripts to remove those based on poor evidence.
# Also filter transcripts to remove those on unknown chromosome to spare us a lot of
# pain down the line when dealing with sequence levels.
transcript_expression_evidence = transcript_expression %>%
    filter(grepl("tsl1", transcript_tsl), !is.na(entrezgene_id), !grepl("CHR_", chromosome_name))
        

transcript_expression_evidence = transcript_expression_evidence %>% 
    filter(!grepl("CHR_", chromosome_name), !is.na(entrezgene_id))
        
# Finally, of all remaining transcript, choose the most expressed one.
most_expressed = transcript_expression_evidence %>% group_by(entrezgene_id) %>% dplyr::slice(which.max(TPM))                         
most_expressed$strand = ifelse(most_expressed$strand == -1, "-", "+")
most_expressed$start = most_expressed$transcript_start
most_expressed$end = most_expressed$transcript_end
most_expressed$chromosome_name = paste0("chr", most_expressed$chromosome_name)

all_genes = GenomicRanges::GRanges(most_expressed)

write.table(as.data.frame(all_genes), file="output/analyses/tss_gene_coordinates.txt", col.names=TRUE, row.names=FALSE, sep="\t")
