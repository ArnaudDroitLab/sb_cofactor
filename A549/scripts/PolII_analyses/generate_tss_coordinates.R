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
ENCODExplorer::downloadEncode(file_acc="ENCFF289BZX")
transcript_expression = read.table("ENCFF289BZX.tsv", sep="\t", header=TRUE)
transcript_expression$ensembl_transcript_id = gsub("\\..*$", "", transcript_expression$transcript_id)

# Annotate individual transcripts.
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")
attribute_set = c('ensembl_transcript_id', 'entrezgene', 'external_gene_name',
                  'chromosome_name', 'start_position', 'end_position', 'strand',
                  'transcript_start', 'transcript_end', 'transcription_start_site', 
                  'transcript_tsl')
gene_coordinates = getBM(attributes=attribute_set, 
                         filters = 'ensembl_transcript_id', 
                         values = transcript_expression$ensembl_transcript_id, 
                         mart = ensembl)
transcript_expression = dplyr::left_join(transcript_expression, gene_coordinates, by="ensembl_transcript_id")
                         
# Filter transcripts to remove those based on poor evidence.
transcript_expression_evidence = transcript_expression %>% filter(grepl("tsl1", transcript_tsl), !is.na(entrezgene))
                         
# Finally, of all remaining transcript, choose the most expressed one.
most_expressed = transcript_expression_evidence %>% group_by(gene_id) %>% dplyr::slice(which.max(TPM))                         
most_expressed$strand = ifelse(most_expressed$strand == -1, "-", "+")
most_expressed$start = most_expressed$transcript_start
most_expressed$end = most_expressed$transcript_end
most_expressed$chromosome_name = paste0("chr", most_expressed$chromosome_name)

all_genes = GenomicRanges::GRanges(most_expressed)

write.table(as.data.frame(all_genes), file="output/analyses/tss_gene_coordinates.txt", col.names=TRUE, row.names=FALSE, sep="\t")
