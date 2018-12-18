library(Rsamtools)
library(knitr)

# countBam(file, index=file, ..., param=ScanBamParam())
# file: The character(1) file name of the `BAM' ('SAM' for asBam) file to be processed.
# index: The character(1) name of the index file of the 'BAM' file being processed; this is given without the '.bai' extension.
# param: An instance of ScanBamParam. This influences what fields and which records are imported.

# ScanBamParam class are provide on its help page;
# several salient points are reiterated here.
# ScanBamParam can contain a field what, specifying the components of the BAM records to be returned.
# Valid values of what are available with scanBamWhat.
# ScanBamParam can contain an argument which that specifies a subset of reads to return.
# This requires that the BAM file be indexed, and that the file be named following samtools convention as .bai.
# ScanBamParam can contain an argument tag to specify which tags will be extracted.

which <- IRangesList(seq1=IRanges(1000, 2000), seq2=IRanges(c(100, 1000), c(1000, 2000)))
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)

bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools")
bam <- scanBam(bamFile, param=param)
count <- countBam(bamFile, param=param)

#### Gather all GR peaks all over the time frame
gr_regions <- load_reddy_gr_binding_consensus()

all_gr_regions <- GRanges()
for (name in names(gr_regions)) {
  all_gr_regions <- c(all_gr_regions, gr_regions[[name]])
  message(length(gr_regions[[name]]))
}

all_gr_regions_reduced <- reduce(all_gr_regions)
summary(width(all_gr_regions_reduced))

##### List of bam GR
all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
gr_bam <- all_chip_bam %>% filter(target == "NR3C1", assembly == "GRCh38", lab == "Tim Reddy, Duke")
report_gr_bam <- gr_bam %>% select(accession, file_accession, submitted_by, file_format, target, treatment, treatment_duration, treatment_duration_unit, biological_replicates, controls)
report_gr_bam$treatment_duration[is.na(report_gr_bam$treatment_duration)] <- 0
report_gr_bam$treatment_duration_unit[is.na(report_gr_bam$treatment_duration_unit)] <- "minute"
report_gr_bam <- report_gr_bam %>% arrange(desc(treatment_duration_unit), treatment_duration, submitted_by, biological_replicates)
kable(report_gr_bam)

bam_path <- paste0("/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/",
                   paste0(report_gr_bam$target, "_", report_gr_bam$treatment_duration, report_gr_bam$treatment_duration_unit,
                          "_rep", report_gr_bam$biological_replicates,
                          "_", report_gr_bam$file_accession, ".bam"))

##### Count reads
which <- all_gr_regions_reduced
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)

bamFile <- bam_path[1]
count <- countBam(bamFile, index = "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/NR3C1_0minute_rep1_ENCFF668EHX", param=param)
