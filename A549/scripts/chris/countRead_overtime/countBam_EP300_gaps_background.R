setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")
source("scripts/chris/countRead_overtime/load_Reddy_EP300_consensus.R")
library(Rsamtools)
library(knitr)
library(dplyr)

####
generate_coordVector <- function(granges) {
  df <- as.data.frame(granges)
  coord_Vector <- paste0(df$seqnames, ":", df$start, "-", df$end)
  return (coord_Vector)
}

#### Gather all EP300 peaks all over the time frame
ep300_regions <- load_reddy_ep300_binding_consensus()

all_ep300_regions <- GRanges()
for (name in names(ep300_regions)) {
  all_ep300_regions <- c(all_ep300_regions, ep300_regions[[name]])
  message(length(ep300_regions[[name]]))
}

all_ep300_regions_reduced <- reduce(all_ep300_regions)
summary(width(all_ep300_regions_reduced))

#### Get the background: all regions where GR are not binding on the genome
gapbackground_ep300_regions_reduced <- gaps(all_ep300_regions_reduced)
summary(width(gapbackground_ep300_regions_reduced))
length(gapbackground_ep300_regions_reduced)

peaks_EP300_coordVector <- generate_coordVector(gapbackground_ep300_regions_reduced)

##### List of bam EP300
all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
ep300_bam <- all_chip_bam %>% filter(target == "EP300", assembly == "GRCh38", lab == "Tim Reddy, Duke")
report_ep300_bam <- ep300_bam %>% select(accession, file_accession, submitted_by, file_format, target, treatment, treatment_duration, treatment_duration_unit, biological_replicates, controls)
report_ep300_bam$treatment_duration[is.na(report_ep300_bam$treatment_duration)] <- 0
report_ep300_bam$treatment_duration_unit[is.na(report_ep300_bam$treatment_duration_unit)] <- "minute"
report_ep300_bam <- report_ep300_bam %>% arrange(desc(treatment_duration_unit), treatment_duration, submitted_by, biological_replicates)
kable(report_ep300_bam)

bam_path <- paste0("/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/",
                   paste0(report_ep300_bam$target, "_", report_ep300_bam$treatment_duration, report_ep300_bam$treatment_duration_unit,
                          "_rep", report_ep300_bam$biological_replicates,
                          "_", report_ep300_bam$file_accession, ".bam"))

##### index BAM
index_bam <- function(bam_path) {
  for (bam in bam_path) {
    bai <- gsub("\\.bam", ".bai", bam)
    cmd_line <- paste("samtools index", bam, bai)
    cat(cmd_line, "\n")
    system(cmd_line)
  }
}

# index_bam(bam_path)


##### Count reads
which <- gapbackground_ep300_regions_reduced
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)

count_total_EP300 <- data.frame(peaks_EP300_coordVector)
for (bam in bam_path) {
  message(bam)
  index_bam <- gsub("\\.bam", "", bam)
  count_bySample <- countBam(bam, index = index_bam, param = param)
  records <- count_bySample$records
  count_total_EP300 <- data.frame(count_total_EP300, records)
}

names(count_total_EP300) <- c("Coordinates", paste0(report_ep300_bam$target, "_", report_ep300_bam$treatment_duration, report_ep300_bam$treatment_duration_unit,
                                                 "_rep", report_ep300_bam$biological_replicates,
                                                 "_", report_ep300_bam$file_accession))

output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"
save(count_total_EP300, file = file.path(output_path, "count_total_EP300_gaps_background.RData"))
write.table(count_total_EP300, file = file.path(output_path, "count_total_EP300_gaps_background.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
