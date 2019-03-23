# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

source("scripts/ckn_utils.R")
source("scripts/load_reddy.R")
library(ef.utils)
library(ChIPseeker)

# Loading peaks
peaks_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NB"
gainNB <- rtracklayer::import(con = file.path(peaks_dir, "NB_specific_DEX.bed")); print(length(gainNB)) # 803
lossNB <- rtracklayer::import(con = file.path(peaks_dir, "NB_specific_CTRL.bed")); print(length(lossNB)) # 5952
commonNB <- rtracklayer::import(con = file.path(peaks_dir, "NB_common.bed")); print(length(commonNB)) # 1680

###
chrom_state_hg19 <- download_chromatin_states("A549", "hg19")
hg19_to_hg38_chain <- rtracklayer::import.chain("input/hg19ToHg38.over.chain")
chrom_state_hg38_unordered <- unlist(rtracklayer::liftOver(chrom_state_hg19, hg19_to_hg38_chain))
chrom_state_hg38 = chrom_state_hg38_unordered[order(chrom_state_hg38_unordered$name)]

###
assign_chrom_state_simplify <- function(regions, chrom_state) {
  unique.states = sort(unique(chrom_state$name))
  state.overlap.indices = findOverlaps(regions, chrom_state, select="first")
  matching.chrom.state = factor(chrom_state$name[state.overlap.indices], levels=unique.states)
  regions_chrom_state <- simplify_chrom_state(matching.chrom.state)
  message(sum(!is.na(regions_chrom_state)), " / ", length(regions))
  
  table_chrom_state <- table(regions_chrom_state, useNA = "ifany")
  print(table_chrom_state)
  
  print(plot_chrom_state(table_chrom_state))
  
  return(regions_chrom_state)
}

assign_chrom_state_whole <- function(regions, chrom_state) {
  unique.states = sort(unique(chrom_state$name))
  state.overlap.indices = findOverlaps(regions, chrom_state, select="first")
  matching.chrom.state = factor(chrom_state$name[state.overlap.indices], levels=unique.states)
  regions_chrom_state <- matching.chrom.state
  message(sum(!is.na(regions_chrom_state)), " / ", length(regions))
  
  table_chrom_state <- table(regions_chrom_state, useNA = "ifany")
  print(table_chrom_state)
  
  print(plot_chrom_state(table_chrom_state))
  
  return(regions_chrom_state)
}

plot_chrom_state <- function(table_chrom_state) {
  data <- as.data.frame(table_chrom_state)
  data$regions_chrom_state <- as.character(data$regions_chrom_state)
  data$regions_chrom_state[is.na(data$regions_chrom_state)] <- "NA"
  p <- plot_ly(data, labels = ~regions_chrom_state, values = ~Freq, type = "pie",
               textinfo = "label+percent",
               insidetextfont = list(color = "#FFFFFF"))
  return(p)
}

###
lossNB_chrom <- assign_chrom_state_simplify(lossNB, chrom_state_hg38)
commonNB_chrom <- assign_chrom_state_simplify(commonNB, chrom_state_hg38)
gainNB_chrom <- assign_chrom_state_simplify(gainNB, chrom_state_hg38)

lossNB_chrom <- assign_chrom_state_whole(lossNB, chrom_state_hg38)
commonNB_chrom <- assign_chrom_state_whole(commonNB, chrom_state_hg38)
gainNB_chrom <- assign_chrom_state_whole(gainNB, chrom_state_hg38)

###
lossNB_anno <- annotatePeaks(lossNB, output = "anno")
commonNB_anno <- annotatePeaks(commonNB, output = "anno")
gainNB_anno <- annotatePeaks(gainNB, output = "anno")

l <- list("NB_spe_CTRL" = lossNB_anno, "NB_common" = commonNB_anno, "NB_spe_DEX" = gainNB_anno)

plotAnnoPie(gainNB_anno)
plotAnnoPie(commonNB_anno)
plotAnnoPie(lossNB_anno)

plotAnnoBar(l)
plotDistToTSS(l)
