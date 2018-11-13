# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(ef.utils)
library(plotly)

nbc_peaks_dir <- file.path("output/chip-pipeline-GRCh38/peak_call/A549_NBC")

###
chrom_state_hg19 <- download_chromatin_states("A549", "hg19")
hg19_to_hg38_chain <- rtracklayer::import.chain("input/hg19ToHg38.over.chain")
chrom_state_hg38_unordered <- unlist(rtracklayer::liftOver(chrom_state_hg19, hg19_to_hg38_chain))
chrom_state_hg38 = chrom_state_hg38_unordered[order(chrom_state_hg38_unordered$name)]

###
nbc_spectrl_path <- file.path(nbc_peaks_dir, "A549_NBC_CTRL_specific.bed")
nbc_common_path <- file.path(nbc_peaks_dir, "A549_NBC_common.bed")
nbc_spedex_path <- file.path(nbc_peaks_dir, "A549_NBC_DEX_specific.bed")

nbc_spectrl <- rtracklayer::import(nbc_spectrl_path) # 4176
nbc_common <- rtracklayer::import(nbc_common_path) # 1376
nbc_spedex <- rtracklayer::import(nbc_spedex_path) # 808

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

nbc_spectrl_chrom <- assign_chrom_state_simplify(nbc_spectrl, chrom_state_hg38)
nbc_common_chrom <- assign_chrom_state_simplify(nbc_common, chrom_state_hg38)
nbc_spedex_chrom <- assign_chrom_state_simplify(nbc_spedex, chrom_state_hg38)

nbc_spectrl_chrom <- assign_chrom_state_whole(nbc_spectrl, chrom_state_hg38)
nbc_common_chrom <- assign_chrom_state_whole(nbc_common, chrom_state_hg38)
nbc_spedex_chrom <- assign_chrom_state_whole(nbc_spedex, chrom_state_hg38)