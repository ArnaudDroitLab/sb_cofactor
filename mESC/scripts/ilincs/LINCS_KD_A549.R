setwd("/home/chris/Bureau/sb_cofactor_hr/mESC")

library(dplyr)
library(httr)
library(jsonlite)
source("scripts/ilincs/lincs.utils.R")

mutated_cofactors <- c("NIPBL", "SMC1A", "SMC3", "RAD21", "HDAC8", "BRD4", "CREBBP", "EP300", "KAT6A", "KAT6B",
                       "SRCAP", "BRPF1", "RAI1", "MBD5", "EHMT1", "HDAC4", "ARID1A", "ARID1B", "SMARCA4", "SMARCB1",
                       "SMARCE1", "ANKRD11", "KMT2A", "SMARCA2", "AFF4", "KMT2D", "KDM6A", "TAF6", "ESCO2", "CHD7",
                       "DDX11", "STAG1", "STAG2", "SETD5", "MED12", "MED13L", "MED17", "MED23", "MED25") # 39

all_signatures <- read.csv("input/LINCS/LINCS_consensus_gene_KD_Signatures_all_37275.xls", sep = "\t")

all_cell_lines <- table(all_signatures$CellLine)

A549_cell_line <- all_signatures %>% filter(CellLine == "A549")
unique(A549_cell_line$Time)

signIds_A549 <- get_signIds(A549_cell_line, cell_line = "A549")
targets_A549 <- get_target(A549_cell_line, cell_line = "A549")
sum(mutated_cofactors %in% target_A549)

###

downloadSignature <- function(signId) {
  www <- "http://www.ilincs.org/api/ilincsR/downloadSignature"
  
  call <- paste(www, "?", "sigID=", signId, sep="")
  request <- GET(call)
  content <- content(request, "text")
  df <- fromJSON(txt = content)
  res <- df$data$signature %>% select(Name_GeneSymbol, Value_LogDiffExp)
  return(res)
}

downloadSignatureInBatch <- function(signIds, targets) {
  l <- list()
  for (i in 1:length(signIds)) {
    sign <- signIds[i]
    target <- targets[i]
    tmp <- downloadSignature(sign)
    tmp_l <- list(tmp)
    names(tmp_l) <- target
    l <- append(l, tmp_l)
  }
  return(l)
}

# signIds <- "LINCSKD_5664"
ten_signIds <- signIds_A549[1:5]
ten_target_A549 <- targets_A549[1:5]

t <- downloadSignatureInBatch(ten_signIds, ten_target_A549)
