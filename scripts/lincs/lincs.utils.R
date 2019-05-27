library(dplyr)
library(httr)
library(jsonlite)

##### Allow to retrieve nth element of each list from a "list of list"
get_nth_element <- function(lst, n) {
  sapply(lst, "[", n)
}

#####
get_mutated_cofactors <- function() {
  # 39 mutated cofactors associated with developmental syndromes / trasncriptomopathies
  mutated_cofactors <- c("NIPBL", "SMC1A", "SMC3", "RAD21", "HDAC8", "BRD4", "CREBBP", "EP300", "KAT6A", "KAT6B",
                         "SRCAP", "BRPF1", "RAI1", "MBD5", "EHMT1", "HDAC4", "ARID1A", "ARID1B", "SMARCA4", "SMARCB1",
                         "ARID2", "SMARCC2", "DFP2",
                         "SMARCE1", "ANKRD11", "KMT2A", "SMARCA2", "AFF4", "KMT2D", "KDM6A", "TAF6", "ESCO2", "CHD7",
                         "DDX11", "STAG1", "STAG2", "SETD5", "MED12", "MED13L", "MED17", "MED23", "MED25")
  return(mutated_cofactors)
}

#####
get_signIds <- function(df, cell_line, time) {
  res <- df %>% filter(CellLine == cell_line, Time == time) %>% pull(SignatureId) %>% as.character
  return(res)
}

#####
get_target <- function(df, cell_line, time) {
  res <- df %>% filter(CellLine == cell_line, Time == time) %>% pull(TargetGene) %>% as.character
  return(res)
}

####
downloadSignature <- function(signId, pval = FALSE) {
  www <- "http://www.ilincs.org/api/ilincsR/downloadSignature"
  
  call <- paste(www, "?", "sigID=", signId, sep="")
  request <- GET(call)
  content <- content(request, "text")
  df <- fromJSON(txt = content)
  
  if (pval == TRUE) {
    res <- df$data$signature %>% select(Name_GeneSymbol, Significance_pvalue)
  } else {
    res <- df$data$signature %>% select(Name_GeneSymbol, Value_LogDiffExp)
  }
  return(res)
}

#####
downloadSignatureInBatch <- function(signIds, targets) {
  message(1, " / ", length(signIds))
  sign1 <- signIds[1]
  target1 <- targets[1]
  tmp1 <- downloadSignature(sign1)
  colnames(tmp1) <- c("Name_GeneSymbol", target1)
  
  signMat <- data.frame(tmp1)
  if (length(signIds) > 1) {
    for (i in 2:length(signIds)) {
      message(i, " / ", length(signIds))
      sign <- signIds[i]
      target <- targets[i]
      tmp <- downloadSignature(sign)
      colnames(tmp) <- c("Name_GeneSymbol", target)
      signMat <- merge(signMat, tmp, by = "Name_GeneSymbol")
      # tmp_l <- list(tmp)
      # names(tmp_l) <- target
      # l <- append(l, tmp_l)
    }
  }
  return(signMat)
}

#### Download all KD signatures from a specific cell line
# Take too much time
downloadSignature_KD_CellLine <- function(CellLine, time, output_path) {
  message("#####\t Cell line : ", CL)
  CL_cell_line <- all_signatures %>% filter(CellLine == CL, Time == time)
  message("Number of KD signatures : ", nrow(CL_cell_line))
  
  signIds_CL <- get_signIds(CL_cell_line, cell_line = CL, time = time)
  targets_CL <- get_target(CL_cell_line, cell_line = CL, time = time)
  message("Including ", sum(MUTATED_COFACTORS %in% targets_CL), " transcriptional coregulators (mutated in dev. syndromes)")
  
  CL_ILINCs_KD_matrix <- downloadSignatureInBatch(signIds_CL, targets_CL)
  
  time_withoutspace <- gsub(" ", "", time)
  output_filename <- file.path(output_path, paste0(CL, "_", time_withoutspace, "_LINCS_KD_matrix.rds"))
  saveRDS(CL_ILINCs_KD_matrix, file = output_filename)
  message(" >>> Matrix has been save in ", output_filename)
}

#### Save signature matrix
saveSignMat <- function(matrix, cellLine, output_dir, output_file) {
  output_filepath <- file.path(output_dir, paste0(output_file, ".txt"))
  write.table(x = matrix, file = output_filepath, sep = "\t", quote = FALSE, row.names = FALSE, dec = ".")
  message(" > Signature matrix saved in ", output_filepath)
}

#### Save heatmap
saveHeatmap <- function(heatmap_obj, output_dir, output_file, width_val = 25, height_val = 22, format = "pdf") {
  output_filepath <- file.path(output_dir, paste0(output_file, ".", format))
  pdf(file = output_filepath, width = width_val, height = height_val)
  print(heatmap_obj)
  dev.off()
  message(" > Heatmap saved in ", output_filepath)
}

#### Get signature matrix of MUTATED_COFACTORS (KD and OE), in one cell line, at one time point
get_signMat_KDOE <- function(cellLine, time, KD_matrix, OE_matrix) {
  # download KD matrix
  dfKD_cLine <- KD_matrix %>% filter(CellLine == cellLine, Time == time, TargetGene %in% MUTATED_COFACTORS)
  if (nrow(dfKD_cLine) != 0) {
  signIds_KD <- get_signIds(dfKD_cLine, cell_line = cellLine, time = time)
  targets_KD <- paste0("KD_", get_target(dfKD_cLine, cell_line = cellLine, time = time))
  message("> Downloading KD...")
  KD_mat <- downloadSignatureInBatch(signIds_KD, targets_KD)
  }
  
  # download OE matrix if available and merge with KD matrix
  dfOE_cLine <- OE_matrix %>% filter(CellLine == cellLine, Time == time, TargetGene %in% MUTATED_COFACTORS)
  if (nrow(dfOE_cLine) != 0) {
    signIds_OE <- get_signIds(dfOE_cLine, cell_line = cellLine, time = time)
    targets_OE <- paste0("OE_", get_target(dfOE_cLine, cell_line = cellLine, time = time))
    message("> Downloading OE...")
    OE_mat <- downloadSignatureInBatch(signIds_OE, targets_OE)
  }
  
  if ((nrow(dfKD_cLine) != 0) && (nrow(dfOE_cLine) != 0)) {
    signMat_wGeneName <- merge(KD_mat, OE_mat, by = "Name_GeneSymbol")
  } else if (nrow(dfOE_cLine) == 0) {
    signMat_wGeneName <- KD_mat
  } else {
    signMat_wGeneName <- OE_mat
  }
  return(signMat_wGeneName)
}