mutated_cofactors <- c("NIPBL", "SMC1A", "SMC3", "RAD21", "HDAC8", "BRD4", "CREBBP", "EP300", "KAT6A", "KAT6B",
                       "SRCAP", "BRPF1", "RAI1", "MBD5", "EHMT1", "HDAC4", "ARID1A", "ARID1B", "SMARCA4", "SMARCB1",
                       "SMARCE1", "ANKRD11", "KMT2A", "SMARCA2", "AFF4", "KMT2D", "KDM6A", "TAF6", "ESCO2", "CHD7",
                       "DDX11", "STAG1", "STAG2", "SETD5", "MED12", "MED13L", "MED17", "MED23", "MED25") # 39

get_signIds <- function(df, cell_line) {
  res <- df %>% filter(CellLine == cell_line) %>% pull(SignatureId) %>% as.character
  return(res)
}

get_target <- function(df, cell_line) {
  res <- df %>% filter(CellLine == cell_line) %>% pull(TargetGene) %>% as.character
  return(res)
}

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
  cat(1, " ")
  sign1 <- signIds[1]
  target1 <- targets[1]
  tmp1 <- downloadSignature(sign1)
  colnames(tmp1) <- c("Name_GeneSymbol", target1)
  
  l <- data.frame(tmp1)
  for (i in 2:length(signIds)) {
    cat(i, " ")
    sign <- signIds[i]
    target <- targets[i]
    tmp <- downloadSignature(sign)
    colnames(tmp) <- c("Name_GeneSymbol", target)
    l <- merge(l, tmp, by = "Name_GeneSymbol")
    # tmp_l <- list(tmp)
    # names(tmp_l) <- target
    # l <- append(l, tmp_l)
  }
  return(l)
}