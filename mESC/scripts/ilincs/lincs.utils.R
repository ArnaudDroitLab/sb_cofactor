
get_signIds <- function(df, cell_line) {
  res <- df %>% filter(CellLine == cell_line) %>% pull(SignatureId) %>% as.character
  return(res)
}

get_target <- function(df, cell_line) {
  res <- df %>% filter(CellLine == cell_line) %>% pull(TargetGene) %>% as.character
  return(res)
}
