setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(DiffBind)
source("scripts/ckn_utils.R")

###
perform_diffbind <- function(pol, cst, effect, peak) {
  message("######\t", pol, " | ", cst, " | ", effect, " effect | ", peak, "Peak")
  filename = paste0("sSheet_", pol, "_", cst, "_", effect, "_effect_", peak, "Peak.csv")
  message("   # >>> ", filename)
  sSheet <- read.table(file.path("scripts/chris/diffbind_pol2/sampleSheet", filename), sep= "," , header = TRUE)
  
  dba <- dba(sampleSheet = sSheet)
  
  count <- dba.count(dba)

  if (effect == "DEX") {category = DBA_TREATMENT} else if (effect == "shNIPBL") {category = DBA_CONDITION}
  message(effect)
  message(category)
  contrast <- dba.contrast(count, categories = category, minMembers = 2)

  analyze <- dba.analyze(contrast)

  report <- dba.report(analyze, bCounts = T)

  annodf <- annotatePeaks2(report, output = "df"); print(nrow(annodf)); print(sum(annodf$Annot == "Promoter")); print(sort(annodf$SYMBOL))

  annodf_filename <- paste0("diffbind_", pol, "_", cst, "_", effect, "_effect_", peak, ".txt")
  output_path <- file.path("output/analyses/diffbind_pol2", annodf_filename)
  write.table(annodf, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE)
}

###
for (pol in c("POL2")) {
  for (peak in c("narrow", "broad")) {
    for (cst in c("CTRL", "DEX")) {
      effect <- "shNIPBL"
      perform_diffbind(pol, cst, effect, peak)
    }
    for (cst in c("shCTRL", "shNIPBL")) {
      effect <- "DEX"
      perform_diffbind(pol, cst, effect, peak)
    }
  }
}

# 
# # Test Annotation
# pol <- "POL2"
# cst <- "shCTRL"
# effect <- "DEX"
# peak <- "narrow"
# 
# message("######\t", pol, " | ", cst, " | ", effect, " effect | ", peak, "Peak")
# filename = paste0("sSheet_", pol, "_", cst, "_", effect, "_effect_", peak, "Peak.csv")
# message("   # >>> ", filename)
# sSheet <- read.table(file.path("scripts/chris/diffbind_pol2/sampleSheet", filename), sep= "," , header = TRUE)
# 
# dba <- dba(sampleSheet = sSheet)
# 
# count <- dba.count(dba)
# 
# if (effect == "DEX") {category = DBA_TREATMENT} else if (effect == "shNIPBL") {category = DBA_CONDITION}
# message(effect)
# message(category)
# contrast <- dba.contrast(count, categories = category, minMembers = 2)
# 
# analyze <- dba.analyze(contrast)
# 
# report <- dba.report(analyze, bCounts = T, bUsePval = TRUE)
# report <- dba.report(analyze, bCounts = T)
# 
# annodf1 <- annotatePeaks(report, output = "df"); print(nrow(annodf1)); print(sum(annodf1$Annot == "Promoter")); print(sort(annodf1$SYMBOL))
# annodf2 <- annotatePeaks2(report, output = "df"); print(nrow(annodf2)); print(sum(annodf2$Annot == "Promoter")); print(sort(annodf2$SYMBOL))
# 
# data.frame(annodf1$Fold, annodf2$Fold)
# annodf1$Fold == annodf2$Fold
# annodf1 <- annodf1 %>% arrange(desc(Fold)) 
# annodf2 <- annodf2 %>% arrange(desc(Fold))
# 
# comparison <- data.frame(annodf1$Fold, annodf1$p.value, annodf1$FDR,
#                          annodf1$Annot, annodf1$distanceToTSS, annodf1$SYMBOL,
#                          annodf2$Annot, annodf2$distanceToTSS, annodf2$SYMBOL,
#                          annodf1$Coordinates)
# 
# comparison$annodf1.SYMBOL <- as.character(comparison$annodf1.SYMBOL)
# comparison$annodf2.SYMBOL <- as.character(comparison$annodf2.SYMBOL)
# 
# comparison %>% filter(annodf1.Annot == "Promoter" | annodf2.Annot == "Promoter")
# comparison[comparison$annodf1.SYMBOL != comparison$annodf2.SYMBOL, ]
