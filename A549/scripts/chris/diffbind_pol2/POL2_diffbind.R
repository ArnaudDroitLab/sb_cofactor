setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(DiffBind)
source("scripts/ckn_utils.R")
source("scripts/reddy_time_series/draw_graph_log2FC_0-12h.R")

###
perform_diffbind <- function(pol, cst, effect, peak, pval) {
  message("######\t", pol, " | ", cst, " | ", effect, " effect | ", peak, "Peak")
  
  basename <- paste0(pol, "_", cst, "_", effect, "_effect_", peak)
  
  filename = paste0("sSheet_", basename, "Peak.csv")
  message("   # >>> ", filename)
  
  sSheet <- read.table(file.path("scripts/chris/diffbind_pol2/sampleSheet", filename), sep= "," , header = TRUE)
  
  dba <- dba(sampleSheet = sSheet)
  
  message("Counting...")
  count <- dba.count(dba)

  if (effect == "DEX") {
    category = DBA_TREATMENT
  } else if (effect == "shNIPBL") {
      category = DBA_CONDITION
      }
  message("effect : ", effect , " | ", category)
  
  contrast <- dba.contrast(count, categories = category, minMembers = 2)

  analyze <- dba.analyze(contrast)
  
  if (pval) {
  report <- dba.report(analyze, bCounts = T, bUsePval = pval)
  } else {
    report <- dba.report(analyze, bCounts = T)
  }

  annodf <- annotatePeaks(report, output = "df", tss = 5000)
  # print(nrow(annodf)); print(sum(annodf$Annot == "Promoter")); print(sort(unique(annodf$SYMBOL)))

  if (pval) {
    annodf_filename <- paste0("diffbind_", basename, "_pval.txt")
  } else {
      annodf_filename <- paste0("diffbind_", basename, "_fdr.txt")
      }
  
  output_path <- file.path("output/analyses/diffbind_pol2", annodf_filename)
  write.table(annodf, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE)
}

###
stats_diffbind <- function(pol, cst, effect, peak, pval, geneList = FALSE) {
  message("######\t", pol, " | ", cst, " | ", effect, " effect | ", peak, "Peak")
  
  basename <- paste0(pol, "_", cst, "_", effect, "_effect_", peak)
  
  if (pval) {
    annodf_filename <- paste0("diffbind_", basename, "_pval.txt")
  } else {
    annodf_filename <- paste0("diffbind_", basename, "_fdr.txt")
    }
  
  output_path <- file.path("output/analyses/diffbind_pol2", annodf_filename)
  message("   # >>> ", output_path)
  
  annodf <- read.delim(output_path, sep = "\t" , header = TRUE)
  annodf$geneId <- as.character(annodf$geneId)
  annodf$SYMBOL <- as.character(annodf$SYMBOL)
  
  annodf_up <- annodf %>% filter(Fold > 0)
  annodf_down <- annodf %>% filter(Fold < 0)
  
  message("    Number of differential regions : ", nrow(annodf))
  message("       Increased signal : ", nrow(annodf_up), " (including ", nrow(annodf_up %>% filter(Annot == "Promoter")), " regions at promoters)")
  message("       Decreased signal : ", nrow(annodf_down), " (including ", nrow(annodf_down %>% filter(Annot == "Promoter")), " regions at promoters)")
  
  if (geneList == "increase") {
    geneList <- annodf_up %>% filter(Annot == "Promoter") %>% pull(geneId) %>% unique
    message("   Length of geneList (increase) : ", length(geneList))
    return(geneList)
  } else if (geneList == "decrease") {
    geneList <- annodf_down %>% filter(Annot == "Promoter") %>% pull(geneId) %>% unique
    message("   Length of geneList (decrease) : ", length(geneList))
    return(geneList)
  }
}

### Differential analysis
# for (pol in c("POL2", "POL2-ser2")) {
#   for (peak in c("narrow", "broad")) {
#     for (cst in c("shCTRL", "shNIPBL")) {
#       perform_diffbind(pol, cst, effect = "DEX", peak, pval = FALSE)
#     }
#     for (cst in c("CTRL", "DEX")) {
#       perform_diffbind(pol, cst, effect = "shNIPBL", peak, pval = FALSE)
#     }
#   }
# }

### Stats
for (pol in c("POL2", "POL2-ser2")) {
  for (peak in c("narrow", "broad")) {
    for (cst in c("shCTRL", "shNIPBL")) {
      stats_diffbind(pol, cst, effect = "DEX", peak, pval = FALSE)
    }
    for (cst in c("CTRL", "DEX")) {
      stats_diffbind(pol, cst, effect = "shNIPBL", peak, pval = FALSE)
    }
  }
}

# Retrieve geneList and plot RNA expression over time with Reddy's data
shCTRL_DEX_increase <- stats_diffbind(pol = "POL2", cst = "shCTRL", effect = "DEX", peak = "narrow", pval = FALSE, geneList = "increase") # fdr: 6 | pval: 110
shCTRL_DEX_decrease <- stats_diffbind(pol = "POL2", cst = "shCTRL", effect = "DEX", peak = "narrow", pval = FALSE, geneList = "decrease") # fdr: 4 | pval: 33
shNIPBL_DEX_increase <- stats_diffbind(pol = "POL2", cst = "shNIPBL", effect = "DEX", peak = "narrow", pval = FALSE, geneList = "increase") # fdr: 29 | pval: 511
shNIPBL_DEX_decrease <- stats_diffbind(pol = "POL2", cst = "shNIPBL", effect = "DEX", peak = "narrow", pval = FALSE, geneList = "decrease") # fdr: 0 | pval: 7
CTRL_shNIPBL_increase <- stats_diffbind(pol = "POL2", cst = "CTRL", effect = "shNIPBL", peak = "narrow", pval = FALSE, geneList = "increase") # fdr: 0 | pval: 29
CTRL_shNIPBL_decrease <- stats_diffbind(pol = "POL2", cst = "CTRL", effect = "shNIPBL", peak = "narrow", pval = FALSE, geneList = "decrease") # fdr: 2 | pval: 53
DEX_shNIPBL_increase <- stats_diffbind(pol = "POL2", cst = "DEX", effect = "shNIPBL", peak = "narrow", pval = FALSE, geneList = "increase") # fdr: 0 | pval: 250
DEX_shNIPBL_decrease <- stats_diffbind(pol = "POL2", cst = "DEX", effect = "shNIPBL", peak = "narrow", pval = FALSE, geneList = "decrease") # fdr: 0 | pval: 35

draw_time_course_FC(shCTRL_DEX_increase)
draw_time_course_FC(shCTRL_DEX_decrease)
draw_time_course_FC(shNIPBL_DEX_increase)
draw_time_course_FC(shNIPBL_DEX_decrease)
draw_time_course_FC(CTRL_shNIPBL_increase)
draw_time_course_FC(CTRL_shNIPBL_decrease)
draw_time_course_FC(DEX_shNIPBL_increase)
draw_time_course_FC(DEX_shNIPBL_decrease)

geneList_group <- list("shCTRL_DEX_increase" = shCTRL_DEX_increase,
                          "shCTRL_DEX_decrease" = shCTRL_DEX_decrease,
                          "shNIPBL_DEX_increase" = shNIPBL_DEX_increase,
                          "shNIPBL_DEX_decrease" = shNIPBL_DEX_decrease,
                          "CTRL_shNIPBL_increase" = CTRL_shNIPBL_increase,
                          "CTRL_shNIPBL_decrease" = CTRL_shNIPBL_decrease,
                          "DEX_shNIPBL_increase" = DEX_shNIPBL_increase,
                          "DEX_shNIPBL_decrease" = DEX_shNIPBL_decrease)

draw_time_course_pergroup_FC(geneList_group)

# Test Annotation
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
# #
# # report <- dba.report(analyze, bCounts = T, bUsePval = TRUE)
# report <- dba.report(analyze, bCounts = T)
# #
# annodf1 <- annotatePeaks(report, output = "df", tss = 5000); print(nrow(annodf1)); print(sum(annodf1$Annot == "Promoter")); print(sort(annodf1$SYMBOL))
# annodf2 <- annotatePeaks2(report, output = "df", tss = 5000); print(nrow(annodf2)); print(sum(annodf2$Annot == "Promoter")); print(sort(annodf2$SYMBOL))
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
