setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/chris/countRead_overtime/countReads.utils.R")
library(knitr)
library(ggplot2)
library(scales)
library(dplyr)

# Get values
GR <- retrieve_sumcount("count_total_GR.txt") %>% assignTime
# GR_WCE <- retrieve_sumcount("count_total_GR_wce.txt") %>% assignTime
# GR_background <- retrieve_sumcount("count_total_GR_background.txt") %>% assignTime
GR_gapbackground <- retrieve_sumcount("count_total_GR_gaps_background.txt") %>% assignTime

EP300 <- retrieve_sumcount("count_total_EP300.txt") %>% assignTime
# EP300_WCE <- retrieve_sumcount("count_total_EP300_wce.txt") %>% assignTime
# EP300_background <- retrieve_sumcount("count_total_EP300_background.txt") %>% assignTime
EP300_gapbackground <- retrieve_sumcount("count_total_EP300_gaps_background.txt") %>% assignTime

CTCF <- retrieve_sumcount("count_total_CTCF.txt") %>% assignTime
CTCF_gapbackground <- retrieve_sumcount("count_total_CTCF_gaps_background.txt") %>% assignTime

RAD21 <- retrieve_sumcount("count_total_RAD21.txt") %>% assignTime
RAD21_gapbackground <- retrieve_sumcount("count_total_RAD21_gaps_background.txt") %>% assignTime

SMC3 <- retrieve_sumcount("count_total_SMC3.txt") %>% assignTime
SMC3_gapbackground <- retrieve_sumcount("count_total_SMC3_gaps_background.txt") %>% assignTime

FOSL2 <- retrieve_sumcount("count_total_FOSL2.txt") %>% assignTime
FOSL2_gapbackground <- retrieve_sumcount("count_total_FOSL2_gaps_background.txt") %>% assignTime

BCL3 <- retrieve_sumcount("count_total_BCL3.txt") %>% assignTime
BCL3_gapbackground <- retrieve_sumcount("count_total_BCL3_gaps_background.txt") %>% assignTime

JUN <- retrieve_sumcount("count_total_JUN.txt") %>% assignTime
JUN_gapbackground <- retrieve_sumcount("count_total_JUN_gaps_background.txt") %>% assignTime

JUNB <- retrieve_sumcount("count_total_JUNB.txt") %>% assignTime
JUNB_gapbackground <- retrieve_sumcount("count_total_JUNB_gaps_background.txt") %>% assignTime

HES2 <- retrieve_sumcount("count_total_HES2.txt") %>% assignTime
HES2_gapbackground <- retrieve_sumcount("count_total_HES2_gaps_background.txt") %>% assignTime

CEBPB <- retrieve_sumcount("count_total_CEBPB.txt") %>% assignTime
CEBPB_gapbackground <- retrieve_sumcount("count_total_CEBPB_gaps_background.txt") %>% assignTime

##### Plot figures
# Raw count values
bigdf_GR_allTargets <- bigdfForPlot(list(GR, EP300, CTCF, RAD21, SMC3, FOSL2, BCL3, JUN, JUNB, HES2, CEBPB),
                               protein = c("GR", "EP300", "CTCF", "RAD21", "SMC3", "FOSL2", "BCL3", "JUN", "JUNB", "HES2", "CEBPB"))
plotReadCount(bigdf_GR_allTargets)

# Raw count values with WCE
bigdf_GR_EP300_withWCE <- bigdfForPlot(list(GR, EP300, GR_WCE, EP300_WCE), protein = c("GR", "EP300", "GR-WCE", "EP300-WCE"))
plotReadCount(bigdf_GR_EP300_withWCE)

# Raw count / background (shift 10000)
bigdf_GR_EP300_background <- bigdfForPlot(list(GR/GR_background, EP300/EP300_background), protein = c("GR", "EP300"))
plotReadCount(bigdf_GR_EP300_background)

# Raw count / gapbackground (all regions that are not binded by protein)
bigdf_GR_allTargets_gapbackground <- bigdfForPlot(list(GR/GR_gapbackground,
                                                       EP300/EP300_gapbackground,
                                                       CTCF/CTCF_gapbackground,
                                                       RAD21/RAD21_gapbackground,
                                                       SMC3/SMC3_gapbackground,
                                                       FOSL2/FOSL2_gapbackground,
                                                       BCL3/BCL3_gapbackground,
                                                       JUN/JUN_gapbackground,
                                                       JUNB/JUNB_gapbackground,
                                                       HES2/HES2_gapbackground,
                                                       CEBPB/CEBPB_gapbackground),
                                                  protein = c("GR", "EP300", "CTCF", "RAD21", "SMC3", "FOSL2", "BCL3", "JUN", "JUNB", "HES2", "CEBPB"))
plotReadCount(bigdf_GR_allTargets_gapbackground)

