# get OS:http://conjugateprior.org/2015/06/identifying-the-os-from-r/
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

output_dir <- "output/analyses/NBC_balance"

source("scripts/ckn_utils.R")
library(ef.utils)
library(GenomicRanges)
library(eulerr)

##########
cofactors_peaks <- load_cofactor_stdchr_peaks(cofactors = c("NIPBL", "BRD4", "CDK9"))

NIPBL_DEX <- cofactors_peaks[["NIPBL_DEX"]]; print(length(NIPBL_DEX)) # 2733
BRD4_DEX <- cofactors_peaks[["BRD4_DEX"]]; print(length(BRD4_DEX)) # 21225
CDK9_DEX <- cofactors_peaks[["CDK9_DEX"]]; print(length(CDK9_DEX)) # 6061

NBC_DEX <- GenomicRangesList("NIPBL_DEX" = NIPBL_DEX, "BRD4_DEX" = BRD4_DEX, "CDK9_DEX" = CDK9_DEX)
plotVenn(NBC_DEX)

##########
gr_regions <- load_reddy_gr_binding_consensus()
names(gr_regions)

##########
gr_5m <- gr_regions[["5 minute"]]; print(length(gr_5m)) # 
NIPBL_DEX_gr_5m <- subsetByOverlaps(NIPBL_DEX, gr_5m); print(length(NIPBL_DEX_gr_5m)) # 
BRD4_DEX_gr_5m <- subsetByOverlaps(BRD4_DEX, gr_5m); print(length(BRD4_DEX_gr_5m)) # 
CDK9_DEX_gr_5m <- subsetByOverlaps(CDK9_DEX, gr_5m); print(length(CDK9_DEX_gr_5m)) # 
NBC_DEX_gr_5m <- GenomicRangesList("NIPBL_DEX" = NIPBL_DEX_gr_5m,
                                   "BRD4_DEX" = BRD4_DEX_gr_5m,
                                   "CDK9_DEX" = CDK9_DEX_gr_5m)
plotVenn(NBC_DEX_gr_5m)

##########
gr_10m <- gr_regions[["10 minute"]]; print(length(gr_10m)) # 
NIPBL_DEX_gr_10m <- subsetByOverlaps(NIPBL_DEX, gr_10m); print(length(NIPBL_DEX_gr_10m)) # 
BRD4_DEX_gr_10m <- subsetByOverlaps(BRD4_DEX, gr_10m); print(length(BRD4_DEX_gr_10m)) # 
CDK9_DEX_gr_10m <- subsetByOverlaps(CDK9_DEX, gr_10m); print(length(CDK9_DEX_gr_10m)) # 
NBC_DEX_gr_10m <- GenomicRangesList("NIPBL_DEX" = NIPBL_DEX_gr_10m,
                                    "BRD4_DEX" = BRD4_DEX_gr_10m,
                                    "CDK9_DEX" = CDK9_DEX_gr_10m)
plotVenn(NBC_DEX_gr_10m)

##########
gr_15m <- gr_regions[["15 minute"]]; print(length(gr_15m)) # 
NIPBL_DEX_gr_15m <- subsetByOverlaps(NIPBL_DEX, gr_15m); print(length(NIPBL_DEX_gr_15m)) # 
BRD4_DEX_gr_15m <- subsetByOverlaps(BRD4_DEX, gr_15m); print(length(BRD4_DEX_gr_15m)) # 
CDK9_DEX_gr_15m <- subsetByOverlaps(CDK9_DEX, gr_15m); print(length(CDK9_DEX_gr_15m)) # 
NBC_DEX_gr_15m <- GenomicRangesList("NIPBL_DEX" = NIPBL_DEX_gr_15m,
                                    "BRD4_DEX" = BRD4_DEX_gr_15m,
                                    "CDK9_DEX" = CDK9_DEX_gr_15m)
plotVenn(NBC_DEX_gr_15m)

##########
gr_20m <- gr_regions[["20 minute"]]; print(length(gr_20m)) # 
NIPBL_DEX_gr_20m <- subsetByOverlaps(NIPBL_DEX, gr_20m); print(length(NIPBL_DEX_gr_20m)) # 
BRD4_DEX_gr_20m <- subsetByOverlaps(BRD4_DEX, gr_20m); print(length(BRD4_DEX_gr_20m)) # 
CDK9_DEX_gr_20m <- subsetByOverlaps(CDK9_DEX, gr_20m); print(length(CDK9_DEX_gr_20m)) # 
NBC_DEX_gr_20m <- GenomicRangesList("NIPBL_DEX" = NIPBL_DEX_gr_20m,
                                    "BRD4_DEX" = BRD4_DEX_gr_20m,
                                    "CDK9_DEX" = CDK9_DEX_gr_20m)
plotVenn(NBC_DEX_gr_20m)

##########
gr_25m <- gr_regions[["25 minute"]]; print(length(gr_25m)) # 
NIPBL_DEX_gr_25m <- subsetByOverlaps(NIPBL_DEX, gr_25m); print(length(NIPBL_DEX_gr_25m)) # 
BRD4_DEX_gr_25m <- subsetByOverlaps(BRD4_DEX, gr_25m); print(length(BRD4_DEX_gr_25m)) # 
CDK9_DEX_gr_25m <- subsetByOverlaps(CDK9_DEX, gr_25m); print(length(CDK9_DEX_gr_25m)) # 
NBC_DEX_gr_25m <- GenomicRangesList("NIPBL_DEX" = NIPBL_DEX_gr_25m,
                                    "BRD4_DEX" = BRD4_DEX_gr_25m,
                                    "CDK9_DEX" = CDK9_DEX_gr_25m)
plotVenn(NBC_DEX_gr_25m)

##########
gr_30m <- gr_regions[["30 minute"]]; print(length(gr_30m)) # 
NIPBL_DEX_gr_30m <- subsetByOverlaps(NIPBL_DEX, gr_30m); print(length(NIPBL_DEX_gr_30m)) # 
BRD4_DEX_gr_30m <- subsetByOverlaps(BRD4_DEX, gr_30m); print(length(BRD4_DEX_gr_30m)) # 
CDK9_DEX_gr_30m <- subsetByOverlaps(CDK9_DEX, gr_30m); print(length(CDK9_DEX_gr_30m)) # 
NBC_DEX_gr_30m <- GenomicRangesList("NIPBL_DEX" = NIPBL_DEX_gr_30m,
                                    "BRD4_DEX" = BRD4_DEX_gr_30m,
                                    "CDK9_DEX" = CDK9_DEX_gr_30m)
plotVenn(NBC_DEX_gr_30m)

##########
gr_1h <- gr_regions[["1 hour"]]; print(length(gr_1h)) # 8116
NIPBL_DEX_gr_1h <- subsetByOverlaps(NIPBL_DEX, gr_1h); print(length(NIPBL_DEX_gr_1h)) # 2244
BRD4_DEX_gr_1h <- subsetByOverlaps(BRD4_DEX, gr_1h); print(length(BRD4_DEX_gr_1h)) # 6442
CDK9_DEX_gr_1h <- subsetByOverlaps(CDK9_DEX, gr_1h); print(length(CDK9_DEX_gr_1h)) # 3266
NBC_DEX_gr_1h <- GenomicRangesList("NIPBL_DEX" = NIPBL_DEX_gr_1h,
                             "BRD4_DEX" = BRD4_DEX_gr_1h,
                             "CDK9_DEX" = CDK9_DEX_gr_1h)
plotVenn(NBC_DEX_gr_1h)

##########
gr_2h <- gr_regions[["2 hour"]]; print(length(gr_2h)) # 7077
NIPBL_DEX_gr_2h <- subsetByOverlaps(NIPBL_DEX, gr_2h); print(length(NIPBL_DEX_gr_2h)) # 2126
BRD4_DEX_gr_2h <- subsetByOverlaps(BRD4_DEX, gr_2h); print(length(BRD4_DEX_gr_2)) # 5779
CDK9_DEX_gr_2h <- subsetByOverlaps(CDK9_DEX, gr_2h); print(length(CDK9_DEX_gr_2h)) # 3053
NBC_DEX_gr_2h <- GenomicRangesList("NIPBL_DEX" = NIPBL_DEX_gr_2h,
                                   "BRD4_DEX" = BRD4_DEX_gr_2h,
                                   "CDK9_DEX" = CDK9_DEX_gr_2h)
plotVenn(NBC_DEX_gr_2h)

##########
plot_grid(plotVenn(NBC_DEX_gr_5m, labels = FALSE), plotVenn(NBC_DEX_gr_10m, labels = FALSE),
          plotVenn(NBC_DEX_gr_15m, labels = FALSE), plotVenn(NBC_DEX_gr_20m, labels = FALSE),
          plotVenn(NBC_DEX_gr_25m, labels = FALSE), plotVenn(NBC_DEX_gr_30m, labels = FALSE),
          plotVenn(NBC_DEX_gr_1h, labels = FALSE), plotVenn(NBC_DEX_gr_2h, labels = FALSE),
          labels = c("5min", "10min", "15min", "20min",
                     "25min", "30min", "1h", "2h"),
          ncol = 4)
