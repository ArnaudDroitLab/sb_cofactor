setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")

### regions to assess
gr_regions <- load_reddy_gr_binding_consensus()

all_gr_regions <- GRanges()
for (name in names(gr_regions)) {
  all_gr_regions <- c(all_gr_regions, gr_regions[[name]])
  message(length(gr_regions[[name]]))
}

all_gr_regions_reduced <- reduce(all_gr_regions)
summary(width(all_gr_regions_reduced))