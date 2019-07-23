# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")


for (i in 1:(ltp-1)) {
  for (j in (i+1):ltp) {
    for (TF in c(TRUE)) {
      tp1 <- timepoint[i]
      tp2 <- timepoint[j]
      report <- open_diffBind("GR", tp1, tp2, pval = TF, output_dir = "output/analyses/GR_diffbind")
    }
  }
}