# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

# Load cofactors differential binding linear annotation
genes_reg_by_cofactors_viaLinear <- readRDS(file = "output/analyses/ecosystem/genes_reg_by_cofactors_viaLinear.rds")

# Load cofactors differential binding 3D annotation
genes_reg_by_cofactors_via3D <- readRDS(file = "output/analyses/ecosystem/genes_reg_by_cofactors_via3D.rds")

##### Load induced and repressed gene categories
deg <- readRDS(file = "output/analyses/deg.rds")
upreg <- deg$gene_list$FC1$upreg
downreg <- deg$gene_list$FC1$downreg

cofactors_cat <- names(genes_reg_by_cofactors_viaLinear)

for (elt in cofactors_cat) {
  message("##### ", elt)
  
  gene_linear <- genes_reg_by_cofactors_viaLinear[[elt]]
  message("   # Number of genes determined with linear annotation : ", length(gene_linear))
  
  gene_3D <- genes_reg_by_cofactors_via3D[[elt]]
  message("   # Number of genes determined with 3D annotation : ", length(gene_3D))
  
  all_genes <- c(gene_linear, gene_3D)
  message("   # Number of total genes : ", length(all_genes))

  all_genes_unique <- unique(all_genes)
  message("   # Number of total unique genes : ", length(all_genes_unique))
  
  # genes_added_by3D <- setdiff(gene_3D, gene_linear)
  # message("   # Number of genes added with 3D annotation : ", length(genes_added_by3D))
  message(" ")
  
  # activated/repressed over 12h
  res <- cbind(all_genes_unique, "upreg" = all_genes_unique %in% upreg, "downreg" = all_genes_unique %in% downreg) %>% as.data.frame
  res$all_genes_unique <- as.character(all_genes_unique)
  res$upreg <- as.logical(res$upreg)
  res$downreg <- as.logical(res$downreg)
  
  nb_upreg <- sum(res$upreg)
  nb_downreg <- sum(res$downreg)
  message("     Associated with ",  nb_upreg, " induced genes > ", round(nb_upreg/nrow(res)*100, 2), " %")
  message("     Associated with ",  nb_downreg, " repressed genes > ", round(nb_downreg/nrow(res)*100, 2), " %")
  message(" ")
}