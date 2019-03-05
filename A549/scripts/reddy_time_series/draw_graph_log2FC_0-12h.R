setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(dplyr)
library(ggplot2)
library(scales)
library(tidyverse)

FC_mat <- read.csv("results/a549_dex_time_points/FC_mat.csv", sep=",", header=TRUE, row.names=1)
colnames(FC_mat) <- c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12)

draw_time_course_FC <- function(genes_list, counts_mat = FC_mat) {
  mat_for_graph <- counts_mat %>% filter(rownames(counts_mat) %in% genes_list)
  # inds <- rownames(counts_mat) %in% genes_list
  # mat_for_graph <- counts_mat[inds,]
  
  inds2 <- genes_list %in% rownames(counts_mat)
  # rownames(mat_for_graph) <- genes_list[inds2]
  
  print (genes_list[!inds2])
  
  # Data
  hours <- as.numeric(colnames(FC_mat))
  
  # select(hours, as.factor(rownames(mat_for_graph))) %>%
  df <- data.frame(hours, t(mat_for_graph)) %>% gather(key = "genes", value = "counts", -hours)
  
  # Calculation of mean
  df_mean <- aggregate(df$counts, by=list(df$hours, df$genes), mean)
  colnames(df_mean) <- c("hours", "genes", "mean")
  
  # Make the plot
  plot <- ggplot(df, aes(x = hours, y = counts)) +
    geom_point(aes(color = genes), size = 1) +
    theme(legend.position='none') +
    geom_line(data=df_mean, aes(x=hours, y=mean, group=genes, color=genes)) +
    scale_x_continuous(name="Hours", labels=c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12), breaks=c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12)) +
    scale_y_continuous(name="log2(FC)", labels=comma)
    
  return (plot)
}

draw_time_course_pergroup_FC <- function(genes_group_list, counts_mat = FC_mat) {
  df <- data.frame(hours = numeric(0),
                   mean = numeric(0),
                   Gene_group = character(0))
  
  hours <- as.numeric(colnames(counts_mat))
  
  for (i in 1:length(genes_group_list)) {
    group_n <- names(genes_group_list)[i]
    group <- genes_group_list[[i]]
    inds <- rownames(counts_mat) %in% group
    mat_for_graph <- counts_mat[inds,]
    means <- colMeans(mat_for_graph)
    mat_for_graph_mean <- data.frame(hours, means)
    df_mean <- aggregate(mat_for_graph_mean$means, by=list(mat_for_graph_mean$hours), mean)
    Gene_group <- rep.int(group_n, nrow(df_mean))
    df_mean <- data.frame(df_mean, Gene_group)
    colnames(df_mean) <- c("hours", "mean", "Gene_group")
    df <- rbind(df, df_mean)
  }
  
  # Make the plot
  plot <- ggplot(df, aes(x = hours, y = mean)) +
    geom_point(aes(color = Gene_group), size = 1) +
    #    theme(legend.position='none') +
    geom_line(data=df, aes(x=hours, y=mean, group=Gene_group, color=Gene_group)) +
    scale_x_continuous(name="Hours", labels=c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12), breaks=c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12)) +
    scale_y_continuous(name="mean(log2(FC))", labels=comma)
  
  return (plot)
}

# geneList1 <- c("ENSG00000128016", "ENSG00000167772")
# geneList2 <- c("ENSG00000095752", "ENSG00000069399")
# gene <- list("test1" = geneList1, "gene2" = geneList2)
# draw_time_course_FC(geneList1)
# draw_time_course_FC(geneList2)
# draw_time_course_pergroup_FC(gene)
