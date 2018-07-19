#setwd("C:/Dev/Projects/sb_cofactor/A549/")
source("scripts/load_reddy.R")

library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(cowplot)

intersect_obj = load_reddy_gr_binding_intersect()
# Subset of promoters
all_combo = apply(intersect_obj$Matrix[grepl("Promoter", intersect_obj$Regions$annotation),]>0, 1, paste, collapse=";")

combo_frequencies = sort(table(all_combo), decreasing=TRUE)
top_combo_frequencies = combo_frequencies[1:50]
top_matrix = matrix(unlist(lapply(as.list(names(top_combo_frequencies)), function(x) { as.logical(unlist(strsplit(x, ";")))})), nrow=50, byrow=T)
colnames(top_matrix) = intersect_obj$Names

top_df = as.data.frame(top_matrix)
top_df$Frequency=as.numeric(top_combo_frequencies)
top_df$Order = 1:50

binding_df = reshape2::melt(top_df, id.vars=c("Order", "Frequency"), variable.name="Time")
binding_plot = ggplot(binding_df, aes(x=Time, y=-Order, fill=value)) + 
    geom_tile(color="black") + 
    scale_fill_manual(values=c("TRUE"="red", "FALSE"="white")) + 
    theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    guides(fill=FALSE)

frequency_plot = ggplot(top_df) + 
    geom_bar(mapping=aes(x=-Order, y=Frequency), stat="identity", fill="red", color="black") + 
    coord_flip() +
    theme(axis.line.y=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank())

pdf("GR binding frequency promoters.pdf", width=14, height=7)
plot_grid(binding_plot, frequency_plot, align='h')
dev.off()
