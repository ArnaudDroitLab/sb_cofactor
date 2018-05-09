# setwd("/home/chris/Bureau/explore_metagene/")

library(plotly)
library(wesanderson)

###############################################################################
# Generate metagene plots 
###############################################################################

# cofactor_list <- c("BRD4", "CDK9", "NIPBL", "SMC1A", "MED1")
cofactor_list <- c("BRD4")

output_dir <- "output/chip-pipeline-GRCh38/metagene/metagene_cofactor"

colorReplicate = wes_palette(n=2, name="Cavalcanti1")

for (cofactor in cofactor_list) {
	metagene_obj_name <- file.path(output_dir, paste0(cofactor, "_GR30_3utr_GR30_5utr_GR30_distalintergenic_GR30_downstream_GR30_exon_GR30_intron_GR30_promoter", "_metagene_obj.RData"))
	load(metagene_obj_name)
	
	mg_df <- mg$get_data_frame()
	mg_df$Condition <- as.factor(ifelse(grepl("DEX", mg_df$group), "DEX", "EtOH"))
	mg_df$Replicate <- as.factor(ifelse(grepl("rep1", mg_df$group), "Replicate 1", "Replicate 2"))
	mg_df$Cofactor <- as.factor(cofactor)
	mg_df$Annot[grepl("3utr", mg_df$group)] <- "3' UTR"
	mg_df$Annot[grepl("5utr", mg_df$group)] <- "5' UTR"
	mg_df$Annot[grepl("distalintergenic", mg_df$group)] <- "Distal Intergenic"
	mg_df$Annot[grepl("downstream", mg_df$group)] <- "Downstream"
	mg_df$Annot[grepl("exon", mg_df$group)] <- "Exon"
	mg_df$Annot[grepl("intron", mg_df$group)] <- "Intron"
	mg_df$Annot[grepl("promoter", mg_df$group)] <- "Promoter"
	
	title <- paste0(cofactor, "_GRBound_30min_annotated")
	
	ggplot(mg_df, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group = Replicate)) +
		geom_ribbon(aes(fill = Replicate), alpha = 0.3) +
		scale_fill_manual(values=colorReplicate) +
		geom_line(aes(color = Replicate), size = 1) + 
		scale_color_manual(values=colorReplicate) +
		theme_bw(base_size = 20) +
		theme(panel.grid.major = element_line(),
			  panel.grid.minor = element_line(),
			  panel.background = element_blank()) +
		theme(legend.position = "bottom",
			  legend.direction = "vertical") +
		scale_x_continuous(breaks = seq(0, 200, 50),
						   labels = c(-300, -150, 0, 150, 300)) +
		theme(axis.text.x = element_text(size=10),
			  axis.text.y = element_text(size=10)) +
		ggtitle(title) +
		xlab("Position (bp)") +
		ylab("Mean coverage (RPM)") +
		facet_grid(Cofactor*Annot ~ Condition)
	
	output_filename <- file.path(output_dir, paste0(cofactor, "_GRBound_30min_w600_annotated_metagene.pdf"))
	
	ggsave(output_filename, width = 10, height = 24, units = "in")
	message("Saved in ", output_filename)
}
