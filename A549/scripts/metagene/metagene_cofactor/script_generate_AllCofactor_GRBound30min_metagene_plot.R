# setwd("/home/chris/Bureau/explore_metagene/")

library(plotly)
library(wesanderson)

###############################################################################
# Generate metagene plots 
###############################################################################

cofactor_list <- c("BRD4", "CDK9", "NIPBL", "SMC1A", "MED1")

output_dir <- "output/chip-pipeline-GRCh38/metagene/metagene_cofactor"

colorReplicate = wes_palette(n=3, name="Cavalcanti1")

for (cofactor in cofactor_list) {
	metagene_obj_name <- file.path(output_dir, paste0(cofactor, "_UpRegulatedBoundTSS_UpRegulatedUnboundTSS_DownRegulatedBoundTSS_DownRegulatedUnboundTSS_UpRegulatedAllGRTSS_DownRegulatedAllGRTSS", "_metagene_obj.RData"))
	load(metagene_obj_name)
	
	mg_df <- mg$get_data_frame()
	mg_df$Condition <- as.factor(ifelse(grepl("DEX", mg_df$group), "DEX", "EtOH"))
	mg_df$Replicate <- as.factor(ifelse(grepl("rep1", mg_df$group), "Replicate 1", "Replicate 2"))
	mg_df$Cofactor <- as.factor(cofactor)
	mg_df$Direction <- as.factor(ifelse(grepl("UpRegulated", mg_df$group), "UpRegulatedTSS", "DownRegulatedTSS"))
	mg_df$GR <- as.factor(ifelse(grepl("Bound", mg_df$group), "GR-Bound",
								 ifelse(grepl("Unbound", mg_df$group), "GR-Unbound", "GR-All")))
	
	title <- paste0(cofactor, "_DEX_2hours_GRBound_30min")
	
	ggplot(mg_df, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group = GR)) +
		geom_ribbon(aes(fill = GR), alpha = 0.3) +
		scale_fill_manual(values=colorReplicate) +
		geom_line(aes(color = GR), size = 1) + 
		scale_color_manual(values=colorReplicate) +
		theme_bw(base_size = 20) +
		theme(panel.grid.major = element_line(),
			  panel.grid.minor = element_line(),
			  panel.background = element_blank()) +
		theme(legend.position = "bottom",
			  legend.direction = "vertical") +
		scale_x_continuous(breaks = seq(0, 200, 50),
						   labels = c(c(-0.5, -0.25), "TSS", c(0.25, 0.5))) +
		theme(axis.text.x = element_text(size=10),
			  axis.text.y = element_text(size=10)) +
		ggtitle(title) +
		xlab("Position from TSS (kb)") +
		ylab("Mean coverage (RPM)") +
		facet_grid(Cofactor*Replicate ~ Direction*Condition)
	
	output_filename <- file.path(output_dir, paste0(cofactor, "_DEX_2hours_GRBound_30min_w1000_metagene_withAllGR.pdf"))
	
	ggsave(output_filename, width = 14, height = 9, units = "in")
	message("Saved in ", output_filename)
}
