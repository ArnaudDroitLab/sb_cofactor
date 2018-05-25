library(plotly)
# library(wesanderson)

###############################################################################
# Generate metagene plots 
###############################################################################

binding_diff_dir <- "output/chip-pipeline-GRCh38/binding_diff"
customColors <- c("#D8B70A", "#02401B", "#94A187", "#C5AFA0", "#E9BCB7")

# cofactor_list <- c("BRD4", "CDK9", "NIPBL", "SMC1A", "MED1")
for (cofactor in cofactor_list) {
	output_dir <- file.path(binding_diff_dir, paste0("A549_", cofactor), "output_filters")
	
	metagene_WCE_obj_name <- file.path(output_dir, paste0("WCE_", cofactor, "_BindUpRegions_BindDownRegions_metagene_obj.RData"))
	load(metagene_WCE_obj_name)
	wce_df <- mg$get_data_frame()

	metagene_obj_name <- file.path(output_dir, paste0(cofactor, "_BindUpRegions_BindDownRegions_metagene_obj.RData"))
	load(metagene_obj_name)
	mg_df <- mg$get_data_frame()
	
	bigdf <- rbind(wce_df, mg_df)
	
	bigdf$Condition <- as.factor(ifelse(grepl("DEX", bigdf$group), "DEX", "EtOH"))
	bigdf$Cofactor <- as.factor(cofactor)
	
	bigdf$Replicate <- ifelse(grepl("rep1", bigdf$group), paste0(cofactor, " Replicate 1"), paste0(cofactor, " Replicate 2"))
	bigdf$Replicate[grepl("WCE_rep1", bigdf$group)] <- "WCE Replicate 1"
	bigdf$Replicate[grepl("WCE_rep2", bigdf$group)] <- "WCE Replicate 2"
	bigdf$Replicate[grepl("WCE_rep3", bigdf$group)] <- "WCE Replicate 3"
	bigdf$Replicate <- as.factor(bigdf$Replicate)
	
	bigdf$Regions <- as.factor(ifelse(grepl("BindUp", bigdf$group), "BindUpRegions", "BindDownRegions"))

	title <- paste0(cofactor, "_Differential_Binding")
	
	ggplot(bigdf, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group = Replicate)) +
		geom_ribbon(aes(fill = Replicate), alpha = 0.3) +
		scale_fill_manual(values=customColors) +
		geom_line(aes(color = Replicate), size = 1) + 
		scale_color_manual(values=customColors) +
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
		facet_grid(Regions ~ Condition)
	
	output_filename <- file.path(output_dir, paste0(cofactor, "_Differential_Binding_metagene.pdf"))
	
	ggsave(output_filename, width = 14, height = 9, units = "in")
	message("Saved in ", output_filename)
}
