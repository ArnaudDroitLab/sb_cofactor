setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(plotly)

###############################################################################
# Generate metagene plots 
###############################################################################

target_list <- c("POL2-ser2", "POL2", "WCE")
sh_list <-  c("shCTRL-1", "shCTRL-2", "shNIPBL-3", "shNIPBL-5")

output_dir <- "output/chip-pipeline-PolII-GRCh38/metagene_polII"
customColors <- c("#D8B70A", "#02401B", "#94A187", "#C5AFA0", "#E9BCB7")
#colorReplicate = wes_palette(n=5, name="Cavalcanti1")

bigdf <- data.frame()

for (target in target_list) {
	for (sh in sh_list) {
		metagene_obj_name <- file.path(output_dir, paste0(target, "_", sh, "_GR_Regions_30m_metagene_obj.RData"))
		load(metagene_obj_name)
		
		mg_df <- mg$get_data_frame()
		bigdf <- rbind(bigdf, mg_df)
	}
}


bigdf$Target <- as.factor(ifelse(grepl("POL2-ser2", bigdf$group), "POL2-ser2",
								 ifelse(grepl("POL2", bigdf$group), "POL2", "WCE")))
bigdf$Condition <- as.factor(ifelse(grepl("Dex", bigdf$group), "DEX", "EtOH"))
bigdf$Sh <- as.factor(ifelse(grepl("shCTRL-1", bigdf$group), "shCTRL-1",
							 ifelse(grepl("shCTRL-2", bigdf$group), "shCTRL-2",
							 	   ifelse(grepl("shNIPBL-3", bigdf$group), "shNIPBL-3", "shNIPBL-5"))))

title <- paste0("POL2_POL2-ser2_WCE_GRBound_30min")

ggplot(bigdf, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group = Condition)) +
	geom_ribbon(aes(fill = Condition), alpha = 0.3) +
	scale_fill_manual(values=customColors) +
	geom_line(aes(color = Condition), size = 1) + 
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
	facet_grid(Sh ~ Target)

output_filename <- file.path(output_dir, paste0("POL2_POL2-ser2_WCE_GRBound_30min_w600_metagene.pdf"))

ggsave(output_filename, width = 7, height = 10, units = "in")
message("Saved in ", output_filename)
