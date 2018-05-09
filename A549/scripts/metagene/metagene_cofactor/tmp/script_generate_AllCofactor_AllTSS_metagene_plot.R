# setwd("/home/chris/Bureau/explore_metagene/")

library(plotly)
library(wesanderson)

###############################################################################
# Generate metagene plots
###############################################################################

cofactor_list <- c("BRD4", "CDK9", "NIPBL", "SMC1A", "MED1")

bigdf <- data.frame(region = character(),
				   design = character(),
				   bin = integer(),
				   value = double(),
				   strand = character(),
				   qinf = double(),
				   qsup = double(),
				   group = factor(),
				   Condition = factor(),
				   Replicate = factor(),
				   Cofactor = factor())

output_dir <- "output/chip-pipeline-GRCh38/metagene_cofactor"

for (cofactor in cofactor_list) {
	metagene_obj_name <- file.path(output_dir, paste0(cofactor, regions_names, "_metagene_obj.RData"))
	load(metagene_obj_name)
	
	mg_df <- mg$get_data_frame()
	mg_df$Condition <- as.factor(ifelse(grepl("DEX", mg_df$group), "DEX", "EtOH"))
	mg_df$Replicate <- as.factor(ifelse(grepl("rep1", mg_df$group), "Replicate 1", "Replicate 2"))
	mg_df$Cofactor <- as.factor(cofactor)
	
	bigdf <- rbind(bigdf, mg_df)
}

colorReplicate = wes_palette(n=2, name="Cavalcanti1")

title <- paste0("AllTSS regions")

ggplot(bigdf, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group=Replicate)) +
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
	scale_x_continuous(breaks = seq(0, 100, 25),
					   labels = c(c(-1.5, -0.75), "TSS", c(0.75, 1.5))) +
	theme(axis.text.x = element_text(size=10),
		  axis.text.y = element_text(size=10)) +
	ggtitle(title) +
	xlab("Position from TSS (kb)") +
	ylab("Mean coverage (RPM)") +
	facet_grid(Cofactor ~ Condition)

output_filename <- file.path(output_dir, "AllCofactor_AllTSS_w3000_metagene.pdf")

ggsave(output_filename, width = 7, height = 14, units = "in")