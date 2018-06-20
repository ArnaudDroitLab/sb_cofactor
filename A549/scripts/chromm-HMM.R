library(ef.utils)
library(ggplot2)

source("scripts/load_reddy.R")

# Load all binding sites for all cofactors.
lab_chip = load_cofactor_binding()

# Get chromatin states
chromatin_states = rtracklayer::import("input/E114_18_core_K27ac_mnemonics_hg38_liftOver.bed", format="BED")

# Get H3K27ac marks.

# Load all ENCODE chips.
h3k27ac_samples = ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format="bed",
                                                    assay="ChIP-seq", target="H3K27ac",
                                                    lab="Tim Reddy, Duke", status="released")

# Keeping the number of replicates at three makes ChIPs more compaable. The latest Reddy
# batch all coems in triplicates. So we'll remove older batches when there are more
h3k27ac_samples = h3k27ac_samples %>% dplyr::filter(!(date_released=="2016-07-21" &is.na(treatment)))
h3k27ac_regions = summarize_GGR_chip(h3k27ac_samples)

# Get induced and depleted regions
cofactor_names <- c("BRD4", "CDK9", "NIPBL", "SMC1A", "MED1")
cofactor_lists <- list(Up=list(), Down=list(), DEX=list(), CTRL=list())
for(cofactor in cofactor_names) {
    binding_diff_dir <- file.path("output/chip-pipeline-GRCh38/binding_diff", paste0("A549_", cofactor), "output_filters")
    BindUpRegions_path <- file.path(binding_diff_dir, paste0("A549_DEX_", cofactor, "_rep1_peaks.narrowPeak_M_above_1.0_biased_peaks.bed"))
    BindDownRegions_path <- file.path(binding_diff_dir, paste0("A549_CTRL_", cofactor, "_rep1_peaks.narrowPeak_M_below_-1.0_biased_peaks.bed"))
    
    # CstRegions <- rtracklayer::import(CstRegions_path)
    cofactor_lists[["Up"]][[cofactor]] = rtracklayer::import(BindUpRegions_path)
    cofactor_lists[["Down"]][[cofactor]] = rtracklayer::import(BindDownRegions_path)
    cofactor_lists[["DEX"]][[cofactor]] = lab_chip[["DEX"]][[cofactor]]
    cofactor_lists[["CTRL"]][[cofactor]] = lab_chip[["CTRL"]][[cofactor]]
}

# Perform enrichment on all groups of cofactor binding sites.
output_dir = "output/analyses/differential_binding_sites"
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
enrichment_results = list()
for(region_list in names(cofactor_lists)) {
    region_dir = file.path(output_dir, paste0("Chromatin states enrichments ", region_list))
    enrichment_results[[region_list]] = multiple_region_enrichment(cofactor_lists[[region_list]],
                                                                   chromatin_states, 
                                                                   file.prefix=region_dir)
}

#### Produce graph of induced/depleted region chromatin state
#### enrichment compared to EtOH condition.

# Calculate enrichments.
up_enrichment = log2(enrichment_results$Up$Data$Proportion[-1,] / enrichment_results$CTRL$Data$Proportion[-1,])
down_enrichment = log2(enrichment_results$Down$Data$Proportion[-1,] / enrichment_results$CTRL$Data$Proportion[-1,])                           
                        
# Get everything into a single long-form data-frame.                        
up_df = reshape2::melt(up_enrichment)
up_df$Direction="Induced"
                           
down_df = reshape2::melt(down_enrichment)
down_df$Direction="Depleted"

all_df = rbind(up_df, down_df)
colnames(all_df) = c("Factor", "ChromatinState", "Enrichment", "Direction")
all_df$ChromatinState = gsub("^([1-9])_", "0\\1_", all_df$ChromatinState)
all_df$ChromatinState = factor(all_df$ChromatinState, levels=sort(unique(all_df$ChromatinState)))
all_df$Direction = factor(all_df$Direction, levels=c("Induced", "Depleted"))

# Replace infinite values with values just outside of the scale.
all_df$Enrichment[all_df$Enrichment==Inf] = 10
all_df$Enrichment[all_df$Enrichment==-Inf] = -10 

# Plot it all.
ggplot(all_df, aes(x=ChromatinState, y=Enrichment, fill=Direction, shape=Direction)) +
    geom_point(color="black", size=3) +
    facet_grid(Factor~.) +
    scale_shape_manual(values=c(Induced=24, Depleted=25)) +
    scale_fill_manual(values=c(Induced="#239b56", Depleted="#c0392b")) +
    labs(y="Enrichment compared to all cofactor binding sites in EtOH") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file="output/analyses/differential_binding_sites/Induced and depleted enrichment vs all EtOH.pdf", width=14, height=7, units="in")    


# Calculate enrichments.
up_enrichment = enrichment_results$Up$Data$Proportion[-1,] - enrichment_results$CTRL$Data$Proportion[-1,]
down_enrichment = enrichment_results$Down$Data$Proportion[-1,] - enrichment_results$CTRL$Data$Proportion[-1,]
                        
# Get everything into a single long-form data-frame.                        
up_df = reshape2::melt(up_enrichment)
up_df$Direction="Induced"
                           
down_df = reshape2::melt(down_enrichment)
down_df$Direction="Depleted"

all_df = rbind(up_df, down_df)
colnames(all_df) = c("Factor", "ChromatinState", "Enrichment", "Direction")
all_df$ChromatinState = gsub("^([1-9])_", "0\\1_", all_df$ChromatinState)
all_df$ChromatinState = factor(all_df$ChromatinState, levels=sort(unique(all_df$ChromatinState)))
all_df$Direction = factor(all_df$Direction, levels=c("Induced", "Depleted"))

# Replace infinite values with values just outside of the scale.
all_df$Enrichment[all_df$Enrichment==Inf] = 10
all_df$Enrichment[all_df$Enrichment==-Inf] = -10 

# Plot it all.
ggplot(all_df, aes(x=ChromatinState, y=Enrichment, fill=Direction, shape=Direction)) +
    geom_point(color="black", size=3) +
    geom_hline(yintercept=0) +
    facet_grid(Factor~.) +
    scale_shape_manual(values=c(Induced=24, Depleted=25)) +
    scale_fill_manual(values=c(Induced="#239b56", Depleted="#c0392b")) +
    labs(y="Difference between proportion in all EtOH binding sites and induced/depleted regions") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file="output/analyses/differential_binding_sites/Induced and depleted enrichment difference with EtOH.pdf", width=14, height=7, units="in")    



#### Produce graph of induced/depleted region chromatin state proportions.
up_proportion = reshape2::melt(enrichment_results$Up$Data$Proportion[-1,])
up_proportion$Group = "Induced"
down_proportion = reshape2::melt(enrichment_results$Down$Data$Proportion[-1,])
down_proportion$Group = "Depleted"
ctrl_proportion = reshape2::melt(enrichment_results$CTRL$Data$Proportion[-1,])
ctrl_proportion$Group = "All_EtOH"

all_df = rbind(up_proportion, down_proportion, ctrl_proportion)
colnames(all_df) = c("Factor", "ChromatinState", "Proportion", "Group")
all_df$ChromatinState = gsub("^([1-9])_", "0\\1_", all_df$ChromatinState)
all_df$ChromatinState = factor(all_df$ChromatinState, levels=sort(unique(all_df$ChromatinState)))
all_df$Group = factor(all_df$Group, levels=c("Induced", "All_EtOH", "Depleted"))

ggplot(all_df, aes(x=ChromatinState, y=Proportion, fill=Group, shape=Group)) +
    geom_point(color="black", size=3) +
    facet_grid(Factor~.) +
    scale_shape_manual(values=c(Induced=24, Depleted=25, All_EtOH=21)) +
    scale_fill_manual(values=c(Induced="#239b56", Depleted="#c0392b", All_EtOH="#3399ff")) +
    labs(y="Proportion of binding sites in chromain state") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
ggsave(file="output/analyses/differential_binding_sites/Induced, depleted and all EtOH proportions.pdf", width=14, height=7, units="in")    
     
#### Perform comparisons with H3K27ac regions to determine if there is
#### an increase.

for(region_list in names(cofactor_lists)) {
    region_dir = file.path(output_dir, paste0("Chromatin states enrichments ", region_list))
    enrichment_results[[region_list]] = multiple_region_enrichment(cofactor_lists[[region_list]],
                                                                   chromatin_states, 
                                                                   file.prefix=region_dir)
}

     
# Overlap with H3K27ac regions.                  
                  