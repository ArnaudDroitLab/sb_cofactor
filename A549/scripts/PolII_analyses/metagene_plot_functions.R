###############################################################################
# Define functions for metagene plotting.
###############################################################################

# Function for subsetting the metagene data-frame and associating
# sample conditions to profiles.
subset_metagene_df <- function(region_name, antibody=NA) {
    metagene.obj = metagenes[[region_name]]
    
    sh_subset = subset(metagene.obj$get_data_frame(), grepl("sh(C..L|NIPBL).._", group))
    sh_subset$group = gsub("_regions", "", sh_subset$group)
    sh_subset$group = gsub("CRTL", "CTRL", sh_subset$group)
    sh_subset$Condition = ifelse(grepl("Dex", sh_subset$group), "Dex", "EtOH")
    sh_subset$Antibody = ifelse(grepl("POL2.ser2", sh_subset$group), "PolII-ser2", "PolII")
    sh_subset$sh_name = gsub("POL.*(sh.*)_.*", "\\1", sh_subset$group)
    sh_subset$Direction = ifelse(grepl("Regulated", region_name), gsub("(.*Regulated).*", "\\1", region_name), NA)
    sh_subset$Bound = ifelse(grepl("ound", region_name), gsub(".*((Unb|B)ound).*", "\\1", region_name), NA)
    
    if(!is.na(antibody)) {
        sh_subset = sh_subset[sh_subset$Antibody==antibody,]
    }
    
    return(sh_subset)
}

determine_bp_per_bin <- function(region_name, bin_number=200) {
    region_widths = width(metagenes[[region_name]]$get_regions()[[1]])
    if(all(region_widths == mean(region_widths))) {
        # All regions are the same length.
        return(region_widths[1] / bin_number)
    } else {
        return(NA)
    }
}

# Function for plotting the sh effect on a group of regions.
plot_sh_effect <- function(region_name) {
    color_palette = c(shCTRL.1="#40C4FF", shCTRL.2="#00B0FF", shNIPBL.3="#FF8A80", shNIPBL.5="#FF5252")
    do_meta_ggplot_single(region_name, "sh_name", "Condition", color_palette, "sh")
}

# Function for plotting the dex effect on a group of region.
plot_dex_effect <- function(region_name) {
    do_meta_ggplot_single(region_name, "Condition", "sh_name", c(EtOH="#40C4FF", Dex="#FF8A80"), "Dex")
}

do_meta_ggplot_single <- function(region_name, group_var, facet_var, color_palette, file_label, antibody=NA) {
    sh_subset = subset_metagene_df(region_name, antibody)
    bp_per_bin = determine_bp_per_bin(region_name) 
    do_meta_ggplot_generic(sh_subset, region_name, group_var, facet_var, color_palette, file_label, FacetVar~Antibody, bp_per_bin=bp_per_bin)
}

do_meta_ggplot_double <- function(region_name1, region_name2, label1, label2, facet_var, color_palette, file_label, facet_formula, group_label, antibody=NA) {
    subset_1 = subset_metagene_df(region_name1, antibody)
    subset_2 = subset_metagene_df(region_name2, antibody)
    subset_1$RegionName = label1
    subset_2$RegionName = label2
    
    sh_subset = rbind(subset_1, subset_2)
    bp_per_bin = determine_bp_per_bin(region_name1)    
    do_meta_ggplot_generic(sh_subset, group_label, "RegionName", facet_var, color_palette, file_label, facet_formula, bp_per_bin=bp_per_bin)
}

do_meta_ggplot_any <- function(name_list, facet_var, color_palette, file_label, facet_formula, group_label, antibody="both", group_var="RegionName") {
    all_subset = data.frame()
    for(name_index in 1:length(name_list)) {
        this_subset = subset_metagene_df(name_list[[name_index]], antibody)
        this_subset$RegionName = names(name_list)[name_index]
        all_subset = rbind(all_subset, this_subset)
    }

    bp_per_bin = determine_bp_per_bin(name_list[[1]])
    do_meta_ggplot_generic(all_subset, group_label, group_var, facet_var, color_palette, file_label, facet_formula, bp_per_bin=bp_per_bin)
}


do_meta_ggplot_generic <- function(sh_subset, region_label, group_var, facet_var, color_palette, file_label, facet_formula, bp_per_bin=4) {
    sh_subset$GroupVar = sh_subset[[group_var]]
    sh_subset$FacetVar = sh_subset[[facet_var]]
    
    if(grepl("TSS", region_label)) {
        sh_subset$bin = (sh_subset$bin - (max(sh_subset$bin) / 2)) * bp_per_bin
        xlabel = "Distance from TSS (bp)"
    } else {
        xlabel = paste0("Relative distance from TSS (0) to TES (", max(sh_subset$bin), ")")
    }
    
    title_str = paste0("Effect of ", file_label, " on ", region_label)
    ggplot(sh_subset, aes(x=bin, y=value, ymin=qinf, ymax=qsup, group=GroupVar, color=GroupVar, fill=GroupVar)) +
        geom_line() +
        geom_ribbon(alpha = 0.1) +
        facet_grid(facet_formula) +
        scale_color_manual(name=group_var, values=color_palette) +
        scale_fill_manual(name=group_var, values=color_palette) +
        ylab("Mean coverage (RPM)") +
        xlab(xlabel) +
        ggtitle(title_str) +
        theme(plot.title = element_text(hjust = 0.5))
    
    file_name = paste0("Metagene - ", title_str, ".pdf")
    ggsave(file.path(output_dir, file_name))
}
