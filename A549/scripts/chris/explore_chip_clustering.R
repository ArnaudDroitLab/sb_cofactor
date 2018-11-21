setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")

# Load all ENCODE chips.
all_chip = ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format="bed", assay="ChIP-seq")
all_chip = all_chip %>% dplyr::filter(status!="revoked")

# The most reliable ChIPs are those from Tim Reddy'S lab: 
# they're the most recent and have the most time points. So we'll
# process those first.
all_reddy = all_chip %>% dplyr::filter(lab=="Tim Reddy, Duke")

# Keeping the number of replicates at three makes ChIPs more compaable. The latest Reddy
# batch all coems in triplicates. So we'll remove older batches when there are more
# than 3 replicates.
# Remove a couple of extraneous JUN replicates
all_reddy = all_reddy %>% dplyr::filter(!(target=="JUN" & date_released=="2016-12-20" &is.na(treatment)))

# Do the same thing for some H3K27ac replicates.
all_reddy = all_reddy %>% dplyr::filter(!(target=="H3K27ac" & date_released=="2016-07-21" &is.na(treatment)))
all_targets = unique(all_reddy$target)

all_chip_regions = list()
for(target_name in all_targets) {
    all_chip_regions[[target_name]] = summarize_GGR_chip(all_reddy %>% dplyr::filter(target==target_name),
                                                         diagnostic_dir="output/analyses/cofactors")
}

# Load our own ChIP.
lab_chip = load_cofactor_binding()
for(target in names(lab_chip[["DEX"]])) {
    all_targets = c(all_targets, target)
    all_chip_regions[[target]] = GRangesList("1 hour"   = lab_chip[["DEX"]][[target]],
                                             "0 minute" = lab_chip[["CTRL"]][[target]])
}

# Now move on to other ChIP assays. Start by those with results in 
# GRCh38.
not_reddy_GRCh38 = all_chip %>% dplyr::filter(lab!="Tim Reddy, Duke", assembly=="GRCh38", !(target %in% names(all_chip_regions)))

# Remove combined peak calls, and keep only those originating from a single replicate.
not_reddy_GRCh38 = not_reddy_GRCh38 %>%
                      dplyr::filter(!grepl(";", biological_replicates)) %>%
                      dplyr::filter(!grepl("pseudoreplicated", output_type))

# REST and SIN3A have two separate datasets. Drop one of each.
not_reddy_GRCh38 = not_reddy_GRCh38 %>%          
                      dplyr::filter(accession!="ENCSR892DRK") %>%
                      dplyr::filter(accession!="ENCSR513XQX")

# Load regions.
not_reddy_targets = unique(not_reddy_GRCh38$target)
for(target_name in not_reddy_targets) {
    all_chip_regions[[target_name]] = summarize_GGR_chip(not_reddy_GRCh38 %>% dplyr::filter(target==target_name),
                                                         diagnostic_dir="output/analyses/cofactors")
}

# Now look at cofactors that are only available in hg19
not_reddy_hg19 = all_chip %>% dplyr::filter(lab!="Tim Reddy, Duke", assembly=="hg19", !(target %in% names(all_chip_regions)))               

# Remove combined peak calls.
not_reddy_hg19 = not_reddy_hg19 %>%
                      dplyr::filter(!grepl(";", biological_replicates)) %>%
                      dplyr::filter(status != "revoked")

# Remove ambiguous peak files which might be duplicates
not_reddy_hg19 = not_reddy_hg19 %>%
                      dplyr::filter(!(target=="MAX" & is.na(biological_replicates))) %>%
                      dplyr::filter(!(target=="USF1" & is.na(biological_replicates))) %>%
                      dplyr::filter(!(target=="USF1" & treatment=="ethanol" & date_released=="2011-07-18")) %>%
                      dplyr::filter(!(target=="BHLHE40" & output_type=="peaks"))
                      
# Load regions.
not_reddy_targets = unique(not_reddy_hg19$target)
for(target_name in not_reddy_targets) {
    all_chip_regions[[target_name]] = summarize_GGR_chip(not_reddy_hg19 %>% dplyr::filter(target==target_name),
                                                         diagnostic_dir="output/analyses/cofactors")
}
                      
                      
all_time_points = unique(unlist(lapply(all_chip_regions, names)))
for(time_point in all_time_points) {
    target_at_time = list()
    for(target in names(all_chip_regions)) {
        if(time_point %in% names(all_chip_regions[[target]])) {
            target_at_time[[target]] = all_chip_regions[[target]][[time_point]]
        }
    }
    target_at_time = GenomicRanges::GRangesList(target_at_time)
    
    time_intersect = ef.utils::build_intersect(target_at_time)
    pdf_width = ifelse(time_point=="0 minute", 14, 7)
    pdf(paste0("output/analyses/cofactors/Clustering of cofactors at time ", time_point, ".pdf"), width=pdf_width)
    plot(hclust(dist(t(time_intersect$Matrix))))
    dev.off()
}
