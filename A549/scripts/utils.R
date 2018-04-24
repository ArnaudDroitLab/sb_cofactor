library(metagene)

cached_metagene <- function(regions, design, label, bin_count, cache_path,
                            new_args=list(assay='chipseq', force_seqlevels=TRUE),
                            produce_table_args=list(normalization="RPM", flip_regions=TRUE),
                            produce_data_frame_args=list()) {
    dir.create(cache_path, recursive=TRUE, showWarnings=FALSE)

    loaded.cache.filename = file.path(cache_path, paste0(region_name, " Loaded.RData"))
    if(!file.exists(loaded.cache.filename)) {
        metagene.obj = do.call(metagene$new, c(list(regions=regions, bam_files=design$Sample), new_args))
        save(metagene.obj, file=loaded.cache.filename)
    } else {
        load(loaded.cache.filename)
    }
    
    matrix.cache.filename = file.path(cache_path, paste0(region_name, " matrix.RData"))
    if(!file.exists(matrix.cache.filename)) {
        do.call(metagene.obj$produce_table, c(list(design = design, bin_count=bin_count), produce_table_args))
        save(metagene.obj, file=matrix.cache.filename)
    } else {
        load(matrix.cache.filename)
    }
    
    df.cache.filename = file.path(cache_path, paste0(region_name, " data-frame.RData"))
    if(!file.exists(df.cache.filename)) {
        do.call(metagene.obj$produce_data_frame, produce_data_frame_args)
        save(metagene.obj, file=df.cache.filename)
    } else {
        load(df.cache.filename)
    }
    
    return(metagene.obj)
}