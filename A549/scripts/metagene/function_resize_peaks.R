
resize_all_peaks <- function(granges.obj, window) {
	ranges(granges.obj) <- resize(ranges(granges.obj), width = window, fix="center")
	return (granges.obj)
}

resize_broad_peaks <- function(granges.obj, window) {
	l <- length(granges.obj)
	message("Number of regions : ", l)
	for (i in seq(1, l)) {
		w <- width(ranges(granges.obj[i]))
		if (w > window) {
			ranges(granges.obj[i]) <- resize(ranges(granges.obj[i]), width = window, fix="center")
			message("Region number ", i, " / ", l, " : ", w, " nt", " --» resized to ", window, " nt")
		}
		else {
			message("Region number ", i, " / ", l, " : ", w, " nt")
		}
	}
	return (granges.obj)
}
