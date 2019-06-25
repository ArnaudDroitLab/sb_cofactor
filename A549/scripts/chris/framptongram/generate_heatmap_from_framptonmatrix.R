setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

### Call python script to generate framptongram matrix
python_script_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/scripts/chris/framptongram/generate_framptongram_matrix.py"
input_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/scripts/chris/framptongram/bedfile_list_Reddy_GR_EP300_JUN_CTCF_SMC3_RAD21_CEBPB_BCL3_FOSL2_HES2_JUNB_onlyH.txt"
matrix_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/scripts/chris/framptongram/framptongram_matrix_Reddy_GR_EP300_JUN_CTCF_SMC3_RAD21_CEBPB_BCL3_FOSL2_HES2_JUNB_onlyH.txt"

call <- paste("python", python_script_path, input_path, matrix_path)
message(call)
system(call)

### Generate heatmap from framptongram matrix
#
get_nth_element <- function(lst, n) {
  sapply(lst, "[", n)
}

# Make proper sample names from list of files
make_sample_names <- function(data) {
  tmp1 <- as.character(data$NAME) %>% basename 
  tmp2 <- gsub("NR3C1", "GR", tmp1) %>% strsplit(split = "_")
  sample_names <- paste(get_nth_element(tmp2, 1), get_nth_element(tmp2, 2), sep = "_")
  return(sample_names)
}

# Open framptongram matrix file and return a formatted matrix, proper to be plot with Heatmap function from ComplexHeatmap package
process_frampton_matrix <- function(matrix_filename) {
  data <- read.table(matrix_filename, header = TRUE)
  sample_names <- make_sample_names(data)
  
  mat <- as.matrix(data[, 4:ncol(data)])
  colnames(mat) <- sample_names
  rownames(mat) <- sample_names
  
  return(mat)
}

## TODO Add annotation on side of the framptongram

# Load framptomgram matrix and process it
mat <- process_frampton_matrix(matrix_path)

# Generate heatmap
customColors = colorRamp2(c(-1, 0, 1), c("#0f4259", "white", "#800020"))
hm2 <- Heatmap(mat, name = "Correlation",
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 11),
              row_dend_side = "right",
              column_names_side = "top",
              column_names_gp = gpar(fontsize = 11),
              column_names_rot = 45,
              column_dend_side = "bottom",
              show_column_dend = FALSE,
              row_dend_width = unit(50, "mm"),
              column_dend_height = unit(50, "mm"),
              col = customColors,
              rect_gp = gpar(col = "white", lwd = 1),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 7))
              })

png("output/analyses/heatmap_framptongram/20190624_framptongram_GR_EP300_JUN_CTCF_SMC3_RAD21_CEBPB_BCL3_FOSL2_HES2_JUNB_onlyH.png",
    width = 1500, height = 1300)
hm2
dev.off()
