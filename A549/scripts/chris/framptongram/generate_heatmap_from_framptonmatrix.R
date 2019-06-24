setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(ComplexHeatmap)
library(circlize)

### Call python script to generate framptongram matrix
python_script_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/scripts/chris/framptongram/generate_framptongram_matrix.py"
input_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/scripts/chris/framptongram/bedfile_list_Reddy_allData.txt"
output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/scripts/chris/framptongram/framptongram_matrix_Reddy_allData.txt"

call <- paste("python", python_script_path, input_path, output_path)
message(call)
system(call)

### Generate heatmap from framptongram matrix
data <- read.table(output_path, header = TRUE)

mat <- as.matrix(data[, 4:ncol(data)])

sample_names_tmp <- as.character(data$NAME)
basename(sample_names_tmp)

h2 <- Heatmap(mat,
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
              col = col_fun,
              rect_gp = gpar(col = "white", lwd = 1))
