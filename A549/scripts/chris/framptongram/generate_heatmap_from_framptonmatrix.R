setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

### Call python script to generate framptongram matrix
python_script_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/scripts/chris/framptongram/generate_framptongram_matrix.py"
input_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/scripts/chris/framptongram/bedfile_list_allReddy_withChromatin_onlyH.txt"
matrix_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/scripts/chris/framptongram/framptongram_matrix_Reddy_allReddy_withChromatin_onlyH.txt"

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

# Load framptomgram matrix and process it
mat <- process_frampton_matrix(matrix_path)

# Add annotation
sample_factors <- rownames(mat) %>% strsplit(split = "_") %>% get_nth_element(1)
colors_factors <- c("GR" = "#7B241C", "EP300" = "#633974", "JUN" = "#1A5276", "JUNB" = "#117864",
                    "FOSL2" = "#F39C12", "HES2" = "#BA4A00", "CTCF" = "#212F3D", "SMC3" = "#4E342E",
                    "RAD21" = "#29B6F6", "BCL3" = "#E91E63", "CEBPB" = "#757575", "H3K27ac" = "#32CD32",
                     "H3K4me1" = "#556B2F", "H3K4me2" = "#000000", "H3K4me3" = "#FFEB3B")
                    #32CD32)
colha <- HeatmapAnnotation(Factors = sample_factors,
                           col = list(Factors = colors_factors))
rowha <- rowAnnotation(Factors = sample_factors,
                      col = list(Factors = colors_factors),
                      show_legend = FALSE)

# ha <- HeatmapAnnotation(Factors = samples_factors,
#                         col = list(Factors = c("KD" = "#16DB93", "OE" = "#EFEA5A")))
# rowha = rowAnnotation(Experiment = col_kd_oe,
#                       col = list(Experiment = c("KD" = "#16DB93", "OE" = "#EFEA5A")),
#                       show_legend = FALSE)

# Generate heatmap
customColors = colorRamp2(c(-1, 0, 1), c("#0f4259", "white", "#800020"))
hm2 <- Heatmap(mat, name = "Correlation",
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 11), # change to 10/11 if necessary
              row_dend_side = "right",
              column_names_side = "top",
              column_names_gp = gpar(fontsize = 11), # change to 10/11 if necessary
              column_names_rot = 45,
              column_dend_side = "bottom",
              show_column_dend = FALSE,
              row_dend_width = unit(50, "mm"),
              column_dend_height = unit(50, "mm"),
              col = customColors,
              rect_gp = gpar(col = "white", lwd = 1), # change to 0.5/1 if necessary
              top_annotation = colha,
              left_annotation = rowha)

# to add correlation value
# cell_fun = function(j, i, x, y, width, height, fill) {
# grid.text(sprintf("%.2f", cor.mutcof[i, j]), x, y, gp = gpar(fontsize = 10))
# }
              

png("output/analyses/heatmap_framptongram/20190625_framptongram_allReddy_withChromatin_onlyH.png",
    width = 4000, height = 3700) # change to w = 3000/2000/1500 and h = 2700/1800/1300
hm2
dev.off()

