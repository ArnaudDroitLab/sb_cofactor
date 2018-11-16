setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/ckn_utils.R")
library(plotly)

count_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_perRegion_perSamples"

### speCTRL
countTable_speCTRL_path <- file.path(count_path, "countTable_speNBC_CTRL.txt")
countTable_speCTRL <- read.table(countTable_speCTRL_path, header = TRUE)

calculateCor <- function(countTable) {
}

cor(countTable_speCTRL$NIPBL_CTRL, countTable_speCTRL$BRD4_CTRL)
cor(countTable_speCTRL$NIPBL_CTRL, countTable_speCTRL$CDK9_CTRL)
cor(countTable_speCTRL$BRD4_CTRL, countTable_speCTRL$CDK9_CTRL)
cor(countTable_speCTRL$NIPBL_DEX, countTable_speCTRL$BRD4_DEX)
cor(countTable_speCTRL$NIPBL_DEX, countTable_speCTRL$CDK9_DEX)
cor(countTable_speCTRL$BRD4_DEX, countTable_speCTRL$CDK9_DEX)

plot_ly(data = countTable_speCTRL, x = ~NIPBL_CTRL, y = ~BRD4_CTRL, type = "scatter", mode = "markers")
plot_ly(data = countTable_speCTRL, x = ~NIPBL_CTRL, y = ~CDK9_CTRL, type = "scatter", mode = "markers")
plot_ly(data = countTable_speCTRL, x = ~BRD4_CTRL, y = ~CDK9_CTRL, type = "scatter", mode = "markers")
plot_ly(data = countTable_speCTRL, x = ~NIPBL_DEX, y = ~BRD4_DEX, type = "scatter", mode = "markers")
plot_ly(data = countTable_speCTRL, x = ~NIPBL_DEX, y = ~CDK9_DEX, type = "scatter", mode = "markers")
plot_ly(data = countTable_speCTRL, x = ~BRD4_DEX, y = ~CDK9_DEX, type = "scatter", mode = "markers")

plot_ly(data = countTable_speCTRL, y = ~NIPBL_CTRL, type = "box") %>% add_trace(y = ~NIPBL_DEX) %>%
  add_trace(y = ~BRD4_CTRL) %>% add_trace(y = ~BRD4_DEX) %>%
  add_trace(y = ~CDK9_CTRL) %>% add_trace(y = ~CDK9_DEX)

t.test(countTable_speCTRL$NIPBL_CTRL, countTable_speCTRL$NIPBL_DEX)
t.test(countTable_speCTRL$BRD4_CTRL, countTable_speCTRL$BRD4_DEX)
t.test(countTable_speCTRL$CDK9_CTRL, countTable_speCTRL$CDK9_DEX)

# correlation entre ces deux listes de chiffres
# faire une PCA