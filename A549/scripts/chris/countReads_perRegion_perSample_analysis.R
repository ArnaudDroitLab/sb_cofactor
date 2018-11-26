setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/ckn_utils.R")
library(plotly)

count_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_perRegion_perSamples"

###
# correlation entre ces deux listes de chiffres
# faire une PCA
### speCTRL
countTable_speCTRL_path <- file.path(count_path, "countTable_speNBC_CTRL.txt")
countTable_speCTRL <- read.table(countTable_speCTRL_path, header = TRUE)

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

### common
countTable_common_path <- file.path(count_path, "countTable_NBC_common.txt")
countTable_common <- read.table(countTable_common_path, header = TRUE)

cor(countTable_common$NIPBL_CTRL, countTable_common$BRD4_CTRL)
cor(countTable_common$NIPBL_CTRL, countTable_common$CDK9_CTRL)
cor(countTable_common$BRD4_CTRL, countTable_common$CDK9_CTRL)
cor(countTable_common$NIPBL_DEX, countTable_common$BRD4_DEX)
cor(countTable_common$NIPBL_DEX, countTable_common$CDK9_DEX)
cor(countTable_common$BRD4_DEX, countTable_common$CDK9_DEX)

plot_ly(data = countTable_common, x = ~NIPBL_CTRL, y = ~BRD4_CTRL, type = "scatter", mode = "markers")
plot_ly(data = countTable_common, x = ~NIPBL_CTRL, y = ~CDK9_CTRL, type = "scatter", mode = "markers")
plot_ly(data = countTable_common, x = ~BRD4_CTRL, y = ~CDK9_CTRL, type = "scatter", mode = "markers")
plot_ly(data = countTable_common, x = ~NIPBL_DEX, y = ~BRD4_DEX, type = "scatter", mode = "markers")
plot_ly(data = countTable_common, x = ~NIPBL_DEX, y = ~CDK9_DEX, type = "scatter", mode = "markers")
plot_ly(data = countTable_common, x = ~BRD4_DEX, y = ~CDK9_DEX, type = "scatter", mode = "markers")

plot_ly(data = countTable_common, y = ~NIPBL_CTRL, type = "box") %>% add_trace(y = ~NIPBL_DEX) %>%
  add_trace(y = ~BRD4_CTRL) %>% add_trace(y = ~BRD4_DEX) %>%
  add_trace(y = ~CDK9_CTRL) %>% add_trace(y = ~CDK9_DEX)

t.test(countTable_common$NIPBL_CTRL, countTable_common$NIPBL_DEX)
t.test(countTable_common$BRD4_CTRL, countTable_common$BRD4_DEX)
t.test(countTable_common$CDK9_CTRL, countTable_speCTRL$CDK9_DEX)

### speDEX
countTable_speDEX_path <- file.path(count_path, "countTable_speNBC_DEX.txt")
countTable_speDEX <- read.table(countTable_speDEX_path, header = TRUE)

cor(countTable_speDEX$NIPBL_CTRL, countTable_speDEX$BRD4_CTRL)
cor(countTable_speDEX$NIPBL_CTRL, countTable_speDEX$CDK9_CTRL)
cor(countTable_speDEX$BRD4_CTRL, countTable_speDEX$CDK9_CTRL)
cor(countTable_speDEX$NIPBL_DEX, countTable_speDEX$BRD4_DEX)
cor(countTable_speDEX$NIPBL_DEX, countTable_speDEX$CDK9_DEX)
cor(countTable_speDEX$BRD4_DEX, countTable_speDEX$CDK9_DEX)

plot_ly(data = countTable_speDEX, x = ~NIPBL_CTRL, y = ~BRD4_CTRL, type = "scatter", mode = "markers")
plot_ly(data = countTable_speDEX, x = ~NIPBL_CTRL, y = ~CDK9_CTRL, type = "scatter", mode = "markers")
plot_ly(data = countTable_speDEX, x = ~BRD4_CTRL, y = ~CDK9_CTRL, type = "scatter", mode = "markers")
plot_ly(data = countTable_speDEX, x = ~NIPBL_DEX, y = ~BRD4_DEX, type = "scatter", mode = "markers")
plot_ly(data = countTable_speDEX, x = ~NIPBL_DEX, y = ~CDK9_DEX, type = "scatter", mode = "markers")
plot_ly(data = countTable_speDEX, x = ~BRD4_DEX, y = ~CDK9_DEX, type = "scatter", mode = "markers")

plot_ly(data = countTable_speDEX, y = ~NIPBL_CTRL, type = "box") %>% add_trace(y = ~NIPBL_DEX) %>%
  add_trace(y = ~BRD4_CTRL) %>% add_trace(y = ~BRD4_DEX) %>%
  add_trace(y = ~CDK9_CTRL) %>% add_trace(y = ~CDK9_DEX)

t.test(countTable_speDEX$NIPBL_CTRL, countTable_speDEX$NIPBL_DEX)
t.test(countTable_speDEX$BRD4_CTRL, countTable_speDEX$BRD4_DEX)
t.test(countTable_speDEX$CDK9_CTRL, countTable_speCTRL$CDK9_DEX)
