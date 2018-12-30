setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/chris/countRead_overtime/countReads.utils.R")
library(knitr)
library(ggplot2)
library(scales)

# Get ratio
EP300 <- retrieve_sumcount("count_total_EP300.txt")
GR <- retrieve_sumcount("count_total_GR2.txt")
GR_ratio <- retrieve_sumcount("count_total_GR_background.txt")
names(GR) <- names(GR_ratio)

GR_scaled <- assignTime(GR)
EP300_scaled <- assignTime(EP300)

df_GR_scaled <- dfForPlot(GR_scaled, "GR")
kable(df_GR_scaled)
df_EP300_scaled <- dfForPlot(EP300_scaled, "EP300")
kable(df_EP300_scaled)

df_GR_EP300 <- rbind(df_GR_scaled, df_EP300_scaled)

df_mean <- aggregate(df_GR_EP300$Value, by=list(df_GR_EP300$Time, df_GR_EP300$Protein), mean)
colnames(df_mean) <- c("time", "protein", "mean")

# Make the plot
ggplot(df_GR_EP300, aes(x = Time, y = Value)) +
  geom_point(aes(color = Protein), size = 1) +
  geom_line(data=df_mean, aes(x=time, y=mean, group=protein, color=protein)) +
  scale_x_continuous(name="Time",
                     labels = c("0h", "", "", "", "", "", "", "1h", "2h", "3h", "4h", "5h", "6h", "7h", "8h", "10h", "12h"),
                     breaks = c(0, 5, 10, 15, 20, 25, 30, 60, 120, 180, 240, 300, 360, 420, 480, 600, 720)) +
  # scale_x_continuous(name="Time", labels=c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12), breaks=c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12)) +
  scale_y_continuous(name="Value", labels=comma)
