# PLOT TIME SERIES OF Sr/Ca RATIOS 
require(ggplot2)
library(readr)

MC_SrCa <- read_csv("M_curvidens_SrCa.csv")

OR <- as.numeric(MC_SrCa$`Otolith radius`) # grab otolith radius variable
otoliths <- as.data.frame(MC_SrCa[2:55]) # select columns containing the chemical ratios

col_names <- c(sprintf("fish_%02d", seq(1,ncol(otoliths)))) # name columns as fish_01, ..., fish_0N
colnames(otoliths) <- col_names

Sr_data <- cbind(OR, otoliths)
Sr_data[-1][Sr_data[-1]>=0.09] <- NA # exclude points that have ablated after the edge
N = ncol(Sr_data)-1
cols <- colnames(Sr_data[,-1]) # get column names

for(i in 1:N) {
  
  plot <- ggplot(data = Sr_data, aes_string(x = OR, y = cols[i])) + 
    geom_line() + theme_classic() +
    geom_hline(yintercept = 0.03, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    xlim(0,800)
  
  ggsave(plot, file = paste0("Fish_", i, ".png"),
         width = 15, height = 8, units = "cm")
  # this will save plots to the current directory
  
}
