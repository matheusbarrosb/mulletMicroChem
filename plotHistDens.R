# PLOT HISTOGRAMS/DENSITIES OF Sr/Ca RATIOS 
require(ggplot2)
library(readr)

# EXECUTE THOSE IN ORDER #######################################################
# BUT CHANGE DIRECTORIES ACCORDINGLY
setwd("/Users/matheusb/Documents/mulletMicroChem")
MC_SrCa <- read_csv("M_curvidens_SrCa.csv")
setwd("/Users/matheusb/Documents/mulletMicroChem/SrCaHistograms")
################################################################################

MC_SrCa <- read_csv("M_curvidens_SrCa.csv")
OR <- as.numeric(MC_SrCa$`Otolith radius`) # grab otolith radius variable
otoliths <- as.data.frame(MC_SrCa[2:55]) # select columns containing the chemical ratios
col_names <- c(sprintf("%02d", seq(1,ncol(otoliths)))) # name columns as 01, ..., 0N
colnames(otoliths) <- col_names
Sr_data <- cbind(OR, otoliths)
Sr_data[-1][Sr_data[-1]>=0.09] <- NA # exclude points that have ablated after the edge
Sr_data <- Sr_data[-c(35)] # exclude wonky fish
# transforming data to long format so fish # is a column
Sr_data2 <- Sr_data %>% 
  pivot_longer(!OR, # ignore first column
               names_to = "Fish_number", values_to = "Sr/Ca")


ggplot(data = Sr_data2, aes(x = `Sr/Ca`)) +
  geom_histogram(fill = 'white',
                 color = 'black',
                 bins = 50) +
  theme_classic() +
  geom_vline(xintercept = 0.03, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  facet_wrap(~Fish_number)


Sr_data2 %>% na.exclude() %>%
  ggplot(aes(x = `Sr/Ca`)) +
  geom_density(fill = 'grey',
                 color = 'black') +
  theme_classic() +
  geom_vline(xintercept = 0.03, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  facet_wrap(~Fish_number)






