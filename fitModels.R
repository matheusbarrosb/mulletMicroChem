list.of.packages <- c("rstan", "bayesplot", "dplyr", "tidyverse", "ggmcmc", "see")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

require(rstan)
require(bayesplot)
require(dplyr)
require(tidyverse)
require(ggmcmc)
require(see)

source("helperFunctions.R")

setwd("/Users/matheusb/Documents/mulletMicroChem")

# DATA WRANGLING ####
MC_SrCa <- read_csv("M_curvidens_SrCa.csv")
OR <- as.numeric(MC_SrCa$`Otolith radius`) # grab otolith radius variable
otoliths <- as.data.frame(MC_SrCa[2:55]) # select columns containing the chemical ratios
col_names <- c(sprintf("Fish_%1d", seq(1,ncol(otoliths)))) # name columns as fish_01, ..., fish_0N
colnames(otoliths) <- col_names
Sr_data <- cbind(OR, otoliths)
Sr_data[-1][Sr_data[-1]>=0.09] <- NA # exclude points that have ablated after the edge
Sr_data <- Sr_data[-c(35)] # exclude wonky fish
# transforming data to long format so fish # is a column
Sr_data2 <- Sr_data %>% 
  pivot_longer(!OR, names_to = "Fish_number", values_to = "Sr/Ca")
# formatting data for stan
Sr_data2 <- na.exclude(Sr_data2)
# grab fish No
selected <- paste0("Fish_", 1:52)
#selected <- paste0("Fish_", 1:length(unique(Sr_data2$Fish_number)))
Sr_data2<-Sr_data2[Sr_data2$Fish_number %in% selected,]
# prepare data for stan
stanMasterDataList <- makeStanData(Sr_data2)

# get initial parameters
inits <- function() {
  inits = list(mu = c(0.03,0.07),
               sigma = c(0.0075, 0.001))
  return(inits)
}

stanc("mixture.stan") # check if model code compiles
# RUNNING ALL MODELS ####
fits <- runModels(masterDataList = stanMasterDataList,
          inits = inits, N = 2, nWarmup = 500,
          nIter = 1000, nChains = 51)

# GRAPHICAL PP CHECKS ####
graphPPchecks(obData = stanMasterDataList,
              fitDataList = fits[[2]], Nfish = 51)







