require(rstan)
require(bayesplot)
require(dplyr)
require(tidyverse)
require(ggmcmc)
require(see)

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
selected <- paste0("Fish_", 1:6)
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
# RUNNING MODELS ####
## FISH 01 ####
fit01 <- rstan::stan(file = "mixture.stan", data = stanMasterDataList[[1]],
                   warmup = 2000, iter = 10000, chains = 2, init = inits)

## FISH 02 ####
fit02 <- rstan::stan(file = "mixture.stan", data = stanMasterDataList[[2]],
                     warmup = 2000, iter = 10000, chains = 2, init = inits)

## FISH 03 ####
fit03 <- rstan::stan(file = "mixture.stan", data = stanMasterDataList[[3]],
                     warmup = 2000, iter = 10000, chains = 2, init = inits)

## FISH 04 ####
fit04 <- rstan::stan(file = "mixture.stan", data = stanMasterDataList[[4]],
                     warmup = 2000, iter = 10000, chains = 2, init = inits)


# GRAPHICAL POSTERIOR PREDICTIVE CHECKS ####
fit01mcmc <- ggs(fit01) %>% filter(grepl('y_hat', Parameter))
fit02mcmc <- ggs(fit02) %>% filter(grepl('y_hat', Parameter))
fit03mcmc <- ggs(fit03) %>% filter(grepl('y_hat', Parameter))
fit04mcmc <- ggs(fit04) %>% filter(grepl('y_hat', Parameter))

fitDataList = list(fit01mcmc, fit02mcmc, fit03mcmc,
                   fit04mcmc)

graphPPchecks(obData = stanMasterDataList,
         fitDataList = fitDataList, Nfish = 4)









