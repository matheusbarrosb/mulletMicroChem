require(rstan)
require(bayesplot)
require(dplyr)
require(tidyverse)
require(ggmcmc)

setwd("/Users/matheusb/Documents/mulletMicroChem")
MC_SrCa <- read_csv("M_curvidens_SrCa.csv")
OR <- as.numeric(MC_SrCa$`Otolith radius`) # grab otolith radius variable
otoliths <- as.data.frame(MC_SrCa[2:55]) # select columns containing the chemical ratios
col_names <- c(sprintf("fish_%02d", seq(1,ncol(otoliths)))) # name columns as fish_01, ..., fish_0N
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
selected <- c("fish_02")

Sr_data2<-Sr_data2[Sr_data2$Fish_number %in% selected,]

# prepare data for stan
Sr_Ca <- Sr_data2$`Sr/Ca`
y <- Sr_Ca/mean(Sr_Ca, na.rm = TRUE) # normalized data
N <- length(y)
K <- 2
Nchains <- 2
stan_data <- list(y = y,
                  N = N,
                  K = K)

# get initial parameters
inits <- getInits(mu = c(0.9,1.3),
                  sigma = c(0.3, 0.01),
                  Nchains = Nchains)

stanc("mixture.stan") # check if model code compiles

# run model
fit <- rstan::stan(file = "mixture.stan", 
                   data = stan_data,
                   warmup = 200,
                   iter = 1000,
                   chains = Nchains,
                   init = inits)

#print main results
print(fit, pars = c("theta", "mu", "sigma"))
plot(fit, pars = c("theta"))





