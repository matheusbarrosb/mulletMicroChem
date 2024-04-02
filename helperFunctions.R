### AUXILIARY FUNCTIONS ####
#------------------------------------------------------------------------------#

makeStanData <- function (inputData) {
# returns nested lists (a list of lists) with input data for Stan model runs
 splitDFlist <- split(x = inputData, f = inputData$Fish_number)
 Nfish <- length(splitDFlist)
  
 masterList <- rep(list(list()),Nfish)
 for (i in 1:Nfish) {
 masterList[[i]] <- list(y = splitDFlist[[i]]$`Sr/Ca`,#/mean(splitDFlist[[i]]$`Sr/Ca`, na.rm = TRUE),
                        N = length(splitDFlist[[i]]$OR), K = 2)
 }
  return(masterList)
}

#------------------------------------------------------------------------------#

runModels <- function (masterDataList, N = 51, 
                       inits, nIter,
                       nChains, nWarmup) {
  
  # runs many models consecutively and stores outputs in merged dataframe
  
  X <- list()
  for (i in 1:N) {
    # store individual fits in a N-dimensional list
    X[[i]]<-stan(file = "mixture.stan", 
                 data = stanMasterDataList[[i]],
                 warmup = nWarmup,
                 iter = nIter,
                 chains = nChains,
                 init = inits)
  }
  
  Y <- list()
  for (i in 1:N) {
  # converts individual fits in previous list to ggs dataframes  
    Y[[i]] <- ggs(X[[i]])
    
  }
  
  Z <- list()
  for (i in 1:N) {
    # include FishID identifier in dataframes
    index <- print(paste("", i, sep = ""))
    Z[[i]] <- Y[[i]] %>% mutate(fishID = paste(index))
    
  }
  
  outputDF <- do.call(rbind, Z) # merge all dataframes
  
  return(outputDF)

}

#------------------------------------------------------------------------------#

graphPPchecks <- function(obDataList, fitDataList, Nfish = 6) {
# performs graphical posterior predictive checks
# object 'fitDataList' has to be a list of dataframes provided by ggs() function from ggmcmc package
  
  #set up error messages
  if (!is.list(obDataList) | !is.list(fitDataList)) {
    stop("Wrong input type. Provide inputs as list of dataframes")
  }
  
observed = list()
for (i in 1:Nfish) {
  observed[[i]] = data.frame(obDataList[[i]]$y,
                             rep(i, length(obDataList[[i]]$y)))
  names(observed[[i]]) = c("value", "Fish_number")
}
observed <- do.call(rbind, observed)
observed <- observed %>% mutate(Source = "Observed")

fits = list()
for (i in 1:Nfish) {
  fits[[i]] = data.frame(log(fitDataList[[i]]$value),
                         rep(i, length(fitDataList[[i]]$value)))
  names(fits[[i]]) = c("value", "Fish_number")
}
fits <- do.call(rbind, fits)
fits <- fits %>% mutate(Source = "Expected")

mergedDF <- rbind(fits, observed)

ggplot(data = mergedDF) + 
  geom_density(aes(x = value, 
                   color = Source,
                   linetype = Source),
               linewidth = 1) +
  facet_wrap(~Fish_number)+ xlim(0,0.08) +
  theme_classic() + ylab("Density") +
  xlab("Sr/Ca")

}

#------------------------------------------------------------------------------#

plotThetas <- function(fitList, Nfish = 6) {
# plot mixture component proportions
  mcmcList = list()
  for (i in 1:Nfish) {
    mcmcList[[i]] = ggs(fitList[[i]]) %>% 
      filter(grepl('theta', Parameter)) %>%
      mutate(Fish_number = i)
  }
  
  merged = do.call(rbind,mcmcList)
  merged = merged %>% mutate(Habitat = ifelse(Parameter == 'theta[1]',
                                              'Brackish', 'Sea'))
  
  ggplot(data = merged, aes(x=as.factor(Fish_number),
                            y=value,
                            color=Habitat)) +
    geom_boxplot() + theme_classic()
}










