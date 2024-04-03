### AUXILIARY FUNCTIONS ####
#------------------------------------------------------------------------------#

makeStanData <- function (inputData, K = NULL) {
# returns nested lists (a list of lists) with input data for Stan model runs
 splitDFlist <- split(x = inputData, f = inputData$Fish_number)
 Nfish <- length(splitDFlist)
  
 masterList <- rep(list(list()),Nfish)
 for (i in 1:Nfish) {
 masterList[[i]] <- list(y = splitDFlist[[i]]$`Sr/Ca`,#/mean(splitDFlist[[i]]$`Sr/Ca`, na.rm = TRUE),
                        N = length(splitDFlist[[i]]$OR), K = K)
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
  
  output <- list(outputDF, Y, X) # merged data frame, list of individual dataframes, and list of stan fits, respectively
  
  return(output)

}

#------------------------------------------------------------------------------#

graphPPchecks <- function(obDataList, fitDataList, Nfish = 51) {
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
               linewidth = 0.6) +
  facet_wrap(~Fish_number, scales = "free_y")+
  scale_linetype_manual(values = c(1,2)) + # 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash
  scale_color_manual(values = c("seagreen", "black")) +
  theme_custom() + ylab("Density") +
  xlab(bquote('Sr:Ca x 10'^{-1})) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_x_continuous(breaks = c(0.01, 0.04, 0.07),
                     limits = c(0,0.09)) +
  scale_y_continuous(n.breaks = 3)

}

#------------------------------------------------------------------------------#

plotThetas <- function(fitList, Nfish = 51) {
# plot mixture component proportions
  mcmcList = list()
  for (i in 1:Nfish) {
    mcmcList[[i]] = ggs(fitList[[i]]) %>% 
      filter(Parameter=='theta[1]') %>%
      mutate(Fish_number = i)
  }
  
  merged = do.call(rbind,mcmcList)
  merged %>%
    mutate(Fish_number = fct_reorder(as.factor(Fish_number),
                                     value, .fun = "median")) %>%
    ggplot(aes(x=reorder(as.factor(Fish_number),value),
                            y=value)) +
    see::geom_violinhalf(scale = "width", color = "seagreen") + coord_flip() + 
    theme_custom()+
    xlab("Fish #") + ylab(bquote(theta[Brackish]))
}

#------------------------------------------------------------------------------#

plotProfiles <- function(data, Nfish = NULL) {
  
  # create new fishnumber column
  data$fishNumber <- as.numeric(gsub(".*?([0-9]+).*", "\\1", as.character(data$Fish_number)))
  
  data <- data %>% filter(fishNumber <= Nfish) 
    ggplot(data, aes(x = OR, y = `Sr/Ca`)) +
      geom_line(color = "seagreen") +
      geom_hline(yintercept = 0.064) +
      facet_wrap(~as.factor(fishNumber), ncol = 2) +
      xlim(0,750) + ylab(bquote('Sr:Ca x 10'^{-1})) +
      xlab('Otolith radius (micromolar)')+
      theme_custom()

}

#------------------------------------------------------------------------------#

getSalinity <- function(Sr) {
  
  Sr = Sr*100
  
  sal = (Sr-2.0595)/0.2385 
  # from Santana et al. 2018
  #Santana, F. M., Morize, E., Labonne, M., Lessa, R., & Clavier, J. (2018).
  #Connectivity between the marine coast and estuary for white mullet (Mugil curema) 
  #in northeastern Brazil revealed by otolith Sr: Ca ratio. 
  #Estuarine, Coastal and Shelf Science, 215, 124-131.
  
  return(sal)
  
}

# cool ggplot theme ------------------------------------------------------------
theme_custom <- function(){ 
  font <- "sans"   #assign font family up front
  
  theme_minimal() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.major = element_line(color = "lightgrey",
                                      linewidth = 0.2,
                                      linetype = 11),
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_blank(),          #strip axis ticks
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 13,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 12),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 10),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 9),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10))
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}










