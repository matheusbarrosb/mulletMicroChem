require(readr)

setwd("/Users/matheusb/Documents/mulletMicroChem/Data")
source("helperFunctions.R")


# DATA WRANGLING ####
bioData = read_csv("McurvidensBioData.csv")
y  = ifelse(as.numeric(bioData$EMG) == 1 & 2, 0, 1)
TL = as.numeric(bioData$`CT (mm)`)/10
df = data.frame(y, TL) 
df = na.exclude(df)
X  = model.matrix(~df$TL)
y  = df$y

plot(y ~ X[,2])

# MATURITY OGIVE
res = MHflatLogit(nIter = 80000, nBurnIn = 1000,
                y = y, X = X, delta = 0.5, plot = TRUE)

# GET L50
# means
alpha = res$posteriorSummary$Mean[1]
beta  = res$posteriorSummary$Mean[2]

# uppers and lowers
upperAlpha = res$posteriorSummary$Mean[1] + res$posteriorSummary$SD[1]
upperBeta  = res$posteriorSummary$Mean[2] + res$posteriorSummary$SD[2]
lowerAlpha = res$posteriorSummary$Mean[1] - res$posteriorSummary$SD[1]
lowerBeta  = res$posteriorSummary$Mean[2] - res$posteriorSummary$SD[2]

L50      = -alpha/beta
L50lower = -upperAlpha/upperBeta 
L50upper = -lowerAlpha/lowerBeta 

# GET PREDICTIONS
XX = seq(10, 40)
p = as.vector((exp(alpha + beta*XX))/(1 + exp(alpha + beta*XX)))
pLower = (exp(lowerAlpha + lowerBeta*XX))/(1 + exp(lowerAlpha + lowerBeta*XX))
pUpper = (exp(upperAlpha + upperBeta*XX))/(1 + exp(upperAlpha + upperBeta*XX))

# PLOTTING
plot(p ~ XX, type = "l")
lines(pUpper ~ XX)
lines(pLower ~ XX)

predDF    = data.frame(p, pUpper, pLower, XX)
alphaDens = as.vector(res$posteriorSamples[,1])
betaDens  = as.vector(res$posteriorSamples[,2])
L50Dens   = -alphaDens/betaDens
densDF    = data.frame(alphaDens, betaDens, L50Dens)

require(ggplot2)

ogivePlot = ggplot(data = predDF, aes(x = XX, y = p)) +
  geom_ribbon(aes(ymax = pUpper, ymin = pLower),
              alpha = 0.55, fill = "lightgreen") +
  geom_line() +
  geom_segment(x = -Inf, xend = 20.55, y = .5, yend = .5,
               color = "red", linetype = "dashed") +
  geom_segment(x = 20.55, xend = 20.55, y = -Inf, yend = .5,
               color = "red", linetype = "dashed") +
  theme_custom() + 
  xlab("TL (cm)") +
  ylab("Proportion mature") +
  geom_label(
    label = "L50 = 22.5 cm",
    x     = 28,
    y     = .5,
    label.padding = unit(0.3, "lines"),
    color = "black",
    fill="grey"
  ) +
  ggtitle("a)")


alphaDensPlot = ggplot(densDF, aes(x = alphaDens)) +
  geom_density(fill = "lightgreen") + theme_custom() +
  xlab("Parameter value") +
  ylab("Probability density") +
  ggtitle("b)")

betaDensPlot = ggplot(densDF, aes(x = betaDens)) +
  geom_density(fill = "lightgreen") + theme_custom() +
  xlab("Parameter value") +
  ylab("Probability density") +
  ggtitle("c)")

L50DensPlot = ggplot(densDF, aes(x = L50Dens)) +
  geom_density(fill = "lightgreen") + theme_custom() +
  xlab("Parameter value") +
  ylab("Probability density") +
  ggtitle("d)") + xlim(15, 26)


require(ggpubr)

ggarrange(ogivePlot, alphaDensPlot, betaDensPlot, L50DensPlot)














  