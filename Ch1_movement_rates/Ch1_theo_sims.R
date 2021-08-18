# Script to run simulations for theoretical populations with three levels of variation in movement rates

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(fitdistrplus)
library(Hmisc)

set.seed(98)

av.mov.rate <- 1
sd.mov.rate <- c(0.3, 0.5, 1)
max.x <- 10

mycols <- c("red", "green", "blue")
#create density plots
curve(dlnorm(x, meanlog=av.mov.rate, sdlog=sd.mov.rate[1]),
      from=0, to=max.x,
      col=mycols[1],
      main = 'Log Normally distributed movement rates for three populations', #add title
      ylab = 'Density', #change y-axis label
)

curve(dlnorm(x, meanlog=av.mov.rate, sdlog=sd.mov.rate[2]),
      from=0, to=max.x,
      col=mycols[2], add=TRUE)

curve(dlnorm(x, meanlog=av.mov.rate, sdlog=sd.mov.rate[3]),
      from=0, to=max.x,
      col=mycols[3], add=TRUE)

#add legend
legend("topright", legend=c("sdlog=.3", "sdlog=.5", "sdlog=1"),
       col=mycols, lty=1, cex=1.2)
