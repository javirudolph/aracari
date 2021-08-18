# Script to run simulations for theoretical populations with three levels of variation in movement rates

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(fitdistrplus)
library(Hmisc)

set.seed(98)

# Three populations ---------------------------------------------------------------------
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

## Simulate movement rates -------------------------------------------------------------
n.individuals <- 30
m_1 <- sort(round(rlnorm(n.individuals, meanlog = av.mov.rate, sdlog = sd.mov.rate[1]),3))
m_2 <- sort(round(rlnorm(n.individuals, meanlog = av.mov.rate, sdlog = sd.mov.rate[2]),3))
m_3 <- sort(round(rlnorm(n.individuals, meanlog = av.mov.rate, sdlog = sd.mov.rate[3]),3))

## Visualize individual differences ----------------------------------------------------

par(mfrow=c(1,3))
mycols <- rainbow(n.individuals, s=1)

range <- c(-0.1,max.x+5)
#create density plots
### Fist population--------------------------------------------------------------------
curve(dexp(x, 1/m_1[1]), #notice rate for exponential is 1/movement rate.
      from=range[1], to=range[2],
      col=mycols[1],
      main = 'Distribution of movement lengths', #add title
      ylab = 'Density', #change y-axis label
      xlab = "x"
)
for(i in 2:5){
  curve(dexp(x, 1/m_1[i]),
        from=range[1], to=range[2],
        col=mycols[i], add=TRUE)
}

#add legend
legend("topright", legend= m_1,
       col=mycols, lty=1, cex=1.2)


### Second population----------------------------------------------------------------
curve(dexp(x, 1/m_2[1]), #notice rate for exponential is 1/movement rate.
      from=range[1], to=range[2],
      col=mycols[1],
      main = 'Distribution of movement lengths', #add title
      ylab = 'Density', #change y-axis label
      xlab = "x"
)
for(i in 2:5){
  curve(dexp(x, 1/m_2[i]),
        from=range[1], to=range[2],
        col=mycols[i], add=TRUE)
}

#add legend
legend("topright", legend= m_2,
       col=mycols, lty=1, cex=1.2)



### Third population ---------------------------------------------------------------
curve(dexp(x, 1/m_3[1]), #notice rate for exponential is 1/movement rate.
      from=range[1], to=range[2],
      col=mycols[1],
      main = 'Distribution of movement lengths', #add title
      ylab = 'Density', #change y-axis label
      xlab = "x"
)
for(i in 2:5){
  curve(dexp(x, 1/m_3[i]),
        from=range[1], to=range[2],
        col=mycols[i], add=TRUE)
}

#add legend
legend("topright", legend= m_3,
       col=mycols, lty=1, cex=1.2)
