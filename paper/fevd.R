
# Fit extreme value distributions to the dispersal kernel data

load("paper/dispersal_kernels.RData")

library(dplyr)
library(ggplot2)
library(cowplot)
library(extRemes)


data("damage", package="extRemes")

plot(damage$Year, damage$Dam, xlab = "Year",
     ylab = "U.S. Hurricane Damage (billion USD)", cex = 1.25,
     cex.lab = 1.25, col = "darkblue", bg = "lightblue", pch = 21)

plot(damage[, "Year"], log(damage[, "Dam"]), xlab = "Year",
     ylab = "ln(Damage)", ylim = c(-10, 5), cex.lab = 1.25,
     col = "darkblue", bg = "lightblue", pch = 21)

qqnorm(log(damage[, "Dam"]), ylim = c(-10, 5), cex.lab = 1.25)
hist(damage$Dam)
threshrange.plot(damage$Dam, r = c(1, 15), nint = 20)

mrlplot(damage$Dam, xlim = c(0, 12))


# My data

null_sample <- sample(null_dispersal$dispersal, 10000)
hist(null_sample)
hist(log(null_sample))
qqnorm(log(null_sample))
qqnorm(null_sample)
# Qq plots don't show straight lines, not consistent with a normal distribution of the data,

threshrange.plot(null_sample, r = c(300, 800), nint = 100)

fit_D <- fevd(null_sample, threshold = 500, type = "GP")
plot(fit_D)
pextRemes(fit_D, c(500, 1000, 1500, 2000, 3000, 3500, 4000), lower.tail = FALSE)

hist(sample(indiv_dispersal$dispersal, 100))
threshrange.plot(sample(indiv_dispersal$dispersal, 100), r = c(1, 700), nint = 100)


