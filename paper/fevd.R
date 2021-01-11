
# Fit extreme value distributions to the dispersal kernel data

load("paper/dispersal_kernels.RData")

library(dplyr)
library(ggplot2)
library(cowplot)
library(extRemes)


# My data

null_sample <- sample(null_dispersal$dispersal, 10000)
hist(null_sample)
#hist(log(null_sample))
#qqnorm(log(null_sample))
qqnorm(null_sample)
# Qq plots don't show straight lines, not consistent with a normal distribution of the data,

a <- threshrange.plot(null_sample, r = c(300, 700), nint = 100)
b <- mrlplot(null_sample, xlim = c(300, 700), nint = 100)

fit_D <- fevd(null_sample, threshold = 500, type = "GP")
fit_D$results
summary(fit_D)
plot(fit_D)
t <- plot(fit_D, type = "qq2")
pextRemes(fit_D, c(500, 1000, 1500, 2000, 3000, 3500, 4000), lower.tail = FALSE)



hist(sample(indiv_dispersal$dispersal, 100))
threshrange.plot(sample(indiv_dispersal$dispersal, 100), r = c(1, 700), nint = 100)


garcia <- read.csv("paper/crap.txt", header = FALSE)
hist(garcia$V1)
threshrange.plot(garcia$V1, r = c(1, 700), nint = 200)
mrlplot(garcia$V1)

fit_garcia <- fevd(garcia$V1, threshold = 170, type = "GP")
summary(fit_garcia)

