
# Fit extreme value distributions to the dispersal kernel data

load("paper/dispersal_kernels.RData")

library(dplyr)
library(ggplot2)
library(cowplot)
library(extRemes)


# Just testing

null_sample <- sample(null_dispersal$dispersal, 10000)
hist(null_sample)
#hist(log(null_sample))
#qqnorm(log(null_sample))
qqnorm(null_sample)
# Qq plots don't show straight lines, not consistent with a normal distribution of the data,

a <- threshrange.plot(null_sample, r = c(300, 700), nint = 200)
b <- mrlplot(null_sample, xlim = c(300, 700), nint = 100)

fit_D <- fevd(null_sample, threshold = 500, type = "GP")
#fit_D$results
summary(fit_D)
plot(fit_D)
t <- plot(fit_D, type = "qq")
pextRemes(fit_D, c(500, 1000, 1500, 2000, 3000, 3500, 4000), lower.tail = FALSE)

scale <- fit_D$results$par[1]
shape <- fit_D$results$par[2]
dens <- revd(10000, scale = scale, shape = shape, type = "GP")
hist(dens)
d <- density(dens)
plot(d)

ggplot(data = data.frame(x = c(0.1, 6)), aes(x)) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = 1, shape = 0)) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = 1, shape = 0.5))
  ylab("") +
  scale_y_continuous(breaks = NULL)


lines <- purrr::map(seq(0, 1, 0.1), function(y) stat_function(fun = devd, args = list(type = "GP", scale = 1, shape = y), color = "grey"))

ggplot(data = data.frame(x = c(0.1, 6)), aes(x)) +
  ylab("") +
  scale_y_continuous(breaks = NULL) -> p

p + lines +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = 1, shape = 0), linetype = "dashed", size = 2) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = 1, shape = 0.5))
