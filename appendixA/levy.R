library(aracari)
library(fitdistrplus)
library(rmutil)
library(dplyr)


data(ptpl)

df <- ptpl %>%
  dplyr::filter(mpm != 0) %>%
  ungroup() %>%
  mutate(Bird_ID = paste0("B", Bird_ID)) %>%
  group_by("id", "burst") %>%
  mutate(cumdt = cumsum(dt)/60) %>%
  ungroup()

mdata <- df %>%
  #dplyr::filter(., Bird_ID == "B1") %>%
  dplyr::filter(., mpm > 50) %>%
  dplyr::select(mpm) %>%
  pull()


dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))
fitgumbel <- fitdist(mdata, "gumbel", start=list(a=10, b=10))
summary(fitgumbel)
plot(fitgumbel)



dlevy <- function(x, m, s) sqrt(s/(2*pi*(x-m)^3))*exp(-s/(2*(x-m)))
plevy <- function(q, m, s) 2*(1-pnorm(1/sqrt((q-m)/s)))
qlevy <- function(p, m, s) s/qnorm(1-p/2)^2+m
rlevy <- function(n, m, s) s/qnorm(1-runif(n)/2)^2+m
fitlevy <- fitdist(mdata, "levy", start=list(m=0, s=0.5))
summary(fitlevy)
plot(fitlevy)


fitnorm <- fitdist(mdata, "norm")
summary(fitnorm)
plot(fitnorm)
