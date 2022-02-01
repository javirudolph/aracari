##########################################
## DATA GENERATION CH3
#########################################

# I think the overall question for this is chapter is:
#         Can we accurately estimate the true proportion of rare events produced by a complex model using a simpler hierarchical model?

# Our complex model is one that considers individual variation in animal movement, with increasing mean and variance for specific groups
# We use a mixture modeling approach, where each component distribution is associated to an individual's state or age.
# We assume that an individual's stage will affect their movement.
# We assume, for example, very young birds might not move as far, as will nesting or molting birds, whereas adults will move more.

# The overall assumption for this model is that individuals vary in how much they move, which we determine based on step lengths only.
# With some groups having an overall larger variance and mean than other groups.
# The proportions, or weights, of the mixture, should follow a bell curve for now
# With the average movers being the majority of the group, and low or high movement individuals being less abundant.

# Libraries -----------------------------------------------
set.seed(20220201)

library(aracari)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(fitdistrplus)
library(Hmisc)
library(plotly)

## Functions ----------------------------------------------------
# Building this function so that we get the parameters for a desired mean and standard deviation.
desired_mean_sd <- function(mean, sd){

  mu <- log(mean*2/(sqrt(mean*2+sd*2)))
  sigsq <- log(1+(sd*2/mean*2))

  return(data.frame(mu = mu, sigsq = sigsq))
}

exp.val <- function(mu, sigsq){
  exp(mu + sqrt(sigsq))
}


# MIXTURE MODEL --------------------------------------------------
## Mixture proportions ----------------------------------------
# We assume four categories of individuals with increasing movement.
pis <- c( 0.15, 0.3, 0.5, 0.05)
histogram(pis)

## Mixture components --------------------------------------
# Since movement lengths are only positive, we use a lognormal distribution to describe them

# CASE1: same variance for all, but change expected values

pars <- desired_mean_sd(mean = c(160, 300, 600, 1000), sd = c(90, 120, 160, 200))

# Special case I am considering where for each category the mean doubles, and then sd=mean*0.6
# pars <- desired_mean_sd(mean = c(150, 300, 600, 1200), sd = c(90, 180, 360, 720))

mus <- pars$mu
sigsqs <- pars$sigsq
dens_cols <- c("black", "blue", "green", "red")

### PLOT the mixture components -------------------------------
lnorm_densities <- purrr::map(1:4, function(y) stat_function(fun = dlnorm,
                                                             args = list(meanlog = mus[y], sdlog = sigsqs[y]),
                              color = dens_cols[y]))
ggplot() +
  lnorm_densities +
  theme_minimal() +
  lims(x = c(0, 150))


# SAMPLING --------------------------------------------------------------










