# 20211207 Data analysis

# Libraries -----------------------------------------------
set.seed(20211207)

library(aracari)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(fitdistrplus)

# Data ----------------------------------------------------
# The script for this data is in the data-raw folder
load("data/ptpl.rda")

ptpl %>%
  filter(., mpm != 0) -> ptpl


dist.used <- c("gamma", "weibull", "lnorm")
param <- c("shape", "rate", "shape", "scale", "meanlog", "sdlog")
dist <- c("gamma", "gamma", "weibull", "weibull", "lnorm", "lnorm")
kpars <- c(2, 2, 2, 2, 2, 2)

# Functions ----------------------------------------------------

build_fits_df <- function(x){

  nm <- deparse(substitute(x))

  x %>%
    map(., `[`, c("estimate", "sd", "loglik", "n")) %>%
    map(as_tibble) %>%
    bind_rows()
}

# Figure 1 ---------------------------------------------
# Distribution of velocities at the different pooling levels

ptpl %>%
  ggplot(., aes((x = mpm))) +
  geom_line(aes(x = mpm, group = id), stat = "density", alpha = 0.4) +
  # geom_line(aes(x = mpm, group = fam_g), stat = "density", alpha = 0.4,
  #           color = "blue") +
  geom_line(stat = "density", color = "#1191d1", lwd = 1) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  labs(y = "Density", x = "Velocity in m/min") +
  theme_classic()

# Complete pooling fits -------------------------------------
pop <- lapply(dist.used, function(x){fitdist(ptpl$mpm, distr = x)})

# Get data for qqplot at population level
qqpop <- qqcomp(pop, plotstyle = "ggplot")$data
names(qqpop) <- c("theoretical", "fdist", "empirical")

# Save the population level parameters

build_fits_df <- function(x){

  nm <- deparse(substitute(x))

  x %>%
    map(., `[`, c("estimate", "sd", "loglik", "n")) %>%
    map(as_tibble) %>%
    bind_rows()
}

prms.df <- build_fits_df(pop) %>%
  mutate(param = param,
         dist = dist,
         kpars = kpars,
         data = "pop")


# NP -------------------------------------------------
# Repeating the process for data at the individual level
inds <- unique(as.character(nu_df$B_id))

# Fitting the four distributions for each individual and saving the parameters in a dataframe
ind.models <- NULL
for(i in 1:length(inds)){
  fit.df <- nu_df %>%
    dplyr::filter(., B_id == inds[i])

  fit <- lapply(dist.used, function(x){fitdist(fit.df$nu, distr = x)})

  ind.models[[i]] <- fit

  assign(paste(inds[i]), fit)

  fit.prms <- build_fits_df(fit) %>%
    mutate(param = param,
           dist = dist,
           kpars = kpars,
           data = inds[i])

  prms.df <- bind_rows(prms.df, fit.prms)
}


# Getting qqplot data for each individual
inds.models.data <- NULL
for(i in 1:length(inds)){
  x <- ind.models[[i]]

  qqdata <- qqcomp(x, plotstyle = "ggplot")$data %>%
    mutate(data = inds[i])
  names(qqdata) <- c("theoretical", "fdist", "empirical", "individual")
  inds.models.data <- bind_rows(inds.models.data, qqdata)
}



