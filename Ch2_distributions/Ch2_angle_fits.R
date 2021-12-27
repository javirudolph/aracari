# 20211207 Data analysis

# Libraries -----------------------------------------------
set.seed(20211207)

library(aracari)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(fitdistrplus)
library(CircStats)

# Data ----------------------------------------------------
# The script for this data is in the data-raw folder
load("data/ptpl.rda")

# Wrapped Cauchy fit

?dwrpcauchy

#vector of angles in radians:

ptpl %>%
  drop_na(rel.angle) -> angle_data


# Estimate for complete pooling

mu0 <- circ.mean(angle_data$rel.angle)
rho0 <- est.rho(angle_data$rel.angle)
a <- wrpcauchy.ml(angle_data$rel.angle, mu0 = mu0, rho0 = rho0)





