# For the simulations, I want to consider four different scenarios

# Libraries -----------------------------------------------
set.seed(20220201)

library(aracari)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(fitdistrplus)

source("Ch3_samplesize/Ch3_functions.R")

# ORIGINAL
scenario <- "_original"
desired_means <- c(28, 32, 40, 50)
desired_sds <- c(49.7, 39.9, 33.3, 31.1)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("Ch3_samplesize/Ch3_simulation.R")

# 1) All groups have the same mean but different variance
scenario <- "_case1"
desired_means <- c(30, 30, 30, 30)
desired_sds <- c(20, 30, 40, 50)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("Ch3_samplesize/Ch3_simulation.R")

# 2) All groups have different means but the same variance
scenario <- "_case2"
desired_means <- c(30, 40, 50, 60)
desired_sds <- c(25, 25, 25, 25)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("Ch3_samplesize/Ch3_simulation.R")

# 3) Increasing mean and variance for groups
scenario <- "_case3"
desired_means <- c(30, 40, 50, 60)
desired_sds <- c(20, 30, 40, 50)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("Ch3_samplesize/Ch3_simulation.R")

# 4) Increasing mean and decreasing variance for groups
scenario <- "_case4"
desired_means <- c(30, 40, 50, 60)
desired_sds <- c(30, 25, 15, 10)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("Ch3_samplesize/Ch3_simulation.R")







