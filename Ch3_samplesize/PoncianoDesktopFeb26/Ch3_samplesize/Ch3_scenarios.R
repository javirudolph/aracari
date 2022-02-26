
# Libraries -----------------------------------------------
set.seed(22227)

library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(fitdistrplus)
library(extRemes)


source("Ch3_samplesize/Ch3_functions.R")


# ORIGINAL SCENARIO --------------------------------------
dir_scenario <- "original"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

points_for_boxplots <- 50
B <- 1000

desired_means <- c(28, 32, 40, 50)
desired_sds <- c(49.7, 39.9, 33.3, 31.1)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("Ch3_samplesize/Ch3_process.R")
