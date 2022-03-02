# Get results for this on the allinone script



# Libraries -----------------------------------------------
set.seed(271367)

library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(fitdistrplus)
library(extRemes)


source("Ch3_samplesize/Ch3_functions.R")

# points_for_boxplots <- 50
# B <- 1000

points_for_boxplots <- 30
B <- 100



# SAME MEAN INCREASING SD ------------------------------
dir_scenario <- "ain1_samemean_upsd"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

desired_means <- c(30, 30, 30, 30)
desired_sds <- c(20, 30, 40, 50)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)

#INCREASING MEANS DECRESING SD --------------------------------------
dir_scenario <- "ain1_upmean_downsd"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

desired_means <- c(30, 35, 40, 50)
desired_sds <- c(50, 40, 35, 30)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)

# SAME SD INCREASING MEAN ------------------------------
dir_scenario <- "ain1_samesd_upmean"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

desired_means <- c(30, 40, 50, 60)
desired_sds <- c(25, 25, 25, 25)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)

# INCREASING SD INCREASING MEAN ------------------------------
dir_scenario <- "ain1_upsd_upmean"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

desired_means <- c(30, 40, 50, 60)
desired_sds <- c(20, 30, 40, 50)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)

# DECREASING BOTH ------------------------------
dir_scenario <- "ain1_downsd_downmean"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

desired_means <- c(60, 50, 40, 30)
desired_sds <- c(50, 40, 30, 20)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)


