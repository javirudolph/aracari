
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

# ORIGINAL SCENARIO --------------------------------------
dir_scenario <- "original"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

desired_means <- c(28, 32, 40, 50)
desired_sds <- c(49.7, 39.9, 33.3, 31.1)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("Ch3_samplesize/Ch3_process.R")



################################################################
### These are only sort of test scenarios
# Doing less simulations so they don't take as long

points_for_boxplots <- 30
B <- 100

# SAME MEAN INCREASING SD ------------------------------
dir_scenario <- "samemean_upsd"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

desired_means <- c(30, 30, 30, 30)
desired_sds <- c(20, 30, 40, 50)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("Ch3_samplesize/Ch3_process.R")

#INCREASING MEANS DECRESING SD --------------------------------------
dir_scenario <- "upmean_downsd"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

desired_means <- c(30, 35, 40, 50)
desired_sds <- c(50, 40, 35, 30)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("Ch3_samplesize/Ch3_process.R")

# SAME SD INCREASING MEAN ------------------------------
dir_scenario <- "samesd_upmean"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

desired_means <- c(30, 40, 50, 60)
desired_sds <- c(25, 25, 25, 25)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("Ch3_samplesize/Ch3_process.R")

# INCREASING SD INCREASING MEAN ------------------------------
dir_scenario <- "upsd_upmean"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

desired_means <- c(30, 40, 50, 60)
desired_sds <- c(20, 30, 40, 50)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("Ch3_samplesize/Ch3_process.R")

# DECREASING BOTH ------------------------------
dir_scenario <- "downsd_downmean"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

desired_means <- c(60, 50, 40, 30)
desired_sds <- c(50, 40, 30, 20)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("Ch3_samplesize/Ch3_process.R")


